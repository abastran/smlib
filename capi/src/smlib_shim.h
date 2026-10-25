// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// What every shim file shares: the handle type, the one library lock, and Call, which runs an
// entry point under that lock, catches every C++ exception, logs one line, records the call in
// flight for crash reports and, when the call fails, dumps its inputs for replay. Internal to the
// library; include/smlib.h is the only public header.

#ifndef SMLIB_SHIM_H
#define SMLIB_SHIM_H

#include "smlib.h"

#include <SmApiGeneral.h>
#include <SmApiImportExport.h>
#include <SmBrep.h>
#include <SmExtent2d.h>
#include <SmMessages.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <functional>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

struct smlib_brep
{
    SmBrep *p;
};

// A surface with its UV domain (src/smlib_ssi.cpp makes them; src/smlib_split.cpp cuts with them).
struct smlib_surface
{
    SmSurface *s = nullptr;
    SmExtent2d domain;
    // The description in the SSI corpus case format, for dumps and the log.
    std::string json;

    // In src/smlib_ssi.cpp, where SmSurface is complete.
    ~smlib_surface();
};

namespace smlib_shim
{

inline std::mutex &lock()
{
    static std::mutex m;
    return m;
}

// Diagnostics settings; changed and read under the lock.
inline smlib_log_fn g_log = nullptr;
inline void *g_log_user = nullptr;
inline std::string g_dump_dir;
inline unsigned long g_dump_seq = 0;

// The call in flight, readable without the lock by a crash handler. g_current_len is 0 while no
// call runs; the text is written before the length is published.
inline constexpr size_t CURRENT_CAP = 512;
inline char g_current[CURRENT_CAP];
inline std::atomic<size_t> g_current_len{0};

// SMLib's status as the C API's: SMLIB_OK for success; an offset whose walls would meet
// (SM_ERR_OFFSET_WALLS_TOUCH, from any offset: a shell, a solid offset, a face move)
// SMLIB_ERR_WALLS_TOUCH; any other status unchanged.
inline int32_t status(SmApiStatus s)
{
    switch (s)
    {
    case SM_SUCCESS:
        return SMLIB_OK;
    case SM_ERR_OFFSET_WALLS_TOUCH:
        return SMLIB_ERR_WALLS_TOUCH;
    default:
        return static_cast<int32_t>(s);
    }
}

inline bool finite(const double *v, int n)
{
    for (int i = 0; i < n; ++i)
    {
        if (!std::isfinite(v[i]))
        {
            return false;
        }
    }
    return true;
}

inline std::string num(double x)
{
    char b[32];
    std::snprintf(b, sizeof b, "%.17g", x);
    return b;
}

inline std::string vec(const double *v, int n)
{
    if (v == nullptr)
    {
        return "null";
    }
    std::string s = "[";
    for (int i = 0; i < n; ++i)
    {
        s += (i ? "," : "") + num(v[i]);
    }
    return s + "]";
}

inline std::string handle(const smlib_brep *b)
{
    char t[32];
    std::snprintf(t, sizeof t, "%p", static_cast<const void *>(b));
    return b == nullptr ? std::string("null") : std::string("brep@") + t;
}

// SMLib text (TCHAR is char here and wchar_t on Windows) as ASCII, at most `cap` characters:
// anything else becomes '?'.
inline std::string ascii(const TCHAR *t, size_t cap)
{
    std::string out;
    for (size_t i = 0; t != nullptr && i < cap && t[i] != 0; ++i)
    {
        const long c = static_cast<long>(t[i]);
        out += c > 0 && c < 128 ? static_cast<char>(c) : '?';
    }
    return out;
}

// SMLib's own reports (smos_ErrorMessage: warnings, messages, and the errors an SER logs), as
// one line of the shim log, tagged with the call in flight: "SMLib <kind> at <file>:<line>:
// <message> (in <call>)". SMLib keeps its error callback per thread (SmThreadLocalStorage), so
// Call::run installs it on the calling thread; it runs under the library lock, like the call.
inline void on_smlib_message(SmStatus error, const TCHAR *const file, ULONG line, const TCHAR *const message)
{
    if (g_log == nullptr)
    {
        return;
    }
    auto text = [](const TCHAR *t) { return ascii(t, CURRENT_CAP); };
    std::string kind = error == SM_ERR_WARNING ? "warning" : error == SM_ERR_MESSAGE ? "message" : "error " + std::to_string(error);
    std::string line_text = "SMLib " + kind + " at " + text(file) + ":" + std::to_string(line);
    std::string what = text(message);
    if (!what.empty())
    {
        line_text += ": " + what;
    }
    size_t n = g_current_len.load(std::memory_order_acquire);
    if (n > 0)
    {
        line_text += " (in " + std::string(g_current, n) + ")";
    }
    g_log(g_log_user, 1, line_text.c_str());
}

// One API call: its line for the log, the dump and crash reports, and its B-rep inputs.
class Call
{
public:
    Call(const char *fn, std::string args) : m_fn(fn), m_line(std::string(fn) + "(" + args + ")")
    {
    }

    void input(const char *name, const smlib_brep *b)
    {
        if (b != nullptr && b->p != nullptr)
        {
            m_inputs.emplace_back(name, b->p);
        }
    }

    // A non-B-rep input, dumped by `write` (given the file path) as `<n>-<function>-<name>.txt`
    // when the call fails.
    void input_text(const char *name, std::function<void(const std::string &path)> write)
    {
        m_texts.emplace_back(name, std::move(write));
    }

    // Runs `f` under the lock, with SMLib's context created, and returns its status.
    template <typename F>
    int32_t run(F &&f)
    {
        int32_t s = SMLIB_ERR_EXCEPTION;
        try
        {
            std::lock_guard<std::mutex> g(lock());
            publish();
            auto t0 = std::chrono::steady_clock::now();
            try
            {
                if (SmApiGetOrCreateContext() == nullptr)
                {
                    s = SMLIB_ERR_EXCEPTION;
                }
                else
                {
                    smos_SetErrorCallback(on_smlib_message);
                    s = f();
                }
            }
            catch (...)
            {
                s = SMLIB_ERR_EXCEPTION;
            }
            double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
            finish(s, ms);
        }
        catch (...)
        {
            // The lock itself or the logging failed; the status stands.
        }
        g_current_len.store(0, std::memory_order_release);
        return s;
    }

private:
    void publish()
    {
        size_t n = std::min(m_line.size(), CURRENT_CAP - 1);
        std::memcpy(g_current, m_line.data(), n);
        g_current[n] = '\0';
        g_current_len.store(n, std::memory_order_release);
    }

    void finish(int32_t s, double ms)
    {
        char tail[64];
        std::snprintf(tail, sizeof tail, " -> %d (%.3f ms)", static_cast<int>(s), ms);
        std::string line = m_line + tail;
        if (g_log != nullptr)
        {
            g_log(g_log_user, s == SMLIB_OK ? 0 : 1, line.c_str());
        }
        if (s != SMLIB_OK && !g_dump_dir.empty())
        {
            dump(line);
        }
    }

    // Writes the call's line and its B-rep inputs into the dump directory. Best effort: a dump
    // that cannot be written never changes the call's status.
    void dump(const std::string &line)
    {
        char stem[64];
        std::snprintf(stem, sizeof stem, "%lu-%s", ++g_dump_seq, m_fn);
        std::string base = g_dump_dir + "/" + stem;
        if (FILE *f = std::fopen((base + ".txt").c_str(), "w"))
        {
            std::fprintf(f, "%s\n", line.c_str());
            for (auto &in : m_inputs)
            {
                std::fprintf(f, "%s: %s-%s.smb\n", in.first, stem, in.first);
            }
            for (auto &in : m_texts)
            {
                std::fprintf(f, "%s: %s-%s.txt\n", in.first, stem, in.first);
            }
            std::fclose(f);
        }
        for (auto &in : m_inputs)
        {
            std::string path = base + "-" + in.first + ".smb";
            (void)SmApiWriteBrepToFile(in.second, path.c_str(), FALSE); // best effort, as above
        }
        for (auto &in : m_texts)
        {
            in.second(base + "-" + in.first + ".txt");
        }
    }

    const char *m_fn;
    std::string m_line;
    std::vector<std::pair<const char *, const SmBrep *>> m_inputs;
    std::vector<std::pair<const char *, std::function<void(const std::string &)>>> m_texts;
};


// Hands `p` to the caller, or deletes it when it is no usable body.
inline int32_t take(SmBrep *p, smlib_brep **out)
{
    if (p == nullptr || p->GetNumFaces() == 0)
    {
        delete p;
        return SMLIB_ERR_EMPTY;
    }
    *out = new smlib_brep{p};
    return SMLIB_OK;
}

inline bool valid(const smlib_brep *b)
{
    return b != nullptr && b->p != nullptr;
}

} // namespace smlib_shim

#endif // SMLIB_SHIM_H
