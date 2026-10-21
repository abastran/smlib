// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// The C API in include/smlib.h. Every entry point checks its arguments, then runs under the
// library lock inside a Call, which catches every C++ exception, logs one line, records the call
// in flight for crash reports and, when the call fails, dumps its inputs for replay.

#include "smlib.h"

#include <SmApiGeneral.h>
#include <SmApiImportExport.h>
#include <SmApiPrimitives.h>
#include <SmApiQueries.h>
#include <SmBrep.h>
#include <SmVector3d.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

struct smlib_brep
{
    SmBrep *p;
};

namespace
{

std::mutex &lock()
{
    static std::mutex m;
    return m;
}

// Diagnostics settings; changed and read under the lock.
smlib_log_fn g_log = nullptr;
void *g_log_user = nullptr;
std::string g_dump_dir;
unsigned long g_dump_seq = 0;

// The call in flight, readable without the lock by a crash handler. g_current_len is 0 while no
// call runs; the text is written before the length is published.
constexpr size_t CURRENT_CAP = 512;
char g_current[CURRENT_CAP];
std::atomic<size_t> g_current_len{0};

int32_t status(SmApiStatus s)
{
    return s == SM_SUCCESS ? SMLIB_OK : static_cast<int32_t>(s);
}

bool finite(const double *v, int n)
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

std::string num(double x)
{
    char b[32];
    std::snprintf(b, sizeof b, "%.17g", x);
    return b;
}

std::string vec(const double *v, int n)
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

std::string handle(const smlib_brep *b)
{
    char t[32];
    std::snprintf(t, sizeof t, "%p", static_cast<const void *>(b));
    return b == nullptr ? std::string("null") : std::string("brep@") + t;
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
                s = SmApiGetOrCreateContext() == nullptr ? SMLIB_ERR_EXCEPTION : f();
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
            std::fclose(f);
        }
        for (auto &in : m_inputs)
        {
            std::string path = base + "-" + in.first + ".smb";
            SmApiWriteBrepToFile(in.second, path.c_str(), FALSE);
        }
    }

    const char *m_fn;
    std::string m_line;
    std::vector<std::pair<const char *, const SmBrep *>> m_inputs;
};

// Hands `p` to the caller, or deletes it when it is no usable body.
int32_t take(SmBrep *p, smlib_brep **out)
{
    if (p == nullptr || p->GetNumFaces() == 0)
    {
        delete p;
        return SMLIB_ERR_EMPTY;
    }
    *out = new smlib_brep{p};
    return SMLIB_OK;
}

bool valid(const smlib_brep *b)
{
    return b != nullptr && b->p != nullptr;
}

} // namespace

extern "C" int32_t smlib_abi_version(void)
{
    return SMLIB_ABI_VERSION;
}

extern "C" int32_t smlib_init(void)
{
    return Call("smlib_init", "").run([] { return SMLIB_OK; });
}

extern "C" void smlib_set_log(smlib_log_fn log, void *user)
{
    try
    {
        std::lock_guard<std::mutex> g(lock());
        g_log = log;
        g_log_user = user;
    }
    catch (...)
    {
    }
}

extern "C" int32_t smlib_set_dump_dir(const char *dir)
{
    try
    {
        std::lock_guard<std::mutex> g(lock());
        g_dump_dir = dir == nullptr ? std::string() : std::string(dir);
        return SMLIB_OK;
    }
    catch (...)
    {
        return SMLIB_ERR_EXCEPTION;
    }
}

extern "C" int64_t smlib_current_call(char *buf, int64_t cap)
{
    size_t n = g_current_len.load(std::memory_order_acquire);
    if (buf == nullptr || cap <= 0)
    {
        return static_cast<int64_t>(n);
    }
    size_t m = std::min(n, static_cast<size_t>(cap - 1));
    std::memcpy(buf, g_current, m);
    buf[m] = '\0';
    return static_cast<int64_t>(m);
}

extern "C" void smlib_brep_free(smlib_brep *brep)
{
    if (brep == nullptr)
    {
        return;
    }
    try
    {
        std::lock_guard<std::mutex> g(lock());
        delete brep->p;
    }
    catch (...)
    {
    }
    delete brep;
}

extern "C" int32_t smlib_brep_read(const char *path, int32_t ascii, smlib_brep **out)
{
    if (path == nullptr || out == nullptr)
    {
        return SMLIB_ERR_ARG;
    }
    return Call("smlib_brep_read", std::string("path=") + path).run(
        [&]
        {
            SmBrep *r = nullptr;
            int32_t s = status(SmApiReadBrepFromFile(path, ascii ? TRUE : FALSE, FALSE, r));
            if (s != SMLIB_OK)
            {
                delete r;
                return s;
            }
            return take(r, out);
        }
    );
}

extern "C" int32_t smlib_brep_write(const smlib_brep *brep, const char *path, int32_t ascii)
{
    if (!valid(brep) || path == nullptr)
    {
        return SMLIB_ERR_ARG;
    }
    Call c("smlib_brep_write", handle(brep) + ", path=" + path);
    c.input("brep", brep);
    return c.run([&] { return status(SmApiWriteBrepToFile(brep->p, path, ascii ? TRUE : FALSE)); });
}

extern "C" int32_t smlib_brep_copy(const smlib_brep *brep, smlib_brep **out)
{
    if (!valid(brep) || out == nullptr)
    {
        return SMLIB_ERR_ARG;
    }
    Call c("smlib_brep_copy", handle(brep));
    c.input("brep", brep);
    return c.run(
        [&]
        {
            SmBrep *r = nullptr;
            int32_t s = status(SmApiBrepCopy(brep->p, r));
            if (s != SMLIB_OK)
            {
                delete r;
                return s;
            }
            return take(r, out);
        }
    );
}

extern "C" int32_t smlib_make_box(const double origin[3], const double size[3], smlib_brep **out)
{
    if (origin == nullptr || size == nullptr || out == nullptr || !finite(origin, 3) || !finite(size, 3) || !(size[0] > 0.0) ||
        !(size[1] > 0.0) || !(size[2] > 0.0))
    {
        return SMLIB_ERR_ARG;
    }
    return Call("smlib_make_box", "origin=" + vec(origin, 3) + ", size=" + vec(size, 3)).run(
        [&]
        {
            SmBrep *r = nullptr;
            int32_t s = status(SmApiCreateBox(SmVector3d(origin[0], origin[1], origin[2]), size[0], size[1], size[2], r));
            if (s != SMLIB_OK)
            {
                delete r;
                return s;
            }
            return take(r, out);
        }
    );
}

extern "C" int32_t smlib_brep_volume(const smlib_brep *brep, double relative_accuracy, double *out)
{
    if (!valid(brep) || out == nullptr || !std::isfinite(relative_accuracy))
    {
        return SMLIB_ERR_ARG;
    }
    Call c("smlib_brep_volume", handle(brep) + ", relative_accuracy=" + num(relative_accuracy));
    c.input("brep", brep);
    return c.run(
        [&]
        {
            double v = 0.0;
            int32_t s = status(SmApiBrepComputeVolume(brep->p, relative_accuracy, v));
            if (s == SMLIB_OK)
            {
                *out = v;
            }
            return s;
        }
    );
}

extern "C" int32_t smlib_brep_counts(const smlib_brep *brep, int64_t out[3])
{
    if (!valid(brep) || out == nullptr)
    {
        return SMLIB_ERR_ARG;
    }
    return Call("smlib_brep_counts", handle(brep)).run(
        [&]
        {
            out[0] = static_cast<int64_t>(brep->p->GetNumFaces());
            out[1] = static_cast<int64_t>(brep->p->GetNumEdges());
            out[2] = static_cast<int64_t>(brep->p->GetNumVertices());
            return SMLIB_OK;
        }
    );
}
