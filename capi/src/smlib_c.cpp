// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// The C API in include/smlib.h. Every entry point checks its arguments, then runs under the
// library lock inside a Call (smlib_shim.h), which catches every C++ exception, logs one line,
// records the call in flight for crash reports and, when the call fails, dumps its inputs for
// replay.

#include "smlib_shim.h"

#include <SmApiPrimitives.h>
#include <SmApiQueries.h>
#include <SmVector3d.h>

using namespace smlib_shim;

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
