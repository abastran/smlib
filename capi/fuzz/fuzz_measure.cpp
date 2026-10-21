// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// libFuzzer target for smlib_brep_mass_properties and smlib_brep_area: arbitrary accuracies and
// NULLs, on primitives and sheets from 10 um to 5 m, at the origin and far from it, some with
// vertices pushed anywhere. Any status is fine; a crash, hang, leak, exception or non-finite result is not.
//
//   cmake -S capi -B build-fuzz -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -DSMLIB_CAPI_FUZZ=ON
//   cmake --build build-fuzz --target smlib_fuzz_measure && build-fuzz/smlib_fuzz_measure -max_total_time=60

#include "../src/smlib_shim.h"

#include <SmApiPrimitives.h>
#include <SmApiQueries.h>
#include <SmVertex.h>

#include <sanitizer/lsan_interface.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>

namespace
{

struct Input
{
    const uint8_t *data;
    size_t size;

    uint8_t byte()
    {
        if (size == 0)
        {
            return 0;
        }
        uint8_t b = *data;
        ++data;
        --size;
        return b;
    }

    double real()
    {
        double x = 0.0;
        size_t n = size < sizeof x ? size : sizeof x;
        std::memcpy(&x, data, n);
        data += n;
        size -= n;
        return x;
    }

    double place()
    {
        double x = real();
        return std::isfinite(x) ? std::fmod(x, 1e6) : 0.0;
    }
};

// The bodies, made once before fuzzing (LLVMFuzzerInitialize): primitive construction is
// SM_API's business, not this target's, and it leaks surface-cache trees and, for some
// degenerate sizes, curves. Each input works on a copy.
const SmTArray<SmBrep *> &bodies()
{
    static SmTArray<SmBrep *> all = []
    {
        SmTArray<SmBrep *> v;
        const SmVector3d at0(0, 0, 0), far(1e5, -2e5, 3e4);
        auto add = [&](SmApiStatus s, SmBrep *p)
        {
            if (s == SM_SUCCESS && p != nullptr)
            {
                v.Add(p);
            }
        };
        for (const SmVector3d &o : {at0, far})
        {
            for (double k : {1e-2, 1.0, 1e3})
            {
                SmBrep *p = nullptr;
                add(SmApiCreateBox(o, 4 * k, 3 * k, 2 * k, p), p);
                p = nullptr;
                add(SmApiCreateCylinder(o, k, 5 * k, p), p);
                p = nullptr;
                add(SmApiCreateSphere(o, k, p), p);
                p = nullptr;
                add(SmApiCreateCone(o, 2 * k, k, 3 * k, p), p);
                p = nullptr;
                add(SmApiCreateCone(o, 2 * k, 0.0, 3 * k, p), p);
                p = nullptr;
                add(SmApiCreateTorus(o, 3 * k, k, p), p);
                p = nullptr;
                add(SmApiCreatePlane(o, 4 * k, 3 * k, p), p);
                p = nullptr;
                add(SmApiCreatePlanarCircle(o, k, p), p);
            }
        }
        return v;
    }();
    return all;
}

} // namespace

extern "C" int LLVMFuzzerInitialize(int *, char ***)
{
    if (smlib_init() == SMLIB_OK)
    {
        // What construction leaks is not this target's finding.
        __lsan_disable();
        bodies();
        __lsan_enable();
    }
    return 0;
}

extern "C" int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size)
{
    Input in{data, size};
    if (smlib_init() != SMLIB_OK)
    {
        return 0;
    }
    const SmTArray<SmBrep *> &all = bodies();
    SmBrep *p = nullptr;
    if (all.GetSize() == 0 || SmApiBrepCopy(all[in.byte() % all.GetSize()], p) != SM_SUCCESS || p == nullptr)
    {
        delete p;
        return 0;
    }
    smlib_brep *brep = new smlib_brep{p};
    SmTArray<SmVertex *> vertices;
    p->GetVertices(vertices);
    for (uint8_t moves = in.byte() % 3; moves > 0 && vertices.GetSize() > 0; --moves)
    {
        SmVertex *v = vertices[in.byte() % vertices.GetSize()];
        SmPoint3d q = v->GetPoint();
        q.x += in.place();
        q.y += in.place();
        q.z += in.place();
        v->SetPoint(q);
    }
    for (int calls = 0; calls < 4 && in.size > 0; ++calls)
    {
        const uint8_t flags = in.byte();
        // Mostly sensible accuracies (fine ones are slow), sometimes anything at all.
        double accuracy = flags & 4 ? in.real() : std::pow(10.0, -1.0 - (in.byte() % 6));
        const smlib_brep *b = flags & 1 ? brep : nullptr;
        if (flags & 8)
        {
            double a = 0.0;
            if (smlib_brep_area(b, accuracy, flags & 2 ? &a : nullptr) == SMLIB_OK && !std::isfinite(a))
            {
                std::abort();
            }
        }
        else
        {
            double m[11];
            if (smlib_brep_mass_properties(b, accuracy, flags & 2 ? m : nullptr) == SMLIB_OK)
            {
                for (double x : m)
                {
                    if (!std::isfinite(x))
                    {
                        std::abort();
                    }
                }
            }
        }
    }
    smlib_brep_free(brep);
    return 0;
}
