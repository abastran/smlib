// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// libFuzzer target for the operations (src/smlib_ops.cpp): booleans, fillet, chamfer and shell
// on boxes, with arbitrary sizes and selection points, NaNs, NULLs and out-of-range values, and
// chained so later operations see the bodies earlier ones made. Half the time a value is pulled
// into a range SMLib can build, so the kernel itself runs, not only the argument checks. A call
// may return any status; it must never crash, hang, leak or let an exception out.
//
//   cmake -S capi -B build-fuzz -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -DSMLIB_CAPI_FUZZ=ON
//   cmake --build build-fuzz --target smlib_fuzz_ops && build-fuzz/smlib_fuzz_ops -max_total_time=60

#include "smlib.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
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

    // Arbitrary half the time; otherwise in [lo, hi).
    double tame(double lo, double hi)
    {
        double x = real();
        if (byte() & 1)
        {
            return x;
        }
        double f = x == x && std::isfinite(x) ? std::fmod(std::fabs(x), 1.0) : 0.5;
        return lo + f * (hi - lo);
    }
};

// A point on the box [0, 20]^3: an edge midpoint, a face centre, or anything.
void point(Input &in, double p[3])
{
    uint8_t k = in.byte();
    switch (k % 3)
    {
    case 0:
    {
        // Edge midpoints: one coordinate 10, the others 0 or 20.
        int axis = (k >> 2) % 3;
        for (int i = 0; i < 3; ++i)
        {
            p[i] = i == axis ? 10.0 : ((k >> (4 + i)) & 1) * 20.0;
        }
        break;
    }
    case 1:
    {
        // Face centres.
        int axis = (k >> 2) % 3;
        for (int i = 0; i < 3; ++i)
        {
            p[i] = i == axis ? ((k >> 4) & 1) * 20.0 : 10.0;
        }
        break;
    }
    default:
        for (int i = 0; i < 3; ++i)
        {
            p[i] = in.real();
        }
    }
}

} // namespace

extern "C" int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size)
{
    Input in{data, size};
    const double o[3] = {0, 0, 0};
    const double s[3] = {20, 20, 20};
    smlib_brep *body = nullptr;
    if (smlib_make_box(o, s, &body) != SMLIB_OK)
    {
        return 0;
    }
    // A few operations per input keep each run short.
    for (int step = 0; step < 6 && in.size > 0; ++step)
    {
        smlib_brep *r = nullptr;
        int64_t failed = -1;
        double pts[3 * 4];
        int64_t n = static_cast<int64_t>(in.byte() % 6) - 1; // -1 .. 4
        for (int64_t i = 0; i < 4; ++i)
        {
            point(in, pts + 3 * i);
        }
        const double *p = in.byte() % 8 == 0 ? nullptr : pts;
        smlib_brep *target = in.byte() % 16 == 0 ? nullptr : body;
        int32_t st = SMLIB_OK;
        switch (in.byte() % 7)
        {
        case 0:
        case 1:
        {
            double org[3] = {in.tame(-5, 25), in.tame(-5, 25), in.tame(-5, 25)};
            double ext[3] = {in.tame(0.5, 30), in.tame(0.5, 30), in.tame(0.5, 30)};
            smlib_brep *tool = nullptr;
            if (smlib_make_box(org, ext, &tool) == SMLIB_OK)
            {
                int32_t which = static_cast<int32_t>(in.byte() % 5);
                // 0..3: a boolean (3 is refused); 4: the box imprinted on the body.
                st = which == 4 ? smlib_imprint(target, tool, &r) : smlib_boolean(target, tool, which, &r);
                smlib_brep_free(tool);
            }
            break;
        }
        case 2:
            st = smlib_fillet(target, p, n, in.tame(0.01, 1), in.tame(0.1, 12), &r, &failed);
            break;
        case 3:
            st = smlib_chamfer(target, p, n, in.tame(0.01, 1), in.tame(0.1, 12), &r, &failed);
            break;
        case 5:
        {
            double sizes[4] = {in.tame(0.1, 12), in.tame(0.1, 12), in.tame(0.1, 12), in.tame(0.1, 12)};
            st = smlib_blend(target, static_cast<int32_t>(in.byte() % 5), p, sizes, n, in.tame(0.01, 1), &r, &failed);
            break;
        }
        case 6:
        {
            const double knots[4] = {0, 0, 1, 1};
            double radii[2] = {in.tame(0, 12), in.tame(0, 12)};
            st = smlib_fillet_variable(target, p, in.tame(0.01, 1), pts + 3, 1, knots, 4, radii, 2, &r);
            break;
        }
        default:
            st = smlib_shell(target, p, n, in.tame(0.01, 1), in.tame(-12, 12), &r, &failed);
            break;
        }
        if (st == SMLIB_OK && r != nullptr)
        {
            int64_t solids = 0;
            smlib_brep_solid_count(r, &solids);
            // Keep the result for the next step half the time.
            if (in.byte() & 1)
            {
                smlib_brep_free(body);
                body = r;
            }
            else
            {
                smlib_brep_free(r);
            }
        }
    }
    smlib_brep_free(body);
    return 0;
}
