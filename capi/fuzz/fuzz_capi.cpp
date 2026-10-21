// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// libFuzzer target for the C API: every shim, fed arbitrary arguments (NaN, infinities, huge
// and tiny values, NULLs). A shim may return any status; it must never crash, hang, leak or let
// an exception out. Each shim function needs a case here.
//
//   cmake -S capi -B build-fuzz -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -DSMLIB_CAPI_FUZZ=ON
//   cmake --build build-fuzz --target smlib_fuzz && build-fuzz/smlib_fuzz -max_total_time=60

#include "smlib.h"

#include <cstddef>
#include <cstdint>
#include <cmath>
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
};

} // namespace

extern "C" int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size)
{
    Input in{data, size};
    smlib_brep *box = nullptr;
    double origin[3] = {in.real(), in.real(), in.real()};
    double extent[3] = {in.real(), in.real(), in.real()};
    // Half the time, a box that is valid whatever the bytes say, so later shims see a real body.
    if (in.byte() & 1)
    {
        for (double &x : extent)
        {
            x = 1.0 + (x == x ? std::fmod(x < 0 ? -x : x, 1000.0) : 0.0);
        }
        for (double &x : origin)
        {
            x = x == x ? std::fmod(x, 1e6) : 0.0;
        }
    }
    smlib_make_box(in.byte() & 1 ? origin : nullptr, extent, &box);
    while (in.size > 0)
    {
        switch (in.byte() % 4)
        {
        case 0:
        {
            double v = 0.0;
            smlib_brep_volume(in.byte() & 1 ? box : nullptr, in.real(), &v);
            break;
        }
        case 1:
        {
            int64_t c[3];
            smlib_brep_counts(in.byte() & 1 ? box : nullptr, c);
            break;
        }
        case 2:
        {
            smlib_brep *copy = nullptr;
            if (smlib_brep_copy(in.byte() & 1 ? box : nullptr, &copy) == SMLIB_OK)
            {
                smlib_brep_free(box);
                box = copy;
            }
            break;
        }
        default:
        {
            char buf[64];
            smlib_current_call(in.byte() & 1 ? buf : nullptr, static_cast<int64_t>(in.byte()) - 8);
            break;
        }
        }
    }
    smlib_brep_free(box);
    return 0;
}
