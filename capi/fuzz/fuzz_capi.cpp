// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// libFuzzer target for the C API: every shim, fed arbitrary arguments (NaN, infinities, huge
// and tiny values, NULLs). A shim may return any status; it must never crash, hang, leak or let
// an exception out. Each shim function needs a case here.
//
//   cmake -S capi -B build-fuzz -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -DSMLIB_CAPI_FUZZ=ON
//   cmake --build build-fuzz --target smlib_fuzz && build-fuzz/smlib_fuzz -max_total_time=60
//
// SMLIB_FUZZ_SMB=<native .smb> adds that solid as a second starting point for the flat-array
// case (a box has only analytic faces; a blend-filleted box adds NURBS faces and UV trims).

#include "smlib.h"

#include <cstddef>
#include <cstdint>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <vector>

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
        switch (in.byte() % 7)
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
        case 4:
        {
            int32_t m = 0;
            smlib_brep_is_manifold_solid(in.byte() & 1 ? box : nullptr, &m);
            break;
        }
        case 5:
        {
            // Export the body, mutate its arrays from the bytes, import the result.
            // With SMLIB_FUZZ_SMB set, half the time start from that solid (NURBS faces, UV trims).
            static smlib_brep *seed = []
            {
                smlib_brep *b = nullptr;
                const char *path = std::getenv("SMLIB_FUZZ_SMB");
                if (path != nullptr)
                {
                    smlib_brep_read(path, 0, &b);
                }
                return b;
            }();
            smlib_flat *flat = nullptr;
            if (smlib_brep_export(seed != nullptr && (in.byte() & 1) ? seed : box, &flat) != SMLIB_OK)
            {
                flat = smlib_flat_new();
            }
            int edits = in.byte() % 8;
            for (int e = 0; e < edits && in.size > 0; ++e)
            {
                int32_t id = static_cast<int32_t>(in.byte() % 64) + ((in.byte() & 1) ? 100 : 0);
                uint8_t how = in.byte();
                const int64_t *pi = nullptr;
                const double *pf = nullptr;
                int64_t n = 0;
                if (smlib_flat_get_i64(flat, id, &pi, &n) == SMLIB_OK)
                {
                    std::vector<int64_t> v(pi, pi + n);
                    if (how % 3 == 0 && !v.empty())
                    {
                        v.resize(static_cast<size_t>(in.byte()) % v.size());
                    }
                    else if (how % 3 == 1)
                    {
                        v.push_back(static_cast<int64_t>(in.real()));
                    }
                    else if (!v.empty())
                    {
                        int64_t x = 0;
                        double r = in.real();
                        std::memcpy(&x, &r, sizeof x);
                        v[static_cast<size_t>(in.byte()) % v.size()] = (how & 8) ? x : x % 16;
                    }
                    smlib_flat_set_i64(flat, id, v.data(), static_cast<int64_t>(v.size()));
                }
                else if (smlib_flat_get_f64(flat, id, &pf, &n) == SMLIB_OK)
                {
                    std::vector<double> v(pf, pf + n);
                    if (how % 3 == 0 && !v.empty())
                    {
                        v.resize(static_cast<size_t>(in.byte()) % v.size());
                    }
                    else if (how % 3 == 1)
                    {
                        v.push_back(in.real());
                    }
                    else if (!v.empty())
                    {
                        double &x = v[static_cast<size_t>(in.byte()) % v.size()];
                        x = (how & 8) ? in.real() : x + static_cast<double>(static_cast<int8_t>(in.byte())) / 16.0;
                    }
                    smlib_flat_set_f64(flat, id, v.data(), static_cast<int64_t>(v.size()));
                }
                else
                {
                    smlib_flat_set_i64(flat, id, nullptr, 0);
                }
            }
            smlib_brep *imported = nullptr;
            if (smlib_brep_import(in.byte() & 7 ? flat : nullptr, &imported) == SMLIB_OK)
            {
                double v = 0.0;
                smlib_brep_volume(imported, 1e-3, &v);
                smlib_brep_free(imported);
            }
            smlib_flat_free(flat);
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
