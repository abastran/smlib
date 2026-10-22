// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// libFuzzer target for the sweeps (src/smlib_sweeps.cpp): curve lists made from the input, half
// of them well-formed loops (rectangles, discs, polygons and NURBS pieces of random sizes, places
// and orientations) so SMLib's sweeps run, the rest arbitrary kinds and numbers; each one through
// smlib_extrude, smlib_revolve, smlib_sweep, smlib_loft, smlib_loft_to_point, smlib_screw_sweep and
// smlib_sweep_stations with arbitrary or near-sane options. Any status is fine; a crash, hang, leak or
// exception is not.
//
//   cmake -S capi -B build-fuzz -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -DSMLIB_CAPI_FUZZ=ON -DSMLIB_CAPI_TESTS=OFF
//   cmake --build build-fuzz --target smlib_fuzz_sweeps && build-fuzz/smlib_fuzz_sweeps -max_total_time=120

#include "smlib.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
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

    // A number in [lo, hi], from two bytes.
    double in(double lo, double hi)
    {
        double f = (byte() * 256.0 + byte()) / 65535.0;
        return lo + (hi - lo) * f;
    }
};

struct List
{
    std::vector<int64_t> kinds;
    std::vector<double> reals;
};

// A closed polygon of `n` corners in the plane z = `z`, around (cx, cy) at radius `r`.
void polygon(List &l, int n, double cx, double cy, double z, double r)
{
    for (int i = 0; i < n; ++i)
    {
        double a0 = 2 * M_PI * i / n, a1 = 2 * M_PI * (i + 1) / n;
        l.kinds.push_back(SMLIB_CURVE_LINE);
        l.reals.insert(l.reals.end(), {cx + r * std::cos(a0), cy + r * std::sin(a0), z, cx + r * std::cos(a1), cy + r * std::sin(a1), z});
    }
}

void disc(List &l, double cx, double cy, double z, double r)
{
    l.kinds.push_back(SMLIB_CURVE_ARC);
    l.reals.insert(l.reals.end(), {cx, cy, z, 0, 0, 1, 1, 0, 0, r, 0, 2 * M_PI});
}

// A loop of one closed cubic NURBS (a rounded square) in the plane z = `z`.
void blob(List &l, double cx, double cy, double z, double r)
{
    l.kinds.insert(l.kinds.end(), {SMLIB_CURVE_NURBS, 3, 7});
    const double p[7][2] = {{1, 0}, {1, 1}, {-1, 1}, {-1, 0}, {-1, -1}, {1, -1}, {1, 0}};
    for (auto &q : p)
    {
        l.reals.insert(l.reals.end(), {cx + r * q[0], cy + r * q[1], z, 1.0});
    }
    l.reals.insert(l.reals.end(), {0, 0, 0, 0, 0.25, 0.5, 0.75, 1, 1, 1, 1});
}

// A well-formed profile in z = `z`, or arbitrary kinds and numbers.
List profile(Input &in, double z)
{
    List l;
    uint8_t pick = in.byte();
    double r = in.in(1e-3, 50.0);
    double cx = in.in(-100.0, 100.0), cy = in.in(-100.0, 100.0);
    switch (pick % 5)
    {
    case 0:
        polygon(l, 3 + in.byte() % 6, cx, cy, z, r);
        break;
    case 1:
        disc(l, cx, cy, z, r);
        break;
    case 2:
        blob(l, cx, cy, z, r);
        break;
    case 3:
        disc(l, cx, cy, z, r);
        polygon(l, 3 + in.byte() % 4, cx, cy, z, r * in.in(0.01, 1.2));
        break;
    default:
    {
        int n = in.byte() % 12;
        for (int i = 0; i < n; ++i)
        {
            l.kinds.push_back(static_cast<int64_t>(static_cast<int8_t>(in.byte())));
        }
        int m = in.byte() % 64;
        for (int i = 0; i < m; ++i)
        {
            l.reals.push_back(in.real());
        }
    }
    }
    // Now and then, an edit that breaks it.
    if (in.byte() % 7 == 0 && !l.reals.empty())
    {
        l.reals[in.byte() % l.reals.size()] = in.real();
    }
    return l;
}

void drop(int32_t s, smlib_brep *b)
{
    if (s == SMLIB_OK)
    {
        smlib_brep_free(b);
    }
}

} // namespace

extern "C" int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size)
{
    static bool ready = smlib_init() == SMLIB_OK;
    if (!ready)
    {
        return 0;
    }
    Input in{data, size};
    List p = profile(in, 0.0);
    const double *pr = p.reals.empty() ? nullptr : p.reals.data();
    const int64_t *pk = p.kinds.empty() ? nullptr : p.kinds.data();
    int64_t nk = static_cast<int64_t>(p.kinds.size()), nr = static_cast<int64_t>(p.reals.size());
    double dir[3] = {in.in(-1, 1), in.in(-1, 1), in.in(-1, 1)};
    if (in.byte() % 2 == 0)
    {
        dir[0] = dir[1] = 0.0;
        dir[2] = 1.0;
    }
    int32_t cap = in.byte() % 3 == 0 ? 0 : 1;
    smlib_brep *b = nullptr;
    switch (in.byte() % 7)
    {
    case 0:
    {
        // 3 is no corner rule: refused
        int32_t corner = static_cast<int32_t>(in.byte() % 4);
        drop(smlib_extrude(pk, nk, pr, nr, dir, in.in(-1.0, 200.0), in.byte() % 2 ? 0.0 : in.in(-1.7, 1.7), corner, cap, &b, nullptr), b);
        break;
    }
    case 1:
    {
        double o[3] = {in.in(-200, 200), in.in(-200, 200), 0.0};
        double ax[3] = {0.0, 1.0, 0.0};
        drop(smlib_revolve(pk, nk, pr, nr, o, ax, in.in(-0.1, 7.0), cap, &b), b);
        break;
    }
    case 2:
    {
        List path;
        double len = in.in(-1.0, 300.0);
        if (in.byte() % 2)
        {
            path.kinds.push_back(SMLIB_CURVE_LINE);
            path.reals.insert(path.reals.end(), {0, 0, 0, 0, in.in(-5, 5), len});
        }
        else
        {
            path.kinds.push_back(SMLIB_CURVE_ARC);
            double rr = in.in(1e-3, 300.0);
            path.reals.insert(path.reals.end(), {-rr, 0, 0, 0, -1, 0, 1, 0, 0, rr, 0, in.in(-0.5, 6.5)});
        }
        drop(smlib_sweep(pk, nk, pr, nr, path.kinds.data(), static_cast<int64_t>(path.kinds.size()), path.reals.data(),
                         static_cast<int64_t>(path.reals.size()), static_cast<int32_t>(in.byte() % 3), cap, &b),
             b);
        break;
    }
    case 3:
    {
        int sections = 2 + in.byte() % 3;
        List all = p;
        std::vector<int64_t> sizes{static_cast<int64_t>(p.kinds.size())};
        for (int i = 1; i < sections; ++i)
        {
            List q = profile(in, in.in(-50, 50));
            sizes.push_back(static_cast<int64_t>(q.kinds.size()));
            all.kinds.insert(all.kinds.end(), q.kinds.begin(), q.kinds.end());
            all.reals.insert(all.reals.end(), q.reals.begin(), q.reals.end());
        }
        drop(smlib_loft(all.kinds.empty() ? nullptr : all.kinds.data(), static_cast<int64_t>(all.kinds.size()),
                        all.reals.empty() ? nullptr : all.reals.data(), static_cast<int64_t>(all.reals.size()), sizes.data(),
                        static_cast<int64_t>(sizes.size()), static_cast<int32_t>(in.byte() % 5), cap, &b),
             b);
        break;
    }
    case 4:
    {
        // Stations up a line, turning, scaling and shearing by arbitrary amounts (some singular),
        // or wholly arbitrary maps.
        int n = 1 + in.byte() % 40;
        bool sane = in.byte() % 2 == 0;
        double twist = in.in(-30.0, 30.0), scale = in.in(-2.0, 50.0), shear = in.in(-3.0, 3.0), rise = in.in(-1.0, 300.0);
        std::vector<double> maps;
        for (int i = 0; i < n; ++i)
        {
            double f = n > 1 ? static_cast<double>(i) / (n - 1) : 0.0;
            if (sane)
            {
                double a = twist * f, k = 1.0 + (scale - 1.0) * f;
                maps.insert(maps.end(), {k * std::cos(a), k * std::sin(a), 0, -k * std::sin(a), k * std::cos(a), 0, shear * f, 0, 1, 0, 0, rise * f});
            }
            else
            {
                for (int j = 0; j < 12; ++j)
                {
                    maps.push_back(in.in(-1e3, 1e3));
                }
            }
        }
        drop(smlib_sweep_stations(pk, nk, pr, nr, maps.data(), n, static_cast<int32_t>(in.byte() % 4), &b), b);
        break;
    }
    case 5:
    {
        double apex[3] = {in.in(-200, 200), in.in(-200, 200), in.byte() % 4 == 0 ? 0.0 : in.in(-200, 200)};
        drop(smlib_loft_to_point(pk, nk, pr, nr, apex, &b), b);
        break;
    }
    default:
    {
        double o[3] = {in.in(-200, 200), in.in(-200, 200), in.in(-200, 200)};
        drop(smlib_screw_sweep(pk, nk, pr, nr, o, dir, in.in(-20.0, 20.0), in.in(0.0, 6.0), in.in(1e-7, 1.0), &b), b);
    }
    }
    return 0;
}
