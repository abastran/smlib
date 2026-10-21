// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// libFuzzer target for the surface constructors and smlib_ssi: two surfaces of any kind, mostly
// sensible (unit-ish sizes, near the origin, so they tend to meet) with fuzzed perturbations,
// sometimes raw doubles; small NURBS (degree <= 3, <= 6 x 6 control points) so one input stays
// cheap. Any status is fine; a crash, hang, leak or exception is not, and neither is a result
// whose accessors disagree with it.
//
//   cmake -S capi -B build-fuzz -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -DSMLIB_CAPI_FUZZ=ON
//   cmake --build build-fuzz --target smlib_fuzz_ssi && build-fuzz/smlib_fuzz_ssi -timeout=10 -max_total_time=120

#include "smlib.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <utility>
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

    double raw()
    {
        double x = 0.0;
        size_t n = size < sizeof x ? size : sizeof x;
        std::memcpy(&x, data, n);
        data += n;
        size -= n;
        return x;
    }

    // Mostly a value in [lo, hi] from one byte; with `wild` set, sometimes any double at all.
    double in(double lo, double hi, bool wild)
    {
        uint8_t b = byte();
        if (wild && b == 255)
        {
            return raw();
        }
        return lo + (hi - lo) * (b / 254.0);
    }
};

void vec3(Input &in, double *v, double lo, double hi, bool wild)
{
    for (int i = 0; i < 3; ++i)
    {
        v[i] = in.in(lo, hi, wild);
    }
}

// Two axes: usually a rotated orthonormal frame, sometimes anything.
void frame(Input &in, double *x, double *y, bool wild)
{
    if (wild && (in.byte() & 7) == 0)
    {
        vec3(in, x, -1, 1, true);
        vec3(in, y, -1, 1, true);
        return;
    }
    const double a = in.in(0, 6.283, false), b = in.in(-1.5, 1.5, false), c = in.in(0, 6.283, false);
    // x = (cos a cos b, sin a cos b, sin b); y perpendicular to it, turned by c about x.
    const double xs[3] = {std::cos(a) * std::cos(b), std::sin(a) * std::cos(b), std::sin(b)};
    const double p[3] = {-std::sin(a), std::cos(a), 0.0};
    const double q[3] = {xs[1] * p[2] - xs[2] * p[1], xs[2] * p[0] - xs[0] * p[2], xs[0] * p[1] - xs[1] * p[0]};
    for (int i = 0; i < 3; ++i)
    {
        x[i] = xs[i];
        y[i] = std::cos(c) * p[i] + std::sin(c) * q[i];
    }
}

void range(Input &in, double *r, double lo, double hi, bool wild)
{
    r[0] = in.in(lo, hi, wild);
    r[1] = in.in(lo, hi, wild);
    if (r[0] > r[1] && (in.byte() & 3) != 0)
    {
        std::swap(r[0], r[1]);
    }
}

int32_t surface(Input &in, smlib_surface **out)
{
    const uint8_t kind = in.byte();
    const bool wild = (kind & 0x80) != 0;
    double o[3], x[3], y[3], r1[2], r2[2];
    vec3(in, o, -1.5, 1.5, wild);
    frame(in, x, y, wild);
    switch (kind % 6)
    {
    case 0:
        range(in, r1, -3, 3, wild);
        range(in, r2, -3, 3, wild);
        return smlib_surface_plane(o, x, y, r1, r2, out);
    case 1:
        range(in, r1, -360, 360, wild);
        range(in, r2, -3, 3, wild);
        return smlib_surface_cylinder(o, x, y, in.in(0, 2, wild), r1, r2, out);
    case 2:
        range(in, r1, 0, 3, wild);
        return smlib_surface_cone(o, x, y, in.in(0, 90, wild), r1, out);
    case 3:
        range(in, r1, -360, 360, wild);
        range(in, r2, -90, 90, wild);
        return smlib_surface_sphere(o, x, y, in.in(0, 2, wild), r1, r2, out);
    case 4:
        range(in, r1, -360, 360, wild);
        range(in, r2, -360, 360, wild);
        return smlib_surface_torus(o, x, y, in.in(0, 3, wild), in.in(0, 1.5, wild), r1, r2, out);
    default:
    {
        const uint8_t shape = in.byte();
        const int32_t p = 1 + shape % 3, q = 1 + (shape / 3) % 3;
        const int32_t nu = p + 1 + (shape / 9) % 3, nv = q + 1 + (shape / 27) % 3;
        // Clamped, mostly uniform knot vectors, with an interior knot moved now and then.
        auto knots = [&](int32_t n, int32_t d)
        {
            std::vector<double> k(n + d + 1);
            for (int32_t i = 0; i < n + d + 1; ++i)
            {
                k[i] = i <= d ? 0.0 : i >= n ? 1.0 : static_cast<double>(i - d) / (n - d);
            }
            if (in.byte() < 32)
            {
                k[in.byte() % k.size()] = in.in(-0.5, 1.5, wild);
            }
            return k;
        };
        std::vector<double> ku = knots(nu, p), kv = knots(nv, q);
        // A bumpy grid over [-1.5, 1.5]^2, tilted.
        std::vector<double> pts(3 * nu * nv), w(nu * nv);
        for (int32_t i = 0; i < nu; ++i)
        {
            for (int32_t j = 0; j < nv; ++j)
            {
                const double s = -1.5 + 3.0 * i / (nu - 1), t = -1.5 + 3.0 * j / (nv - 1), h = in.in(-1, 1, wild);
                double *P = &pts[3 * (i * nv + j)];
                for (int c = 0; c < 3; ++c)
                {
                    P[c] = o[c] + s * x[c] + t * y[c] + h * (x[(c + 1) % 3] * y[(c + 2) % 3] - x[(c + 2) % 3] * y[(c + 1) % 3]);
                }
                w[i * nv + j] = in.in(0.2, 3, wild);
            }
        }
        const bool rational = (shape & 0x80) != 0;
        int64_t nku = static_cast<int64_t>(ku.size()), nkv = static_cast<int64_t>(kv.size());
        if (wild && in.byte() == 0)
        {
            nku -= 1;
        }
        return smlib_surface_nurbs(p, q, nu, nv, ku.data(), nku, kv.data(), nkv, pts.data(), rational ? w.data() : nullptr, out);
    }
    }
}

} // namespace

extern "C" int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size)
{
    Input in{data, size};
    smlib_surface *a = nullptr, *b = nullptr;
    surface(in, &a);
    surface(in, &b);
    double dom[4], p[3];
    if (a != nullptr && smlib_surface_domain(a, dom) == SMLIB_OK)
    {
        smlib_surface_evaluate(a, in.in(dom[0], dom[1], true), in.in(dom[2], dom[3], true), p);
    }
    const uint8_t flags = in.byte();
    const int32_t method = (flags & 0x40) ? static_cast<int32_t>(in.byte()) - 1 : flags % 3;
    const double tol = (flags & 0x20) ? in.raw() : std::pow(10.0, -2.0 - (flags >> 2) % 4);
    const int32_t samples = (flags & 0x80) ? static_cast<int32_t>(in.byte()) - 1 : 2 + in.byte() % 32;
    smlib_ssi_result *r = nullptr;
    if (smlib_ssi(a, (flags & 0x10) ? a : b, method, tol, samples, &r) == SMLIB_OK)
    {
        int32_t n = 0;
        if (smlib_ssi_curve_count(r, &n) != SMLIB_OK || n < 0)
        {
            std::abort();
        }
        for (int32_t i = 0; i < n; ++i)
        {
            int32_t type = 0, fl = 0, ns = 0;
            double dev, ivl[2];
            if (smlib_ssi_curve_info(r, i, &type, &fl, &ns, &dev, ivl) != SMLIB_OK || ns != samples || type < SMLIB_SSI_NONE ||
                type > SMLIB_SSI_REGION_BOUNDARY)
            {
                std::abort();
            }
            std::vector<double> xyz(3 * ns), ua(2 * ns), ub(2 * ns), sa(3 * ns), sb(3 * ns);
            if (smlib_ssi_curve_samples(r, i, ns, xyz.data(), ua.data(), ub.data(), sa.data(), sb.data()) != SMLIB_OK)
            {
                std::abort();
            }
        }
        int32_t np = 0;
        smlib_ssi_point_count(r, &np);
        for (int32_t i = 0; i < np; ++i)
        {
            double q[3], u1[2], u2[2];
            if (smlib_ssi_point(r, i, q, u1, u2) != SMLIB_OK)
            {
                std::abort();
            }
        }
    }
    smlib_ssi_result_free(r);
    smlib_surface_free(a);
    smlib_surface_free(b);
    return 0;
}
