// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// smlib_surface_* and smlib_ssi on cases from the SSI corpus (plane x cylinder ellipse, plane x
// sphere circle, two NURBS paraboloids, a sphere touching a plane), checking curve counts, types
// and that each curve, evaluated from its exact B-spline, lies on both surfaces; then hostile
// inputs, which must be rejected.

#include "smlib.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <vector>

static int failures = 0;

#define CHECK(cond)                                                                                \
    do                                                                                             \
    {                                                                                              \
        if (!(cond))                                                                               \
        {                                                                                          \
            std::fprintf(stderr, "%s:%d: CHECK(%s) failed\n", __FILE__, __LINE__, #cond);          \
            failures++;                                                                            \
        }                                                                                          \
    } while (0)

static const double NaN = std::numeric_limits<double>::quiet_NaN();
static const double INF = std::numeric_limits<double>::infinity();
static const double X[3] = {1, 0, 0}, Y[3] = {0, 1, 0}, O[3] = {0, 0, 0};

// A (rational) B-spline as smlib_ssi_curve_nurbs gives it.
struct Nurbs
{
    int32_t status = SMLIB_ERR_EMPTY, degree = 0;
    std::vector<double> points, weights, knots;
};

static Nurbs nurbs(const smlib_ssi_result *r, int32_t i, int32_t which, int dim)
{
    Nurbs c;
    int64_t np = 0, nk = 0;
    c.status = smlib_ssi_curve_nurbs(r, i, which, &c.degree, &np, &nk, 0, nullptr, nullptr, 0, nullptr);
    if (c.status != SMLIB_OK)
    {
        return c;
    }
    CHECK(np > c.degree && nk == np + c.degree + 1);
    c.points.resize(dim * np);
    c.weights.resize(np);
    c.knots.resize(nk);
    CHECK(smlib_ssi_curve_nurbs(r, i, which, nullptr, nullptr, nullptr, np - 1, c.points.data(), nullptr, 0, nullptr) == SMLIB_ERR_BUFFER);
    c.status = smlib_ssi_curve_nurbs(r, i, which, nullptr, nullptr, nullptr, np, c.points.data(), c.weights.data(), nk, c.knots.data());
    CHECK(c.status == SMLIB_OK);
    return c;
}

// The curve's point at `t` (de Boor on the homogeneous points), `dim` coordinates into `out`.
static void point(const Nurbs &c, int dim, double t, double *out)
{
    const int p = c.degree;
    const int n = static_cast<int>(c.weights.size());
    int k = p;
    while (k + 1 < n && c.knots[k + 1] <= t)
    {
        ++k;
    }
    std::vector<double> d((p + 1) * 4);
    for (int j = 0; j <= p; ++j)
    {
        const int q = k - p + j;
        for (int e = 0; e < dim; ++e)
        {
            d[4 * j + e] = c.points[dim * q + e] * c.weights[q];
        }
        d[4 * j + 3] = c.weights[q];
    }
    for (int rr = 1; rr <= p; ++rr)
    {
        for (int j = p; j >= rr; --j)
        {
            const double lo = c.knots[k - p + j], hi = c.knots[k + 1 + j - rr];
            const double a = hi > lo ? (t - lo) / (hi - lo) : 0.0;
            for (int e = 0; e < 4; ++e)
            {
                d[4 * j + e] = (1.0 - a) * d[4 * (j - 1) + e] + a * d[4 * j + e];
            }
        }
    }
    for (int e = 0; e < dim; ++e)
    {
        out[e] = d[4 * p + e] / d[4 * p + 3];
    }
}

// Curve `i` evaluated at `n` parameters spread evenly over each B-spline's knot range: the 3D
// curve, and where the method returned them the UV curves with the surfaces at their points.
struct Samples
{
    int32_t type = -2, flags = 0, n = 0;
    std::vector<double> xyz, uv_a, uv_b, a, b;
};

static void sample(const Nurbs &c, int dim, int32_t n, std::vector<double> &out)
{
    out.assign(static_cast<size_t>(dim) * n, NaN);
    if (c.status != SMLIB_OK)
    {
        return;
    }
    const double t0 = c.knots[c.degree], t1 = c.knots[c.weights.size()];
    for (int32_t k = 0; k < n; ++k)
    {
        point(c, dim, t0 + (t1 - t0) * k / (n - 1), &out[static_cast<size_t>(dim) * k]);
    }
}

static void at(const smlib_surface *s, const std::vector<double> &uv, std::vector<double> &out)
{
    out.assign(uv.size() / 2 * 3, NaN);
    for (size_t k = 0; 2 * k < uv.size(); ++k)
    {
        if (std::isfinite(uv[2 * k]))
        {
            CHECK(smlib_surface_evaluate(s, uv[2 * k], uv[2 * k + 1], &out[3 * k]) == SMLIB_OK);
        }
    }
}

static Samples curve(const smlib_ssi_result *r, int32_t i, const smlib_surface *a = nullptr, const smlib_surface *b = nullptr, int32_t n = 201)
{
    Samples s;
    s.n = n;
    CHECK(smlib_ssi_curve_info(r, i, &s.type, &s.flags, nullptr, nullptr) == SMLIB_OK);
    const Nurbs c3 = nurbs(r, i, SMLIB_SSI_CURVE_3D, 3), ua = nurbs(r, i, SMLIB_SSI_CURVE_UV_A, 2), ub = nurbs(r, i, SMLIB_SSI_CURVE_UV_B, 2);
    CHECK(c3.status == SMLIB_OK);
    CHECK((ua.status == SMLIB_OK) == ((s.flags & SMLIB_SSI_HAS_UV_A) != 0) && (ua.status == SMLIB_OK || ua.status == SMLIB_ERR_EMPTY));
    CHECK((ub.status == SMLIB_OK) == ((s.flags & SMLIB_SSI_HAS_UV_B) != 0) && (ub.status == SMLIB_OK || ub.status == SMLIB_ERR_EMPTY));
    sample(c3, 3, n, s.xyz);
    sample(ua, 2, n, s.uv_a);
    sample(ub, 2, n, s.uv_b);
    if (a != nullptr && b != nullptr)
    {
        at(a, s.uv_a, s.a);
        at(b, s.uv_b, s.b);
    }
    return s;
}

static double dist(const double *p, const double *q)
{
    return std::sqrt((p[0] - q[0]) * (p[0] - q[0]) + (p[1] - q[1]) * (p[1] - q[1]) + (p[2] - q[2]) * (p[2] - q[2]));
}

// Every sample of every curve: `on` (the implicit equations of both surfaces, as one residual)
// within `tol`; where the method returned UV curves, the surfaces at those UVs within `tol` of
// the 3D sample too.
template <typename F>
static void check_on_both(const smlib_ssi_result *r, const smlib_surface *a, const smlib_surface *b, double tol, F on)
{
    int32_t n = 0;
    CHECK(smlib_ssi_curve_count(r, &n) == SMLIB_OK);
    double worst = 0.0;
    for (int32_t i = 0; i < n; ++i)
    {
        Samples s = curve(r, i, a, b);
        for (int32_t k = 0; k < s.n; ++k)
        {
            const double *p = &s.xyz[3 * k];
            CHECK(std::isfinite(p[0]) && std::isfinite(p[1]) && std::isfinite(p[2]));
            worst = std::max(worst, on(p));
            if (s.flags & SMLIB_SSI_HAS_UV_A)
            {
                worst = std::max(worst, dist(p, &s.a[3 * k]));
            }
            else
            {
                CHECK(std::isnan(s.uv_a[2 * k]) && std::isnan(s.a[3 * k]));
            }
            if (s.flags & SMLIB_SSI_HAS_UV_B)
            {
                worst = std::max(worst, dist(p, &s.b[3 * k]));
            }
        }
    }
    CHECK(worst < tol);
}

static int32_t count(const smlib_ssi_result *r)
{
    int32_t n = -1;
    CHECK(smlib_ssi_curve_count(r, &n) == SMLIB_OK);
    return n;
}

static int32_t type(const smlib_ssi_result *r, int32_t i)
{
    int32_t t = -2;
    CHECK(smlib_ssi_curve_info(r, i, &t, nullptr, nullptr, nullptr) == SMLIB_OK);
    return t;
}

// c1_plane_cyl_ellipse: a tilted plane through a cylinder's axis, one closed ellipse.
static void plane_cylinder()
{
    const double py[3] = {0.0, 0.8944271909999159, 0.4472135954999579};
    const double u[2] = {-3, 3}, v[2] = {-3, 3}, angle[2] = {0, 360}, height[2] = {-2, 2};
    smlib_surface *plane = nullptr, *cyl = nullptr;
    CHECK(smlib_surface_plane(O, X, py, u, v, &plane) == SMLIB_OK);
    CHECK(smlib_surface_cylinder(O, X, Y, 1.0, angle, height, &cyl) == SMLIB_OK);
    double p[3];
    CHECK(smlib_surface_evaluate(plane, 1.0, 2.0, p) == SMLIB_OK);
    CHECK(std::fabs(p[0] - 1.0) < 1e-12 && std::fabs(p[1] - 2 * py[1]) < 1e-12 && std::fabs(p[2] - 2 * py[2]) < 1e-12);
    double dom[4];
    CHECK(smlib_surface_domain(plane, dom) == SMLIB_OK && dom[0] == -3 && dom[1] == 3 && dom[2] == -3 && dom[3] == 3);
    // On the plane: n . p = 0 with n = X x py; on the cylinder: x^2 + y^2 = 1.
    const double nrm[3] = {0, -py[2], py[1]};
    auto on = [&](const double *q)
    {
        return std::max(std::fabs(q[0] * nrm[0] + q[1] * nrm[1] + q[2] * nrm[2]), std::fabs(std::hypot(q[0], q[1]) - 1.0));
    };
    for (int32_t method : {SMLIB_SSI_GLOBAL, SMLIB_SSI_ADVANCED})
    {
        smlib_ssi_result *r = nullptr;
        CHECK(smlib_ssi(plane, cyl, method, 1e-6, &r) == SMLIB_OK);
        if (r == nullptr)
        {
            continue;
        }
        CHECK(count(r) == 1);
        CHECK(type(r, 0) == SMLIB_SSI_CROSSING);
        int32_t np = -1;
        CHECK(smlib_ssi_point_count(r, &np) == SMLIB_OK && np == 0);
        check_on_both(r, plane, cyl, 1e-5, on);
        if (method == SMLIB_SSI_ADVANCED)
        {
            Samples s = curve(r, 0);
            CHECK(s.flags == (SMLIB_SSI_HAS_UV_A | SMLIB_SSI_HAS_UV_B));
        }
        else
        {
            // The analytic special case returns no UV curves.
            Samples s = curve(r, 0);
            CHECK(s.flags == 0);
        }
        smlib_ssi_result_free(r);
    }
    smlib_surface_free(plane);
    smlib_surface_free(cyl);
}

// c1_plane_sphere_circle: z = 1 through a sphere of radius 2, a circle of radius sqrt(3).
static void plane_sphere()
{
    const double o[3] = {0, 0, 1}, u[2] = {-3, 3}, lon[2] = {0, 360}, lat[2] = {-90, 90};
    smlib_surface *plane = nullptr, *sphere = nullptr;
    CHECK(smlib_surface_plane(o, X, Y, u, u, &plane) == SMLIB_OK);
    CHECK(smlib_surface_sphere(O, X, Y, 2.0, lon, lat, &sphere) == SMLIB_OK);
    smlib_ssi_result *r = nullptr;
    CHECK(smlib_ssi(plane, sphere, SMLIB_SSI_GLOBAL, 0.0, &r) == SMLIB_OK);
    if (r != nullptr)
    {
        CHECK(count(r) == 1);
        CHECK(type(r, 0) == SMLIB_SSI_CROSSING);
        check_on_both(r, plane, sphere, 1e-6, [](const double *q) { return std::max(std::fabs(q[2] - 1.0), std::fabs(dist(q, O) - 2.0)); });
        smlib_ssi_result_free(r);
    }
    smlib_surface_free(plane);
    smlib_surface_free(sphere);
}

// c1n_parab_parab_circle: z = x^2 + y^2 against z = 2 - x^2 - y^2 over [-2, 2]^2, both as
// biquadratic Beziers; they cross in the unit circle at z = 1.
static void nurbs_paraboloids()
{
    const double k[6] = {0, 0, 0, 1, 1, 1};
    const double up[27] = {-2, -2, 8, -2, 0, 0, -2, 2, 8, 0, -2, 0, 0, 0, -8, 0, 2, 0, 2, -2, 8, 2, 0, 0, 2, 2, 8};
    const double down[27] = {-2, -2, -6, -2, 0, 2, -2, 2, -6, 0, -2, 2, 0, 0, 10, 0, 2, 2, 2, -2, -6, 2, 0, 2, 2, 2, -6};
    const double w[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
    smlib_surface *a = nullptr, *b = nullptr;
    CHECK(smlib_surface_nurbs(2, 2, 3, 3, k, 6, k, 6, up, nullptr, &a) == SMLIB_OK);
    CHECK(smlib_surface_nurbs(2, 2, 3, 3, k, 6, k, 6, down, w, &b) == SMLIB_OK);
    double p[3];
    CHECK(smlib_surface_evaluate(a, 0.5, 0.5, p) == SMLIB_OK && dist(p, O) < 1e-12);
    for (int32_t method : {SMLIB_SSI_GLOBAL, SMLIB_SSI_ADVANCED, SMLIB_SSI_BASIC})
    {
        smlib_ssi_result *r = nullptr;
        CHECK(smlib_ssi(a, b, method, 1e-6, &r) == SMLIB_OK);
        if (r == nullptr)
        {
            continue;
        }
        CHECK(count(r) == 1);
        CHECK(type(r, 0) == SMLIB_SSI_CROSSING);
        check_on_both(r, a, b, 1e-5, [](const double *q) { return std::max(std::fabs(q[2] - 1.0), std::fabs(std::hypot(q[0], q[1]) - 1.0)); });
        Samples s = curve(r, 0);
        CHECK(s.flags == (SMLIB_SSI_HAS_UV_A | SMLIB_SSI_HAS_UV_B));
        smlib_ssi_result_free(r);
    }
    smlib_surface_free(a);
    smlib_surface_free(b);
}

// c3_sphere_on_plane: a unit sphere resting on z = 0 touches it at the origin.
static void sphere_touching_plane()
{
    const double c[3] = {0, 0, 1}, u[2] = {-3, 3}, lon[2] = {0, 360}, lat[2] = {-90, 90};
    smlib_surface *sphere = nullptr, *plane = nullptr;
    CHECK(smlib_surface_sphere(c, X, Y, 1.0, lon, lat, &sphere) == SMLIB_OK);
    CHECK(smlib_surface_plane(O, X, Y, u, u, &plane) == SMLIB_OK);
    smlib_ssi_result *r = nullptr;
    CHECK(smlib_ssi(sphere, plane, SMLIB_SSI_GLOBAL, 1e-6, &r) == SMLIB_OK);
    if (r != nullptr)
    {
        CHECK(count(r) == 1);
        CHECK(type(r, 0) == SMLIB_SSI_TOUCHING);
        int32_t np = 0;
        CHECK(smlib_ssi_point_count(r, &np) == SMLIB_OK && np == 1);
        double p[3], ua[2], ub[2];
        CHECK(smlib_ssi_point(r, 0, p, ua, ub) == SMLIB_OK && dist(p, O) < 1e-9);
        CHECK(smlib_ssi_point(r, 1, p, ua, ub) == SMLIB_ERR_ARG);
        smlib_ssi_result_free(r);
    }
    smlib_surface_free(sphere);
    smlib_surface_free(plane);
}

// A plane z = 0 against z = f(x) over x in [0, 4], y in [0, 1], with f = 0 for x <= 1 and
// f = s (x - 1)^2 (x - 3)^2 beyond: the surfaces coincide over x in [0, 1] and touch along the
// line x = 3 without crossing. f is an exact quartic B-spline in x (knot 1 of multiplicity 3,
// so C1 there; control x at the Greville abscissae, so x(u) = u) and linear in y. `mirror` maps
// x to 4 - x, so the coincident strip and the tangent line swap ends. The tangent line must
// come out whichever contact the tangent-curve search meets first.
static void partly_coincident_with_tangent_line(bool mirror)
{
    const double s = 0.25;
    const double ku[13] = {0, 0, 0, 0, 0, 1, 1, 1, 4, 4, 4, 4, 4}, kv[4] = {0, 0, 1, 1};
    const double gx[8] = {0, 0.25, 0.5, 0.75, 1.75, 2.5, 3.25, 4};
    const double gz[8] = {0, 0, 0, 0, 0, 6 * s, -9 * s, 9 * s};
    double pts[48];
    for (int i = 0; i < 8; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            double *p = &pts[3 * (i * 2 + j)];
            p[0] = mirror ? 4.0 - gx[i] : gx[i];
            p[1] = j;
            p[2] = gz[i];
        }
    }
    const double xt = mirror ? 1.0 : 3.0;
    smlib_surface *plane = nullptr, *bump = nullptr;
    const double pu[2] = {-1, 5}, pv[2] = {-1, 2};
    CHECK(smlib_surface_plane(O, X, Y, pu, pv, &plane) == SMLIB_OK);
    CHECK(smlib_surface_nurbs(4, 1, 8, 2, ku, 13, kv, 4, pts, nullptr, &bump) == SMLIB_OK);
    double q[3];
    CHECK(smlib_surface_evaluate(bump, 2.0, 0.5, q) == SMLIB_OK && std::fabs(q[0] - 2.0) < 1e-12 && std::fabs(q[2] - s) < 1e-12);
    CHECK(smlib_surface_evaluate(bump, 3.0, 0.5, q) == SMLIB_OK && std::fabs(q[0] - xt) < 1e-12 && std::fabs(q[2]) < 1e-12);
    CHECK(smlib_surface_evaluate(bump, 0.5, 0.5, q) == SMLIB_OK && std::fabs(q[2]) < 1e-15);
    // Not SMLIB_SSI_BASIC: the basic intersector has no seeds inside the domains for the tangent
    // search, and finds no tangent line in the middle of a pair even without the coincident strip.
    for (int32_t method : {SMLIB_SSI_GLOBAL, SMLIB_SSI_ADVANCED})
    {
        smlib_ssi_result *r = nullptr;
        CHECK(smlib_ssi(plane, bump, method, 1e-6, &r) == SMLIB_OK);
        if (r == nullptr)
        {
            continue;
        }
        // The tangent line: a curve all along x = xt, z = 0, from y = 0 to y = 1.
        bool found = false;
        for (int32_t i = 0; i < count(r) && !found; ++i)
        {
            Samples c = curve(r, i);
            double ylo = 1e300, yhi = -1e300, off = 0.0;
            for (int32_t k = 0; k < c.n; ++k)
            {
                const double *p = &c.xyz[3 * k];
                off = std::max(off, std::max(std::fabs(p[0] - xt), std::fabs(p[2])));
                ylo = std::min(ylo, p[1]);
                yhi = std::max(yhi, p[1]);
            }
            found = off < 1e-4 && ylo < 1e-4 && yhi > 1.0 - 1e-4;
        }
        if (!found)
        {
            std::fprintf(stderr, "partly_coincident_with_tangent_line(mirror=%d), method %d: no curve along x = %g\n", mirror ? 1 : 0,
                         method, xt);
        }
        CHECK(found);
        smlib_ssi_result_free(r);
    }
    smlib_surface_free(plane);
    smlib_surface_free(bump);
}

static void hostile()
{
    smlib_surface *s = nullptr;
    const double u[2] = {-1, 1}, bad_u[2] = {1, -1}, nan3[3] = {NaN, 0, 0}, inf3[3] = {0, INF, 0}, zero[3] = {0, 0, 0};
    const double huge3[3] = {1e300, 0, 0}, ang[2] = {0, 360}, wide[2] = {0, 400}, h[2] = {0, 1}, rev_h[2] = {1, 0};
    const double lat_bad[2] = {-100, 90}, neg_h[2] = {-1, 1};
    CHECK(smlib_surface_plane(nan3, X, Y, u, u, &s) == SMLIB_ERR_ARG);
    CHECK(smlib_surface_plane(O, inf3, Y, u, u, &s) == SMLIB_ERR_ARG);
    CHECK(smlib_surface_plane(O, zero, Y, u, u, &s) == SMLIB_ERR_ARG);
    CHECK(smlib_surface_plane(O, X, X, u, u, &s) == SMLIB_ERR_ARG);
    CHECK(smlib_surface_plane(huge3, X, Y, u, u, &s) == SMLIB_ERR_ARG);
    CHECK(smlib_surface_plane(O, X, Y, bad_u, u, &s) == SMLIB_ERR_ARG);
    CHECK(smlib_surface_plane(nullptr, X, Y, u, u, &s) == SMLIB_ERR_ARG);
    CHECK(smlib_surface_plane(O, X, Y, u, u, nullptr) == SMLIB_ERR_ARG);
    CHECK(smlib_surface_cylinder(O, X, Y, 0.0, ang, h, &s) == SMLIB_ERR_ARG);
    CHECK(smlib_surface_cylinder(O, X, Y, -1.0, ang, h, &s) == SMLIB_ERR_ARG);
    CHECK(smlib_surface_cylinder(O, X, Y, NaN, ang, h, &s) == SMLIB_ERR_ARG);
    CHECK(smlib_surface_cylinder(O, X, Y, 1.0, wide, h, &s) == SMLIB_ERR_ARG);
    CHECK(smlib_surface_cylinder(O, X, Y, 1.0, ang, rev_h, &s) == SMLIB_ERR_ARG);
    CHECK(smlib_surface_cone(O, X, Y, 0.0, h, &s) == SMLIB_ERR_ARG);
    CHECK(smlib_surface_cone(O, X, Y, 90.0, h, &s) == SMLIB_ERR_ARG);
    CHECK(smlib_surface_cone(O, X, Y, 30.0, neg_h, &s) == SMLIB_ERR_ARG);
    CHECK(smlib_surface_sphere(O, X, Y, 1.0, ang, lat_bad, &s) == SMLIB_ERR_ARG);
    CHECK(smlib_surface_sphere(O, X, zero, 1.0, ang, h, &s) == SMLIB_ERR_ARG);
    CHECK(smlib_surface_torus(O, X, Y, 1.0, 2.0, ang, ang, &s) == SMLIB_ERR_ARG);
    CHECK(smlib_surface_torus(O, X, Y, 3.0, 0.0, ang, ang, &s) == SMLIB_ERR_ARG);
    CHECK(s == nullptr);

    const double k[4] = {0, 0, 1, 1}, k_desc[4] = {0, 1, 0, 1}, k_flat[4] = {0, 0, 0, 0}, k_nan[4] = {0, 0, NaN, 1};
    const double pts[12] = {0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0};
    const double pts_nan[12] = {0, 0, 0, 0, NaN, 0, 1, 0, 0, 1, 1, 0};
    const double w[4] = {1, 1, 1, 1}, w_zero[4] = {1, 0, 1, 1}, w_neg[4] = {1, -1, 1, 1};
    CHECK(smlib_surface_nurbs(1, 1, 2, 2, k, 4, k, 4, pts, w, &s) == SMLIB_OK);
    smlib_surface_free(s);
    s = nullptr;
    CHECK(smlib_surface_nurbs(1, 1, 2, 2, k_desc, 4, k, 4, pts, w, &s) == SMLIB_ERR_ARG);
    CHECK(smlib_surface_nurbs(1, 1, 2, 2, k_flat, 4, k, 4, pts, w, &s) == SMLIB_ERR_ARG);
    CHECK(smlib_surface_nurbs(1, 1, 2, 2, k_nan, 4, k, 4, pts, w, &s) == SMLIB_ERR_ARG);
    CHECK(smlib_surface_nurbs(1, 1, 2, 2, k, 3, k, 4, pts, w, &s) == SMLIB_ERR_ARG);
    CHECK(smlib_surface_nurbs(1, 1, 2, 2, k, 4, k, 4, pts_nan, w, &s) == SMLIB_ERR_ARG);
    CHECK(smlib_surface_nurbs(1, 1, 2, 2, k, 4, k, 4, pts, w_zero, &s) == SMLIB_ERR_ARG);
    CHECK(smlib_surface_nurbs(1, 1, 2, 2, k, 4, k, 4, pts, w_neg, &s) == SMLIB_ERR_ARG);
    CHECK(smlib_surface_nurbs(0, 1, 2, 2, k, 4, k, 4, pts, w, &s) == SMLIB_ERR_ARG);
    CHECK(smlib_surface_nurbs(2, 1, 2, 2, k, 4, k, 4, pts, w, &s) == SMLIB_ERR_ARG);
    CHECK(smlib_surface_nurbs(1, 1, 2, 2, nullptr, 4, k, 4, pts, w, &s) == SMLIB_ERR_ARG);
    CHECK(smlib_surface_nurbs(1, 1, 2, 2, k, 4, k, 4, nullptr, w, &s) == SMLIB_ERR_ARG);
    // Huge counts are rejected before any array is read.
    CHECK(smlib_surface_nurbs(1, 1, 2000000, 2, k, 2000003, k, 4, pts, w, &s) == SMLIB_ERR_ARG);
    CHECK(smlib_surface_nurbs(1, 1, 100000, 100000, k, 100003, k, 100003, pts, w, &s) == SMLIB_ERR_ARG);
    CHECK(smlib_surface_nurbs(40, 1, 41, 2, k, 82, k, 4, pts, w, &s) == SMLIB_ERR_ARG);
    CHECK(s == nullptr);

    smlib_surface *a = nullptr;
    CHECK(smlib_surface_plane(O, X, Y, u, u, &a) == SMLIB_OK);
    smlib_ssi_result *r = nullptr;
    CHECK(smlib_ssi(a, nullptr, 0, 0.0, &r) == SMLIB_ERR_ARG);
    CHECK(smlib_ssi(a, a, 3, 0.0, &r) == SMLIB_ERR_ARG);
    CHECK(smlib_ssi(a, a, -1, 0.0, &r) == SMLIB_ERR_ARG);
    CHECK(smlib_ssi(a, a, 0, NaN, &r) == SMLIB_ERR_ARG);
    CHECK(smlib_ssi(a, a, 0, 0.0, nullptr) == SMLIB_ERR_ARG);
    CHECK(r == nullptr);
    CHECK(smlib_ssi_curve_count(nullptr, nullptr) == SMLIB_ERR_ARG);
    CHECK(smlib_surface_evaluate(a, NaN, 0.0, nullptr) == SMLIB_ERR_ARG);
    smlib_surface_free(a);
    smlib_surface_free(nullptr);
    smlib_ssi_result_free(nullptr);
}

int main()
{
    plane_cylinder();
    plane_sphere();
    nurbs_paraboloids();
    sphere_touching_plane();
    partly_coincident_with_tangent_line(false);
    partly_coincident_with_tangent_line(true);
    hostile();
    if (failures != 0)
    {
        std::fprintf(stderr, "%d check(s) failed\n", failures);
        return 1;
    }
    std::printf("ssi: all checks passed\n");
    return 0;
}
