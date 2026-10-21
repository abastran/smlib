// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// Surface/surface intersection: plain-C surface descriptions (built exactly as the SSI corpus
// runner, ssi-corpus/runner/ssi_runner.cpp, builds them from a case file), SMLib's three SSI
// entry points, and the curves sampled the way that runner samples them.

#include "smlib_shim.h"

#include <SmAdvSurfaceIntersector.h>
#include <SmAxis2Placement.h>
#include <SmBSplineSurface.h>
#include <SmCone.h>
#include <SmCurve.h>
#include <SmCylinder.h>
#include <NL_Geometry.h>
#include <NL_Nurbsdef.h>
#include <SmNurbs.h>
#include <SmNurbsSrf.h>
#include <SmPlane.h>
#include <SmSphere.h>
#include <SmSurface.h>
#include <SmSurfaceIntersector.h>
#include <SmTorus.h>

#include <limits>
#include <memory>
#include <set>

struct smlib_surface
{
    SmSurface *s = nullptr;
    SmExtent2d domain;
    // The description in the SSI corpus case format, for dumps and the log.
    std::string json;

    ~smlib_surface()
    {
        delete s;
    }
};

namespace
{

struct Curve
{
    int32_t type;
    int32_t flags;
    double deviation;
    double interval[2];
    // Per sample: 3, 2, 2, 3, 3 values.
    std::vector<double> xyz, uv_a, uv_b, a_at_uv, b_at_uv;
};

struct Point
{
    double xyz[3], uv_a[2], uv_b[2];
};

} // namespace

struct smlib_ssi_result
{
    int32_t samples = 0;
    std::vector<Curve> curves;
    std::vector<Point> points;
};

using namespace smlib_shim;

namespace
{

constexpr double NaN = std::numeric_limits<double>::quiet_NaN();

bool bounded(const double *v, int n, double limit = SMLIB_SURFACE_MAX_COORD)
{
    if (v == nullptr || !finite(v, n))
    {
        return false;
    }
    for (int i = 0; i < n; ++i)
    {
        if (std::fabs(v[i]) > limit)
        {
            return false;
        }
    }
    return true;
}

// Two finite, non-zero axes that are not (nearly) parallel.
bool axes(const double *x, const double *y)
{
    if (!bounded(x, 3) || !bounded(y, 3))
    {
        return false;
    }
    const double lx = std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
    const double ly = std::sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]);
    if (!(lx > 1e-12) || !(ly > 1e-12))
    {
        return false;
    }
    const double c[3] = {x[1] * y[2] - x[2] * y[1], x[2] * y[0] - x[0] * y[2], x[0] * y[1] - x[1] * y[0]};
    return std::sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]) > 1e-9 * lx * ly;
}

// An increasing interval [r[0], r[1]] inside [lo, hi] spanning at most `span`.
bool range(const double *r, double lo, double hi, double span)
{
    return bounded(r, 2) && r[0] < r[1] && r[0] >= lo && r[1] <= hi && r[1] - r[0] <= span;
}

bool length(double x)
{
    return std::isfinite(x) && x > 1e-9 && x <= SMLIB_SURFACE_MAX_COORD;
}

SmVector3d v3(const double *v)
{
    return SmVector3d(v[0], v[1], v[2]);
}

// The runner's unit_cross.
SmVector3d unit_cross(const SmVector3d &a, const SmVector3d &b)
{
    SmVector3d c(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
    double l = std::sqrt(c.x * c.x + c.y * c.y + c.z * c.z);
    return SmVector3d(c.x / l, c.y / l, c.z / l);
}

std::string field(const char *name, const std::string &value)
{
    return std::string(",\"") + name + "\":" + value;
}

// Runs one constructor: `make` builds the SMLib surface under the lock.
template <typename F>
int32_t construct(const char *fn, std::string json, smlib_surface **out, F &&make)
{
    Call c(fn, json);
    c.input_text("surface",
                 [&](const std::string &path)
                 {
                     if (FILE *f = std::fopen(path.c_str(), "w"))
                     {
                         std::fprintf(f, "%s\n", json.c_str());
                         std::fclose(f);
                     }
                 });
    return c.run(
        [&]
        {
            const SmContext &ctx = *SmApiGetOrCreateContext();
            std::unique_ptr<smlib_surface> s(new smlib_surface);
            s->json = json;
            int32_t st = make(ctx, s->s);
            if (st != SMLIB_OK)
            {
                return st;
            }
            if (s->s == nullptr)
            {
                return static_cast<int32_t>(SM_ERR);
            }
            s->domain = s->s->GetNaturalUVDomain();
            const double d[4] = {s->domain.GetMin().x, s->domain.GetMax().x, s->domain.GetMin().y, s->domain.GetMax().y};
            if (!finite(d, 4) || !(d[0] < d[1]) || !(d[2] < d[3]))
            {
                return static_cast<int32_t>(SM_ERR_INVALID_INPUT);
            }
            *out = s.release();
            return SMLIB_OK;
        }
    );
}

// The runner's sample_curve, except that a failed evaluation is NaN instead of being dropped.
void sample(const SmCurve *c, int32_t n, int dim, std::vector<double> &out, double *interval)
{
    out.assign(static_cast<size_t>(n) * dim, NaN);
    if (c == nullptr)
    {
        return;
    }
    SmExtent1d ivl = c->GetNaturalInterval();
    if (interval != nullptr)
    {
        interval[0] = ivl.GetMin();
        interval[1] = ivl.GetMax();
    }
    for (int32_t i = 0; i < n; ++i)
    {
        double t = ivl.GetMin() + (ivl.GetMax() - ivl.GetMin()) * (static_cast<double>(i) / (n - 1));
        SmPoint3d p;
        if (c->EvaluatePoint(t, p) == SM_SUCCESS)
        {
            const double q[3] = {p.x, p.y, p.z};
            std::memcpy(&out[static_cast<size_t>(i) * dim], q, dim * sizeof(double));
        }
    }
}

void evaluate(const SmSurface *s, const std::vector<double> &uv, std::vector<double> &out)
{
    const size_t n = uv.size() / 2;
    out.assign(n * 3, NaN);
    for (size_t i = 0; i < n; ++i)
    {
        if (!std::isfinite(uv[2 * i]) || !std::isfinite(uv[2 * i + 1]))
        {
            continue;
        }
        SmPoint3d q;
        if (s->EvaluatePoint(SmPoint2d(uv[2 * i], uv[2 * i + 1]), q) == SM_SUCCESS)
        {
            out[3 * i] = q.x;
            out[3 * i + 1] = q.y;
            out[3 * i + 2] = q.z;
        }
    }
}

// Owns what an intersector hands back.
struct Curves
{
    SmTArray<SmCurve *> c3, uv1, uv2;

    ~Curves()
    {
        std::set<SmCurve *> all;
        for (SmTArray<SmCurve *> *a : {&c3, &uv1, &uv2})
        {
            for (ULONG i = 0; i < a->GetSize(); ++i)
            {
                all.insert((*a)[i]);
            }
        }
        for (SmCurve *c : all)
        {
            delete c;
        }
    }
};

constexpr int64_t MAX_TOTAL_SAMPLES = 10000000;

} // namespace

extern "C" int32_t smlib_surface_plane(const double origin[3], const double xaxis[3], const double yaxis[3], const double u[2],
                                       const double v[2], smlib_surface **out)
{
    if (out == nullptr || !bounded(origin, 3) || !axes(xaxis, yaxis) || !range(u, -SMLIB_SURFACE_MAX_COORD, SMLIB_SURFACE_MAX_COORD, 2 * SMLIB_SURFACE_MAX_COORD) ||
        !range(v, -SMLIB_SURFACE_MAX_COORD, SMLIB_SURFACE_MAX_COORD, 2 * SMLIB_SURFACE_MAX_COORD))
    {
        return SMLIB_ERR_ARG;
    }
    std::string j = "{\"type\":\"plane\"" + field("origin", vec(origin, 3)) + field("xaxis", vec(xaxis, 3)) +
                    field("yaxis", vec(yaxis, 3)) + field("u", vec(u, 2)) + field("v", vec(v, 2)) + "}";
    return construct("smlib_surface_plane", j, out,
                     [&](const SmContext &ctx, SmSurface *&s)
                     {
                         SmExtent2d dom(u[0], v[0], u[1], v[1]);
                         s = new (ctx) SmPlane(v3(origin), v3(xaxis), v3(yaxis), SmVector2d(1.0, 1.0), dom, &ctx);
                         return SMLIB_OK;
                     });
}

extern "C" int32_t smlib_surface_cylinder(const double origin[3], const double xaxis[3], const double yaxis[3], double radius,
                                          const double angle_deg[2], const double height[2], smlib_surface **out)
{
    if (out == nullptr || !bounded(origin, 3) || !axes(xaxis, yaxis) || !length(radius) || !range(angle_deg, -360.0, 360.0, 360.0) ||
        !range(height, -SMLIB_SURFACE_MAX_COORD, SMLIB_SURFACE_MAX_COORD, SMLIB_SURFACE_MAX_COORD))
    {
        return SMLIB_ERR_ARG;
    }
    std::string j = "{\"type\":\"cylinder\"" + field("origin", vec(origin, 3)) + field("xaxis", vec(xaxis, 3)) +
                    field("yaxis", vec(yaxis, 3)) + field("radius", num(radius)) + field("angle", vec(angle_deg, 2)) +
                    field("height", vec(height, 2)) + "}";
    return construct("smlib_surface_cylinder", j, out,
                     [&](const SmContext &ctx, SmSurface *&s)
                     {
                         SmVector3d x = v3(xaxis), y = v3(yaxis);
                         SmVector3d z = unit_cross(x, y);
                         SmPoint3d o = v3(origin) + z * height[0];
                         s = new (ctx) SmCylinder(o, x, y, radius, angle_deg[0], angle_deg[1], height[1] - height[0], FALSE, FALSE,
                                                  FALSE, &ctx);
                         return SMLIB_OK;
                     });
}

extern "C" int32_t smlib_surface_cone(const double apex[3], const double xaxis[3], const double yaxis[3], double half_angle_deg,
                                      const double height[2], smlib_surface **out)
{
    if (out == nullptr || !bounded(apex, 3) || !axes(xaxis, yaxis) || !std::isfinite(half_angle_deg) || !(half_angle_deg >= 0.01) ||
        !(half_angle_deg <= 89.0) || !range(height, 0.0, SMLIB_SURFACE_MAX_COORD, SMLIB_SURFACE_MAX_COORD) ||
        !(height[1] - height[0] > 1e-9))
    {
        return SMLIB_ERR_ARG;
    }
    std::string j = "{\"type\":\"cone\"" + field("apex", vec(apex, 3)) + field("xaxis", vec(xaxis, 3)) + field("yaxis", vec(yaxis, 3)) +
                    field("half_angle_deg", num(half_angle_deg)) + field("height", vec(height, 2)) + "}";
    return construct("smlib_surface_cone", j, out,
                     [&](const SmContext &ctx, SmSurface *&s)
                     {
                         SmVector3d x = v3(xaxis), y = v3(yaxis);
                         SmVector3d z = unit_cross(x, y);
                         const double h0 = height[0], h1 = height[1];
                         const double t = std::tan(half_angle_deg * M_PI / 180.0);
                         SmPoint3d o = v3(apex) + z * h0;
                         SmAxis2Placement pl(o, x, y);
                         SmCone *cone = nullptr;
                         int32_t st = status(SmCone::CreateCanonical(ctx, pl, h0 * t, h1 * t, h1 - h0, cone));
                         if (st != SMLIB_OK)
                         {
                             delete cone;
                             return st;
                         }
                         s = cone;
                         return SMLIB_OK;
                     });
}

extern "C" int32_t smlib_surface_sphere(const double center[3], const double xaxis[3], const double yaxis[3], double radius,
                                        const double lon_deg[2], const double lat_deg[2], smlib_surface **out)
{
    if (out == nullptr || !bounded(center, 3) || !axes(xaxis, yaxis) || !length(radius) || !range(lon_deg, -360.0, 360.0, 360.0) ||
        !range(lat_deg, -90.0, 90.0, 180.0))
    {
        return SMLIB_ERR_ARG;
    }
    std::string j = "{\"type\":\"sphere\"" + field("center", vec(center, 3)) + field("xaxis", vec(xaxis, 3)) +
                    field("yaxis", vec(yaxis, 3)) + field("radius", num(radius)) + field("lon", vec(lon_deg, 2)) +
                    field("lat", vec(lat_deg, 2)) + "}";
    return construct("smlib_surface_sphere", j, out,
                     [&](const SmContext &ctx, SmSurface *&s)
                     {
                         SmExtent2d anal(lon_deg[0], lat_deg[0], lon_deg[1], lat_deg[1]);
                         s = new (ctx) SmSphere(v3(center), v3(xaxis), v3(yaxis), anal, radius, FALSE, FALSE, NULL, &ctx);
                         return SMLIB_OK;
                     });
}

extern "C" int32_t smlib_surface_torus(const double center[3], const double xaxis[3], const double yaxis[3], double major_radius,
                                       double minor_radius, const double major_deg[2], const double minor_deg[2], smlib_surface **out)
{
    if (out == nullptr || !bounded(center, 3) || !axes(xaxis, yaxis) || !length(major_radius) || !length(minor_radius) ||
        !(minor_radius < major_radius) || !range(major_deg, -360.0, 360.0, 360.0) || !range(minor_deg, -360.0, 360.0, 360.0))
    {
        return SMLIB_ERR_ARG;
    }
    std::string j = "{\"type\":\"torus\"" + field("center", vec(center, 3)) + field("xaxis", vec(xaxis, 3)) +
                    field("yaxis", vec(yaxis, 3)) + field("R", num(major_radius)) + field("r", num(minor_radius)) +
                    field("major", vec(major_deg, 2)) + field("minor", vec(minor_deg, 2)) + "}";
    return construct("smlib_surface_torus", j, out,
                     [&](const SmContext &ctx, SmSurface *&s)
                     {
                         SmExtent2d anal(major_deg[0], minor_deg[0], major_deg[1], minor_deg[1]);
                         s = new (ctx) SmTorus(v3(center), v3(xaxis), v3(yaxis), anal, major_radius, minor_radius, FALSE, FALSE, NULL,
                                               &ctx);
                         return SMLIB_OK;
                     });
}

namespace
{

// A knot vector of n + p + 1 finite, non-decreasing values, with a non-empty domain
// [k[p], k[n]], no interior knot repeated more than p times and no end knot more than p + 1.
bool knots(const double *k, int64_t count, int32_t n, int32_t p)
{
    if (k == nullptr || count != static_cast<int64_t>(n) + p + 1 || !bounded(k, static_cast<int>(count)))
    {
        return false;
    }
    int run = 1;
    for (int64_t i = 1; i < count; ++i)
    {
        if (k[i] < k[i - 1])
        {
            return false;
        }
        run = k[i] == k[i - 1] ? run + 1 : 1;
        if (run > p + 1)
        {
            return false;
        }
    }
    // Interior multiplicity: a knot strictly inside (k[p], k[n]) at most p times.
    run = 1;
    for (int64_t i = 1; i < count; ++i)
    {
        run = k[i] == k[i - 1] ? run + 1 : 1;
        if (run > p && k[i] > k[p] && k[i] < k[n])
        {
            return false;
        }
    }
    return k[p] < k[n];
}

} // namespace

extern "C" int32_t smlib_surface_nurbs(int32_t degree_u, int32_t degree_v, int32_t nu, int32_t nv, const double *knots_u,
                                       int64_t n_knots_u, const double *knots_v, int64_t n_knots_v, const double *points,
                                       const double *weights, smlib_surface **out)
{
    if (out == nullptr || degree_u < 1 || degree_v < 1 || degree_u > SMLIB_SURFACE_MAX_DEGREE || degree_v > SMLIB_SURFACE_MAX_DEGREE ||
        nu < degree_u + 1 || nv < degree_v + 1 || nu > SMLIB_SURFACE_MAX_POLES || nv > SMLIB_SURFACE_MAX_POLES ||
        static_cast<int64_t>(nu) * nv > SMLIB_SURFACE_MAX_POLES || !knots(knots_u, n_knots_u, nu, degree_u) ||
        !knots(knots_v, n_knots_v, nv, degree_v))
    {
        return SMLIB_ERR_ARG;
    }
    const int n = nu * nv;
    if (!bounded(points, 3 * n))
    {
        return SMLIB_ERR_ARG;
    }
    if (weights != nullptr)
    {
        for (int i = 0; i < n; ++i)
        {
            if (!std::isfinite(weights[i]) || !(weights[i] >= 1e-6) || !(weights[i] <= 1e6))
            {
                return SMLIB_ERR_ARG;
            }
        }
    }
    // The description in the case format: points as [x,y,z] triples, weights always present.
    std::string pts = "[", wts = "[";
    for (int i = 0; i < n; ++i)
    {
        pts += (i ? "," : "") + vec(points + 3 * i, 3);
        wts += (i ? "," : "") + num(weights ? weights[i] : 1.0);
    }
    pts += "]";
    wts += "]";
    const double deg[2] = {static_cast<double>(degree_u), static_cast<double>(degree_v)};
    std::string j = "{\"type\":\"nurbs\"" + field("degree", vec(deg, 2)) + field("knots_u", vec(knots_u, static_cast<int>(n_knots_u))) +
                    field("knots_v", vec(knots_v, static_cast<int>(n_knots_v))) + field("nu", std::to_string(nu)) +
                    field("nv", std::to_string(nv)) + field("points", pts) + field("weights", wts) + "}";
    // The log line gets the sizes only; the description goes to the dump.
    std::string line = "degree=[" + std::to_string(degree_u) + "," + std::to_string(degree_v) + "], nu=" + std::to_string(nu) +
                       ", nv=" + std::to_string(nv);
    Call c("smlib_surface_nurbs", line);
    c.input_text("surface",
                 [&](const std::string &path)
                 {
                     if (FILE *f = std::fopen(path.c_str(), "w"))
                     {
                         std::fprintf(f, "%s\n", j.c_str());
                         std::fclose(f);
                     }
                 });
    return c.run(
        [&]
        {
            const SmContext &ctx = *SmApiGetOrCreateContext();
            // What the runner does: allocate NLib's surface and fill it directly (homogeneous
            // control points), sidestepping SmBSplineSurface's component constructor.
            gw_SURFACE *gs = sm_AllocateNurbSurface(nu - 1, nv - 1, static_cast<gw_DEGREE>(degree_u), static_cast<gw_DEGREE>(degree_v),
                                                    static_cast<gw_INDEX>(n_knots_u - 1), static_cast<gw_INDEX>(n_knots_v - 1));
            if (gs == nullptr)
            {
                return static_cast<int32_t>(SM_ERR);
            }
            for (int64_t i = 0; i < n_knots_u; ++i)
            {
                gs->knu->U[i] = knots_u[i];
            }
            for (int64_t i = 0; i < n_knots_v; ++i)
            {
                gs->knv->U[i] = knots_v[i];
            }
            for (int32_t i = 0; i < nu; ++i)
            {
                for (int32_t k = 0; k < nv; ++k)
                {
                    const double *P = points + 3 * (static_cast<size_t>(i) * nv + k);
                    const double w = weights ? weights[static_cast<size_t>(i) * nv + k] : 1.0;
                    gs->net->Pw[i][k].x = P[0] * w;
                    gs->net->Pw[i][k].y = P[1] * w;
                    gs->net->Pw[i][k].z = P[2] * w;
                    gs->net->Pw[i][k].w = w;
                }
            }
            std::unique_ptr<smlib_surface> s(new smlib_surface);
            s->json = j;
            s->s = new (ctx) SmBSplineSurface(gs, FALSE, &ctx);
            s->domain = s->s->GetNaturalUVDomain();
            const double d[4] = {s->domain.GetMin().x, s->domain.GetMax().x, s->domain.GetMin().y, s->domain.GetMax().y};
            if (!finite(d, 4) || !(d[0] < d[1]) || !(d[2] < d[3]))
            {
                return static_cast<int32_t>(SM_ERR_INVALID_INPUT);
            }
            *out = s.release();
            return SMLIB_OK;
        }
    );
}

extern "C" void smlib_surface_free(smlib_surface *s)
{
    if (s == nullptr)
    {
        return;
    }
    try
    {
        std::lock_guard<std::mutex> g(lock());
        delete s;
    }
    catch (...)
    {
    }
}

extern "C" int32_t smlib_surface_domain(const smlib_surface *s, double out[4])
{
    if (s == nullptr || out == nullptr)
    {
        return SMLIB_ERR_ARG;
    }
    out[0] = s->domain.GetMin().x;
    out[1] = s->domain.GetMax().x;
    out[2] = s->domain.GetMin().y;
    out[3] = s->domain.GetMax().y;
    return SMLIB_OK;
}

extern "C" int32_t smlib_surface_evaluate(const smlib_surface *s, double u, double v, double out[3])
{
    if (s == nullptr || out == nullptr || !std::isfinite(u) || !std::isfinite(v))
    {
        return SMLIB_ERR_ARG;
    }
    Call c("smlib_surface_evaluate", "u=" + num(u) + ", v=" + num(v));
    return c.run(
        [&]
        {
            SmPoint3d q;
            int32_t st = status(s->s->EvaluatePoint(SmPoint2d(u, v), q));
            if (st != SMLIB_OK)
            {
                return st;
            }
            const double r[3] = {q.x, q.y, q.z};
            if (!finite(r, 3))
            {
                return static_cast<int32_t>(SM_ERR);
            }
            std::memcpy(out, r, sizeof r);
            return SMLIB_OK;
        }
    );
}

extern "C" int32_t smlib_ssi(const smlib_surface *a, const smlib_surface *b, int32_t method, double approx_tol, int32_t samples_per_curve,
                             smlib_ssi_result **out)
{
    if (a == nullptr || b == nullptr || out == nullptr || method < SMLIB_SSI_GLOBAL || method > SMLIB_SSI_BASIC ||
        !std::isfinite(approx_tol) || samples_per_curve < 2 || samples_per_curve > SMLIB_SSI_MAX_SAMPLES)
    {
        return SMLIB_ERR_ARG;
    }
    Call c("smlib_ssi", "method=" + std::to_string(method) + ", approx_tol=" + num(approx_tol) +
                            ", samples_per_curve=" + std::to_string(samples_per_curve));
    // A case file the corpus runner replays: ssi_runner <file> gsi|adv|basic <tol> <out.json>.
    c.input_text("case",
                 [&](const std::string &path)
                 {
                     if (FILE *f = std::fopen(path.c_str(), "w"))
                     {
                         std::fprintf(f, "{\"id\":\"smlib_ssi\",\"surfaces\":[%s,%s]}\n", a->json.c_str(), b->json.c_str());
                         std::fclose(f);
                     }
                 });
    return c.run(
        [&]
        {
            const SmContext &ctx = *SmApiGetOrCreateContext();
            SmApproxTol3d tol(approx_tol > 0.0 ? approx_tol : 0.0);
            const SmApproxTol3d *ptol = approx_tol > 0.0 ? &tol : nullptr;
            Curves cs;
            SmTArray<SmTsectCurveType> types;
            SmTArray<double> devs;
            SmStatus stat = SM_SUCCESS;
            if (method == SMLIB_SSI_GLOBAL)
            {
                SmBoolean edges[2] = {TRUE, TRUE};
                stat = a->s->GlobalSurfaceIntersect(ctx, a->domain, *b->s, b->domain, edges, ptol, NULL, &cs.c3, &cs.uv1, &cs.uv2, &types,
                                                    &devs);
            }
            else if (method == SMLIB_SSI_ADVANCED)
            {
                SmAdvSurfaceIntersector ssi(*a->s, a->domain, *b->s, b->domain);
                ssi.SetUseSurfaceEdges(TRUE, TRUE);
                stat = ssi.DoIntersection(ctx, ptol, NULL, &cs.c3, &cs.uv1, &cs.uv2, &types, &devs);
            }
            else
            {
                SmSurfaceIntersector ssi(*a->s, a->domain, *b->s, b->domain);
                ssi.SetUseSurfaceEdges(TRUE, TRUE);
                stat = ssi.DoIntersection(ctx, ptol, NULL, &cs.c3, &cs.uv1, &cs.uv2, &types, &devs);
            }
            int32_t st = status(stat);
            if (st != SMLIB_OK)
            {
                return st;
            }
            const ULONG n = cs.c3.GetSize();
            if (static_cast<int64_t>(n) * samples_per_curve > MAX_TOTAL_SAMPLES)
            {
                return static_cast<int32_t>(SMLIB_ERR_ARG);
            }
            std::unique_ptr<smlib_ssi_result> r(new smlib_ssi_result);
            r->samples = samples_per_curve;
            r->curves.resize(n);
            for (ULONG i = 0; i < n; ++i)
            {
                Curve &k = r->curves[i];
                k.type = i < types.GetSize() ? static_cast<int32_t>(types[i]) : -1;
                k.deviation = i < devs.GetSize() ? devs[i] : NaN;
                const SmCurve *ua = i < cs.uv1.GetSize() ? cs.uv1[i] : nullptr;
                const SmCurve *ub = i < cs.uv2.GetSize() ? cs.uv2[i] : nullptr;
                k.flags = (ua != nullptr ? SMLIB_SSI_HAS_UV_A : 0) | (ub != nullptr ? SMLIB_SSI_HAS_UV_B : 0);
                k.interval[0] = k.interval[1] = NaN;
                sample(cs.c3[i], samples_per_curve, 3, k.xyz, k.interval);
                sample(ua, samples_per_curve, 2, k.uv_a, nullptr);
                sample(ub, samples_per_curve, 2, k.uv_b, nullptr);
                evaluate(a->s, k.uv_a, k.a_at_uv);
                evaluate(b->s, k.uv_b, k.b_at_uv);
                if (k.type == SMLIB_SSI_TOUCHING)
                {
                    Point p;
                    std::memcpy(p.xyz, k.xyz.data(), sizeof p.xyz);
                    std::memcpy(p.uv_a, k.uv_a.data(), sizeof p.uv_a);
                    std::memcpy(p.uv_b, k.uv_b.data(), sizeof p.uv_b);
                    r->points.push_back(p);
                }
            }
            *out = r.release();
            return SMLIB_OK;
        }
    );
}

extern "C" void smlib_ssi_result_free(smlib_ssi_result *r)
{
    delete r;
}

extern "C" int32_t smlib_ssi_curve_count(const smlib_ssi_result *r, int32_t *out)
{
    if (r == nullptr || out == nullptr)
    {
        return SMLIB_ERR_ARG;
    }
    *out = static_cast<int32_t>(r->curves.size());
    return SMLIB_OK;
}

extern "C" int32_t smlib_ssi_curve_info(const smlib_ssi_result *r, int32_t i, int32_t *type, int32_t *flags, int32_t *samples,
                                        double *deviation, double interval[2])
{
    if (r == nullptr || i < 0 || static_cast<size_t>(i) >= r->curves.size())
    {
        return SMLIB_ERR_ARG;
    }
    const Curve &k = r->curves[i];
    if (type != nullptr)
    {
        *type = k.type;
    }
    if (flags != nullptr)
    {
        *flags = k.flags;
    }
    if (samples != nullptr)
    {
        *samples = r->samples;
    }
    if (deviation != nullptr)
    {
        *deviation = k.deviation;
    }
    if (interval != nullptr)
    {
        interval[0] = k.interval[0];
        interval[1] = k.interval[1];
    }
    return SMLIB_OK;
}

extern "C" int32_t smlib_ssi_curve_samples(const smlib_ssi_result *r, int32_t i, int64_t cap, double *xyz, double *uv_a, double *uv_b,
                                           double *a_at_uv, double *b_at_uv)
{
    if (r == nullptr || i < 0 || static_cast<size_t>(i) >= r->curves.size())
    {
        return SMLIB_ERR_ARG;
    }
    if (cap < r->samples)
    {
        return SMLIB_ERR_BUFFER;
    }
    const Curve &k = r->curves[i];
    auto copy = [](double *dst, const std::vector<double> &src)
    {
        if (dst != nullptr && !src.empty())
        {
            std::memcpy(dst, src.data(), src.size() * sizeof(double));
        }
    };
    copy(xyz, k.xyz);
    copy(uv_a, k.uv_a);
    copy(uv_b, k.uv_b);
    copy(a_at_uv, k.a_at_uv);
    copy(b_at_uv, k.b_at_uv);
    return SMLIB_OK;
}

extern "C" int32_t smlib_ssi_point_count(const smlib_ssi_result *r, int32_t *out)
{
    if (r == nullptr || out == nullptr)
    {
        return SMLIB_ERR_ARG;
    }
    *out = static_cast<int32_t>(r->points.size());
    return SMLIB_OK;
}

extern "C" int32_t smlib_ssi_point(const smlib_ssi_result *r, int32_t i, double xyz[3], double uv_a[2], double uv_b[2])
{
    if (r == nullptr || i < 0 || static_cast<size_t>(i) >= r->points.size())
    {
        return SMLIB_ERR_ARG;
    }
    const Point &p = r->points[i];
    if (xyz != nullptr)
    {
        std::memcpy(xyz, p.xyz, sizeof p.xyz);
    }
    if (uv_a != nullptr)
    {
        std::memcpy(uv_a, p.uv_a, sizeof p.uv_a);
    }
    if (uv_b != nullptr)
    {
        std::memcpy(uv_b, p.uv_b, sizeof p.uv_b);
    }
    return SMLIB_OK;
}
