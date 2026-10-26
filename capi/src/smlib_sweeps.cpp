// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// The sweeps of include/smlib.h: extrusion (straight or with draft), revolution, sweep along a
// path, loft through sections, the screw sweep of a thread, and the sweep through stations of a
// profile that twists, scales or follows rails. Profiles, paths and sections
// cross as curve lists (lines, circular arcs and NURBS, include/smlib.h); each call builds its
// own SMLib curves from them, so nothing the caller owns is touched.
//
// Every entry point checks its numbers before SMLib sees them, runs SM_API's sweep, and hands
// the result over only when it is what was asked for: a closed manifold solid when the ends are
// capped (SMLIB_ERR_NOT_SOLID otherwise), and never a body without faces (SMLIB_ERR_EMPTY).

#include "smlib_shim.h"
#include "smlib_pick.h"

#include <SmApiBrep.h>
#include <SmApiCurves.h>
#include <SmApiPrimitives.h>
#include <SmApiTrimmedSurfaces.h>
#include <SmBSplineCurve.h>
#include <SmBSplineSurface.h>
#include <NL_Geometry.h>
#include <NL_Nurbsdef.h>
#include <SmCircle.h>
#include <SmCylinder.h>
#include <SmFace.h>
#include <SmLine.h>
#include <SmPrimitiveCreation.h>

#include <memory>
#include <vector>

using namespace smlib_shim;

namespace
{

constexpr double TWO_PI = 6.283185307179586;

// A coordinate no larger than this (mm, 1000 km), so that SMLib's tolerances stay meaningful.
constexpr double MAX_COORD = 1e9;

// Owns curves made here; frees what is left of them.
struct Curves
{
    std::vector<SmCurve *> v;

    Curves() = default;
    Curves(const Curves &) = delete;
    Curves &operator=(const Curves &) = delete;
    ~Curves()
    {
        for (SmCurve *c : v)
        {
            delete c;
        }
    }

    SmTArray<SmCurve *> array(size_t from = 0, size_t to = SIZE_MAX) const
    {
        SmTArray<SmCurve *> a;
        for (size_t i = from; i < v.size() && i < to; ++i)
        {
            a.Add(v[i]);
        }
        return a;
    }
};

bool coord_ok(const double *p, int n)
{
    for (int i = 0; i < n; ++i)
    {
        if (!std::isfinite(p[i]) || std::fabs(p[i]) > MAX_COORD)
        {
            return false;
        }
    }
    return true;
}

double length(const double *v)
{
    return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

// A unit vector from `v`, or false when it has no usable length.
bool unit(const double *v, SmVector3d &out)
{
    if (!finite(v, 3))
    {
        return false;
    }
    double n = length(v);
    if (!(n > 1e-12) || !std::isfinite(n))
    {
        return false;
    }
    out = SmVector3d(v[0] / n, v[1] / n, v[2] / n);
    return true;
}

// The number of reals curve kind `kind` takes (NURBS: with `degree` and `count`), or -1.
int64_t reals_of(int64_t kind, int64_t degree, int64_t count)
{
    switch (kind)
    {
    case SMLIB_CURVE_LINE:
        return 6;
    case SMLIB_CURVE_ARC:
        return 12;
    case SMLIB_CURVE_NURBS:
        return 4 * count + count + degree + 1;
    default:
        return -1;
    }
}

int32_t line(const SmContext &ctx, const double *r, SmCurve *&out)
{
    if (!coord_ok(r, 6))
    {
        return SMLIB_ERR_ARG;
    }
    SmPoint3d a(r[0], r[1], r[2]);
    SmPoint3d b(r[3], r[4], r[5]);
    if (!(a.DistanceBetween(b) > 1e-9))
    {
        return SMLIB_ERR_ARG;
    }
    SmLine *l = nullptr;
    SmStatus s = SmLine::CreateLineSegment(ctx, 3, a, b, l);
    if (s != SM_SUCCESS || l == nullptr)
    {
        delete l;
        return s != SM_SUCCESS ? status(s) : SM_ERR;
    }
    out = l;
    return SMLIB_OK;
}

// A circular arc: center, unit axis, unit reference direction (perpendicular to the axis),
// radius, start and end angles (radians, counter-clockwise about the axis from the reference).
int32_t arc(const SmContext &ctx, const double *r, SmCurve *&out)
{
    if (!coord_ok(r, 3) || !finite(r + 3, 9))
    {
        return SMLIB_ERR_ARG;
    }
    SmVector3d axis, ref;
    if (!unit(r + 3, axis) || !unit(r + 6, ref) || std::fabs(axis.Dot(ref)) > 1e-9)
    {
        return SMLIB_ERR_ARG;
    }
    double radius = r[9];
    double start = r[10];
    double sweep = r[11] - r[10];
    if (!(radius > 1e-9 && radius < MAX_COORD) || std::fabs(start) > 1e6 || !(sweep > 1e-9 && sweep <= TWO_PI * (1.0 + 1e-12)))
    {
        return SMLIB_ERR_ARG;
    }
    // Start the circle's own parameter at the arc's start: SMLib takes intervals in degrees
    // within [-360, 360].
    SmVector3d y = axis * ref;
    SmVector3d x0 = ref * std::cos(start) + y * std::sin(start);
    SmVector3d y0 = axis * x0;
    SmAxis2Placement frame(SmPoint3d(r[0], r[1], r[2]), x0, y0);
    SmExtent1d degrees(0.0, std::min(360.0, sweep * 180.0 / SM_PI));
    SmCircle *c = nullptr;
    SmStatus s = SmCircle::CreateCanonical(ctx, frame, radius, c, &degrees);
    if (s != SM_SUCCESS || c == nullptr)
    {
        delete c;
        return s != SM_SUCCESS ? status(s) : SM_ERR;
    }
    out = c;
    return SMLIB_OK;
}

// A NURBS curve: `count` control points (x, y, z, weight), then count + degree + 1 knots,
// clamped (end multiplicity degree + 1, inner ones at most degree), non-decreasing.
int32_t nurbs(const SmContext &ctx, int64_t degree, int64_t count, const double *r, SmCurve *&out)
{
    if (degree < 1 || degree > SMLIB_CURVE_MAX_DEGREE || count < degree + 1 || count > SMLIB_CURVE_MAX_POINTS)
    {
        return SMLIB_ERR_ARG;
    }
    SmTArray<SmPoint3d> points;
    SmTArray<double> weights;
    bool rational = false;
    for (int64_t i = 0; i < count; ++i)
    {
        const double *p = r + 4 * i;
        if (!coord_ok(p, 3) || !std::isfinite(p[3]) || !(p[3] > 1e-12 && p[3] < 1e12))
        {
            return SMLIB_ERR_ARG;
        }
        points.Add(SmPoint3d(p[0], p[1], p[2]));
        weights.Add(p[3]);
        rational = rational || p[3] != 1.0;
    }
    const double *k = r + 4 * count;
    int64_t nk = count + degree + 1;
    if (!finite(k, static_cast<int>(nk)))
    {
        return SMLIB_ERR_ARG;
    }
    SmTArray<double> unique;
    SmTArray<ULONG> mult;
    for (int64_t i = 0; i < nk; ++i)
    {
        if (i > 0 && k[i] < k[i - 1])
        {
            return SMLIB_ERR_ARG;
        }
        if (i > 0 && k[i] == k[i - 1])
        {
            mult[mult.GetSize() - 1] += 1;
        }
        else
        {
            unique.Add(k[i]);
            mult.Add(1);
        }
    }
    ULONG d = static_cast<ULONG>(degree);
    ULONG n = unique.GetSize();
    if (n < 2 || mult[0] != d + 1 || mult[n - 1] != d + 1 || !(unique[n - 1] - unique[0] > 1e-12))
    {
        return SMLIB_ERR_ARG;
    }
    for (ULONG i = 1; i + 1 < n; ++i)
    {
        if (mult[i] > d)
        {
            return SMLIB_ERR_ARG;
        }
    }
    SmBSplineCurve *c = nullptr;
    SmStatus s = SmBSplineCurve::CreateCanonical(ctx, 3, d, points, SM_CF_UNSPECIFIED, mult, unique, SM_KT_UNSPECIFIED,
                                                 rational ? &weights : nullptr, nullptr, c);
    if (s != SM_SUCCESS || c == nullptr)
    {
        delete c;
        return s != SM_SUCCESS ? status(s) : SM_ERR;
    }
    out = c;
    return SMLIB_OK;
}

// The curves of a curve list (include/smlib.h). Both arrays must be used up exactly.
int32_t parse(const int64_t *kinds, int64_t n_kinds, const double *reals, int64_t n_reals, Curves &out)
{
    if (n_kinds < 1 || n_kinds > 3 * SMLIB_CURVE_MAX_CURVES || kinds == nullptr || n_reals < 0 || n_reals > SMLIB_FLAT_MAX_LEN ||
        (n_reals > 0 && reals == nullptr))
    {
        return SMLIB_ERR_ARG;
    }
    const SmContext &ctx = *SmApiGetOrCreateContext();
    int64_t i = 0;
    int64_t at = 0;
    while (i < n_kinds)
    {
        if (static_cast<int64_t>(out.v.size()) >= SMLIB_CURVE_MAX_CURVES)
        {
            return SMLIB_ERR_ARG;
        }
        int64_t kind = kinds[i++];
        int64_t degree = 0;
        int64_t count = 0;
        if (kind == SMLIB_CURVE_NURBS)
        {
            if (i + 2 > n_kinds)
            {
                return SMLIB_ERR_ARG;
            }
            degree = kinds[i++];
            count = kinds[i++];
            if (degree < 1 || degree > SMLIB_CURVE_MAX_DEGREE || count < 2 || count > SMLIB_CURVE_MAX_POINTS)
            {
                return SMLIB_ERR_ARG;
            }
        }
        int64_t need = reals_of(kind, degree, count);
        if (need < 0 || at + need > n_reals)
        {
            return SMLIB_ERR_ARG;
        }
        SmCurve *c = nullptr;
        int32_t s = kind == SMLIB_CURVE_LINE  ? line(ctx, reals + at, c)
                    : kind == SMLIB_CURVE_ARC ? arc(ctx, reals + at, c)
                                              : nurbs(ctx, degree, count, reals + at, c);
        if (s != SMLIB_OK)
        {
            return s;
        }
        out.v.push_back(c);
        at += need;
    }
    return at == n_reals ? SMLIB_OK : SMLIB_ERR_ARG;
}

// A curve list for the log line: its size.
std::string list(int64_t n_kinds, int64_t n_reals)
{
    return "curves[kinds=" + std::to_string(n_kinds) + ", reals=" + std::to_string(n_reals) + "]";
}

// Dumps a curve list as text (kinds, then reals) for replay.
void dump_curves(Call &c, const char *name, const int64_t *kinds, int64_t n_kinds, const double *reals, int64_t n_reals)
{
    c.input_text(name,
                 [=](const std::string &path)
                 {
                     if (FILE *f = std::fopen(path.c_str(), "w"))
                     {
                         std::fprintf(f, "kinds %lld\n", static_cast<long long>(n_kinds));
                         for (int64_t i = 0; kinds != nullptr && i < n_kinds; ++i)
                         {
                             std::fprintf(f, "%lld\n", static_cast<long long>(kinds[i]));
                         }
                         std::fprintf(f, "reals %lld\n", static_cast<long long>(n_reals));
                         for (int64_t i = 0; reals != nullptr && i < n_reals; ++i)
                         {
                             std::fprintf(f, "%.17g\n", reals[i]);
                         }
                         std::fclose(f);
                     }
                 });
}

// Hands `r` to the caller when it is what was asked for: a closed manifold solid when `solid`.
int32_t finish(SmBrep *r, bool solid, smlib_brep **out)
{
    std::unique_ptr<SmBrep> owned(r);
    if (r == nullptr || r->GetNumFaces() == 0)
    {
        return SMLIB_ERR_EMPTY;
    }
    if (solid && !r->IsManifoldSolid())
    {
        return SMLIB_ERR_NOT_SOLID;
    }
    *out = new smlib_brep{owned.release()};
    return SMLIB_OK;
}

// Whether the curves, in order, run end to end round one closed loop (within `tol`).
bool one_loop(const Curves &c, double tol)
{
    SmPoint3d first, prev;
    for (size_t i = 0; i < c.v.size(); ++i)
    {
        SmPoint3d a, b;
        c.v[i]->GetEnds(a, b);
        if (i == 0)
        {
            first = a;
        }
        else if (a.DistanceBetween(prev) > tol)
        {
            return false;
        }
        prev = b;
    }
    return !c.v.empty() && prev.DistanceBetween(first) <= tol;
}

// Whether the path is lines end to end along one direction: then `dir` is that unit direction
// and `distance` the path's length.
bool straight(const Curves &path, SmVector3d &dir, double &distance)
{
    distance = 0.0;
    SmPoint3d prev;
    for (size_t i = 0; i < path.v.size(); ++i)
    {
        if (SM_CAST_PTR(SmLine, path.v[i]) == nullptr)
        {
            return false;
        }
        SmPoint3d a, b;
        path.v[i]->GetEnds(a, b);
        if (i > 0 && a.DistanceBetween(prev) > 1e-9)
        {
            return false;
        }
        prev = b;
        SmVector3d d = b - a;
        double l = d.Length();
        if (!(l > 1e-12))
        {
            return false;
        }
        d = d / l;
        if (i == 0)
        {
            dir = d;
        }
        else if (d.Dot(dir) < 1.0 - 1e-12)
        {
            return false;
        }
        distance += l;
    }
    return !path.v.empty();
}

// Whether the path (parsed as `path` from its curve list) is arcs end to end about one axis: then
// `centre`, the unit `axis` they turn about counter-clockwise and the total `angle` (radians).
bool circular(const Curves &path, const int64_t *kinds, int64_t n_kinds, const double *reals, SmPoint3d &centre, SmVector3d &axis,
              double &angle)
{
    angle = 0.0;
    if (kinds == nullptr || reals == nullptr || n_kinds != static_cast<int64_t>(path.v.size()) || path.v.empty())
    {
        return false;
    }
    SmPoint3d prev;
    double radius = 0.0;
    for (int64_t i = 0; i < n_kinds; ++i)
    {
        if (kinds[i] != SMLIB_CURVE_ARC)
        {
            return false;
        }
        const double *r = reals + 12 * i;
        SmPoint3d c(r[0], r[1], r[2]);
        SmVector3d k(r[3], r[4], r[5]);
        k = k / k.Length();
        SmPoint3d a, b;
        path.v[static_cast<size_t>(i)]->GetEnds(a, b);
        if (i == 0)
        {
            centre = c;
            axis = k;
            radius = r[9];
        }
        else if (c.DistanceBetween(centre) > 1e-9 || k.Dot(axis) < 1.0 - 1e-12 || std::fabs(r[9] - radius) > 1e-9 ||
                 a.DistanceBetween(prev) > 1e-9)
        {
            return false;
        }
        prev = b;
        angle += r[11] - r[10];
    }
    return angle <= TWO_PI * (1.0 + 1e-12);
}

// A path of end-to-end curves as one B-spline (CreateCurveSweep takes a single curve).
int32_t joined(const Curves &path, std::unique_ptr<SmBSplineCurve> &out)
{
    const SmContext &ctx = *SmApiGetOrCreateContext();
    SmTArray<SmBSplineCurve *> parts;
    std::vector<std::unique_ptr<SmBSplineCurve>> keep;
    for (SmCurve *c : path.v)
    {
        SmBSplineCurve *b = SM_CAST_PTR(SmBSplineCurve, c);
        if (b == nullptr)
        {
            return SMLIB_ERR_ARG;
        }
        // A plain B-spline copy: SMLib's lines and circles keep an analytic parameterisation.
        keep.emplace_back(new (ctx) SmBSplineCurve(*b));
        parts.Add(keep.back().get());
    }
    if (parts.GetSize() == 1)
    {
        out = std::move(keep.front());
        return SMLIB_OK;
    }
    SmBSplineCurve *j = nullptr;
    SmStatus s = SmBSplineCurve::CreateByJoining(ctx, parts, nullptr, j);
    if (s != SM_SUCCESS || j == nullptr)
    {
        delete j;
        return s != SM_SUCCESS ? status(s) : SM_ERR;
    }
    out.reset(j);
    return SMLIB_OK;
}

} // namespace

extern "C" int32_t smlib_extrude(const int64_t *kinds, int64_t n_kinds, const double *reals, int64_t n_reals, const double direction[3],
                                 double distance, double draft, int32_t corner, int32_t cap, smlib_brep **out, double out_corner[3])
{
    if (out == nullptr || direction == nullptr || !std::isfinite(distance) || !(distance > 1e-9 && distance < MAX_COORD) ||
        !std::isfinite(draft) || std::fabs(draft) >= SM_PI / 2.0 - 1e-6 || (cap != 0 && cap != 1) ||
        (corner != SMLIB_DRAFT_CORNER_CURVATURE && corner != SMLIB_DRAFT_CORNER_TANGENT && corner != SMLIB_DRAFT_CORNER_ROUND))
    {
        return SMLIB_ERR_ARG;
    }
    Call c("smlib_extrude", list(n_kinds, n_reals) + ", direction=" + vec(direction, 3) + ", distance=" + num(distance) +
                                ", draft=" + num(draft) + ", corner=" + std::to_string(corner) + ", cap=" + std::to_string(cap));
    dump_curves(c, "profile", kinds, n_kinds, reals, n_reals);
    return c.run(
        [&]
        {
            SmVector3d dir;
            if (!unit(direction, dir))
            {
                return SMLIB_ERR_ARG;
            }
            Curves profile;
            int32_t s = parse(kinds, n_kinds, reals, n_reals, profile);
            if (s != SMLIB_OK)
            {
                return s;
            }
            SmBrep *r = nullptr;
            if (draft == 0.0)
            {
                s = status(SmApiCreateLinearSweep(profile.array(), dir, distance, cap ? TRUE : FALSE, r));
            }
            else
            {
                // SMLib takes the draft as the walls' angle to the profile's normal, in degrees:
                // positive leaning in, negative leaning out (SmApiBrep.h documents 90 to 180 for
                // outward, but its offset only reads the angle's sign and tangent). It finds the
                // profile's normal itself; the height's sign picks the side, so it is matched to
                // `direction` here.
                // SMLib's draft extrusion takes a single chain of curves (it asserts one).
                if (!one_loop(profile, 1e-6))
                {
                    return SMLIB_ERR_ARG;
                }
                SmTArray<SmCurve *> a = profile.array();
                SmVector3d normal;
                if (SmCurve::ComputeCrvsNormal(a, normal) != SM_SUCCESS)
                {
                    return SMLIB_ERR_ARG;
                }
                double along = normal.Dot(dir);
                if (std::fabs(along) < 1.0 - 1e-9)
                {
                    // A draft extrusion runs square to the profile.
                    return SMLIB_ERR_ARG;
                }
                double height = along > 0.0 ? distance : -distance;
                double angle = draft * 180.0 / SM_PI;
                SmOffsetCornerType rule = corner == SMLIB_DRAFT_CORNER_TANGENT ? SM_OC_LINEAR_EXTENSION
                                          : corner == SMLIB_DRAFT_CORNER_ROUND ? SM_OC_FILLET_CORNER
                                                                               : SM_OC_CURVE_EXTENSION;
                SmTaperCornerRefusal refusal = SM_TCR_NONE;
                SmPoint3d at(0.0, 0.0, 0.0);
                s = status(SmApiCreateTaperExtrude(a, height, angle, cap ? 3 : 0, rule, r, &refusal, &at));
                if (refusal != SM_TCR_NONE)
                {
                    delete r;
                    if (out_corner != nullptr)
                    {
                        out_corner[0] = at.x;
                        out_corner[1] = at.y;
                        out_corner[2] = at.z;
                    }
                    return refusal == SM_TCR_APART       ? SMLIB_ERR_CORNER_APART
                           : refusal == SM_TCR_WALL_CLOSES ? SMLIB_ERR_WALL_CLOSES
                                                           : SMLIB_ERR_NOT_EXACT;
                }
            }
            if (s != SMLIB_OK)
            {
                delete r;
                return s;
            }
            return finish(r, cap != 0, out);
        }
    );
}

extern "C" int32_t smlib_revolve(const int64_t *kinds, int64_t n_kinds, const double *reals, int64_t n_reals, const double origin[3],
                                 const double axis[3], double angle, int32_t cap, smlib_brep **out)
{
    if (out == nullptr || origin == nullptr || axis == nullptr || !coord_ok(origin, 3) || !std::isfinite(angle) ||
        !(angle > 1e-9 && angle <= TWO_PI * (1.0 + 1e-12)) || (cap != 0 && cap != 1))
    {
        return SMLIB_ERR_ARG;
    }
    Call c("smlib_revolve", list(n_kinds, n_reals) + ", origin=" + vec(origin, 3) + ", axis=" + vec(axis, 3) + ", angle=" + num(angle) +
                                ", cap=" + std::to_string(cap));
    dump_curves(c, "profile", kinds, n_kinds, reals, n_reals);
    return c.run(
        [&]
        {
            SmVector3d dir;
            if (!unit(axis, dir))
            {
                return SMLIB_ERR_ARG;
            }
            Curves profile;
            int32_t s = parse(kinds, n_kinds, reals, n_reals, profile);
            if (s != SMLIB_OK)
            {
                return s;
            }
            double degrees = std::min(360.0, angle * 180.0 / SM_PI);
            SmBrep *r = nullptr;
            s = status(SmApiCreateRotationalSweep(profile.array(), SmPoint3d(origin[0], origin[1], origin[2]), dir, degrees, cap ? TRUE : FALSE, r));
            if (s != SMLIB_OK)
            {
                delete r;
                return s;
            }
            return finish(r, cap != 0, out);
        }
    );
}

namespace
{

// Writes the motion `k` (0 first, 1 middle, 2 last) of smlib_sweep_placed: the columns, then the
// translation.
void put_motion(double *placed, int k, const SmVector3d &x, const SmVector3d &y, const SmVector3d &z, const SmPoint3d &t)
{
    double *m = placed + 12 * k;
    const SmVector3d c[3] = {x, y, z};
    for (int i = 0; i < 3; ++i)
    {
        m[3 * i] = c[i].x;
        m[3 * i + 1] = c[i].y;
        m[3 * i + 2] = c[i].z;
    }
    m[9] = t.x;
    m[10] = t.y;
    m[11] = t.z;
}

// `v` turned `angle` right-handed about the unit `k` (Rodrigues).
SmVector3d turned(const SmVector3d &v, const SmVector3d &k, double angle)
{
    return v * std::cos(angle) + (k * v) * std::sin(angle) + k * (k.Dot(v) * (1.0 - std::cos(angle)));
}

} // namespace

extern "C" int32_t smlib_sweep(const int64_t *kinds, int64_t n_kinds, const double *reals, int64_t n_reals, const int64_t *path_kinds,
                               int64_t n_path_kinds, const double *path_reals, int64_t n_path_reals, int32_t frame, int32_t cap,
                               smlib_brep **out)
{
    return smlib_sweep_placed(kinds, n_kinds, reals, n_reals, path_kinds, n_path_kinds, path_reals, n_path_reals, frame, cap, out, nullptr);
}

extern "C" int32_t smlib_sweep_placed(const int64_t *kinds, int64_t n_kinds, const double *reals, int64_t n_reals,
                                      const int64_t *path_kinds, int64_t n_path_kinds, const double *path_reals, int64_t n_path_reals,
                                      int32_t frame, int32_t cap, smlib_brep **out, double *placed)
{
    if (out == nullptr || (frame != SMLIB_SWEEP_FOLLOW && frame != SMLIB_SWEEP_TRANSLATE) || (cap != 0 && cap != 1))
    {
        return SMLIB_ERR_ARG;
    }
    Call c("smlib_sweep", list(n_kinds, n_reals) + ", path " + list(n_path_kinds, n_path_reals) + ", frame=" + std::to_string(frame) +
                              ", cap=" + std::to_string(cap));
    dump_curves(c, "profile", kinds, n_kinds, reals, n_reals);
    dump_curves(c, "path", path_kinds, n_path_kinds, path_reals, n_path_reals);
    return c.run(
        [&]
        {
            Curves profile;
            Curves path;
            int32_t s = parse(kinds, n_kinds, reals, n_reals, profile);
            if (s == SMLIB_OK)
            {
                s = parse(path_kinds, n_path_kinds, path_reals, n_path_reals, path);
            }
            if (s != SMLIB_OK)
            {
                return s;
            }
            // Along a straight path the sweep is an extrusion: SM_API's linear sweep builds the walls
            // exactly (planes, cylinders) where the curve sweep fits B-splines to them.
            SmVector3d dir;
            double distance = 0.0;
            if (straight(path, dir, distance))
            {
                SmBrep *r = nullptr;
                s = status(SmApiCreateLinearSweep(profile.array(), dir, distance, cap ? TRUE : FALSE, r));
                if (s != SMLIB_OK)
                {
                    delete r;
                    return s;
                }
                if (placed != nullptr)
                {
                    for (int k = 0; k < 3; ++k)
                    {
                        const SmVector3d t = dir * (distance * 0.5 * k);
                        put_motion(placed, k, SmVector3d(1, 0, 0), SmVector3d(0, 1, 0), SmVector3d(0, 0, 1), SmPoint3d(t.x, t.y, t.z));
                    }
                }
                return finish(r, cap != 0, out);
            }
            // Turning with a path of arcs about one axis, the profile turns about it: SM_API's
            // rotational sweep builds the walls exactly (tori, cones, planes, ...).
            SmPoint3d centre;
            SmVector3d axis;
            double angle = 0.0;
            if (frame == SMLIB_SWEEP_FOLLOW && circular(path, path_kinds, n_path_kinds, path_reals, centre, axis, angle))
            {
                SmBrep *r = nullptr;
                double degrees = std::min(360.0, angle * 180.0 / SM_PI);
                s = status(SmApiCreateRotationalSweep(profile.array(), centre, axis, degrees, cap ? TRUE : FALSE, r));
                if (s != SMLIB_OK)
                {
                    delete r;
                    return s;
                }
                if (placed != nullptr)
                {
                    SmVector3d k = axis;
                    if (k.Unitize() != SM_SUCCESS)
                    {
                        delete r;
                        return SMLIB_ERR_ARG;
                    }
                    const double turn = degrees * SM_PI / 180.0;
                    for (int j = 0; j < 3; ++j)
                    {
                        const double a = turn * 0.5 * j;
                        const SmVector3d c(centre.x, centre.y, centre.z);
                        const SmVector3d t = c - turned(c, k, a);
                        put_motion(placed, j, turned(SmVector3d(1, 0, 0), k, a), turned(SmVector3d(0, 1, 0), k, a), turned(SmVector3d(0, 0, 1), k, a),
                                   SmPoint3d(t.x, t.y, t.z));
                    }
                }
                return finish(r, cap != 0, out);
            }
            std::unique_ptr<SmBSplineCurve> along;
            s = joined(path, along);
            if (s != SMLIB_OK)
            {
                return s;
            }
            // Round a closed path the profile comes back onto itself: the walls close on each
            // other there, where caps would lie on top of each other inside the solid.
            SmPoint3d first, last, ignored;
            path.v.front()->GetEnds(first, ignored);
            path.v.back()->GetEnds(ignored, last);
            const bool closed = first.DistanceBetween(last) <= 1e-6;
            SmSweepOptions options;
            options.bCapEndsArg = cap && !closed ? TRUE : FALSE;
            options.bTranslationalSweep = frame == SMLIB_SWEEP_TRANSLATE ? TRUE : FALSE;
            options.bProfileAtPathStart = TRUE;
            SmTArray<SmCurve *> a = profile.array();
            SmBrep *r = nullptr;
            SmAxis2Placement motions[3];
            s = status(SmApiCreateCurveSweep(a, along.get(), nullptr, nullptr, &options, r, nullptr, nullptr, nullptr, placed != nullptr ? motions : nullptr));
            if (s == SMLIB_OK && placed != nullptr)
            {
                for (int k = 0; k < 3; ++k)
                {
                    put_motion(placed, k, motions[k].GetXAxis(), motions[k].GetYAxis(), motions[k].GetZAxis(), motions[k].GetOriginRef());
                }
            }
            if (s == SMLIB_OK && closed)
            {
                ULONG stitched = 0, lamina = 0;
                double vertex_gap = 0.0, edge_gap = 0.0;
                s = status(r->StitchFaces(r->GetTolerance(), stitched, lamina, vertex_gap, edge_gap));
                if (s == SMLIB_OK && lamina != 0)
                {
                    s = SMLIB_ERR_NOT_SOLID;
                }
            }
            if (s != SMLIB_OK)
            {
                delete r;
                return s;
            }
            return finish(r, cap != 0, out);
        }
    );
}

extern "C" int32_t smlib_loft(const int64_t *kinds, int64_t n_kinds, const double *reals, int64_t n_reals, const int64_t *section_sizes,
                              int64_t n_sections, int32_t degree, int32_t cap, smlib_brep **out)
{
    if (out == nullptr || section_sizes == nullptr || n_sections < 2 || n_sections > SMLIB_CURVE_MAX_CURVES || degree < 1 || degree > 3 ||
        (cap != 0 && cap != 1))
    {
        return SMLIB_ERR_ARG;
    }
    Call c("smlib_loft", list(n_kinds, n_reals) + ", sections=" + std::to_string(n_sections) + ", degree=" + std::to_string(degree) +
                             ", cap=" + std::to_string(cap));
    dump_curves(c, "sections", kinds, n_kinds, reals, n_reals);
    return c.run(
        [&]
        {
            Curves all;
            int32_t s = parse(kinds, n_kinds, reals, n_reals, all);
            if (s != SMLIB_OK)
            {
                return s;
            }
            SmTArray<ULONG> counts;
            int64_t total = 0;
            for (int64_t i = 0; i < n_sections; ++i)
            {
                if (section_sizes[i] < 1 || section_sizes[i] > SMLIB_CURVE_MAX_CURVES)
                {
                    return SMLIB_ERR_ARG;
                }
                total += section_sizes[i];
                counts.Add(static_cast<ULONG>(section_sizes[i]));
            }
            if (total != static_cast<int64_t>(all.v.size()))
            {
                return SMLIB_ERR_ARG;
            }
            SmTArray<const SmCurve *> curves;
            for (SmCurve *x : all.v)
            {
                curves.Add(x);
            }
            SmBrep *r = nullptr;
            s = status(SmApiCreateSkinPrimitive(counts, curves, static_cast<ULONG>(degree), cap ? TRUE : FALSE, r, nullptr, nullptr, nullptr));
            if (s != SMLIB_OK)
            {
                delete r;
                return s;
            }
            return finish(r, cap != 0, out);
        }
    );
}

extern "C" int32_t smlib_loft_to_point(const int64_t *kinds, int64_t n_kinds, const double *reals, int64_t n_reals, const double apex[3],
                                       smlib_brep **out)
{
    if (out == nullptr || apex == nullptr || !coord_ok(apex, 3))
    {
        return SMLIB_ERR_ARG;
    }
    Call c("smlib_loft_to_point", list(n_kinds, n_reals) + ", apex=" + vec(apex, 3));
    dump_curves(c, "section", kinds, n_kinds, reals, n_reals);
    return c.run(
        [&]
        {
            Curves section;
            int32_t s = parse(kinds, n_kinds, reals, n_reals, section);
            if (s != SMLIB_OK)
            {
                return s;
            }
            if (section.v.size() < 3 || !one_loop(section, 1e-6))
            {
                return SMLIB_ERR_ARG;
            }
            // The base's corners, each line's start.
            std::vector<SmPoint3d> corners;
            for (SmCurve *x : section.v)
            {
                if (SM_CAST_PTR(SmLine, x) == nullptr)
                {
                    return SMLIB_ERR_ARG;
                }
                SmPoint3d a, b;
                x->GetEnds(a, b);
                corners.push_back(a);
            }
            SmTArray<SmCurve *> loop = section.array();
            SmVector3d normal;
            if (SmCurve::ComputeCrvsNormal(loop, normal) != SM_SUCCESS)
            {
                return SMLIB_ERR_ARG;
            }
            SmPoint3d top(apex[0], apex[1], apex[2]);
            double size = 0.0;
            for (const SmPoint3d &q : corners)
            {
                size = std::max(size, q.DistanceBetween(corners.front()));
            }
            SmVector3d up = top - corners.front();
            if (!(std::fabs(up.Dot(normal)) > 1e-9 * std::max(size, 1.0)))
            {
                return SMLIB_ERR_ARG;
            }
            const SmContext &ctx = *SmApiGetOrCreateContext();
            std::unique_ptr<SmBrep> brep(new (ctx) SmBrep());
            SmRegion *region = brep->GetInfiniteRegion();
            const double tol = 1e-6;
            // Each face is made of new lines it takes over.
            auto face = [&](const std::vector<SmPoint3d> &pts) -> int32_t
            {
                SmTArray<SmCurve *> edges;
                for (size_t i = 0; i < pts.size(); ++i)
                {
                    SmLine *l = nullptr;
                    SmStatus ls = SmLine::CreateLineSegment(ctx, 3, pts[i], pts[(i + 1) % pts.size()], l);
                    if (ls != SM_SUCCESS || l == nullptr)
                    {
                        delete l;
                        for (ULONG k = 0; k < edges.GetSize(); ++k)
                        {
                            delete edges[k];
                        }
                        return ls != SM_SUCCESS ? status(ls) : static_cast<int32_t>(SM_ERR);
                    }
                    edges.Add(l);
                }
                SmFace *f = nullptr;
                return status(brep->CreatePlanarFaceWith3DCurves(region, edges, tol, f));
            };
            s = face(corners);
            for (size_t i = 0; s == SMLIB_OK && i < corners.size(); ++i)
            {
                s = face({corners[i], corners[(i + 1) % corners.size()], top});
            }
            if (s == SMLIB_OK)
            {
                s = status(brep->StitchAndOrient());
            }
            if (s != SMLIB_OK)
            {
                return s;
            }
            return finish(brep.release(), true, out);
        }
    );
}

namespace
{

// The screw motion of smlib_screw_sweep: about the axis through `o` along unit `k`, rising
// `pitch` per turn, the section carried by angle a (radians, right-handed about k).
struct Screw
{
    SmPoint3d o;
    SmVector3d k;
    double pitch = 0.0;

    // Rodrigues: v turned right-handed about k by angle a.
    SmVector3d turn(const SmVector3d &v, double a) const
    {
        return v * std::cos(a) + (k * v) * std::sin(a) + k * (k.Dot(v) * (1.0 - std::cos(a)));
    }
    double lift(double a) const
    {
        return pitch * a / TWO_PI;
    }
    // Where p goes at angle a.
    SmPoint3d at(const SmPoint3d &p, double a) const
    {
        return o + turn(p - o, a) + k * lift(a);
    }
    // How far p lies from the axis.
    double reach(const SmPoint3d &p) const
    {
        SmVector3d w = p - o;
        return (w - k * w.Dot(k)).Length();
    }
    // The motion as a placement of the section: p -> o + R (p - o) + k lift, its origin o - R o + k lift.
    SmAxis2Placement placement(double a) const
    {
        SmPoint3d origin = o - turn(SmVector3d(o.x, o.y, o.z), a) + k * lift(a);
        return SmAxis2Placement(origin, turn(SmVector3d(1, 0, 0), a), turn(SmVector3d(0, 1, 0), a));
    }
};

// The stations of a piece of the turns: n spans of h radians each, from t[0] to t[n]. Each span
// carries the profile round as an exact rational quadratic arc does a point (its middle control
// point at the span's middle angle, out by 1 / cos(h / 2), weighted cos(h / 2)), and lifts it by
// the same rational quadratic of the lifts at the span's ends and middle. So every point of a
// wall lies exactly as far from the axis as its profile point and exactly at the angle the
// arc turns it to, only lifted off the screw surface there, the same for the whole profile: by
// at most |pitch| / (2 pi) * lift_error(h).
struct Stations
{
    std::vector<double> t;
    double half = 0.0; // h / 2

    Stations(double from, double to, size_t n)
    {
        for (size_t i = 0; i <= n; ++i)
        {
            t.push_back(i == n ? to : from + (to - from) * static_cast<double>(i) / static_cast<double>(n));
        }
        half = 0.5 * (to - from) / static_cast<double>(n);
    }
    size_t spans() const
    {
        return t.size() - 1;
    }
    // The knots across the stations, as unique values and multiplicities: double inside (the
    // arcs meet tangent), triple at the ends.
    void knots(SmTArray<double> &unique, SmTArray<ULONG> &mult) const
    {
        for (size_t i = 0; i < t.size(); ++i)
        {
            unique.Add(t[i]);
            mult.Add(i == 0 || i + 1 == t.size() ? 3 : 2);
        }
    }
    // The weight of control point j (of 2 * spans + 1) across the stations.
    double weight(size_t j) const
    {
        return j % 2 == 0 ? 1.0 : std::cos(half);
    }
};

// How far, in radians of turn, the rational quadratic lift of a span of half-angle `half` leads
// or lags the true lift, at most. A point at s in [-1, 1] across the span turns
// 2 atan(tan(half / 2) s) from the span's middle and is lifted as by half * s / (1 - (1 - s^2)
// (1 - cos half) / 2) radians; their difference is smooth, odd in s and nought at s = 0 and +-1
// (about 0.008 h^3), so a fine sampling finds its largest.
double lift_error(double half)
{
    double worst = 0.0;
    for (int i = 0; i <= 2048; ++i)
    {
        double s = static_cast<double>(i) / 2048.0;
        double lift = half * s / (1.0 - (1.0 - s * s) * (1.0 - std::cos(half)) / 2.0);
        worst = std::max(worst, std::fabs(lift - 2.0 * std::atan(std::tan(half / 2.0) * s)));
    }
    return worst;
}

// The control points of p's path under the screw motion over `st` (2 * spans + 1 of them).
std::vector<SmPoint3d> path(const Screw &m, const SmPoint3d &p, const Stations &st)
{
    SmVector3d w = p - m.o;
    double axial = w.Dot(m.k);
    SmVector3d radial = w - m.k * axial;
    std::vector<SmPoint3d> D;
    for (size_t i = 0; i <= st.spans(); ++i)
    {
        D.push_back(m.at(p, st.t[i]));
        if (i < st.spans())
        {
            double mid = 0.5 * (st.t[i] + st.t[i + 1]);
            D.push_back(m.o + m.turn(radial, mid) * (1.0 / std::cos(st.half)) + m.k * (axial + m.lift(mid)));
        }
    }
    return D;
}

// A profile curve as the NURBS its control points pair up in.
struct Nurbs
{
    ULONG degree = 0;
    SmTArray<SmPoint3d> points;
    SmTArray<ULONG> mult;
    SmTArray<double> unique;
    SmTArray<double> weights;
};

int32_t nurbs_of(const SmCurve *x, Nurbs &out)
{
    const SmBSplineCurve *b = SM_CAST_PTR(SmBSplineCurve, x);
    ULONG dim = 0;
    SmBSplineCurveForm form = SM_CF_UNSPECIFIED;
    SmKnotType kt = SM_KT_UNSPECIFIED;
    if (b == nullptr || b->GetCanonical(dim, out.degree, out.points, form, out.mult, out.unique, kt, out.weights) != SM_SUCCESS || dim != 3)
    {
        return static_cast<int32_t>(SM_ERR);
    }
    return SMLIB_OK;
}

// The path of p over `st` as a curve: the edge a profile vertex traces, exactly as the walls
// on either side of it have it.
int32_t path_curve(const SmContext &ctx, const Screw &m, const SmPoint3d &p, const Stations &st, SmCurve *&out)
{
    std::vector<SmPoint3d> D = path(m, p, st);
    SmTArray<SmPoint3d> points;
    SmTArray<double> weights;
    for (size_t j = 0; j < D.size(); ++j)
    {
        points.Add(D[j]);
        weights.Add(st.weight(j));
    }
    SmTArray<ULONG> mult;
    SmTArray<double> knots;
    st.knots(knots, mult);
    SmBSplineCurve *c = nullptr;
    SmStatus status_ = SmBSplineCurve::CreateCanonical(ctx, 3, 2, points, SM_CF_UNSPECIFIED, mult, knots, SM_KT_UNSPECIFIED, &weights, nullptr, c);
    if (status_ != SM_SUCCESS || c == nullptr)
    {
        delete c;
        return status_ != SM_SUCCESS ? status(status_) : static_cast<int32_t>(SM_ERR);
    }
    out = c;
    return SMLIB_OK;
}

// A profile line parallel to the axis sweeps the cylinder its points turn on (no point of it
// leaves that cylinder), between the helices its ends trace and the line where the motion starts
// and ends. `radius` is that cylinder's when `line` runs within a thousandth of `tolerance` of
// parallel to the axis, at least 1e-6 off it; 0 otherwise.
double axial_radius(const Screw &m, const SmCurve *line, double tolerance)
{
    if (SM_CAST_PTR(SmLine, line) == nullptr)
    {
        return 0.0;
    }
    SmPoint3d a, b;
    line->GetEnds(a, b);
    // How far a point of the line strays from the axis lies between the ends' and the nearest
    // point's (the distance is convex along the line).
    SmVector3d w = a - m.o, d = b - a;
    SmVector3d wp = w - m.k * w.Dot(m.k), dp = d - m.k * d.Dot(m.k);
    double t = dp.Dot(dp) > 0.0 ? std::clamp(-wp.Dot(dp) / dp.Dot(dp), 0.0, 1.0) : 0.0;
    double ra = m.reach(a), rb = m.reach(b), nearest = (wp + dp * t).Length();
    double far = std::max(ra, rb), near = std::min(nearest, std::min(ra, rb));
    if (!(near >= 1e-6 && far - near <= 1e-3 * tolerance))
    {
        return 0.0;
    }
    return 0.5 * (ra + rb);
}

// The piece of the cylinder of `radius` the profile line from `a` to `b` sweeps over the stations
// of `st` (at most 202.5 degrees of turn): an exact SmCylinder about the axis, bounded by the
// line placed at both ends of the piece and its ends' paths. The cylinder's own turn starts a
// quarter turn before the piece and its height runs past it, so the piece lies inside its domain.
int32_t cylinder_piece(SmBrep &brep, const SmContext &ctx, const Screw &m, const SmPoint3d &a, const SmPoint3d &b, double radius,
                       const Stations &st)
{
    const double a0 = st.t.front(), a1 = st.t.back();
    SmVector3d w = a - m.o;
    SmVector3d e = w - m.k * w.Dot(m.k);
    e = e * (1.0 / e.Length());
    SmVector3d x = m.turn(e, a0 - 0.25 * TWO_PI);
    SmVector3d y = m.k * x;
    double z0 = std::min((a - m.o).Dot(m.k), (b - m.o).Dot(m.k)) + std::min(m.lift(a0), m.lift(a1));
    double z1 = std::max((a - m.o).Dot(m.k), (b - m.o).Dot(m.k)) + std::max(m.lift(a0), m.lift(a1));
    double margin = 0.25 * (z1 - z0) + 1e-3;
    SmPoint3d base = m.o + m.k * (z0 - margin);
    SmCylinder *cylinder = new (ctx) SmCylinder(base, x, y, radius, 0.0, 360.0, z1 - z0 + 2.0 * margin, FALSE, FALSE, FALSE, &ctx);
    if (cylinder == nullptr)
    {
        return status(SM_ERR_OUT_OF_MEMORY);
    }
    Curves loop;
    auto line = [&](double at) -> int32_t
    {
        SmLine *l = nullptr;
        SmStatus st = SmLine::CreateLineSegment(ctx, 3, m.at(a, at), m.at(b, at), l);
        if (st != SM_SUCCESS || l == nullptr)
        {
            delete l;
            return st != SM_SUCCESS ? status(st) : static_cast<int32_t>(SM_ERR);
        }
        loop.v.push_back(l);
        return SMLIB_OK;
    };
    SmCurve *pa = nullptr, *pb = nullptr;
    int32_t s = line(a0);
    if (s == SMLIB_OK && (s = path_curve(ctx, m, b, st, pb)) == SMLIB_OK)
    {
        loop.v.push_back(pb);
        s = line(a1);
    }
    if (s == SMLIB_OK && (s = path_curve(ctx, m, a, st, pa)) == SMLIB_OK)
    {
        loop.v.push_back(pa);
    }
    if (s != SMLIB_OK)
    {
        delete cylinder;
        return s;
    }
    SmTArray<ULONG> counts;
    counts.Add(4);
    SmTArray<SmCurve *> curves = loop.array();
    SmTArray<SmOrientType> orient;
    orient.Add(SM_OT_SAME);
    orient.Add(SM_OT_SAME);
    orient.Add(SM_OT_OPPOSITE);
    orient.Add(SM_OT_OPPOSITE);
    SmShell *shell = nullptr;
    SmTArray<SmFace *> faces;
    s = status(brep.CreateFacesWithCurves(brep.GetInfiniteRegion(), counts, &curves, orient, cylinder, shell, faces));
    if (s == SMLIB_OK)
    {
        loop.v.clear(); // the face took the curves and the cylinder over
        if (faces.GetSize() != 1)
        {
            return static_cast<int32_t>(SM_ERR);
        }
    }
    return s;
}

// The wall one profile curve sweeps over `st`: across the profile curve its NURBS, across the
// turns each control point's path, the weights the curve's times the arcs'. A wall's point is
// the curve's weighted mean of its control points' paths at the same parameter, which all turn
// and lift alike, so the wall is the profile curve turned and lifted as one.
int32_t screw_wall(const SmContext &ctx, const Screw &m, const Nurbs &c, const Stations &st, SmBSplineSurface *&out)
{
    ULONG count = c.points.GetSize();
    ULONG across = static_cast<ULONG>(2 * st.spans() + 1);
    SmTArray<SmPoint3d> net;
    SmTArray<double> weights;
    net.SetSize(count * across);
    weights.SetSize(count * across);
    for (ULONG l = 0; l < count; ++l)
    {
        std::vector<SmPoint3d> D = path(m, c.points[l], st);
        for (ULONG j = 0; j < across; ++j)
        {
            net[l * across + j] = D[j];
            weights[l * across + j] = (c.weights.GetSize() > 0 ? c.weights[l] : 1.0) * st.weight(j);
        }
    }
    SmTArray<ULONG> vmult;
    SmTArray<double> vknots;
    st.knots(vknots, vmult);
    SmBSplineSurface *s = nullptr;
    SmStatus status_ = SmBSplineSurface::CreateCanonical(ctx, c.degree, 2, net, SM_SF_UNSPECIFIED, c.mult, vmult, c.unique, vknots,
                                                         SM_KT_UNSPECIFIED, &weights, nullptr, s);
    if (status_ != SM_SUCCESS || s == nullptr)
    {
        delete s;
        return status_ != SM_SUCCESS ? status(status_) : static_cast<int32_t>(SM_ERR);
    }
    out = s;
    return SMLIB_OK;
}

} // namespace

extern "C" int32_t smlib_screw_sweep(const int64_t *kinds, int64_t n_kinds, const double *reals, int64_t n_reals, const double origin[3],
                                     const double axis[3], double pitch, double turns, double tolerance, smlib_brep **out)
{
    if (out == nullptr || origin == nullptr || axis == nullptr || !coord_ok(origin, 3) || !std::isfinite(pitch) ||
        !(std::fabs(pitch) > 1e-6 && std::fabs(pitch) < MAX_COORD) || !std::isfinite(turns) || !(turns >= 1e-3 && turns <= SMLIB_SCREW_MAX_TURNS) ||
        !std::isfinite(tolerance) || !(tolerance >= 1e-6 && tolerance < 1e3))
    {
        return SMLIB_ERR_ARG;
    }
    Call c("smlib_screw_sweep", list(n_kinds, n_reals) + ", origin=" + vec(origin, 3) + ", axis=" + vec(axis, 3) + ", pitch=" + num(pitch) +
                                    ", turns=" + num(turns) + ", tolerance=" + num(tolerance));
    dump_curves(c, "profile", kinds, n_kinds, reals, n_reals);
    return c.run(
        [&]
        {
            Screw m;
            if (!unit(axis, m.k))
            {
                return SMLIB_ERR_ARG;
            }
            m.o = SmPoint3d(origin[0], origin[1], origin[2]);
            m.pitch = pitch;
            Curves profile;
            int32_t s = parse(kinds, n_kinds, reals, n_reals, profile);
            if (s != SMLIB_OK)
            {
                return s;
            }
            // SM_API's curve sweep carries the profile on its own frame along a path, which along a
            // helix turns the profile away from the screw motion. So each profile curve's wall is
            // built from the motion itself (Stations, screw_wall): turned exactly, lifted within
            // |pitch| / (2 pi) * lift_error(h) of the screw surface, the spans h no wider than
            // that keeps within `tolerance`, nor than an eighth of a turn.
            std::vector<Nurbs> walls(profile.v.size());
            double reach = 0.0;
            for (size_t i = 0; i < profile.v.size(); ++i)
            {
                s = nurbs_of(profile.v[i], walls[i]);
                if (s != SMLIB_OK)
                {
                    return s;
                }
                for (ULONG l = 0; l < walls[i].points.GetSize(); ++l)
                {
                    reach = std::max(reach, m.reach(walls[i].points[l]));
                }
            }
            if (!(reach > 1e-6))
            {
                return SMLIB_ERR_ARG;
            }
            // Walls of lines parallel to the axis are exact cylinders, which cannot turn on
            // themselves: then every wall comes in pieces of half a turn (the last one up to an
            // eighth longer), each its own face, so the walls on both sides of an edge share its
            // vertices.
            std::vector<double> radius(profile.v.size());
            bool cylinders = false;
            for (size_t i = 0; i < profile.v.size(); ++i)
            {
                radius[i] = axial_radius(m, profile.v[i], tolerance);
                cylinders = cylinders || radius[i] > 0.0;
            }
            const double total = TWO_PI * turns;
            std::vector<double> breaks = {0.0};
            for (double half = 0.5 * TWO_PI; cylinders && half < total - TWO_PI / 16.0; half += 0.5 * TWO_PI)
            {
                breaks.push_back(half);
            }
            breaks.push_back(total);
            const double rise = std::fabs(pitch) / TWO_PI;
            double h = std::min(TWO_PI / 8.0, std::cbrt(tolerance / (0.008 * rise)));
            while (rise * lift_error(0.5 * h) > tolerance)
            {
                h *= 0.99;
            }
            std::vector<Stations> pieces;
            double spans = 0.0;
            for (size_t j = 0; j + 1 < breaks.size(); ++j)
            {
                double n = std::ceil((breaks[j + 1] - breaks[j]) / h);
                spans += n;
                if (!(spans <= SMLIB_SCREW_MAX_STATIONS))
                {
                    return SMLIB_ERR_ARG;
                }
                pieces.emplace_back(breaks[j], breaks[j + 1], static_cast<size_t>(n));
            }
            const SmContext &ctx = *SmApiGetOrCreateContext();
            std::unique_ptr<SmBrep> brep(new (ctx) SmBrep());
            for (const Stations &st : pieces)
            {
                for (size_t i = 0; i < walls.size(); ++i)
                {
                    if (radius[i] > 0.0)
                    {
                        SmPoint3d a, b;
                        profile.v[i]->GetEnds(a, b);
                        s = cylinder_piece(*brep, ctx, m, a, b, radius[i], st);
                        if (s != SMLIB_OK)
                        {
                            return s;
                        }
                        continue;
                    }
                    SmBSplineSurface *surface = nullptr;
                    s = screw_wall(ctx, m, walls[i], st, surface);
                    if (s != SMLIB_OK)
                    {
                        return s;
                    }
                    SmFace *face = nullptr;
                    s = status(brep->CreateFaceFromSurface(surface, surface->GetNaturalUVDomain(), face));
                    if (s != SMLIB_OK)
                    {
                        return s;
                    }
                }
            }
            // The caps: the profile where the motion starts and where it ends, placed exactly.
            for (double a : {0.0, total})
            {
                Curves cap;
                for (SmCurve *src : profile.v)
                {
                    SmCurve *copy = nullptr;
                    if (src->Copy(ctx, copy) != SM_SUCCESS || copy == nullptr)
                    {
                        delete copy;
                        return static_cast<int32_t>(SM_ERR);
                    }
                    cap.v.push_back(copy);
                    if (a != 0.0 && copy->Transform(m.placement(a)) != SM_SUCCESS)
                    {
                        return static_cast<int32_t>(SM_ERR);
                    }
                }
                SmTArray<SmFace *> faces;
                s = status(brep->CreatePlanarFacesWith3DCurves(brep->GetInfiniteRegion(), cap.array(), brep->GetTolerance(), faces));
                if (s != SMLIB_OK)
                {
                    return s;
                }
                cap.v.clear(); // the faces took the curves over
            }
            s = status(brep->StitchAndOrient());
            if (s != SMLIB_OK)
            {
                return s;
            }
            return finish(brep.release(), true, out);
        }
    );
}

namespace
{

// The station maps of smlib_sweep_stations: `n` maps of 12 numbers, each a 3 x 3 matrix by columns
// then a translation, every number finite, each matrix of bounded size and far from singular (a
// section squashed to a line or a point is no section), each translation within MAX_COORD.
bool maps_ok(const double *maps, int64_t n)
{
    for (int64_t i = 0; i < n; ++i)
    {
        const double *m = maps + 12 * i;
        if (!finite(m, 12) || !coord_ok(m + 9, 3))
        {
            return false;
        }
        double c0 = length(m), c1 = length(m + 3), c2 = length(m + 6);
        if (!(c0 < 1e6 && c1 < 1e6 && c2 < 1e6))
        {
            return false;
        }
        double det = m[0] * (m[4] * m[8] - m[7] * m[5]) - m[3] * (m[1] * m[8] - m[7] * m[2]) + m[6] * (m[1] * m[5] - m[4] * m[2]);
        if (!(std::fabs(det) > 1e-9 * c0 * c1 * c2))
        {
            return false;
        }
    }
    return true;
}

// A copy of `c` as a plain B-spline with every control point carried by the affine map `m`
// (x -> M x + t, as maps_ok reads it). An affine map moves a NURBS curve exactly by its control
// points, the weights unchanged.
int32_t placed(const SmContext &ctx, const SmBSplineCurve &c, const double *m, SmBSplineCurve *&out)
{
    ULONG dim = 0, degree = 0;
    SmTArray<SmPoint3d> points;
    SmBSplineCurveForm form = SM_CF_UNSPECIFIED;
    SmTArray<ULONG> mult;
    SmTArray<double> unique;
    SmKnotType kt = SM_KT_UNSPECIFIED;
    SmTArray<double> weights;
    if (c.GetCanonical(dim, degree, points, form, mult, unique, kt, weights) != SM_SUCCESS || dim != 3)
    {
        return static_cast<int32_t>(SM_ERR);
    }
    for (ULONG i = 0; i < points.GetSize(); ++i)
    {
        const SmPoint3d p = points[i];
        points[i] = SmPoint3d(m[0] * p.x + m[3] * p.y + m[6] * p.z + m[9], m[1] * p.x + m[4] * p.y + m[7] * p.z + m[10],
                              m[2] * p.x + m[5] * p.y + m[8] * p.z + m[11]);
    }
    SmBSplineCurve *r = nullptr;
    SmStatus s = SmBSplineCurve::CreateCanonical(ctx, 3, degree, points, SM_CF_UNSPECIFIED, mult, unique, SM_KT_UNSPECIFIED,
                                                 weights.GetSize() > 0 ? &weights : nullptr, nullptr, r);
    if (s != SM_SUCCESS || r == nullptr)
    {
        delete r;
        return s != SM_SUCCESS ? status(s) : SM_ERR;
    }
    out = r;
    return SMLIB_OK;
}

} // namespace

extern "C" int32_t smlib_sweep_stations(const int64_t *kinds, int64_t n_kinds, const double *reals, int64_t n_reals, const double *maps,
                                        int64_t n_stations, int32_t ends, smlib_brep **out)
{
    if (out == nullptr || maps == nullptr || n_stations < 2 || n_stations > SMLIB_SWEEP_MAX_STATIONS ||
        (ends != SMLIB_ENDS_OPEN && ends != SMLIB_ENDS_CAPPED && ends != SMLIB_ENDS_JOINED))
    {
        return SMLIB_ERR_ARG;
    }
    Call c("smlib_sweep_stations", list(n_kinds, n_reals) + ", stations=" + std::to_string(n_stations) + ", ends=" + std::to_string(ends));
    dump_curves(c, "profile", kinds, n_kinds, reals, n_reals);
    c.input_text("maps",
                 [=](const std::string &path)
                 {
                     if (FILE *f = std::fopen(path.c_str(), "w"))
                     {
                         for (int64_t i = 0; i < 12 * n_stations; ++i)
                         {
                             std::fprintf(f, "%.17g%c", maps[i], i % 12 == 11 ? '\n' : ' ');
                         }
                         std::fclose(f);
                     }
                 });
    return c.run(
        [&]
        {
            if (!maps_ok(maps, n_stations))
            {
                return SMLIB_ERR_ARG;
            }
            Curves profile;
            int32_t s = parse(kinds, n_kinds, reals, n_reals, profile);
            if (s != SMLIB_OK)
            {
                return s;
            }
            if (static_cast<double>(profile.v.size()) * static_cast<double>(n_stations) > SMLIB_SWEEP_MAX_STATION_CURVES)
            {
                return SMLIB_ERR_ARG;
            }
            // Every section is the profile carried by its station's map: the sections share the
            // profile's degrees and knots, which the skin interpolates control point by control
            // point (cubic across the stations, or as high as their number allows).
            const SmContext &ctx = *SmApiGetOrCreateContext();
            Curves sections;
            SmTArray<ULONG> counts;
            SmTArray<const SmCurve *> all;
            for (int64_t i = 0; i < n_stations; ++i)
            {
                for (SmCurve *x : profile.v)
                {
                    const SmBSplineCurve *b = SM_CAST_PTR(SmBSplineCurve, x);
                    if (b == nullptr)
                    {
                        return SMLIB_ERR_ARG;
                    }
                    SmBSplineCurve *p = nullptr;
                    s = placed(ctx, *b, maps + 12 * i, p);
                    if (s != SMLIB_OK)
                    {
                        return s;
                    }
                    sections.v.push_back(p);
                    all.Add(p);
                }
                counts.Add(static_cast<ULONG>(profile.v.size()));
            }
            ULONG degree = static_cast<ULONG>(std::min<int64_t>(3, n_stations - 1));
            SmBrep *r = nullptr;
            s = status(SmApiCreateSkinPrimitive(counts, all, degree, ends == SMLIB_ENDS_CAPPED ? TRUE : FALSE, r, nullptr, nullptr, nullptr));
            if (s == SMLIB_OK && ends == SMLIB_ENDS_JOINED)
            {
                // The last section lies on the first: the walls' ends are glued to their starts.
                SmBoolean solid = FALSE;
                ULONG stitched = 0;
                double vertex_gap = 0.0, edge_gap = 0.0;
                s = status(SmApiStitchIntoSolid(r, solid, stitched, vertex_gap, edge_gap));
                if (s == SMLIB_OK && !solid)
                {
                    s = SMLIB_ERR_NOT_SOLID;
                }
            }
            if (s != SMLIB_OK)
            {
                delete r;
                return s;
            }
            return finish(r, ends != SMLIB_ENDS_OPEN, out);
        }
    );
}

extern "C" int32_t smlib_planar_faces(const int64_t *kinds, int64_t n_kinds, const double *reals, int64_t n_reals, smlib_brep **out)
{
    if (out == nullptr)
    {
        return SMLIB_ERR_ARG;
    }
    Call c("smlib_planar_faces", list(n_kinds, n_reals));
    dump_curves(c, "loops", kinds, n_kinds, reals, n_reals);
    return c.run(
        [&]
        {
            Curves loops;
            int32_t s = parse(kinds, n_kinds, reals, n_reals, loops);
            if (s != SMLIB_OK)
            {
                return s;
            }
            SmBrep *r = nullptr;
            // On success the faces' edges own the curves; on failure they stay ours.
            s = status(SmApiCreatePlanarFaces(loops.array(), r));
            if (s != SMLIB_OK)
            {
                delete r;
                return s;
            }
            loops.v.clear();
            return finish(r, false, out);
        }
    );
}

// ---- Curves projected onto a face -------------------------------------------------------------

// The curves a projection made, as a curve list (NURBS only) with each curve's interval.
struct smlib_curves
{
    std::vector<int64_t> kinds;
    std::vector<double> reals;
    std::vector<double> intervals;
};

namespace
{

// `c` appended to `out` as the NURBS it is (SMLib holds every curve here as its B-spline):
// control points (x, y, z, weight), all knots, and its interval. SMLIB_ERR_UNSUPPORTED for a
// curve that is no B-spline or whose knots are not clamped (the curve list takes clamped ones).
int32_t put_curve(const SmCurve *c, smlib_curves &out)
{
    SmBSplineCurve *b = c == nullptr ? nullptr : c->GetRootCurve();
    gw_CURVE *n = b == nullptr ? nullptr : b->GetOrCreateGwNurbPointer();
    if (n == nullptr || n->knt == nullptr || n->pol == nullptr)
    {
        return SMLIB_ERR_UNSUPPORTED;
    }
    const int64_t degree = n->p;
    const int64_t count = n->pol->n + 1;
    const int64_t knots = n->knt->m + 1;
    if (degree < 1 || degree > SMLIB_CURVE_MAX_DEGREE || count < degree + 1 || count > SMLIB_CURVE_MAX_POINTS || knots != count + degree + 1)
    {
        return SMLIB_ERR_UNSUPPORTED;
    }
    const double *u = n->knt->U;
    for (int64_t i = 1; i <= degree; ++i)
    {
        if (u[i] != u[0] || u[knots - 1 - i] != u[knots - 1])
        {
            return SMLIB_ERR_UNSUPPORTED;
        }
    }
    out.kinds.insert(out.kinds.end(), {SMLIB_CURVE_NURBS, degree, count});
    for (int64_t i = 0; i < count; ++i)
    {
        const NL_CPOINT &q = n->pol->Pw[i];
        const double w = q.w != NL_NOW ? q.w : 1.0;
        out.reals.insert(out.reals.end(), {q.x / w, q.y / w, q.z / w, w});
    }
    out.reals.insert(out.reals.end(), u, u + knots);
    SmExtent1d ivl = c->GetNaturalInterval();
    out.intervals.insert(out.intervals.end(), {ivl.GetMin(), ivl.GetMax()});
    return SMLIB_OK;
}

} // namespace

extern "C" int32_t smlib_project_curves(const smlib_brep *brep, const double face_point[3], double pick_tolerance, const int64_t *kinds,
                                        int64_t n_kinds, const double *reals, int64_t n_reals, const double direction[3], double tolerance,
                                        smlib_curves **out)
{
    if (brep == nullptr || brep->p == nullptr || face_point == nullptr || !finite(face_point, 3) || direction == nullptr ||
        !finite(direction, 3) || out == nullptr || !std::isfinite(pick_tolerance) || !(pick_tolerance > 0.0 && pick_tolerance < 1e6) ||
        !std::isfinite(tolerance) || !(tolerance >= 1e-9 && tolerance < 1e3))
    {
        return SMLIB_ERR_ARG;
    }
    const double len = length(direction);
    if (!(len > 1e-12))
    {
        return SMLIB_ERR_ARG;
    }
    Call c("smlib_project_curves", handle(brep) + ", face_point=" + vec(face_point, 3) + ", pick_tolerance=" + num(pick_tolerance) + ", " +
                                       list(n_kinds, n_reals) + ", direction=" + vec(direction, 3) + ", tolerance=" + num(tolerance));
    c.input("brep", brep);
    dump_curves(c, "curves", kinds, n_kinds, reals, n_reals);
    return c.run(
        [&]
        {
            Curves curves;
            int32_t s = parse(kinds, n_kinds, reals, n_reals, curves);
            if (s != SMLIB_OK)
            {
                return s;
            }
            std::vector<SmFace *> picked;
            s = smlib_pick::pick_faces(brep->p, face_point, 1, pick_tolerance, picked, nullptr);
            if (s != SMLIB_OK)
            {
                return s;
            }
            SmTArray<SmFace *> faces;
            faces.Add(picked.front());
            const SmContext &ctx = *SmApiGetOrCreateContext();
            const SmVector3d dir(direction[0] / len, direction[1] / len, direction[2] / len);
            const SmApproxTol3d tol = SM_CAST_APPROXTOL3D(tolerance);
            auto result = std::make_unique<smlib_curves>();
            for (SmCurve *curve : curves.v)
            {
                SmBSplineCurve *b = curve->GetRootCurve();
                if (b == nullptr)
                {
                    return SMLIB_ERR_ARG;
                }
                SmTArray<SmCurve *> made;
                SmObjsDelete<SmCurve *> clean(&made);
                s = status(brep->p->CreateParallelProjectionCurves(ctx, *b, dir, nullptr, &tol, nullptr, &made, nullptr, nullptr, &faces));
                if (s != SMLIB_OK)
                {
                    return s;
                }
                for (ULONG i = 0; i < made.GetSize(); ++i)
                {
                    s = put_curve(made[i], *result);
                    if (s != SMLIB_OK)
                    {
                        return s;
                    }
                }
            }
            if (result->intervals.empty())
            {
                return SMLIB_ERR_EMPTY;
            }
            *out = result.release();
            return SMLIB_OK;
        }
    );
}

extern "C" int32_t smlib_face_section(const smlib_brep *brep, const double face_point[3], double pick_tolerance, const double plane_point[3],
                                      const double plane_normal[3], double tolerance, smlib_curves **out)
{
    if (brep == nullptr || brep->p == nullptr || face_point == nullptr || !finite(face_point, 3) || plane_point == nullptr ||
        !finite(plane_point, 3) || plane_normal == nullptr || !finite(plane_normal, 3) || out == nullptr || !std::isfinite(pick_tolerance) ||
        !(pick_tolerance > 0.0 && pick_tolerance < 1e6) || !std::isfinite(tolerance) || !(tolerance >= 1e-9 && tolerance < 1e3))
    {
        return SMLIB_ERR_ARG;
    }
    SmVector3d n(plane_normal[0], plane_normal[1], plane_normal[2]);
    if (!(n.Length() > 1e-12) || n.Unitize() != SM_SUCCESS)
    {
        return SMLIB_ERR_ARG;
    }
    Call c("smlib_face_section", handle(brep) + ", face_point=" + vec(face_point, 3) + ", pick_tolerance=" + num(pick_tolerance) +
                                     ", plane_point=" + vec(plane_point, 3) + ", plane_normal=" + vec(plane_normal, 3) + ", tolerance=" + num(tolerance));
    c.input("brep", brep);
    return c.run(
        [&]
        {
            std::vector<SmFace *> picked;
            int32_t s = smlib_pick::pick_faces(brep->p, face_point, 1, pick_tolerance, picked, nullptr);
            if (s != SMLIB_OK)
            {
                return s;
            }
            const SmContext &ctx = *SmApiGetOrCreateContext();
            const SmApproxTol3d tol = SM_CAST_APPROXTOL3D(tolerance);
            SmTArray<SmCurve *> made, uv;
            SmTArray<SmFace *> faces;
            SmObjsDelete<SmCurve *> clean(&made);
            SmObjsDelete<SmCurve *> clean_uv(&uv);
            s = status(brep->p->CreatePlanarSectionCurves(ctx, SmPoint3d(plane_point[0], plane_point[1], plane_point[2]), n, &tol, nullptr, &made,
                                                          &uv, &faces));
            if (s != SMLIB_OK)
            {
                return s;
            }
            if (faces.GetSize() != made.GetSize())
            {
                return SMLIB_ERR_UNSUPPORTED;
            }
            auto result = std::make_unique<smlib_curves>();
            for (ULONG i = 0; i < made.GetSize(); ++i)
            {
                if (faces[i] == picked.front())
                {
                    s = put_curve(made[i], *result);
                    if (s != SMLIB_OK)
                    {
                        return s;
                    }
                }
            }
            if (result->intervals.empty())
            {
                return SMLIB_ERR_EMPTY;
            }
            *out = result.release();
            return SMLIB_OK;
        }
    );
}

extern "C" int32_t smlib_curves_size(const smlib_curves *curves, int64_t *n_kinds, int64_t *n_reals, int64_t *n_curves)
{
    if (curves == nullptr || n_kinds == nullptr || n_reals == nullptr || n_curves == nullptr)
    {
        return SMLIB_ERR_ARG;
    }
    *n_kinds = static_cast<int64_t>(curves->kinds.size());
    *n_reals = static_cast<int64_t>(curves->reals.size());
    *n_curves = static_cast<int64_t>(curves->intervals.size() / 2);
    return SMLIB_OK;
}

extern "C" int32_t smlib_curves_read(const smlib_curves *curves, int64_t *kinds, int64_t n_kinds, double *reals, int64_t n_reals,
                                     double *intervals, int64_t n_intervals)
{
    if (curves == nullptr || kinds == nullptr || reals == nullptr || intervals == nullptr ||
        n_kinds != static_cast<int64_t>(curves->kinds.size()) || n_reals != static_cast<int64_t>(curves->reals.size()) ||
        n_intervals != static_cast<int64_t>(curves->intervals.size()))
    {
        return SMLIB_ERR_ARG;
    }
    std::copy(curves->kinds.begin(), curves->kinds.end(), kinds);
    std::copy(curves->reals.begin(), curves->reals.end(), reals);
    std::copy(curves->intervals.begin(), curves->intervals.end(), intervals);
    return SMLIB_OK;
}

extern "C" void smlib_curves_free(smlib_curves *curves)
{
    delete curves;
}
