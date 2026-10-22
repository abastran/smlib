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
#include <NL_Geometry.h>
#include <NL_Nurbsdef.h>
#include <SmCircle.h>
#include <SmFace.h>
#include <SmLine.h>
#include <SmPrimitiveCreation.h>

#include <memory>

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
            SmVector3d k;
            if (!unit(axis, k))
            {
                return SMLIB_ERR_ARG;
            }
            Curves profile;
            int32_t s = parse(kinds, n_kinds, reals, n_reals, profile);
            if (s != SMLIB_OK)
            {
                return s;
            }
            // SM_API's curve sweep carries the profile on its own frame along a path, which along a
            // helix turns the profile away from the screw motion. So the profile is placed by the
            // screw motion itself at stations along the turns, and the walls are skinned through
            // them (cubic across the stations, which they interpolate).
            SmPoint3d o(origin[0], origin[1], origin[2]);
            double reach = 0.0;
            for (SmCurve *x : profile.v)
            {
                SmExtent1d ivl = x->GetNaturalInterval();
                for (int j = 0; j <= 16; ++j)
                {
                    SmPoint3d p;
                    if (x->EvaluatePoint(ivl.Evaluate(j / 16.0), p) != SM_SUCCESS)
                    {
                        return static_cast<int32_t>(SM_ERR);
                    }
                    SmVector3d w = p - o;
                    reach = std::max(reach, (w - k * w.Dot(k)).Length());
                }
            }
            if (!(reach > 1e-6))
            {
                return SMLIB_ERR_ARG;
            }
            // A cubic through stations h radians apart strays about reach * h^4 / 384 from the
            // circle each point turns on.
            double h = std::pow(384.0 * tolerance / reach, 0.25);
            double per_turn = std::max(12.0, std::ceil(TWO_PI / std::min(h, TWO_PI)));
            double n = std::ceil(turns * per_turn);
            if (!(n >= 2.0 && n <= SMLIB_SCREW_MAX_STATIONS))
            {
                return SMLIB_ERR_ARG;
            }
            ULONG stations = static_cast<ULONG>(n) + 1;
            const SmContext &ctx = *SmApiGetOrCreateContext();
            Curves placed;
            SmTArray<ULONG> counts;
            SmTArray<const SmCurve *> all;
            // Rodrigues: v turned right-handed about k by angle a.
            auto turn = [&](const SmVector3d &v, double a) { return v * std::cos(a) + (k * v) * std::sin(a) + k * (k.Dot(v) * (1.0 - std::cos(a))); };
            for (ULONG i = 0; i < stations; ++i)
            {
                double f = static_cast<double>(i) / static_cast<double>(stations - 1);
                double a = TWO_PI * turns * f;
                SmVector3d x = turn(SmVector3d(1, 0, 0), a);
                SmVector3d y = turn(SmVector3d(0, 1, 0), a);
                // p -> o + R (p - o) + k * lift, as a placement: origin o - R o + k * lift.
                SmPoint3d at = o - turn(o, a) + k * (pitch * turns * f);
                SmAxis2Placement motion(at, x, y);
                for (SmCurve *src : profile.v)
                {
                    SmCurve *copy = nullptr;
                    if (src->Copy(ctx, copy) != SM_SUCCESS || copy == nullptr)
                    {
                        delete copy;
                        return static_cast<int32_t>(SM_ERR);
                    }
                    placed.v.push_back(copy);
                    if (copy->Transform(motion) != SM_SUCCESS)
                    {
                        return static_cast<int32_t>(SM_ERR);
                    }
                    all.Add(copy);
                }
                counts.Add(static_cast<ULONG>(profile.v.size()));
            }
            SmBrep *r = nullptr;
            s = status(SmApiCreateSkinPrimitive(counts, all, 3, TRUE, r, nullptr, nullptr, nullptr));
            if (s != SMLIB_OK)
            {
                delete r;
                return s;
            }
            return finish(r, true, out);
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
