// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// The sweeps of include/smlib.h: extrusion (straight or with draft), revolution, sweep along a
// path, loft through sections and the screw sweep of a thread. Profiles, paths and sections
// cross as curve lists (lines, circular arcs and NURBS, include/smlib.h); each call builds its
// own SMLib curves from them, so nothing the caller owns is touched.
//
// Every entry point checks its numbers before SMLib sees them, runs SM_API's sweep, and hands
// the result over only when it is what was asked for: a closed manifold solid when the ends are
// capped (SMLIB_ERR_NOT_SOLID otherwise), and never a body without faces (SMLIB_ERR_EMPTY).

#include "smlib_shim.h"

#include <SmApiBrep.h>
#include <SmApiCurves.h>
#include <SmApiPrimitives.h>
#include <SmBSplineCurve.h>
#include <SmCircle.h>
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
    settle_tolerances(r);
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
                                 double distance, double draft, int32_t cap, smlib_brep **out)
{
    if (out == nullptr || direction == nullptr || !std::isfinite(distance) || !(distance > 1e-9 && distance < MAX_COORD) ||
        !std::isfinite(draft) || std::fabs(draft) >= SM_PI / 2.0 - 1e-6 || (cap != 0 && cap != 1))
    {
        return SMLIB_ERR_ARG;
    }
    Call c("smlib_extrude", list(n_kinds, n_reals) + ", direction=" + vec(direction, 3) + ", distance=" + num(distance) +
                                ", draft=" + num(draft) + ", cap=" + std::to_string(cap));
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
                s = status(SmApiCreateTaperExtrude(a, height, angle, cap ? 3 : 0, FALSE, r));
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

extern "C" int32_t smlib_sweep(const int64_t *kinds, int64_t n_kinds, const double *reals, int64_t n_reals, const int64_t *path_kinds,
                               int64_t n_path_kinds, const double *path_reals, int64_t n_path_reals, int32_t frame, int32_t cap,
                               smlib_brep **out)
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
            std::unique_ptr<SmBSplineCurve> along;
            s = joined(path, along);
            if (s != SMLIB_OK)
            {
                return s;
            }
            SmSweepOptions options;
            options.bCapEndsArg = cap ? TRUE : FALSE;
            options.bTranslationalSweep = frame == SMLIB_SWEEP_TRANSLATE ? TRUE : FALSE;
            options.bProfileAtPathStart = TRUE;
            SmTArray<SmCurve *> a = profile.array();
            SmBrep *r = nullptr;
            s = status(SmApiCreateCurveSweep(a, along.get(), nullptr, nullptr, &options, r, nullptr, nullptr, nullptr));
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
