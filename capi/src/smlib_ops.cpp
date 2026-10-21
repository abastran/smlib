// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// The operations of include/smlib.h: booleans, fillet, chamfer and shell. Each one copies its
// inputs and works on the copies only (SM_API edits operand A in place and deletes operand B),
// checks that a reported success is a usable body, and frees every copy on failure.
//
// Edges and faces are chosen by points: each point picks the edge (face) nearest to it, which
// must lie within the caller's tolerance and be nearer than any other by more than the body's
// modelling tolerance. The API thereby never depends on SMLib's internal topology order.

#include "smlib_shim.h"

#include <SmApiBrep.h>
#include <SmApiQueries.h>
#include <SmEdge.h>
#include <SmEdgeuse.h>
#include <SmFace.h>
#include <SmFilletExecutive.h>
#include <SmFilletStandardSolver.h>
#include <SmVertex.h>

#include <memory>
#include <unordered_set>

using namespace smlib_shim;

namespace
{

// Upper bound on the points of one call: far more than any selection, small enough that the
// nearest-entity search (points x entities) stays cheap.
constexpr int64_t MAX_POINTS = 4096;

// Two edges continue each other when their directions at the shared vertex differ by less
// than SMLib's own angle for parallel vectors (SM_CONTINUITY_ANGLE, 0.5 degrees).
const double TANGENT_COS = std::cos(SM_CONTINUITY_ANGLE * SM_PI / 180.0);

// Owns an SmBrep until it is handed on.
using Owned = std::unique_ptr<SmBrep>;

// A deep copy, or null.
Owned copy_of(const smlib_brep *b, int32_t &s)
{
    SmBrep *r = nullptr;
    s = status(SmApiBrepCopy(b->p, r));
    if (s != SMLIB_OK)
    {
        delete r;
        return nullptr;
    }
    if (r == nullptr)
    {
        s = SMLIB_ERR_EMPTY;
    }
    return Owned(r);
}

bool points_ok(const double *points, int64_t n, bool allow_none)
{
    if (n < 0 || n > MAX_POINTS || (n == 0 && !allow_none))
    {
        return false;
    }
    return n == 0 || (points != nullptr && finite(points, static_cast<int>(n * 3)));
}

bool tolerance_ok(double t)
{
    return std::isfinite(t) && t > 0.0 && t < 1e6;
}

// `points` (3 per point), `n` of them, for the log line.
std::string pts(const double *points, int64_t n)
{
    return n == 0 ? std::string("[]") : vec(points, static_cast<int>(n * 3));
}

// The entity nearest each point under the picking rule. On failure, the index of the point
// that matched none or several goes to `failed`.
template <typename T, typename Dist>
int32_t pick(const SmTArray<T *> &all, const double *points, int64_t n, double tol, double model_tol, Dist dist,
             std::vector<T *> &picked, int64_t *failed)
{
    for (int64_t i = 0; i < n; ++i)
    {
        SmPoint3d p(points[3 * i], points[3 * i + 1], points[3 * i + 2]);
        T *best = nullptr;
        double d0 = HUGE_VAL;
        double d1 = HUGE_VAL;
        for (ULONG k = 0; k < all.GetSize(); ++k)
        {
            double d = 0.0;
            if (!dist(all[k], p, d) || !std::isfinite(d))
            {
                continue;
            }
            if (d < d0)
            {
                d1 = d0;
                d0 = d;
                best = all[k];
            }
            else if (d < d1)
            {
                d1 = d;
            }
        }
        int32_t s = SMLIB_OK;
        if (best == nullptr || d0 > tol)
        {
            s = SMLIB_ERR_PICK_NONE;
        }
        else if (d1 <= tol && d1 - d0 <= model_tol)
        {
            s = SMLIB_ERR_PICK_AMBIGUOUS;
        }
        if (s != SMLIB_OK)
        {
            if (failed != nullptr)
            {
                *failed = i;
            }
            return s;
        }
        if (std::find(picked.begin(), picked.end(), best) == picked.end())
        {
            picked.push_back(best);
        }
    }
    return SMLIB_OK;
}

bool edge_distance(SmEdge *e, const SmPoint3d &p, double &d)
{
    SmPoint3d q;
    double t = 0.0;
    return SmApiEdgeClosestPoint(e, p, q, t, d) == SM_SUCCESS;
}

bool face_distance(SmFace *f, const SmPoint3d &p, double &d)
{
    SmPoint3d q;
    return SmApiFaceClosestPoint(f, p, q, d) == SM_SUCCESS;
}

int32_t pick_edges(SmBrep *b, const double *points, int64_t n, double tol, std::vector<SmEdge *> &out, int64_t *failed)
{
    SmTArray<SmEdge *> all;
    b->GetEdges(all);
    return pick(all, points, n, tol, static_cast<double>(b->GetTolerance()), edge_distance, out, failed);
}

int32_t pick_faces(SmBrep *b, const double *points, int64_t n, double tol, std::vector<SmFace *> &out, int64_t *failed)
{
    SmTArray<SmFace *> all;
    b->GetFaces(all);
    return pick(all, points, n, tol, static_cast<double>(b->GetTolerance()), face_distance, out, failed);
}

// The unit direction pointing into `e` from its end at `v`; false for a closed edge (both ends
// at `v`) or an edge whose tangent cannot be evaluated.
bool direction_from(const SmEdge *e, const SmVertex *v, SmVector3d &dir)
{
    SmCurve *c = e->GetCurve();
    if (c == nullptr || e->IsClosed())
    {
        return false;
    }
    SmExtent1d ivl = e->GetInterval();
    SmVector3d lo[2];
    SmVector3d hi[2];
    if (c->Evaluate(ivl.GetMin(), 1, TRUE, lo) != SM_SUCCESS || c->Evaluate(ivl.GetMax(), 1, TRUE, hi) != SM_SUCCESS)
    {
        return false;
    }
    SmVector3d at = v->GetPoint();
    bool start = (lo[0] - at).Length() <= (hi[0] - at).Length();
    dir = start ? lo[1] : -hi[1];
    double len = dir.Length();
    if (!(len > 0.0) || !std::isfinite(len))
    {
        return false;
    }
    dir = dir / len;
    return true;
}

// Grows `edges` by every filletable edge that continues one of them tangentially at a shared
// vertex, until no more do: a pick on one edge of a smooth chain rounds the whole chain.
void tangent_chains(std::vector<SmEdge *> &edges)
{
    std::unordered_set<SmEdge *> in(edges.begin(), edges.end());
    for (size_t i = 0; i < edges.size(); ++i)
    {
        SmEdge *e = edges[i];
        SmVertex *ends[2] = {nullptr, nullptr};
        e->GetVertices(ends[0], ends[1]);
        for (SmVertex *v : ends)
        {
            SmVector3d de;
            if (v == nullptr || !direction_from(e, v, de))
            {
                continue;
            }
            SmTArray<SmEdge *> around;
            v->GetEdges(around);
            for (ULONG k = 0; k < around.GetSize(); ++k)
            {
                SmEdge *f = around[k];
                SmVector3d df;
                if (f == nullptr || in.count(f) != 0 || f->GetBlendEdgeuse() == nullptr || !direction_from(f, v, df))
                {
                    continue;
                }
                // Continuing means leaving the vertex the opposite way from `e`.
                if (de.Dot(df) < -TANGENT_COS)
                {
                    in.insert(f);
                    edges.push_back(f);
                }
            }
        }
    }
}

// The radius of the ball whose tangent points on the two faces of `e` lie `distance` back from
// the edge, at the edge's midpoint: distance * tan(a / 2), with a the angle between b0 and b1,
// the unit directions into the faces perpendicular to the edge (distance itself on a right-
// angled edge). A constant-radius solver with a linear cross-section then cuts the chamfer.
bool chamfer_radius(SmEdge *e, double distance, double &radius)
{
    SmEdgeuse *u0 = e->GetBlendEdgeuse();
    SmEdgeuse *u1 = u0 == nullptr ? nullptr : u0->GetRadial();
    if (u1 == nullptr)
    {
        return false;
    }
    double t = e->GetInterval().Evaluate(0.5);
    SmPoint3d p0;
    SmPoint3d p1;
    SmVector3d b0;
    SmVector3d b1;
    if (u0->EvaluateBinormal(t, FALSE, p0, b0) != SM_SUCCESS || u1->EvaluateBinormal(t, FALSE, p1, b1) != SM_SUCCESS)
    {
        return false;
    }
    double c = std::max(-1.0, std::min(1.0, b0.Dot(b1)));
    radius = distance * std::sqrt((1.0 - c) / (1.0 + c));
    return std::isfinite(radius) && radius > 0.0;
}

// Chamfers `edges` of `b` in place: a constant-radius solver with a linear cross-section per
// edge (SMLib's constant-distance solver, which SmApiChamferFillet uses, leaves a hole rim
// open), bevelled corners, then stitching as SmApiFilletEdges does. The executive refers to `b`
// until it is destroyed, so it lives here alone and is gone before the caller hands `b` on or
// frees it.
int32_t chamfer_edges(SmBrep *b, const std::vector<SmEdge *> &edges, double distance)
{
    SmTArray<SmEdge *> es;
    SmTArray<SmFilletSurfaceGeneratorType> kinds;
    SmTArray<SmFilletSolverType> solvers;
    SmTArray<double> radii;
    for (SmEdge *e : edges)
    {
        double radius = 0.0;
        if (!chamfer_radius(e, distance, radius))
        {
            return static_cast<int32_t>(SM_ERR_INVALID_INPUT);
        }
        es.Add(e);
        kinds.Add(SM_FSG_LINEAR);
        solvers.Add(SM_FS_CONST_RADIUS);
        radii.Add(radius);
    }
    SmFilletExecutive exec(*SmApiGetOrCreateContext(), b);
    int32_t s = status(exec.SetFilletParameters(b, es, kinds, solvers, radii, nullptr, -1.0, -1.0, -1.0, -1.0));
    if (s != SMLIB_OK)
    {
        return s;
    }
    // SetFilletParameters chose the executive's options as SmApiFilletEdges and
    // SmApiSetBevelCorners run it (no global merge, which leaves a hole rim open).
    s = status(exec.CreateFilletCorners());
    if (s != SMLIB_OK)
    {
        return s;
    }
    SmTArray<SmFilletCorner *> &corners = exec.GetFilletCorners();
    for (ULONG j = 0; j < corners.GetSize(); ++j)
    {
        corners[j]->SetBevel(TRUE);
    }
    s = status(exec.DoFilleting());
    if (s != SMLIB_OK)
    {
        return s;
    }
    return status(b->StitchAndOrient());
}

// Hands a finished operation's body to the caller: it must have faces, and a solid input must
// give a closed manifold solid.
int32_t finish(Owned r, bool solid_in, smlib_brep **out)
{
    if (r == nullptr || r->GetNumFaces() == 0)
    {
        return SMLIB_ERR_EMPTY;
    }
    if (solid_in && !r->IsManifoldSolid())
    {
        return SMLIB_ERR_NOT_SOLID;
    }
    return take(r.release(), out);
}

} // namespace

extern "C" int32_t smlib_boolean(const smlib_brep *a, const smlib_brep *b, int32_t op, smlib_brep **out)
{
    SmBooleanOperationType kind;
    switch (op)
    {
    case SMLIB_BOOLEAN_UNION:
        kind = SM_BO_UNION;
        break;
    case SMLIB_BOOLEAN_SUBTRACT:
        kind = SM_BO_DIFFERENCE;
        break;
    case SMLIB_BOOLEAN_INTERSECT:
        kind = SM_BO_INTERSECTION;
        break;
    default:
        return SMLIB_ERR_ARG;
    }
    if (!valid(a) || !valid(b) || out == nullptr)
    {
        return SMLIB_ERR_ARG;
    }
    Call c("smlib_boolean", handle(a) + ", " + handle(b) + ", op=" + std::to_string(op));
    c.input("a", a);
    c.input("b", b);
    return c.run(
        [&]
        {
            int32_t s = SMLIB_OK;
            Owned ca = copy_of(a, s);
            if (s != SMLIB_OK)
            {
                return s;
            }
            Owned cb = copy_of(b, s);
            if (s != SMLIB_OK)
            {
                return s;
            }
            bool solid_in = ca->IsManifoldSolid() && cb->IsManifoldSolid();
            SmBrep *r = nullptr;
            s = status(SmApiBoolean(ca.get(), cb.get(), kind, r));
            if (s != SMLIB_OK)
            {
                // Both operands are still allocated (possibly modified); the owners free them.
                return s;
            }
            // Success: B is deleted and the result is A, edited in place.
            cb.release();
            if (r != nullptr && r != ca.get())
            {
                // Not what SM_API documents; SMLib owns A now, so leave it be.
                ca.release();
                return finish(Owned(r), solid_in, out);
            }
            return finish(std::move(ca), solid_in, out);
        }
    );
}

extern "C" int32_t smlib_brep_solid_count(const smlib_brep *brep, int64_t *out)
{
    if (!valid(brep) || out == nullptr)
    {
        return SMLIB_ERR_ARG;
    }
    Call c("smlib_brep_solid_count", handle(brep));
    c.input("brep", brep);
    return c.run(
        [&]
        {
            long solids = 0;
            long voids = 0;
            int32_t s = status(SmApiBrepMaterialCensus(brep->p, solids, voids));
            if (s == SMLIB_OK)
            {
                *out = static_cast<int64_t>(solids);
            }
            return s;
        }
    );
}

extern "C" int32_t smlib_fillet(const smlib_brep *brep, const double *points, int64_t n_points, double tolerance, double radius,
                                smlib_brep **out, int64_t *out_failed_point)
{
    if (!valid(brep) || out == nullptr || !points_ok(points, n_points, false) || !tolerance_ok(tolerance) || !std::isfinite(radius) ||
        !(radius > 0.0) || radius > 1e6)
    {
        return SMLIB_ERR_ARG;
    }
    Call c("smlib_fillet", handle(brep) + ", points=" + pts(points, n_points) + ", tolerance=" + num(tolerance) + ", radius=" + num(radius));
    c.input("brep", brep);
    return c.run(
        [&]
        {
            int32_t s = SMLIB_OK;
            Owned w = copy_of(brep, s);
            if (s != SMLIB_OK)
            {
                return s;
            }
            std::vector<SmEdge *> picked;
            s = pick_edges(w.get(), points, n_points, tolerance, picked, out_failed_point);
            if (s != SMLIB_OK)
            {
                return s;
            }
            tangent_chains(picked);
            bool solid_in = w->IsManifoldSolid();
            SmTArray<SmEdge *> edges;
            for (SmEdge *e : picked)
            {
                edges.Add(e);
            }
            s = status(SmApiFilletEdges(w.get(), edges, radius, 1, 1, 1.0));
            if (s != SMLIB_OK)
            {
                return s;
            }
            return finish(std::move(w), solid_in, out);
        }
    );
}

extern "C" int32_t smlib_chamfer(const smlib_brep *brep, const double *points, int64_t n_points, double tolerance, double distance,
                                 smlib_brep **out, int64_t *out_failed_point)
{
    if (!valid(brep) || out == nullptr || !points_ok(points, n_points, false) || !tolerance_ok(tolerance) || !std::isfinite(distance) ||
        !(distance > 0.0) || distance > 1e6)
    {
        return SMLIB_ERR_ARG;
    }
    Call c("smlib_chamfer", handle(brep) + ", points=" + pts(points, n_points) + ", tolerance=" + num(tolerance) + ", distance=" + num(distance));
    c.input("brep", brep);
    return c.run(
        [&]
        {
            int32_t s = SMLIB_OK;
            Owned w = copy_of(brep, s);
            if (s != SMLIB_OK)
            {
                return s;
            }
            std::vector<SmEdge *> picked;
            s = pick_edges(w.get(), points, n_points, tolerance, picked, out_failed_point);
            if (s != SMLIB_OK)
            {
                return s;
            }
            tangent_chains(picked);
            bool solid_in = w->IsManifoldSolid();
            s = chamfer_edges(w.get(), picked, distance);
            if (s != SMLIB_OK)
            {
                return s;
            }
            return finish(std::move(w), solid_in, out);
        }
    );
}

extern "C" int32_t smlib_shell(const smlib_brep *brep, const double *points, int64_t n_points, double tolerance, double thickness,
                               smlib_brep **out, int64_t *out_failed_point)
{
    if (!valid(brep) || out == nullptr || !points_ok(points, n_points, true) || !tolerance_ok(tolerance) || !std::isfinite(thickness) ||
        thickness == 0.0 || std::fabs(thickness) > 1e6)
    {
        return SMLIB_ERR_ARG;
    }
    Call c("smlib_shell", handle(brep) + ", points=" + pts(points, n_points) + ", tolerance=" + num(tolerance) + ", thickness=" + num(thickness));
    c.input("brep", brep);
    return c.run(
        [&]
        {
            int32_t s = SMLIB_OK;
            Owned w = copy_of(brep, s);
            if (s != SMLIB_OK)
            {
                return s;
            }
            std::vector<SmFace *> picked;
            s = pick_faces(w.get(), points, n_points, tolerance, picked, out_failed_point);
            if (s != SMLIB_OK)
            {
                return s;
            }
            bool solid_in = w->IsManifoldSolid();
            // An inward wall at least as thick as the body's inradius leaves no cavity, and
            // SMLib's offset of such a wall corrupts its working body (use-after-free under ASan
            // in GlueCoincidentBrepVertices and DeleteEdge). Half the smallest extent of the
            // (loose) bounding box is at least the inradius, so a wall that thick can only fail:
            // refuse it before SMLib sees it.
            if (thickness < 0.0)
            {
                SmPoint3d lo;
                SmPoint3d hi;
                s = status(SmApiBrepBoundingBox(w.get(), FALSE, lo, hi));
                if (s != SMLIB_OK)
                {
                    return s;
                }
                double least = std::min(hi[0] - lo[0], std::min(hi[1] - lo[1], hi[2] - lo[2]));
                if (-thickness >= 0.5 * least)
                {
                    return SMLIB_ERR_ARG;
                }
            }
            SmTArray<const SmFace *> open;
            for (SmFace *f : picked)
            {
                open.Add(f);
            }
            SmBrep *r = nullptr;
            s = status(SmApiShellBrepFull(w.get(), thickness, TRUE, TRUE, TRUE, open, r));
            if (s != SMLIB_OK)
            {
                delete r;
                return s;
            }
            return finish(Owned(r), solid_in, out);
        }
    );
}
