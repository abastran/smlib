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
#include <SmBSplineCurve.h>
#include <SmEdge.h>
#include <SmEdgeuse.h>
#include <SmFace.h>
#include <SmFaceuse.h>
#include <SmFilletExecutive.h>
#include <SmFilletStandardSolver.h>
#include <SmLoopuse.h>
#include <SmOffsetGeometryCreation.h>
#include <SmPlane.h>
#include <SmRegion.h>
#include <SmSSIData.h>
#include <SmShell.h>
#include <SmSurface.h>
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

// The entity nearest each point under the picking rule, each once; with `first`, the index of
// the point that first picked each. On failure, the index of the point that matched none or
// several goes to `failed`.
template <typename T, typename Dist>
int32_t pick(const SmTArray<T *> &all, const double *points, int64_t n, double tol, double model_tol, Dist dist,
             std::vector<T *> &picked, int64_t *failed, std::vector<int64_t> *first = nullptr)
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
            if (first != nullptr)
            {
                first->push_back(i);
            }
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

int32_t pick_edges(SmBrep *b, const double *points, int64_t n, double tol, std::vector<SmEdge *> &out, int64_t *failed,
                   std::vector<int64_t> *first = nullptr)
{
    SmTArray<SmEdge *> all;
    b->GetEdges(all);
    return pick(all, points, n, tol, static_cast<double>(b->GetTolerance()), edge_distance, out, failed, first);
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
// vertex, until no more do: a pick on one edge of a smooth chain rounds the whole chain. Each
// edge added takes the size (`sizes`, one per edge, grown alongside) of the edge it continues.
void tangent_chains(std::vector<SmEdge *> &edges, std::vector<double> &sizes)
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
                    sizes.push_back(sizes[i]);
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

// Why the last blend in this call failed, as the fillet executive noted it: dumped with the
// call's inputs (Call::input_text) so a failing blend can be read without a debugger.
std::string g_blend_why;

void note_blend_failure(SmFilletExecutive &exec, int32_t s)
{
    SmFilletErrorInfo *e = exec.GetFilletErrorInfo();
    if (e == nullptr)
    {
        return;
    }
    char head[256];
    std::snprintf(head, sizeof head, "status %d, fillet error %d, corner type %d, edge type %d, at (%.9g, %.9g, %.9g)\n", static_cast<int>(s),
                  static_cast<int>(e->GetErrorCode()), static_cast<int>(e->m_eCornerType), static_cast<int>(e->m_eEdgeType), e->m_vLocation[0],
                  e->m_vLocation[1], e->m_vLocation[2]);
    g_blend_why = head;
    for (const TCHAR *c : {e->m_cComment, e->m_cEntityName})
    {
        for (size_t i = 0; i < SM_TBLOCK_SIZE && c[i] != 0; ++i)
        {
            g_blend_why += static_cast<char>(c[i] < 128 ? c[i] : '?');
        }
        g_blend_why += "\n";
    }
}

// Chamfers `edges` of `b` in place, edge i with its first edgeuse's surface offset `first[i]`
// and the other's `second[i]` (equal for an equal-distance chamfer): a constant-radius solver
// with a linear cross-section per edge (SMLib's constant-distance solver, which
// SmApiChamferFillet uses, holds the chord across the chamfer, not the setbacks), bevelled
// corners, then stitching as SmApiFilletEdges does. The executive refers to `b` until it is destroyed, so it lives here
// alone and is gone before the caller hands `b` on or frees it.
int32_t chamfer_offsets(SmBrep *b, const std::vector<SmEdge *> &edges, const std::vector<double> &first, const std::vector<double> &second)
{
    SmTArray<SmEdge *> es;
    SmTArray<SmFilletSurfaceGeneratorType> kinds;
    SmTArray<SmFilletSolverType> solvers;
    SmTArray<double> r0;
    SmTArray<double> r1;
    for (size_t i = 0; i < edges.size() && i < first.size() && i < second.size(); ++i)
    {
        es.Add(edges[i]);
        kinds.Add(SM_FSG_LINEAR);
        solvers.Add(SM_FS_CONST_RADIUS);
        r0.Add(first[i]);
        r1.Add(second[i]);
    }
    SmFilletExecutive exec(*SmApiGetOrCreateContext(), b);
    int32_t s = status(exec.SetFilletParameters(b, es, kinds, solvers, r0, nullptr, -1.0, -1.0, -1.0, -1.0, 1, 1.0, &r1));
    if (s != SMLIB_OK)
    {
        note_blend_failure(exec, s);
        return s;
    }
    // SetFilletParameters chose the executive's options as SmApiFilletEdges and
    // SmApiSetBevelCorners run it (no global merge, which leaves a hole rim open).
    s = status(exec.CreateFilletCorners());
    if (s != SMLIB_OK)
    {
        note_blend_failure(exec, s);
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
        note_blend_failure(exec, s);
        return s;
    }
    return status(b->StitchAndOrient());
}

// Chamfers `edges` of `b` in place, edge i `distances[i]` back along each face.
int32_t chamfer_edges(SmBrep *b, const std::vector<SmEdge *> &edges, const std::vector<double> &distances)
{
    std::vector<double> radii;
    for (size_t i = 0; i < edges.size() && i < distances.size(); ++i)
    {
        double radius = 0.0;
        if (!chamfer_radius(edges[i], distances[i], radius))
        {
            return static_cast<int32_t>(SM_ERR_INVALID_INPUT);
        }
        radii.push_back(radius);
    }
    return chamfer_offsets(b, edges, radii, radii);
}

SmFaceuse *outer_faceuse(const SmBrep *b, const SmFace *f);

// The outward unit normal, at `u`'s edge's midpoint, of the face the edgeuse `u` bounds, and the
// unit direction into that face square to the edge.
bool face_frame(const SmBrep *b, SmEdgeuse *u, SmVector3d &normal, SmVector3d &into)
{
    SmEdge *e = u->GetEdge();
    SmLoopuse *lu = u->GetLoopuse();
    SmFace *f = lu == nullptr || lu->GetFaceuse() == nullptr ? nullptr : lu->GetFaceuse()->GetFace();
    SmFaceuse *fu = f == nullptr ? nullptr : outer_faceuse(b, f);
    SmSurface *surf = f == nullptr ? nullptr : f->GetSurface();
    if (e == nullptr || fu == nullptr || surf == nullptr)
    {
        return false;
    }
    SmPoint3d p;
    SmPoint2d uv;
    if (u->EvaluateBinormal(e->GetInterval().Evaluate(0.5), FALSE, p, into, nullptr, nullptr, &uv) != SM_SUCCESS ||
        surf->EvaluateNormal(uv, TRUE, TRUE, normal) != SM_SUCCESS)
    {
        return false;
    }
    double len = normal.Length();
    if (!(len > 0.0) || !std::isfinite(len))
    {
        return false;
    }
    normal = normal * ((fu->GetOrientation() == SM_OT_SAME ? 1.0 : -1.0) / len);
    return true;
}

// The surface offsets of `e` (its blend edgeuse's surface first) at which chamfer_offsets cuts
// a chamfer `d1` back along the face whose outward normal at the edge's midpoint leans most along
// `up`, and `second` back along the other (`angle`: the angle in radians the chamfer makes with
// the first face). In the section square to the edge, with b0, b1 the directions into the faces
// and theta the angle between them, setbacks s0, s1 lie at offsets
// r0 = (s1 - s0 cos theta) / sin theta and r1 = (s0 - s1 cos theta) / sin theta (s tan(theta/2)
// both, for equal setbacks); an angle phi to the first face gives s1 = s0 sin phi / sin(theta + phi).
bool sides_offsets(const SmBrep *b, SmEdge *e, double d1, double second, bool angle, const SmVector3d &up, double &r0, double &r1)
{
    SmEdgeuse *u0 = e->GetBlendEdgeuse();
    SmEdgeuse *u1 = u0 == nullptr ? nullptr : u0->GetRadial();
    SmVector3d n0;
    SmVector3d n1;
    SmVector3d b0;
    SmVector3d b1;
    if (u1 == nullptr || !face_frame(b, u0, n0, b0) || !face_frame(b, u1, n1, b1))
    {
        return false;
    }
    double c = std::max(-1.0, std::min(1.0, b0.Dot(b1)));
    double sn = std::sqrt(std::max(0.0, 1.0 - c * c));
    if (!(sn > 1e-6))
    {
        return false;
    }
    double s0 = d1;
    double s1 = second;
    if (angle)
    {
        double d = std::sin(std::acos(c) + second);
        if (!(d > 1e-9))
        {
            return false;
        }
        s1 = d1 * std::sin(second) / d;
    }
    if (n1.Dot(up) > n0.Dot(up))
    {
        std::swap(s0, s1);
    }
    r0 = (s1 - s0 * c) / sn;
    r1 = (s0 - s1 * c) / sn;
    return std::isfinite(r0) && std::isfinite(r1) && r0 > 0.0 && r1 > 0.0;
}

// Whether `r` is closed though not manifold: solids touching each other along an edge or at a
// vertex. Bounded regions, no wire edges or shell vertices, and every edge shared by its faces
// in pairs (four edgeuses per pair of face sides: none bounds a face on one side only).
bool closed_nonmanifold(const SmBrep &r)
{
    SmTArray<SmRegion *> regions;
    r.GetRegions(regions);
    if (regions.GetSize() < 2)
    {
        return false;
    }
    SmTArray<SmVertex *> points;
    r.GetShellVertices(points);
    if (points.GetSize() > 0)
    {
        return false;
    }
    SmTArray<SmEdge *> edges;
    r.GetEdges(edges);
    for (ULONG i = 0; i < edges.GetSize(); ++i)
    {
        SmEdge *e = edges[i];
        if (e == nullptr || e->IsWire())
        {
            return false;
        }
        SmTArray<SmEdgeuse *> uses;
        e->GetEdgeuses(uses);
        if (uses.GetSize() < 4 || uses.GetSize() % 4 != 0)
        {
            return false;
        }
    }
    return true;
}

// Rounds `edges` of `b` in place, edge i at `sizes[i]`, as SmApiFilletEdgesPerEdge does
// (SetFilletParameters, DoFillet, stitching) with one cross-section and solver for all: SMLIB_BLEND_ROUND a
// circular section at constant radius, SMLIB_BLEND_CURVATURE a G2 blend-curve section between
// the same contact lines, SMLIB_BLEND_CHORD a circular section at constant width (SMLib's
// constant-distance solver). `laws` (one per edge, or empty) makes each edge's radius vary.
int32_t round_edges(SmBrep *b, const std::vector<SmEdge *> &edges, const std::vector<double> &sizes, int32_t kind,
                    const std::vector<SmFilletLaw *> &laws)
{
    SmTArray<SmEdge *> es;
    SmTArray<double> values;
    SmTArray<SmFilletLaw *> law_array;
    for (size_t i = 0; i < edges.size() && i < sizes.size(); ++i)
    {
        es.Add(edges[i]);
        values.Add(sizes[i]);
    }
    for (SmFilletLaw *l : laws)
    {
        law_array.Add(l);
    }
    SmTArray<SmFilletSurfaceGeneratorType> sections;
    SmTArray<SmFilletSolverType> solvers;
    ULONG continuity = 1;
    switch (kind)
    {
    case SMLIB_BLEND_ROUND:
        sections.Add(SM_FSG_CIRCULAR);
        solvers.Add(laws.empty() ? SM_FS_CONST_RADIUS : SM_FS_VARIABLE_RADIUS);
        break;
    case SMLIB_BLEND_CURVATURE:
        sections.Add(SM_FSG_BLEND_CURVE);
        solvers.Add(SM_FS_CONST_RADIUS);
        continuity = 2;
        break;
    case SMLIB_BLEND_CHORD:
        sections.Add(SM_FSG_CIRCULAR);
        solvers.Add(SM_FS_CONST_DIST);
        break;
    default:
        return SMLIB_ERR_ARG;
    }
    if (!laws.empty() && (laws.size() != edges.size() || kind != SMLIB_BLEND_ROUND))
    {
        return SMLIB_ERR_ARG;
    }
    SmFilletExecutive exec(*SmApiGetOrCreateContext(), b);
    // The curvature-continuous section is the quintic whose inner control points lie a third
    // and two thirds of the way from each contact to the corner its tangents make.
    exec.SetBlendCurveAtThirds(kind == SMLIB_BLEND_CURVATURE ? TRUE : FALSE);
    int32_t s = status(exec.SetFilletParameters(b, es, sections, solvers, values, laws.empty() ? nullptr : &law_array, -1.0, -1.0,
                                                -1.0, -1.0, continuity, 1.0));
    if (s != SMLIB_OK)
    {
        note_blend_failure(exec, s);
        return s;
    }
    // SetFilletParameters chose the executive's options as SmApiFilletEdges runs it (its
    // SetDoGlobalMerge comes before, so SetFilletParameters resets it).
    s = status(exec.DoFillet());
    if (s != SMLIB_OK)
    {
        note_blend_failure(exec, s);
        return s;
    }
    return status(b->StitchAndOrient());
}

// Hands a finished operation's body to the caller: it must have faces, and a solid input must
// give a closed manifold solid, or (`nonmanifold`) a closed one whose solids touch along an edge
// or at a vertex.
int32_t finish(Owned r, bool solid_in, smlib_brep **out, bool nonmanifold = false)
{
    if (r == nullptr || r->GetNumFaces() == 0)
    {
        return SMLIB_ERR_EMPTY;
    }
    if (solid_in && !r->IsManifoldSolid() && !(nonmanifold && closed_nonmanifold(*r)))
    {
        return SMLIB_ERR_NOT_SOLID;
    }
    settle_tolerances(r.get());
    return take(r.release(), out);
}

// The faceuse of `f` on the outside of its body: the one bounding the infinite region (for a
// sheet, both do; then the upward one, along the surface normal).
SmFaceuse *outer_faceuse(const SmBrep *b, const SmFace *f)
{
    SmFaceuse *fu = f->GetUpwardFaceuse();
    if (fu == nullptr)
    {
        return nullptr;
    }
    if (fu->GetShell() == nullptr || fu->GetShell()->GetRegion() != b->GetInfiniteRegion())
    {
        fu = fu->GetMate();
    }
    return fu;
}

// Faces with their loose bounding boxes, for nearest-face queries.
struct FaceSet
{
    std::vector<SmFace *> faces;
    std::vector<std::pair<SmPoint3d, SmPoint3d>> boxes;

    bool add(SmFace *f)
    {
        SmPoint3d lo;
        SmPoint3d hi;
        if (SmApiFaceBoundingBox(f, FALSE, lo, hi) != SM_SUCCESS)
        {
            return false;
        }
        faces.push_back(f);
        boxes.emplace_back(lo, hi);
        return true;
    }

    // The distance from `p` to the nearest face, with the nearest point in `at`; HUGE_VAL when
    // none answers. A face whose box is farther than the nearest face so far is skipped.
    double nearest(const SmPoint3d &p, SmPoint3d *at = nullptr) const
    {
        double best = HUGE_VAL;
        for (size_t i = 0; i < faces.size(); ++i)
        {
            const SmPoint3d &lo = boxes[i].first;
            const SmPoint3d &hi = boxes[i].second;
            double box = 0.0;
            for (int k = 0; k < 3; ++k)
            {
                double out = std::max(lo[k] - p[k], std::max(0.0, p[k] - hi[k]));
                box += out * out;
            }
            if (std::sqrt(box) >= best)
            {
                continue;
            }
            SmPoint3d q;
            double d = HUGE_VAL;
            if (SmApiFaceClosestPoint(faces[i], p, q, d) == SM_SUCCESS && std::isfinite(d) && d < best)
            {
                best = d;
                if (at != nullptr)
                {
                    *at = q;
                }
            }
        }
        return best;
    }
};

// Whether an inward offset of `wall` (> 0) from every face of the solid `b` but `open` leaves a
// cavity: a point inside `b` farther than `wall`, and a margin, from every face that is offset.
// The answer is certified, never guessed: a point counts only as the centre of a ball inside
// `b` (one no face reaches into), reached from a face sample by stepping inward along the
// normal, then from ball to ball. SMLib's offset corrupts memory where no cavity remains
// (tangency-labs/solidcore#380), so a wall this finds no room for is refused before SMLib sees
// it. It may refuse a cavity thinner than the margin; it never accepts a wall at or beyond the
// body's inradius (for a closed shell) or the faces' clearance (for an open one).
bool leaves_cavity(SmBrep *b, const std::vector<SmFace *> &open, double wall)
{
    SmTArray<SmFace *> list;
    b->GetFaces(list);
    FaceSet all;
    FaceSet walls;
    for (ULONG i = 0; i < list.GetSize(); ++i)
    {
        if (!all.add(list[i]))
        {
            return false;
        }
        if (std::find(open.begin(), open.end(), list[i]) == open.end() && !walls.add(list[i]))
        {
            return false;
        }
    }
    if (walls.faces.empty())
    {
        return true;
    }
    double tol = static_cast<double>(b->GetTolerance());
    double margin = std::max(10.0 * tol, 1e-3 * wall);
    double need = wall + margin;
    // Face samples: a grid over each offset face's domain, the points inside its trims, each
    // with the inward unit normal there.
    constexpr int GRID = 5;
    constexpr size_t MAX_SAMPLES = 2048;
    std::vector<std::pair<SmPoint3d, SmVector3d>> samples;
    for (SmFace *f : walls.faces)
    {
        SmFaceuse *fu = outer_faceuse(b, f);
        SmSurface *surf = f->GetSurface();
        if (fu == nullptr || surf == nullptr)
        {
            continue;
        }
        double inward = fu->GetOrientation() == SM_OT_SAME ? -1.0 : 1.0;
        SmExtent2d dom = f->GetUVDomain();
        for (int i = 0; i < GRID && samples.size() < MAX_SAMPLES; ++i)
        {
            for (int j = 0; j < GRID && samples.size() < MAX_SAMPLES; ++j)
            {
                SmPoint2d uv = dom.Evaluate((i + 0.5) / GRID, (j + 0.5) / GRID);
                int cls = 0;
                if (SmApiFaceClassifyUV(f, uv, cls) != SM_SUCCESS || cls != SM_PC_FACE)
                {
                    continue;
                }
                SmVector3d at[1];
                SmVector3d n;
                if (surf->Evaluate(uv, 0, 0, TRUE, TRUE, FALSE, at) != SM_SUCCESS || surf->EvaluateNormal(uv, TRUE, TRUE, n) != SM_SUCCESS)
                {
                    continue;
                }
                double len = n.Length();
                if (!(len > 0.0) || !std::isfinite(len))
                {
                    continue;
                }
                samples.emplace_back(SmPoint3d(at[0]), n * (inward / len));
            }
        }
    }
    // A step of `h` inward from a sample ends at a ball's centre when no face is nearer to it
    // than `h` (the step crossed none, and the ball touches only the sample's face).
    auto ball = [&](const SmPoint3d &s, const SmVector3d &n, double h, SmPoint3d &c, double &r) {
        c = s + n * h;
        r = all.nearest(c);
        return std::isfinite(r) && r >= h * (1.0 - 1e-9);
    };
    // First the direct step of `need` from each sample.
    for (const auto &[s, n] : samples)
    {
        SmPoint3d c;
        double r = 0.0;
        if (ball(s, n, need, c, r) && walls.nearest(c) >= need)
        {
            return true;
        }
    }
    // Then climbing, from a small ball at each sample: step away from the nearest offset face
    // by most of the ball's radius (staying inside the ball, so inside the body) while the
    // distance from the offset faces grows.
    constexpr int STEPS = 48;
    for (const auto &[s, n] : samples)
    {
        SmPoint3d c;
        double r = 0.0;
        if (!ball(s, n, std::max(margin, 1e3 * tol), c, r))
        {
            continue;
        }
        for (int k = 0; k < STEPS; ++k)
        {
            SmPoint3d q;
            double d = walls.nearest(c, &q);
            if (d >= need)
            {
                return true;
            }
            if (!(d > 0.0) || !std::isfinite(d))
            {
                break;
            }
            SmPoint3d next = c + (c - q) * (0.9 * r / d);
            double rn = all.nearest(next);
            if (!(rn > 0.0) || !std::isfinite(rn) || walls.nearest(next) <= d * (1.0 + 1e-9))
            {
                break;
            }
            c = next;
            r = rn;
        }
    }
    return false;
}

// Moves the picked faces of `b` (in place) `distance` along their outward normals: each gets
// its surface's offset, and SMLib's local operation extends or trims the faces around them to
// meet it. The local operation takes the surfaces over, also when it fails.
int32_t move_faces(SmBrep *b, const std::vector<SmFace *> &picked, double distance)
{
    double tol = static_cast<double>(b->GetTolerance());
    SmOffsetGeometryCreation make(std::fabs(distance), tol);
    SmTArray<SmFace *> faces;
    SmTArray<SmSurface *> surfaces;
    for (SmFace *f : picked)
    {
        SmFaceuse *fu = outer_faceuse(b, f);
        if (fu == nullptr)
        {
            return static_cast<int32_t>(SM_ERR_INVALID_INPUT);
        }
        // The offset goes along the faceuse's normal: outward from the outer one, inward from
        // its mate.
        if (distance < 0.0)
        {
            fu = fu->GetMate();
        }
        SmSurface *s = nullptr;
        SmSSIData crossings;
        int32_t st = status(make.FaceuseOffset(*SmApiGetOrCreateContext(), fu, TRUE, s, crossings));
        if (st != SMLIB_OK || s == nullptr)
        {
            delete s;
            for (ULONG k = 0; k < surfaces.GetSize(); ++k)
            {
                delete surfaces[k];
            }
            return st != SMLIB_OK ? st : static_cast<int32_t>(SM_ERR_INVALID_INPUT);
        }
        faces.Add(f);
        surfaces.Add(s);
    }
    return status(b->LocalOperation(faces, surfaces));
}

// Tilts each face of `picked` (planes only) by `angle` radians about its line on the neutral
// plane (through `origin`, normal `normal`), so that it leans in going along `pull`, then lets
// SMLib's local operation extend or trim the faces around it to meet the tilted planes. A face
// must be parallel to the pull (its normal square to it); the hinge is where a line from a point
// inside the face, along the pull, meets the neutral plane.
int32_t draft_faces(SmBrep *b, const std::vector<SmFace *> &picked, const SmPoint3d &origin, const SmVector3d &normal, const SmVector3d &pull,
                    double angle)
{
    const SmContext &ctx = *SmApiGetOrCreateContext();
    SmTArray<SmFace *> faces;
    SmTArray<SmSurface *> surfaces;
    auto give_up = [&](int32_t st)
    {
        for (ULONG k = 0; k < surfaces.GetSize(); ++k)
        {
            delete surfaces[k];
        }
        return st;
    };
    for (SmFace *f : picked)
    {
        SmSurface *surf = f->GetSurface();
        SmFaceuse *fu = outer_faceuse(b, f);
        if (surf == nullptr || fu == nullptr)
        {
            return give_up(static_cast<int32_t>(SM_ERR_INVALID_INPUT));
        }
        if (!surf->IsKindOf(SmPlane_TYPE))
        {
            return give_up(SMLIB_ERR_UNSUPPORTED);
        }
        // The face's outward normal: the plane's, turned with the outer faceuse.
        SmVector3d n = static_cast<const SmPlane *>(surf)->GetPosition().GetZAxis();
        if (fu->GetOrientation() != SM_OT_SAME)
        {
            n = -n;
        }
        if (std::fabs(n.Dot(pull)) > 1e-6)
        {
            return give_up(SMLIB_ERR_UNSUPPORTED);
        }
        SmPoint3d inside;
        if (SmApiFaceInternalPoint(f, inside) != SM_SUCCESS)
        {
            return give_up(static_cast<int32_t>(SM_ERR_INVALID_INPUT));
        }
        double along = pull.Dot(normal);
        if (std::fabs(along) < 1e-12)
        {
            return give_up(SMLIB_ERR_ARG);
        }
        SmPoint3d hinge = inside + pull * ((origin - inside).Dot(normal) / along);
        // Turning n about n x pull by `angle` gives n cos(angle) + pull sin(angle).
        SmVector3d axis = n * pull; // the cross product
        double len = axis.Length();
        if (!(len > 0.5))
        {
            return give_up(SMLIB_ERR_UNSUPPORTED);
        }
        axis = axis / len;
        SmSurface *tilted = nullptr;
        if (surf->Copy(ctx, tilted) != SM_SUCCESS || tilted == nullptr)
        {
            delete tilted;
            return give_up(static_cast<int32_t>(SM_ERR_INVALID_INPUT));
        }
        SmAxis2Placement turn;
        turn.RotateAboutAxisAtPoint(angle, hinge, axis);
        if (tilted->Transform(turn) != SM_SUCCESS)
        {
            delete tilted;
            return give_up(static_cast<int32_t>(SM_ERR_INVALID_INPUT));
        }
        faces.Add(f);
        surfaces.Add(tilted);
    }
    return status(b->LocalOperation(faces, surfaces));
}

} // namespace

extern "C" int32_t smlib_boolean(const smlib_brep *a, const smlib_brep *b, int32_t op, smlib_brep **out)
{
    return smlib_boolean_with(a, b, op, SMLIB_RESULT_NONMANIFOLD, out);
}

extern "C" int32_t smlib_boolean_with(const smlib_brep *a, const smlib_brep *b, int32_t op, int32_t flags, smlib_brep **out)
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
    if (!valid(a) || !valid(b) || out == nullptr || (flags & ~SMLIB_RESULT_NONMANIFOLD) != 0)
    {
        return SMLIB_ERR_ARG;
    }
    const bool nonmanifold = (flags & SMLIB_RESULT_NONMANIFOLD) != 0;
    Call c("smlib_boolean", handle(a) + ", " + handle(b) + ", op=" + std::to_string(op) + ", flags=" + std::to_string(flags));
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
                return finish(Owned(r), solid_in, out, nonmanifold);
            }
            return finish(std::move(ca), solid_in, out, nonmanifold);
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

namespace
{

bool size_ok(double x)
{
    return std::isfinite(x) && x > 0.0 && x <= 1e6;
}

// smlib_blend's work, logged as `fn`: `sizes` holds one size per point, or one for all.
int32_t blend_call(const char *fn, const smlib_brep *brep, int32_t kind, const double *points, const double *sizes, bool one_size,
                   int64_t n_points, double tolerance, smlib_brep **out, int64_t *out_failed_point)
{
    bool kind_ok = kind == SMLIB_BLEND_ROUND || kind == SMLIB_BLEND_CURVATURE || kind == SMLIB_BLEND_CHORD || kind == SMLIB_BLEND_CHAMFER;
    if (!valid(brep) || out == nullptr || !kind_ok || !points_ok(points, n_points, false) || !tolerance_ok(tolerance) || sizes == nullptr)
    {
        return SMLIB_ERR_ARG;
    }
    int64_t n_sizes = one_size ? 1 : n_points;
    for (int64_t i = 0; i < n_sizes; ++i)
    {
        if (!size_ok(sizes[i]))
        {
            return SMLIB_ERR_ARG;
        }
    }
    Call c(fn, handle(brep) + ", kind=" + std::to_string(kind) + ", points=" + pts(points, n_points) + ", sizes=" +
                   vec(sizes, static_cast<int>(n_sizes)) + ", tolerance=" + num(tolerance));
    c.input("brep", brep);
    c.input_text("why", [](const std::string &path) {
        if (FILE *f = std::fopen(path.c_str(), "w"))
        {
            std::fputs(g_blend_why.c_str(), f);
            std::fclose(f);
        }
    });
    return c.run(
        [&]
        {
            g_blend_why.clear();
            int32_t s = SMLIB_OK;
            Owned w = copy_of(brep, s);
            if (s != SMLIB_OK)
            {
                return s;
            }
            std::vector<SmEdge *> picked;
            std::vector<int64_t> first;
            s = pick_edges(w.get(), points, n_points, tolerance, picked, out_failed_point, &first);
            if (s != SMLIB_OK)
            {
                return s;
            }
            std::vector<double> each;
            for (int64_t k : first)
            {
                each.push_back(sizes[one_size ? 0 : k]);
            }
            tangent_chains(picked, each);
            bool solid_in = w->IsManifoldSolid();
            s = kind == SMLIB_BLEND_CHAMFER ? chamfer_edges(w.get(), picked, each) : round_edges(w.get(), picked, each, kind, {});
            if (s != SMLIB_OK)
            {
                return s;
            }
            return finish(std::move(w), solid_in, out);
        }
    );
}

// Whether `knots` (n_knots) and `radii` (n_radii) are a clamped B-spline of `degree` from 0
// to 1: degree 1 to 3, 2 to 100 control radii, knots rising, each end repeated degree + 1 times,
// radii finite and non-negative.
bool law_ok(int32_t degree, const double *knots, int64_t n_knots, const double *radii, int64_t n_radii)
{
    if (degree < 1 || degree > 3 || n_radii < 2 || n_radii > 100 || n_radii <= degree || n_knots != n_radii + degree + 1 ||
        knots == nullptr || radii == nullptr || !finite(knots, static_cast<int>(n_knots)) || !finite(radii, static_cast<int>(n_radii)))
    {
        return false;
    }
    for (int64_t i = 0; i < n_knots; ++i)
    {
        bool low_end = i <= degree;
        bool high_end = i >= n_radii;
        if ((low_end && knots[i] != 0.0) || (high_end && knots[i] != 1.0) || (i > 0 && knots[i] < knots[i - 1]))
        {
            return false;
        }
    }
    for (int64_t i = 0; i < n_radii; ++i)
    {
        if (radii[i] < 0.0 || radii[i] > 1e6)
        {
            return false;
        }
    }
    return true;
}

// Whether the law has a corner: a knot inside (0, 1) repeated `degree` times (where the spline
// is only C0) with the slopes of the control legs on either side apart. SMLib's section fit
// runs smoothly across such a corner and misses the law, so it is refused, never built wrong.
bool law_has_corner(int32_t degree, const double *knots, const double *radii, int64_t n_radii)
{
    for (int64_t i = degree + 1; i < n_radii;)
    {
        int64_t m = 1;
        while (i + m < n_radii && knots[i + m] == knots[i])
        {
            ++m;
        }
        if (m > degree)
        {
            return true; // the radius itself jumps
        }
        // Repeated `degree` times, the curve passes through control radius i - 1 at the knot.
        if (m == degree && i >= 2)
        {
            double hl = knots[i] - knots[i - 1];
            double hr = knots[i + m] - knots[i];
            if (hl > 0.0 && hr > 0.0)
            {
                double left = (radii[i - 1] - radii[i - 2]) / hl;
                double right = (radii[i] - radii[i - 1]) / hr;
                if (std::fabs(left - right) > 1e-9 * (1.0 + std::fabs(left) + std::fabs(right)))
                {
                    return true;
                }
            }
        }
        i += m;
    }
    return false;
}

// The radius law along `e`: the B-spline (radius as its x) over the edge's parameter interval,
// reversed when `start` is nearer the interval's end than its start. Null if it can't be built.
std::unique_ptr<SmBSplineFilletLaw> law_along(SmEdge *e, const double start[3], int32_t degree, const double *knots, int64_t n_knots,
                                              const double *radii, int64_t n_radii)
{
    SmCurve *curve = e->GetCurve();
    if (curve == nullptr)
    {
        return nullptr;
    }
    SmExtent1d ivl = e->GetInterval();
    SmVector3d lo[1];
    SmVector3d hi[1];
    if (curve->Evaluate(ivl.GetMin(), 0, TRUE, lo) != SM_SUCCESS || curve->Evaluate(ivl.GetMax(), 0, TRUE, hi) != SM_SUCCESS)
    {
        return nullptr;
    }
    SmVector3d at(start[0], start[1], start[2]);
    bool reverse = (hi[0] - at).Length() < (lo[0] - at).Length();
    const SmContext &ctx = *SmApiGetOrCreateContext();
    SmTArray<SmPoint3d> cps;
    for (int64_t i = 0; i < n_radii; ++i)
    {
        cps.Add(SmPoint3d(radii[reverse ? n_radii - 1 - i : i], 0.0, 0.0));
    }
    // The knots as SMLib takes them: each distinct value once, with its multiplicity.
    SmTArray<double> unique;
    SmTArray<ULONG> mults;
    for (int64_t i = 0; i < n_knots; ++i)
    {
        double k = reverse ? 1.0 - knots[n_knots - 1 - i] : knots[i];
        if (unique.GetSize() > 0 && unique[unique.GetSize() - 1] == k)
        {
            mults[mults.GetSize() - 1] += 1;
        }
        else
        {
            unique.Add(k);
            mults.Add(1);
        }
    }
    SmBSplineCurve *c = nullptr;
    if (SmBSplineCurve::CreateCanonical(ctx, 3, static_cast<ULONG>(degree), cps, SM_CF_UNSPECIFIED, mults, unique, SM_KT_UNSPECIFIED,
                                        nullptr, nullptr, c) != SM_SUCCESS ||
        c == nullptr)
    {
        delete c;
        return nullptr;
    }
    // The law takes the curve; its map is the whole curve, 0 to 1.
    return std::unique_ptr<SmBSplineFilletLaw>(new (ctx) SmBSplineFilletLaw(c, SmExtent1d(0.0, 1.0)));
}

} // namespace

extern "C" int32_t smlib_blend(const smlib_brep *brep, int32_t kind, const double *points, const double *sizes, int64_t n_points,
                               double tolerance, smlib_brep **out, int64_t *out_failed_point)
{
    return blend_call("smlib_blend", brep, kind, points, sizes, false, n_points, tolerance, out, out_failed_point);
}

extern "C" int32_t smlib_fillet(const smlib_brep *brep, const double *points, int64_t n_points, double tolerance, double radius,
                                smlib_brep **out, int64_t *out_failed_point)
{
    return blend_call("smlib_fillet", brep, SMLIB_BLEND_ROUND, points, &radius, true, n_points, tolerance, out, out_failed_point);
}

extern "C" int32_t smlib_chamfer(const smlib_brep *brep, const double *points, int64_t n_points, double tolerance, double distance,
                                 smlib_brep **out, int64_t *out_failed_point)
{
    return blend_call("smlib_chamfer", brep, SMLIB_BLEND_CHAMFER, points, &distance, true, n_points, tolerance, out, out_failed_point);
}

extern "C" int32_t smlib_fillet_variable(const smlib_brep *brep, const double *point, double tolerance, const double *start,
                                         int32_t degree, const double *knots, int64_t n_knots, const double *radii, int64_t n_radii,
                                         smlib_brep **out)
{
    if (!valid(brep) || out == nullptr || !points_ok(point, 1, false) || start == nullptr || !finite(start, 3) ||
        !tolerance_ok(tolerance) || !law_ok(degree, knots, n_knots, radii, n_radii))
    {
        return SMLIB_ERR_ARG;
    }
    if (law_has_corner(degree, knots, radii, n_radii))
    {
        return SMLIB_ERR_UNSUPPORTED;
    }
    Call c("smlib_fillet_variable", handle(brep) + ", point=" + vec(point, 3) + ", tolerance=" + num(tolerance) + ", start=" + vec(start, 3) +
                                        ", degree=" + std::to_string(degree) + ", knots=" + vec(knots, static_cast<int>(n_knots)) +
                                        ", radii=" + vec(radii, static_cast<int>(n_radii)));
    c.input("brep", brep);
    c.input_text("why", [](const std::string &path) {
        if (FILE *f = std::fopen(path.c_str(), "w"))
        {
            std::fputs(g_blend_why.c_str(), f);
            std::fclose(f);
        }
    });
    return c.run(
        [&]
        {
            g_blend_why.clear();
            int32_t s = SMLIB_OK;
            Owned w = copy_of(brep, s);
            if (s != SMLIB_OK)
            {
                return s;
            }
            std::vector<SmEdge *> picked;
            s = pick_edges(w.get(), point, 1, tolerance, picked, nullptr);
            if (s != SMLIB_OK || picked.size() != 1)
            {
                return s != SMLIB_OK ? s : SMLIB_ERR_PICK_NONE;
            }
            std::unique_ptr<SmBSplineFilletLaw> law = law_along(picked[0], start, degree, knots, n_knots, radii, n_radii);
            if (law == nullptr)
            {
                return static_cast<int32_t>(SM_ERR_INVALID_INPUT);
            }
            bool solid_in = w->IsManifoldSolid();
            // The solver evaluates the law (it keeps a reference), so the law outlives the
            // executive inside round_edges. The size is the start radius.
            std::vector<double> size{std::max(radii[0], radii[n_radii - 1])};
            s = round_edges(w.get(), picked, size, SMLIB_BLEND_ROUND, {law.get()});
            if (s != SMLIB_OK)
            {
                return s;
            }
            return finish(std::move(w), solid_in, out);
        }
    );
}

extern "C" int32_t smlib_chamfer_sides(const smlib_brep *brep, const double *points, int64_t n_points, double tolerance, double distance,
                                       double second, int32_t second_is_angle, const double *first, smlib_brep **out, int64_t *out_failed_point)
{
    bool angle = second_is_angle == 1;
    bool second_ok = angle ? std::isfinite(second) && second > 1e-6 && second < SM_PI - 1e-6 : size_ok(second);
    if (!valid(brep) || out == nullptr || !points_ok(points, n_points, false) || !tolerance_ok(tolerance) || !size_ok(distance) ||
        (second_is_angle != 0 && second_is_angle != 1) || !second_ok || first == nullptr || !finite(first, 3))
    {
        return SMLIB_ERR_ARG;
    }
    SmVector3d up(first[0], first[1], first[2]);
    if (!(up.Length() > 0.0))
    {
        return SMLIB_ERR_ARG;
    }
    Call c("smlib_chamfer_sides", handle(brep) + ", points=" + pts(points, n_points) + ", tolerance=" + num(tolerance) + ", distance=" + num(distance) +
                                      ", second=" + num(second) + ", angle=" + std::to_string(second_is_angle) + ", first=" + vec(first, 3));
    c.input("brep", brep);
    c.input_text("why", [](const std::string &path) {
        if (FILE *f = std::fopen(path.c_str(), "w"))
        {
            std::fputs(g_blend_why.c_str(), f);
            std::fclose(f);
        }
    });
    return c.run(
        [&]
        {
            g_blend_why.clear();
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
            std::vector<double> unused(picked.size(), distance);
            tangent_chains(picked, unused);
            bool solid_in = w->IsManifoldSolid();
            std::vector<double> r0;
            std::vector<double> r1;
            for (SmEdge *e : picked)
            {
                double a = 0.0;
                double b = 0.0;
                if (!sides_offsets(w.get(), e, distance, second, angle, up, a, b))
                {
                    return static_cast<int32_t>(SM_ERR_INVALID_INPUT);
                }
                r0.push_back(a);
                r1.push_back(b);
            }
            s = chamfer_offsets(w.get(), picked, r0, r1);
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
            // An inward wall with no room left inside it corrupts SMLib's working body (use-after-
            // free under ASan in GlueCoincidentBrepVertices and DeleteEdge,
            // tangency-labs/solidcore#380): refuse it before SMLib sees it.
            if (thickness < 0.0 && solid_in && !leaves_cavity(w.get(), picked, -thickness))
            {
                return SMLIB_ERR_NO_ROOM;
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

extern "C" int32_t smlib_offset_faces(const smlib_brep *brep, const double *points, int64_t n_points, double tolerance, double distance,
                                      smlib_brep **out, int64_t *out_failed_point)
{
    if (!valid(brep) || out == nullptr || !points_ok(points, n_points, false) || !tolerance_ok(tolerance) || !std::isfinite(distance) ||
        distance == 0.0 || std::fabs(distance) > 1e6)
    {
        return SMLIB_ERR_ARG;
    }
    Call c("smlib_offset_faces",
           handle(brep) + ", points=" + pts(points, n_points) + ", tolerance=" + num(tolerance) + ", distance=" + num(distance));
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
            if (!solid_in)
            {
                return SMLIB_ERR_ARG;
            }
            // Faces moved in so far that no material is left behind them reach the same offset
            // code as a shell with no cavity: the body must reach deeper than the move behind
            // the moved faces.
            if (distance < 0.0)
            {
                SmTArray<SmFace *> all;
                w->GetFaces(all);
                std::vector<SmFace *> others;
                for (ULONG i = 0; i < all.GetSize(); ++i)
                {
                    if (std::find(picked.begin(), picked.end(), all[i]) == picked.end())
                    {
                        others.push_back(all[i]);
                    }
                }
                if (!leaves_cavity(w.get(), others, -distance))
                {
                    return SMLIB_ERR_NO_ROOM;
                }
            }
            s = move_faces(w.get(), picked, distance);
            if (s != SMLIB_OK)
            {
                return s;
            }
            return finish(std::move(w), solid_in, out);
        }
    );
}

extern "C" int32_t smlib_draft(const smlib_brep *brep, const double *points, int64_t n_points, double tolerance, const double neutral_origin[3],
                               const double neutral_normal[3], const double pull[3], double angle, smlib_brep **out, int64_t *out_failed_point)
{
    auto unit = [](const double *v, SmVector3d &u)
    {
        if (v == nullptr || !finite(v, 3))
        {
            return false;
        }
        u = SmVector3d(v[0], v[1], v[2]);
        double len = u.Length();
        if (!(len > 1e-12) || !std::isfinite(len))
        {
            return false;
        }
        u = u / len;
        return true;
    };
    SmVector3d normal;
    SmVector3d along;
    if (!valid(brep) || out == nullptr || !points_ok(points, n_points, false) || !tolerance_ok(tolerance) || neutral_origin == nullptr ||
        !finite(neutral_origin, 3) || std::fabs(neutral_origin[0]) > 1e9 || std::fabs(neutral_origin[1]) > 1e9 || std::fabs(neutral_origin[2]) > 1e9 ||
        !unit(neutral_normal, normal) || !unit(pull, along) || !std::isfinite(angle) || angle == 0.0 || std::fabs(angle) >= 1.5697963267948966)
    {
        return SMLIB_ERR_ARG;
    }
    Call c("smlib_draft", handle(brep) + ", points=" + pts(points, n_points) + ", tolerance=" + num(tolerance) + ", neutral=" + vec(neutral_origin, 3) +
                              "/" + vec(neutral_normal, 3) + ", pull=" + vec(pull, 3) + ", angle=" + num(angle));
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
            if (!w->IsManifoldSolid())
            {
                return SMLIB_ERR_ARG;
            }
            std::vector<SmFace *> picked;
            s = pick_faces(w.get(), points, n_points, tolerance, picked, out_failed_point);
            if (s != SMLIB_OK)
            {
                return s;
            }
            s = draft_faces(w.get(), picked, SmPoint3d(neutral_origin[0], neutral_origin[1], neutral_origin[2]), normal, along, angle);
            if (s != SMLIB_OK)
            {
                return s;
            }
            return finish(std::move(w), true, out);
        }
    );
}

extern "C" int32_t smlib_thicken(const smlib_brep *brep, double thickness, int32_t symmetric, smlib_brep **out)
{
    if (!valid(brep) || out == nullptr || !std::isfinite(thickness) || thickness == 0.0 || std::fabs(thickness) > 1e6)
    {
        return SMLIB_ERR_ARG;
    }
    Call c("smlib_thicken", handle(brep) + ", thickness=" + num(thickness) + ", symmetric=" + std::to_string(symmetric));
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
            if (w->IsManifoldSolid())
            {
                return SMLIB_ERR_ARG;
            }
            // Symmetric: the sheet moved back by half first, then thickened by the whole.
            if (symmetric != 0)
            {
                SmBrep *moved = nullptr;
                s = status(SmApiOffsetBrepFull(w.get(), -0.5 * thickness, TRUE, TRUE, moved));
                if (s != SMLIB_OK || moved == nullptr)
                {
                    delete moved;
                    return s != SMLIB_OK ? s : SMLIB_ERR_EMPTY;
                }
                // OffsetBrep leaves its input allocated (possibly changed) beside the result.
                if (moved != w.get())
                {
                    w.reset(moved);
                }
            }
            SmTArray<const SmFace *> none;
            SmBrep *r = nullptr;
            s = status(SmApiShellBrepFull(w.get(), thickness, TRUE, TRUE, TRUE, none, r));
            if (s != SMLIB_OK)
            {
                delete r;
                return s;
            }
            return finish(Owned(r), true, out);
        }
    );
}
