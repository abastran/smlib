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

#include <vector>
#include "smlib_pick.h"

#include <SmApiBrep.h>
#include <SmPrimitiveCreation.h>
#include <SmApiQueries.h>
#include <SmAttribute.h>
#include <SmBSplineCurve.h>
#include <SmEdge.h>
#include <SmEdgeuse.h>
#include <SmFace.h>
#include <SmFaceuse.h>
#include <SmFilletExecutive.h>
#include <SmFilletStandardSolver.h>
#include <SmLoopuse.h>
#include <SmOffsetFaces.h>
#include <SmOffsetGeometryCreation.h>
#include <SmPlane.h>
#include <SmRegion.h>
#include <SmSSIData.h>
#include <SmShell.h>
#include <SmSurface.h>
#include <SmVertex.h>

#include <memory>
#include <unordered_map>
#include <unordered_set>

using namespace smlib_shim;
using namespace smlib_pick;

namespace
{

// Upper bound on the points of one call: far more than any selection, small enough that the
// nearest-entity search (points x entities) stays cheap.
constexpr int64_t MAX_POINTS = 4096;

// Owns an SmBrep until it is handed on.
using Owned = std::unique_ptr<SmBrep>;

// A mark the shim puts on entities before SMLib runs: on every face of the working body, or on
// each surface it hands SMLib to move a face onto. SMLib hands an attribute on through every
// split, merge and copy it makes (SmAObject::Notify, as for the history tags of
// smlib_history.cpp) and builds new topology without any, so afterwards a face without a mark is
// new to the operation, and a face whose surface has one lies where the shim asked a face to go.
// (A shell's wall is no such test: SMLib makes it of copies of the body's faces and surfaces.)
class Before : public SmAttribute
{
public:
    // Above SMLib's own ids, beside the history tag's (SM_AI_USER_1 + 7).
    static constexpr ULONG ID = SM_AI_USER_1 + 8;
    static constexpr SM_TYPE TYPE = APPLICATION_BASE_TYPE + 4712;

    Before() : SmAttribute(ID, SM_AB_COPY)
    {
    }
    SmAttribute *MakeCopy(const SmContext &c) const override
    {
        return new (c) Before();
    }
    SM_TYPE GetType() const override
    {
        return TYPE;
    }
    SmBoolean IsKindOf(SM_TYPE t) const override
    {
        return t == TYPE ? TRUE : SmAttribute::IsKindOf(t);
    }

    static bool on(const SmAObject *x)
    {
        SmAttribute *a = x == nullptr ? nullptr : x->FindAttribute(ID);
        return a != nullptr && a->GetType() == TYPE;
    }
    static void mark(SmAObject *x)
    {
        if (x != nullptr && !on(x))
        {
            x->AddAttribute(new (*x->GetContext()) Before());
        }
    }
    static void clear(SmAObject *x)
    {
        if (x == nullptr)
        {
            return;
        }
        if (SmAttribute *a = x->FindAttribute(ID); a != nullptr && a->GetType() == TYPE)
        {
            x->RemoveAttribute(a);
        }
    }
};

// Marks every face of `b` as there before the operation.
void mark_faces(SmBrep *b)
{
    SmTArray<SmFace *> faces;
    b->GetFaces(faces);
    for (ULONG i = 0; i < faces.GetSize(); ++i)
    {
        Before::mark(faces[i]);
    }
}

// The faces of `b` an operation made (without a Before mark), and the faces on a marked surface
// (`on_marked_surface`); every mark on `b`'s faces and surfaces is removed.
void read_marks(SmBrep *b, size_t &made, size_t &on_marked_surface)
{
    made = 0;
    on_marked_surface = 0;
    SmTArray<SmFace *> faces;
    b->GetFaces(faces);
    for (ULONG i = 0; i < faces.GetSize(); ++i)
    {
        if (!Before::on(faces[i]))
        {
            ++made;
        }
        if (Before::on(faces[i]->GetSurface()))
        {
            ++on_marked_surface;
        }
    }
    for (ULONG i = 0; i < faces.GetSize(); ++i)
    {
        Before::clear(faces[i]);
        Before::clear(faces[i]->GetSurface());
    }
}

// The status of SMLib's local operation on `b`, which was handed marked surfaces: when it
// reports success, some face of `b` must lie on one of them, or it moved none
// (SMLIB_ERR_NO_EFFECT). Removes the marks.
int32_t moved_faces(SmBrep *b, int32_t s)
{
    size_t made = 0;
    size_t moved = 0;
    read_marks(b, made, moved);
    return s != SMLIB_OK || moved > 0 ? s : SMLIB_ERR_NO_EFFECT;
}

// Faces made by an operation on `b` whose faces were marked by mark_faces: none means SMLib
// reported success without doing anything (SMLIB_ERR_NO_EFFECT). Removes the marks.
int32_t made_faces(SmBrep *b)
{
    size_t made = 0;
    size_t still = 0;
    read_marks(b, made, still);
    return made > 0 ? SMLIB_OK : SMLIB_ERR_NO_EFFECT;
}

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

// The status a blend that failed with status `s` reports, its reason noted (g_blend_why): a
// blend the executive found larger than its faces allow (SM_FILERR_NO_ROOM: a concave blend whose
// rail leaves the face it rolls on, a round whose rail runs on into the body's material, a cliff
// overflow with no room on its faces' continuations) is SMLIB_ERR_NO_ROOM; a convex blend overflowing a
// cliff that needs a continuation SMLib cannot make exactly (SM_FILERR_NOT_EXACT)
// SMLIB_ERR_NOT_EXACT; any other failure `s`.
int32_t blend_failure(SmFilletExecutive &exec, int32_t s)
{
    SmFilletErrorInfo *e = exec.GetFilletErrorInfo();
    if (e == nullptr)
    {
        return s;
    }
    char head[256];
    std::snprintf(head, sizeof head, "status %d, fillet error %d, corner type %d, edge type %d, at (%.9g, %.9g, %.9g)\n", static_cast<int>(s),
                  static_cast<int>(e->GetErrorCode()), static_cast<int>(e->m_eCornerType), static_cast<int>(e->m_eEdgeType), e->m_vLocation[0],
                  e->m_vLocation[1], e->m_vLocation[2]);
    g_blend_why = head;
    for (const TCHAR *c : {e->m_cComment, e->m_cEntityName})
    {
        g_blend_why += ascii(c, SM_TBLOCK_SIZE) + "\n";
    }
    switch (e->GetErrorCode())
    {
    case SM_FILERR_NO_ROOM:
        return SMLIB_ERR_NO_ROOM;
    case SM_FILERR_NOT_EXACT:
        return SMLIB_ERR_NOT_EXACT;
    default:
        return s;
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
    mark_faces(b);
    SmFilletExecutive exec(*SmApiGetOrCreateContext(), b);
    int32_t s = status(exec.SetFilletParameters(b, es, kinds, solvers, r0, nullptr, -1.0, -1.0, -1.0, -1.0, 1, 1.0, &r1));
    if (s != SMLIB_OK)
    {
        return blend_failure(exec, s);
    }
    // SetFilletParameters chose the executive's options as SmApiFilletEdges and
    // SmApiSetBevelCorners run it (no global merge, which leaves a hole rim open).
    s = status(exec.CreateFilletCorners());
    if (s != SMLIB_OK)
    {
        return blend_failure(exec, s);
    }
    SmTArray<SmFilletCorner *> &corners = exec.GetFilletCorners();
    for (ULONG j = 0; j < corners.GetSize(); ++j)
    {
        corners[j]->SetBevel(TRUE);
    }
    s = status(exec.DoFilleting());
    if (s != SMLIB_OK)
    {
        return blend_failure(exec, s);
    }
    s = status(b->StitchAndOrient());
    // A chamfer that made no face did nothing, whatever SMLib's status said.
    int32_t made = made_faces(b);
    return s != SMLIB_OK ? s : made;
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
        if (e->GetEdgeuses(uses) != SM_SUCCESS)
        {
            return false;
        }
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
// `setbacks` holds corner vertices of `b` with the distance along each of their edges at which
// the corner's patch starts (SmFilletCorner's setback): each must be a corner where every edge
// is rounded, all of one convexity (SMLib's same-convexity N x N corners, the only ones that
// build a setback). With `patches`, the
// faces of each setback corner's patch in `b` afterwards (one list per setback, in order): the
// faces on the surfaces that corner made.
int32_t round_edges(SmBrep *b, const std::vector<SmEdge *> &edges, const std::vector<double> &sizes, int32_t kind,
                    const std::vector<SmFilletLaw *> &laws, const std::vector<std::pair<SmVertex *, double>> &setbacks = {},
                    std::vector<std::vector<SmFace *>> *patches = nullptr)
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
    mark_faces(b);
    SmFilletExecutive exec(*SmApiGetOrCreateContext(), b);
    // The curvature-continuous section is the quintic whose inner control points lie a third
    // and two thirds of the way from each contact to the corner its tangents make.
    exec.SetBlendCurveAtThirds(kind == SMLIB_BLEND_CURVATURE ? TRUE : FALSE);
    int32_t s = status(exec.SetFilletParameters(b, es, sections, solvers, values, laws.empty() ? nullptr : &law_array, -1.0, -1.0,
                                                -1.0, -1.0, continuity, 1.0));
    if (s != SMLIB_OK)
    {
        return blend_failure(exec, s);
    }
    // SetFilletParameters chose the executive's options as SmApiFilletEdges runs it (its
    // SetDoGlobalMerge comes before, so SetFilletParameters resets it). DoFillet's steps, with
    // each setback set on its corner between making the corners and filleting.
    s = status(exec.CreateFilletCorners());
    if (s != SMLIB_OK)
    {
        return blend_failure(exec, s);
    }
    std::vector<SmFilletCorner *> set_back;
    for (const auto &sb : setbacks)
    {
        SmFilletCorner *corner = exec.GetFilletCornerOfVertex(sb.first);
        if (corner == nullptr)
        {
            return SMLIB_ERR_ARG;
        }
        // SMLib sets back an N x N corner whose edges all have one convexity (its
        // SmFilletConvexNxNCorner, convex or concave alike). A mixed one (SmFilletConcaveNxNCorner)
        // hands its setback to that same topology, which takes the rails of neighbouring rounds to
        // meet within their edges; at a mixed corner those on the face between the two rounds of
        // one convexity do not, and its sections there come out kinked (a coincident pair of
        // control points) and the result no solid.
        SmFilletCornerType type = corner->GetCornerType();
        if (type != SM_FCR_N_x_N_CONVEX)
        {
            return SMLIB_ERR_UNSUPPORTED;
        }
        corner->SetSetBackDist(sb.second);
        set_back.push_back(corner);
    }
    exec.SetMarkCornerPatches(patches != nullptr ? TRUE : FALSE);
    s = status(exec.DoFilleting());
    if (s == SMLIB_OK)
    {
        // A fillet surface turns where its radius has a corner: its face split there.
        s = status(exec.SplitFilletFacesAtRadiusCorners());
    }
    if (s != SMLIB_OK)
    {
        return blend_failure(exec, s);
    }
    s = status(b->StitchAndOrient());
    // A round that made no face did nothing, whatever SMLib's status said.
    int32_t made = made_faces(b);
    if (s != SMLIB_OK || made != SMLIB_OK)
    {
        return s != SMLIB_OK ? s : made;
    }
    if (patches != nullptr)
    {
        // The merge gives `b` copies of the faces filleting made, with their attributes: a
        // corner's patch is the faces marked with its index (SetMarkCornerPatches).
        SmTArray<SmFace *> faces;
        b->GetFaces(faces);
        for (SmFilletCorner *corner : set_back)
        {
            long index = static_cast<long>(corner->FindIndexInExec());
            std::vector<SmFace *> patch;
            for (ULONG i = 0; i < faces.GetSize(); ++i)
            {
                SmAttribute *mark = faces[i]->FindAttribute(SM_AI_FILLET_CORNER_FACE);
                if (mark != nullptr && static_cast<SmTagAttribute *>(mark)->m_lPrimaryID == index)
                {
                    patch.push_back(faces[i]);
                }
            }
            if (patch.empty())
            {
                return SMLIB_ERR_UNSUPPORTED;
            }
            patches->push_back(std::move(patch));
        }
        s = status(b->RemoveAttributeFromTopology(SM_AI_FILLET_CORNER_FACE));
        if (s != SMLIB_OK)
        {
            return s;
        }
    }
    return SMLIB_OK;
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

// Why the last face move in this call was refused, as SmOffsetFacesExecutive noted it: dumped
// with the call's inputs (Call::input_text).
std::string g_offset_why;

// Moves the picked faces of `b` (in place) `distance` along their outward normals: each gets
// its surface's offset, and SMLib's local operation extends or trims the faces around them to
// meet it; the tangent blends on them follow (SmOffsetFacesExecutive: taken off, the faces
// moved, put back at their own sizes). A blend that no longer fits is SMLIB_ERR_NO_ROOM, one
// SMLib cannot put back as it was SMLIB_ERR_UNSUPPORTED, one that would need a continuation
// SMLib cannot make exactly SMLIB_ERR_NOT_EXACT, a move that moved nothing SMLIB_ERR_NO_EFFECT.
int32_t move_faces(SmBrep *b, const std::vector<SmFace *> &picked, double distance)
{
    SmTArray<SmFace *> faces;
    for (SmFace *f : picked)
    {
        faces.Add(f);
    }
    SmOffsetFacesExecutive exec(*SmApiGetOrCreateContext(), b);
    int32_t s = status(exec.OffsetFaces(faces, distance));
    g_offset_why = ascii(exec.GetReason(), SM_TBLOCK_SIZE);
    if (s == SMLIB_OK)
    {
        return s;
    }
    switch (exec.GetRefusal())
    {
    case SM_OFR_NO_ROOM:
        return SMLIB_ERR_NO_ROOM;
    case SM_OFR_UNSUPPORTED:
        return SMLIB_ERR_UNSUPPORTED;
    case SM_OFR_NOT_EXACT:
        return SMLIB_ERR_NOT_EXACT;
    case SM_OFR_NO_EFFECT:
        return SMLIB_ERR_NO_EFFECT;
    default:
        return s;
    }
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
    for (ULONG k = 0; k < surfaces.GetSize(); ++k)
    {
        Before::mark(surfaces[k]);
    }
    return moved_faces(b, status(b->LocalOperation(faces, surfaces)));
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

extern "C" int32_t smlib_boolean2d(const smlib_brep *a, const smlib_brep *b, int32_t op, smlib_brep **out)
{
    Sm2DBooleanOperationType kind;
    switch (op)
    {
    case SMLIB_BOOLEAN_UNION:
        kind = SM_2D_UNION;
        break;
    case SMLIB_BOOLEAN_SUBTRACT:
        kind = SM_2D_DIFFERENCE;
        break;
    case SMLIB_BOOLEAN_INTERSECT:
        kind = SM_2D_INTERSECTION;
        break;
    default:
        return SMLIB_ERR_ARG;
    }
    if (!valid(a) || !valid(b) || out == nullptr)
    {
        return SMLIB_ERR_ARG;
    }
    Call c("smlib_boolean2d", handle(a) + ", " + handle(b) + ", op=" + std::to_string(op));
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
            // SmApiBoolean2d returns success without a result when an input has no faces.
            if (ca->GetNumFaces() == 0 || cb->GetNumFaces() == 0)
            {
                return SMLIB_ERR_EMPTY;
            }
            SmBrep *r = nullptr;
            s = status(SmApiBoolean2d(ca.get(), cb.get(), kind, r));
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
                return finish(Owned(r), false, out);
            }
            return finish(std::move(ca), false, out);
        }
    );
}

extern "C" int32_t smlib_boolean_union_all(const smlib_brep *const *breps, int64_t n, smlib_brep **out)
{
    if (breps == nullptr || out == nullptr || n < 1 || n > SMLIB_BOOLEAN_MAX_BREPS)
    {
        return SMLIB_ERR_ARG;
    }
    for (int64_t i = 0; i < n; ++i)
    {
        if (!valid(breps[i]))
        {
            return SMLIB_ERR_ARG;
        }
    }
    Call c("smlib_boolean_union_all", "n=" + std::to_string(n));
    for (int64_t i = 0; i < n && i < 2; ++i)
    {
        c.input(i == 0 ? "first" : "second", breps[i]);
    }
    return c.run(
        [&]
        {
            // Copies, as smlib_boolean's: SM_API edits the first operand in place and deletes the
            // second.
            std::vector<Owned> level;
            level.reserve(static_cast<size_t>(n));
            bool solid_in = true;
            for (int64_t i = 0; i < n; ++i)
            {
                int32_t s = SMLIB_OK;
                Owned copy = copy_of(breps[i], s);
                if (s != SMLIB_OK)
                {
                    return s;
                }
                solid_in = solid_in && copy->IsManifoldSolid();
                level.push_back(std::move(copy));
            }
            // Pairwise, level by level: the bodies of a level are unioned two by two, an odd one
            // carried to the next level as it is.
            while (level.size() > 1)
            {
                std::vector<Owned> next;
                next.reserve((level.size() + 1) / 2);
                for (size_t i = 0; i + 1 < level.size(); i += 2)
                {
                    Owned a = std::move(level[i]);
                    Owned b = std::move(level[i + 1]);
                    SmBrep *r = nullptr;
                    int32_t s = status(SmApiBoolean(a.get(), b.get(), SM_BO_UNION, r));
                    if (s != SMLIB_OK)
                    {
                        // Both operands are still allocated (possibly modified); the owners free them.
                        return s;
                    }
                    // Success: b is deleted and the result is a, edited in place.
                    b.release();
                    if (r != nullptr && r != a.get())
                    {
                        // Not what SM_API documents; SMLib owns a now, so leave it be.
                        a.release();
                        a.reset(r);
                    }
                    next.push_back(std::move(a));
                }
                if (level.size() % 2 == 1)
                {
                    next.push_back(std::move(level.back()));
                }
                level = std::move(next);
            }
            return finish(std::move(level.front()), solid_in, out, true);
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

// Whether each setback fits on the edges of its corner: shorter than every edge there, and
// where an edge's other end is a setback corner too, the two together shorter than the edge.
bool setbacks_fit(const std::vector<std::pair<SmVertex *, double>> &setbacks)
{
    for (const auto &sb : setbacks)
    {
        SmTArray<SmEdge *> edges;
        sb.first->GetEdges(edges);
        for (ULONG k = 0; k < edges.GetSize(); ++k)
        {
            double length = 0.0;
            if (SmApiEdgeComputeLength(edges[k], 1e-6, length) != SM_SUCCESS)
            {
                return false;
            }
            double other = 0.0;
            SmVertex *far = edges[k]->GetOtherVertex(sb.first);
            for (const auto &o : setbacks)
            {
                if (o.first == far && far != sb.first)
                {
                    other = o.second;
                }
            }
            if (sb.second + other >= length)
            {
                return false;
            }
        }
    }
    return true;
}

// One edge round an opened corner's hole: the face beside it, and whether that face's boundary,
// run counter-clockwise about its outward normal (the solid's), runs the edge from its start
// vertex to its end.
struct HoleEdge
{
    SmEdge *edge;
    SmFace *face;
    bool forward;
};

// The faces of each patch in `patches` removed from the solid `b`, with the edges and vertices
// only they used: `holes[j]` the edges round patch j's hole, as the faces beside it run them,
// each edge's start the previous one's end. Read on the solid before anything is removed, where
// each face's outward side is known. Each edge round a hole must lie between a face of its patch
// and one other face, outside every patch, and the edges must close one loop (a patch that is a
// disc), else SMLIB_ERR_UNSUPPORTED.
int32_t open_corners(SmBrep *b, const std::vector<std::vector<SmFace *>> &patches, std::vector<std::vector<HoleEdge>> &holes)
{
    std::unordered_set<SmFace *> any;
    for (const auto &patch : patches)
    {
        any.insert(patch.begin(), patch.end());
    }
    for (const auto &patch : patches)
    {
        std::unordered_set<SmFace *> mine(patch.begin(), patch.end());
        std::unordered_set<SmEdge *> seen;
        std::vector<HoleEdge> rim;
        for (SmFace *f : patch)
        {
            SmTArray<SmEdge *> edges;
            f->GetEdges(edges);
            for (ULONG k = 0; k < edges.GetSize(); ++k)
            {
                SmEdge *e = edges[k];
                if (!seen.insert(e).second)
                {
                    continue;
                }
                SmTArray<SmFace *> beside;
                e->GetFaces(beside);
                SmFace *other = nullptr;
                ULONG others = 0;
                for (ULONG i = 0; i < beside.GetSize(); ++i)
                {
                    if (mine.count(beside[i]) == 0)
                    {
                        other = beside[i];
                        ++others;
                    }
                }
                if (others == 0)
                {
                    continue;
                }
                if (others != 1 || beside.GetSize() != 2 || any.count(other) != 0)
                {
                    return SMLIB_ERR_UNSUPPORTED;
                }
                SmFaceuse *fu = outer_faceuse(b, other);
                SmTArray<SmEdgeuse *> uses;
                if (fu == nullptr || e->GetEdgeuses(uses) != SM_SUCCESS)
                {
                    return SMLIB_ERR_UNSUPPORTED;
                }
                SmEdgeuse *use = nullptr;
                for (ULONG i = 0; i < uses.GetSize(); ++i)
                {
                    if (uses[i]->GetFaceuse() == fu)
                    {
                        use = uses[i];
                    }
                }
                if (use == nullptr)
                {
                    return SMLIB_ERR_UNSUPPORTED;
                }
                rim.push_back({e, other, use->GetOrientation() == SM_OT_SAME});
            }
        }
        auto start = [](const HoleEdge &h) { return h.forward ? h.edge->GetStartVertex() : h.edge->GetEndVertex(); };
        auto end = [](const HoleEdge &h) { return h.forward ? h.edge->GetEndVertex() : h.edge->GetStartVertex(); };
        if (rim.empty())
        {
            return SMLIB_ERR_UNSUPPORTED;
        }
        std::vector<HoleEdge> loop{rim.front()};
        std::vector<bool> used(rim.size(), false);
        used[0] = true;
        while (end(loop.back()) != start(loop.front()))
        {
            size_t next = rim.size();
            for (size_t i = 0; i < rim.size(); ++i)
            {
                if (!used[i] && start(rim[i]) == end(loop.back()))
                {
                    if (next != rim.size())
                    {
                        return SMLIB_ERR_UNSUPPORTED;
                    }
                    next = i;
                }
            }
            if (next == rim.size())
            {
                return SMLIB_ERR_UNSUPPORTED;
            }
            used[next] = true;
            loop.push_back(rim[next]);
        }
        if (loop.size() != rim.size())
        {
            return SMLIB_ERR_UNSUPPORTED;
        }
        holes.push_back(std::move(loop));
    }
    // A face is removed only with the Brep open for editing (SmBrep::DeleteFace refuses
    // otherwise).
    SmTemporaryChangeValue<SmBoolean> editing(b->m_bEditingEnabled, TRUE);
    for (SmFace *f : any)
    {
        int32_t s = status(b->DeleteFace(f, TRUE));
        if (s != SMLIB_OK)
        {
            return s;
        }
    }
    return SMLIB_OK;
}

// Where smlib_blend_setback_open writes its holes (as smlib.h describes).
struct HoleOut
{
    int64_t *buf;
    int64_t cap;
    int64_t *out_len;
};

// smlib_blend's work, logged as `fn`: `sizes` holds one size per point, or one for all;
// `corners` (n_corners points) and `setbacks` (one each) the corner setbacks, as for
// smlib_blend_setback.
int32_t blend_call(const char *fn, const smlib_brep *brep, int32_t kind, const double *points, const double *sizes, bool one_size,
                   int64_t n_points, double tolerance, smlib_brep **out, int64_t *out_failed_point, const double *corners = nullptr,
                   const double *setbacks = nullptr, int64_t n_corners = 0, int64_t *out_failed_corner = nullptr,
                   const HoleOut *holes = nullptr)
{
    bool kind_ok = kind == SMLIB_BLEND_ROUND || kind == SMLIB_BLEND_CURVATURE || kind == SMLIB_BLEND_CHORD || kind == SMLIB_BLEND_CHAMFER;
    if (!valid(brep) || out == nullptr || !kind_ok || !points_ok(points, n_points, false) || !tolerance_ok(tolerance) || sizes == nullptr ||
        !points_ok(corners, n_corners, true) || (n_corners > 0 && setbacks == nullptr) ||
        (holes != nullptr && (holes->out_len == nullptr || holes->cap < 0 || (holes->cap > 0 && holes->buf == nullptr))))
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
    for (int64_t j = 0; j < n_corners; ++j)
    {
        if (!size_ok(setbacks[j]))
        {
            return SMLIB_ERR_ARG;
        }
    }
    // A chamfer's corners are bevelled, which builds no setback.
    if (n_corners > 0 && kind == SMLIB_BLEND_CHAMFER)
    {
        return SMLIB_ERR_UNSUPPORTED;
    }
    std::string args = handle(brep) + ", kind=" + std::to_string(kind) + ", points=" + pts(points, n_points) + ", sizes=" +
                       vec(sizes, static_cast<int>(n_sizes)) + ", tolerance=" + num(tolerance);
    if (n_corners > 0)
    {
        args += ", corners=" + pts(corners, n_corners) + ", setbacks=" + vec(setbacks, static_cast<int>(n_corners));
    }
    Call c(fn, args);
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
            // Each corner point picks a vertex, as a point picks an edge; a vertex picked twice
            // takes the first point's setback.
            std::vector<SmVertex *> at;
            std::vector<int64_t> first_corner;
            s = pick_vertices(w.get(), corners, n_corners, tolerance, at, out_failed_corner, &first_corner);
            if (s != SMLIB_OK)
            {
                return s;
            }
            std::vector<std::pair<SmVertex *, double>> back;
            for (size_t k = 0; k < at.size(); ++k)
            {
                back.emplace_back(at[k], setbacks[first_corner[k]]);
            }
            if (!setbacks_fit(back))
            {
                return SMLIB_ERR_NO_ROOM;
            }
            // The corners' places, read now: filleting deletes the vertices.
            std::vector<SmPoint3d> at_points;
            for (SmVertex *v : at)
            {
                at_points.push_back(v->GetPoint());
            }
            bool solid_in = w->IsManifoldSolid();
            std::vector<std::vector<SmFace *>> patches;
            s = kind == SMLIB_BLEND_CHAMFER ? chamfer_edges(w.get(), picked, each)
                                            : round_edges(w.get(), picked, each, kind, {}, back, holes == nullptr ? nullptr : &patches);
            if (s != SMLIB_OK || holes == nullptr)
            {
                return s != SMLIB_OK ? s : finish(std::move(w), solid_in, out);
            }
            // Each corner opened: its patch removed, the edges round its hole in the result's
            // orders, once per corner point (a vertex picked twice gives its hole twice).
            if (!w->IsManifoldSolid())
            {
                return SMLIB_ERR_NOT_SOLID;
            }
            std::vector<std::vector<HoleEdge>> rims;
            s = open_corners(w.get(), patches, rims);
            if (s != SMLIB_OK)
            {
                return s;
            }
            SmTArray<SmEdge *> all_edges;
            SmTArray<SmFace *> all_faces;
            w->GetEdges(all_edges);
            w->GetFaces(all_faces);
            std::unordered_map<SmEdge *, int64_t> edge_index;
            std::unordered_map<SmFace *, int64_t> face_index;
            for (ULONG k = 0; k < all_edges.GetSize(); ++k)
            {
                edge_index[all_edges[k]] = static_cast<int64_t>(k);
            }
            for (ULONG k = 0; k < all_faces.GetSize(); ++k)
            {
                face_index[all_faces[k]] = static_cast<int64_t>(k);
            }
            std::vector<int64_t> values;
            for (int64_t j = 0; j < n_corners; ++j)
            {
                SmPoint3d p(corners[3 * j], corners[3 * j + 1], corners[3 * j + 2]);
                size_t which = at_points.size();
                double best = HUGE_VAL;
                for (size_t k = 0; k < at_points.size(); ++k)
                {
                    double d = (at_points[k] - p).Length();
                    if (d < best)
                    {
                        best = d;
                        which = k;
                    }
                }
                if (which >= rims.size())
                {
                    return SMLIB_ERR_ARG;
                }
                values.push_back(static_cast<int64_t>(rims[which].size()));
                for (const HoleEdge &h : rims[which])
                {
                    auto e = edge_index.find(h.edge);
                    auto f = face_index.find(h.face);
                    if (e == edge_index.end() || f == face_index.end())
                    {
                        return SMLIB_ERR_UNSUPPORTED;
                    }
                    values.push_back(e->second);
                    values.push_back(f->second);
                    values.push_back(h.forward ? 1 : 0);
                }
            }
            *holes->out_len = static_cast<int64_t>(values.size());
            if (static_cast<int64_t>(values.size()) > holes->cap)
            {
                return SMLIB_ERR_BUFFER;
            }
            std::copy(values.begin(), values.end(), holes->buf);
            return finish(std::move(w), false, out);
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

// Whether the law's radius is zero somewhere inside (0, 1), where the round would pinch to a
// point (its ends may be zero). With control radii >= 0 it is zero at a parameter exactly when
// every control radius whose basis function is nonzero there is: radii i - degree .. i inside the
// span [knots[i], knots[i + 1]), and radii s + m - 1 - degree .. s - 1 at a knot first found at
// index s and repeated m times.
bool law_pinches(int32_t degree, const double *knots, const double *radii, int64_t n_radii)
{
    auto all_zero = [&](int64_t from, int64_t to) {
        for (int64_t j = std::max<int64_t>(from, 0); j <= to && j < n_radii; ++j)
        {
            if (radii[j] != 0.0)
            {
                return false;
            }
        }
        return from <= to;
    };
    int64_t n_knots = n_radii + degree + 1;
    for (int64_t i = degree; i + 1 < n_knots && i < n_radii; ++i)
    {
        if (knots[i] < knots[i + 1] && all_zero(i - degree, i))
        {
            return true;
        }
    }
    for (int64_t s = degree + 1; s < n_radii;)
    {
        int64_t m = 1;
        while (s + m < n_radii && knots[s + m] == knots[s])
        {
            ++m;
        }
        if (all_zero(s + m - 1 - degree, s - 1))
        {
            return true;
        }
        s += m;
    }
    return false;
}

// The radius law along `e`: the B-spline (radius as its x, its parameter as its y) over the edge's parameter interval,
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
    // The knots in the law's direction.
    std::vector<double> kv(static_cast<size_t>(n_knots));
    for (int64_t i = 0; i < n_knots; ++i)
    {
        kv[static_cast<size_t>(i)] = reverse ? 1.0 - knots[n_knots - 1 - i] : knots[i];
    }
    // Radius as x, and as y the parameter itself, as SMLib's own laws keep something monotonic
    // there: control points at the Greville abscissae (the mean of the `degree` knots after
    // each) make y(t) = t exactly. A law curve with y = 0 has no speed where the radius is
    // constant, and SMLib cannot extend it past the edge's ends.
    SmTArray<SmPoint3d> cps;
    for (int64_t i = 0; i < n_radii; ++i)
    {
        double g = 0.0;
        for (int64_t j = 1; j <= degree; ++j)
        {
            g += kv[static_cast<size_t>(i + j)];
        }
        cps.Add(SmPoint3d(radii[reverse ? n_radii - 1 - i : i], g / static_cast<double>(degree), 0.0));
    }
    // The knots as SMLib takes them: each distinct value once, with its multiplicity.
    SmTArray<double> unique;
    SmTArray<ULONG> mults;
    for (int64_t i = 0; i < n_knots; ++i)
    {
        double k = kv[static_cast<size_t>(i)];
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

extern "C" int32_t smlib_blend_setback(const smlib_brep *brep, int32_t kind, const double *points, const double *sizes, int64_t n_points,
                                       const double *corners, const double *setbacks, int64_t n_corners, double tolerance, smlib_brep **out,
                                       int64_t *out_failed_point, int64_t *out_failed_corner)
{
    if (n_corners < 1)
    {
        return SMLIB_ERR_ARG;
    }
    return blend_call("smlib_blend_setback", brep, kind, points, sizes, false, n_points, tolerance, out, out_failed_point, corners, setbacks,
                      n_corners, out_failed_corner);
}

extern "C" int32_t smlib_blend_setback_open(const smlib_brep *brep, int32_t kind, const double *points, const double *sizes, int64_t n_points,
                                            const double *corners, const double *setbacks, int64_t n_corners, double tolerance,
                                            smlib_brep **out, int64_t *holes, int64_t cap, int64_t *out_len, int64_t *out_failed_point,
                                            int64_t *out_failed_corner)
{
    if (n_corners < 1)
    {
        return SMLIB_ERR_ARG;
    }
    HoleOut where{holes, cap, out_len};
    return blend_call("smlib_blend_setback_open", brep, kind, points, sizes, false, n_points, tolerance, out, out_failed_point, corners,
                      setbacks, n_corners, out_failed_corner, &where);
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
    if (law_pinches(degree, knots, radii, n_radii))
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
    c.input_text("why", [](const std::string &path) {
        if (FILE *f = std::fopen(path.c_str(), "w"))
        {
            std::fputs(g_offset_why.c_str(), f);
            std::fclose(f);
        }
    });
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
            g_offset_why.clear();
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

extern "C" int32_t smlib_merge_faces(const smlib_brep *brep, const double *keep, int64_t n_keep, double tolerance, smlib_brep **out,
                                     int64_t *out_failed_point)
{
    if (!valid(brep) || out == nullptr || !points_ok(keep, n_keep, true) || !tolerance_ok(tolerance))
    {
        return SMLIB_ERR_ARG;
    }
    Call c("smlib_merge_faces", handle(brep) + ", keep=" + pts(keep, n_keep) + ", tolerance=" + num(tolerance));
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
            std::vector<SmEdge *> kept;
            if (n_keep > 0)
            {
                s = pick_edges(w.get(), keep, n_keep, tolerance, kept, out_failed_point);
                if (s != SMLIB_OK)
                {
                    return s;
                }
            }
            bool solid_in = w->IsManifoldSolid();
            SmTArray<SmEdge *> keep_edges;
            for (SmEdge *e : kept)
            {
                keep_edges.Add(e);
            }
            SmStatus r = w->RemoveTopologicalEdgesAndVertices(kept.empty() ? nullptr : &keep_edges, nullptr);
            if (r != SM_SUCCESS)
            {
                return static_cast<int32_t>(r);
            }
            return finish(std::move(w), solid_in, out);
        }
    );
}
