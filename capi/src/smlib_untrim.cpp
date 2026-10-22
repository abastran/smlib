// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// smlib_untrim_face: a face restored towards its surface. Picked inner loops (holes) are taken
// away with SmFace::DeleteInnerLoop on a copy where they are lamina; the outer boundary goes
// back to a parameter box of the surface (SmBrep::CreateFaceFromSurface, a face bounded by the
// box's own iso-parameter edges), and the holes kept are rebuilt in that face from their own
// edges' curves and intervals (SmBrep::MakeVertexLoop and MakeEdgeInFace, as
// SmBrep::CreateFacesWithCurves inserts loops: closing a loop splits off the face inside it,
// which is deleted). The new face and the rebuilt edges take the attributes of the face and
// edges they stand for (a copy's attributes, SmAttribute::MakeCopy), so the history tags carry
// through.

#include "smlib_shim.h"

#include <SmApiGeneral.h>
#include <SmApiQueries.h>
#include <SmAttribute.h>
#include <SmBSplineCurve.h>
#include <SmBSplineSurface.h>
#include <SmBrep.h>
#include <SmEdge.h>
#include <SmEdgeuse.h>
#include <SmFace.h>
#include <SmFaceuse.h>
#include <SmLoop.h>
#include <SmLoopuse.h>
#include <SmVertex.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

using namespace smlib_shim;

namespace
{

using Owned = std::unique_ptr<SmBrep>;

constexpr int64_t MAX_HOLE_POINTS = 10000;

double edge_distance(SmEdge *e, const SmPoint3d &p)
{
    SmPoint3d q;
    double t = 0.0;
    double d = HUGE_VAL;
    if (SmApiEdgeClosestPoint(e, p, q, t, d) != SM_SUCCESS || !std::isfinite(d))
    {
        return HUGE_VAL;
    }
    return d;
}

double face_distance(SmFace *f, const SmPoint3d &p)
{
    SmPoint3d q;
    double d = HUGE_VAL;
    if (SmApiFaceClosestPoint(f, p, q, d) != SM_SUCCESS || !std::isfinite(d))
    {
        return HUGE_VAL;
    }
    return d;
}

// The face of `b` nearest `p`, within `tol`.
SmFace *face_at(SmBrep *b, const SmPoint3d &p, double tol)
{
    SmTArray<SmFace *> faces;
    b->GetFaces(faces);
    SmFace *best = nullptr;
    double d0 = HUGE_VAL;
    for (ULONG i = 0; i < faces.GetSize(); ++i)
    {
        double d = face_distance(faces[i], p);
        if (d < d0)
        {
            d0 = d;
            best = faces[i];
        }
    }
    return d0 <= tol ? best : nullptr;
}

double loop_distance(SmLoop *l, const SmPoint3d &p)
{
    SmTArray<SmEdge *> edges;
    l->GetEdges(edges);
    double best = HUGE_VAL;
    for (ULONG i = 0; i < edges.GetSize(); ++i)
    {
        best = std::min(best, edge_distance(edges[i], p));
    }
    return best;
}

// `to` takes a copy of each attribute of `from`.
void copy_attributes(const SmAObject *from, SmAObject *to, const SmContext &ctx)
{
    SmTArray<SmAttribute *> attrs;
    from->GetAttributes(attrs);
    for (ULONG i = 0; i < attrs.GetSize(); ++i)
    {
        SmAttribute *a = attrs[i] != nullptr ? attrs[i]->MakeCopy(ctx) : nullptr;
        if (a != nullptr)
        {
            to->AddAttribute(a);
        }
    }
}

// Whether the surface is one of SMLib's analytic kinds (each is held as a B-spline subclass):
// it has no parameter box of its own, so a face goes back to the box of its own span.
bool analytic(const SmSurface *s)
{
    return s->GetType() != SmBSplineSurface_TYPE;
}

// The box of the face's outer loop in the surface's parameters: the union of the tight boxes of
// its edgeuses' UV trim curves (made where the face has none attached yet).
bool span_of(SmLoop *outer, SmExtent2d &box)
{
    SmLoopuse *a = nullptr;
    SmLoopuse *b = nullptr;
    outer->GetLoopuses(a, b);
    if (a == nullptr)
    {
        return false;
    }
    SmTArray<SmEdgeuse *> eus;
    a->GetEdgeuses(eus);
    bool any = false;
    double lo[2] = {HUGE_VAL, HUGE_VAL};
    double hi[2] = {-HUGE_VAL, -HUGE_VAL};
    for (ULONG i = 0; i < eus.GetSize(); ++i)
    {
        SmBSplineCurve *uv = nullptr;
        if (eus[i]->GetOrCreateUVTrimCurve(uv) != SM_SUCCESS || uv == nullptr)
        {
            return false;
        }
        SmExtent3d e;
        if (uv->CalculateTightBoundingBox(uv->GetNaturalInterval(), &e) != SM_SUCCESS)
        {
            return false;
        }
        lo[0] = std::min(lo[0], e.GetMin().x);
        lo[1] = std::min(lo[1], e.GetMin().y);
        hi[0] = std::max(hi[0], e.GetMax().x);
        hi[1] = std::max(hi[1], e.GetMax().y);
        any = true;
    }
    if (!any || !(hi[0] > lo[0]) || !(hi[1] > lo[1]) || !std::isfinite(lo[0] + lo[1] + hi[0] + hi[1]))
    {
        return false;
    }
    box.SetMinMax(SmPoint2d(lo[0], lo[1]), SmPoint2d(hi[0], hi[1]));
    return true;
}

// One edge of a hole to rebuild: its edge, how the loop runs along it, and its ends as the loop
// runs.
struct Piece
{
    SmEdge *edge;
    SmOrientType orient;
    SmPoint3d from;
    SmPoint3d to;
};

bool pieces_of(SmLoop *l, std::vector<Piece> &out)
{
    SmLoopuse *a = nullptr;
    SmLoopuse *b = nullptr;
    l->GetLoopuses(a, b);
    if (a == nullptr)
    {
        return false;
    }
    SmTArray<SmEdgeuse *> eus;
    a->GetEdgeuses(eus);
    for (ULONG i = 0; i < eus.GetSize(); ++i)
    {
        SmEdge *e = eus[i]->GetEdge();
        if (e == nullptr || e->GetCurve() == nullptr)
        {
            return false;
        }
        SmVertex *v0 = nullptr;
        SmVertex *v1 = nullptr;
        e->GetVertices(v0, v1);
        if (v0 == nullptr || v1 == nullptr)
        {
            return false;
        }
        bool same = eus[i]->GetOrientation() == SM_OT_SAME;
        SmPoint3d p0 = v0->GetPoint();
        SmPoint3d p1 = v1->GetPoint();
        out.push_back({e, same ? SM_OT_SAME : SM_OT_OPPOSITE, same ? p0 : p1, same ? p1 : p0});
    }
    return !out.empty();
}

// Whether `f` is bounded by one loop of exactly the edges `made`.
bool is_disc(SmFace *f, const std::vector<SmEdge *> &made)
{
    SmTArray<SmLoop *> loops;
    f->GetLoops(loops);
    if (loops.GetSize() != 1)
    {
        return false;
    }
    SmTArray<SmEdge *> edges;
    loops[0]->GetEdges(edges);
    if (edges.GetSize() != made.size())
    {
        return false;
    }
    for (ULONG i = 0; i < edges.GetSize(); ++i)
    {
        if (std::find(made.begin(), made.end(), edges[i]) == made.end())
        {
            return false;
        }
    }
    return true;
}

// The hole `pieces` rebuilt in `face` of `r`: vertices and edges made in the face, the face the
// closed loop splits off inside it deleted. `face` is the face that keeps the rest. Making edges
// in a face is an edit of the brep, which SMLib allows only in edit mode (SmBrep.h,
// m_bEditingEnabled; MakeEdgeInFace refuses otherwise).
int32_t rebuild_hole(SmBrep *r, SmFace *&face, const std::vector<Piece> &pieces, const SmContext &ctx)
{
    SmTemporaryChangeValue<SmBoolean> editing(r->m_bEditingEnabled, TRUE);
    double tol = static_cast<double>(r->GetTolerance());
    SmLoop *dummy = nullptr;
    SmVertex *first = nullptr;
    int32_t s = status(r->MakeVertexLoop(face, pieces.front().from, dummy, first));
    if (s != SMLIB_OK || first == nullptr)
    {
        return s != SMLIB_OK ? s : SMLIB_ERR_EMPTY;
    }
    SmVertex *at = first;
    std::vector<SmEdge *> made;
    SmFace *split = nullptr;
    for (size_t i = 0; i < pieces.size(); ++i)
    {
        const Piece &p = pieces[i];
        bool last = i + 1 == pieces.size();
        SmVertex *next = first;
        if (!last)
        {
            s = status(r->MakeVertexLoop(face, p.to, dummy, next));
            if (s != SMLIB_OK || next == nullptr)
            {
                return s != SMLIB_OK ? s : SMLIB_ERR_EMPTY;
            }
        }
        SmCurve *curve = nullptr;
        s = status(p.edge->GetCurve()->Copy(ctx, curve));
        if (s != SMLIB_OK || curve == nullptr)
        {
            delete curve;
            return s != SMLIB_OK ? s : SMLIB_ERR_EMPTY;
        }
        SmEdge *e = nullptr;
        SmLoop *nl = nullptr;
        SmFace *nf = nullptr;
        s = status(r->MakeEdgeInFace(face, at, next, curve, nullptr, p.edge->GetInterval(), p.orient, tol, e, nl, nf));
        if (s != SMLIB_OK || e == nullptr)
        {
            return s != SMLIB_OK ? s : SMLIB_ERR_EMPTY;
        }
        copy_attributes(p.edge, e, ctx);
        made.push_back(e);
        if (nf != nullptr)
        {
            if (!last || split != nullptr)
            {
                return SMLIB_ERR_EXCEPTION;
            }
            split = nf;
        }
        at = next;
    }
    if (split == nullptr)
    {
        return SMLIB_ERR_EXCEPTION;
    }
    // One of the two faces is the disc inside the hole: it goes, its edges stay on the other.
    SmFace *disc = is_disc(split, made) ? split : is_disc(face, made) ? face : nullptr;
    if (disc == nullptr)
    {
        return SMLIB_ERR_EXCEPTION;
    }
    SmFace *rest = disc == split ? face : split;
    s = status(r->DeleteFace(disc, FALSE));
    if (s != SMLIB_OK)
    {
        return s;
    }
    face = rest;
    return SMLIB_OK;
}

} // namespace

extern "C" int32_t smlib_untrim_face(const smlib_brep *body, const double *point, double tol, const double *holes, int64_t n_holes,
                                     int32_t outer, int32_t in_place, smlib_brep **out)
{
    if (!valid(body) || point == nullptr || out == nullptr || !finite(point, 3) || !std::isfinite(tol) || !(tol > 0.0) ||
        n_holes > MAX_HOLE_POINTS || (n_holes > 0 && (holes == nullptr || !finite(holes, static_cast<int>(3 * n_holes)))) ||
        (outer != 0 && outer != 1) || (in_place != 0 && in_place != 1) || (n_holes == 0 && outer == 0))
    {
        return SMLIB_ERR_ARG;
    }
    Call c("smlib_untrim_face", handle(body) + ", " + vec(point, 3) + ", " + num(tol) + ", holes=" + std::to_string(n_holes) +
                                    ", outer=" + std::to_string(outer) + ", in_place=" + std::to_string(in_place));
    c.input("body", body);
    return c.run(
        [&]
        {
            const SmContext &ctx = *SmApiGetOrCreateContext();
            SmPoint3d at(point[0], point[1], point[2]);
            // The face to untrim, in a body of its own (its hole edges lamina there) or in a copy
            // of the whole body.
            Owned w;
            SmFace *f = nullptr;
            if (in_place != 0)
            {
                SmBrep *r = nullptr;
                int32_t s = status(SmApiBrepCopy(body->p, r));
                w.reset(r);
                if (s != SMLIB_OK || r == nullptr)
                {
                    return s != SMLIB_OK ? s : SMLIB_ERR_EMPTY;
                }
                f = face_at(r, at, tol);
            }
            else
            {
                SmFace *src = face_at(body->p, at, tol);
                if (src == nullptr)
                {
                    return SMLIB_ERR_PICK_NONE;
                }
                w.reset(new (ctx) SmBrep());
                w->SetTolerance(body->p->GetTolerance());
                SmTArray<SmFace *> one;
                one.Add(src);
                SmTArray<SmFace *> copied;
                // Through SmMerge, which copies the topology with its attributes (history tags);
                // the plain copy keeps only the face's.
                int32_t s = status(body->p->CopyFaces(one, w.get(), &copied, TRUE));
                if (s != SMLIB_OK || copied.GetSize() != 1)
                {
                    return s != SMLIB_OK ? s : SMLIB_ERR_EMPTY;
                }
                f = copied[0];
            }
            if (f == nullptr)
            {
                return SMLIB_ERR_PICK_NONE;
            }
            SmTArray<SmLoop *> loops;
            f->GetLoops(loops);
            if (loops.GetSize() == 0)
            {
                return SMLIB_ERR_UNSUPPORTED;
            }
            // The holes taken away: all of them, or the one nearest each point.
            std::vector<SmLoop *> removed;
            if (n_holes < 0)
            {
                for (ULONG i = 1; i < loops.GetSize(); ++i)
                {
                    removed.push_back(loops[i]);
                }
            }
            for (int64_t k = 0; k < n_holes; ++k)
            {
                SmPoint3d p(holes[3 * k], holes[3 * k + 1], holes[3 * k + 2]);
                SmLoop *best = nullptr;
                double d0 = HUGE_VAL;
                for (ULONG i = 1; i < loops.GetSize(); ++i)
                {
                    double d = loop_distance(loops[i], p);
                    if (d < d0)
                    {
                        d0 = d;
                        best = loops[i];
                    }
                }
                if (best == nullptr || d0 > tol)
                {
                    return SMLIB_ERR_PICK_NONE;
                }
                if (std::find(removed.begin(), removed.end(), best) == removed.end())
                {
                    removed.push_back(best);
                }
            }
            if (outer == 0)
            {
                for (SmLoop *l : removed)
                {
                    // A hole whose edges another face of the body uses is no hole of a sheet.
                    if (f->DeleteInnerLoop(l) != SM_SUCCESS)
                    {
                        return SMLIB_ERR_UNSUPPORTED;
                    }
                }
                return take(w.release(), out);
            }
            // The outer boundary: the surface's parameter box (a free-form surface's own, an
            // analytic one's the span of the face).
            SmSurface *surf = f->GetSurface();
            if (surf == nullptr)
            {
                return SMLIB_ERR_UNSUPPORTED;
            }
            SmExtent2d box = surf->GetNaturalUVDomain();
            if (analytic(surf) && !span_of(loops[0], box))
            {
                return SMLIB_ERR_UNSUPPORTED;
            }
            std::vector<std::vector<Piece>> kept;
            for (ULONG i = 1; i < loops.GetSize(); ++i)
            {
                if (std::find(removed.begin(), removed.end(), loops[i]) != removed.end())
                {
                    continue;
                }
                std::vector<Piece> p;
                if (!pieces_of(loops[i], p))
                {
                    return SMLIB_ERR_UNSUPPORTED;
                }
                kept.push_back(std::move(p));
            }
            SmSurface *copy = nullptr;
            int32_t s = status(surf->Copy(ctx, copy));
            if (s != SMLIB_OK || copy == nullptr)
            {
                delete copy;
                return s != SMLIB_OK ? s : SMLIB_ERR_EMPTY;
            }
            Owned r;
            SmBrep *target = w.get();
            if (in_place == 0)
            {
                r.reset(new (ctx) SmBrep());
                r->SetTolerance(w->GetTolerance());
                target = r.get();
            }
            SmFace *nf = nullptr;
            s = status(target->CreateFaceFromSurface(copy, box, nf));
            if (s != SMLIB_OK || nf == nullptr)
            {
                return s != SMLIB_OK ? s : SMLIB_ERR_EMPTY;
            }
            for (const std::vector<Piece> &hole : kept)
            {
                s = rebuild_hole(target, nf, hole, ctx);
                if (s != SMLIB_OK)
                {
                    return s;
                }
            }
            copy_attributes(f, nf, ctx);
            // The face faces the way the old one did.
            SmFaceuse *was = f->GetUpwardFaceuse();
            SmFaceuse *now = nf->GetUpwardFaceuse();
            if (was != nullptr && now != nullptr && was->GetOrientation() != now->GetOrientation())
            {
                s = status(nf->FlipFaceOrientation());
                if (s != SMLIB_OK)
                {
                    return s;
                }
            }
            if (in_place != 0)
            {
                // The old face goes, with the edges and vertices only it used.
                s = status(w->DeleteFace(f, TRUE));
                if (s != SMLIB_OK)
                {
                    return s;
                }
                return take(w.release(), out);
            }
            return take(r.release(), out);
        }
    );
}
