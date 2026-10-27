// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// smlib_extend_sheet: a one-face sheet bounded by its surface's parameter lines, continued past
// one of those sides by SMLib's extension of the surface (SmBSplineSurface::CreateExtendedSurface:
// its end spans continued, or a G1 or G2 run-off). SMLib holds every surface, its analytic ones
// included, as a B-spline (SmPlane, SmCylinder, SmCone, SmSphere and SmTorus derive from
// SmBSplineSurface), so one path serves all: the extended B-spline continues a plane, a circle or
// any conic section on itself, and CopyAndAddAnalytics gives it back its analytic kind where it
// is one. The new face, and each of its sides, takes the attributes of the face and side it
// continues (history tags).

#include "smlib_shim.h"

#include <SmApiGeneral.h>
#include <SmApiQueries.h>
#include <SmAttribute.h>
#include <SmBSplineSurface.h>
#include <SmBrep.h>
#include <SmEdge.h>
#include <SmFace.h>
#include <SmFaceuse.h>
#include <SmLoop.h>

#include <cmath>
#include <memory>

using namespace smlib_shim;

namespace
{

struct Side
{
    SmSurfParamType which;
    bool u;    // the side is a u = const line (the face grows in u)
    bool high; // at the box's max
};

constexpr Side SIDES[4] = {
    {SM_SP_UMIN, true, false},
    {SM_SP_UMAX, true, true},
    {SM_SP_VMIN, false, false},
    {SM_SP_VMAX, false, true},
};

// The point of side `s` of `dom` at `t` in [0, 1] along it, in UV.
SmPoint2d on_side(const SmExtent2d &dom, const Side &s, double t)
{
    SmPoint2d lo = dom.GetMin();
    SmPoint2d hi = dom.GetMax();
    if (s.u)
    {
        return SmPoint2d(s.high ? hi.x : lo.x, lo.y + t * (hi.y - lo.y));
    }
    return SmPoint2d(lo.x + t * (hi.x - lo.x), s.high ? hi.y : lo.y);
}

bool eval(const SmSurface *surf, const SmPoint2d &uv, SmVector3d &p, SmVector3d *du, SmVector3d *dv)
{
    // 1d layout for (1, 1): [D, Dv, Du, Duv].
    SmVector3d d[4];
    if (surf->Evaluate(uv, 1, 1, TRUE, TRUE, TRUE, d) != SM_SUCCESS)
    {
        return false;
    }
    p = d[0];
    if (du != nullptr)
    {
        *du = d[2];
    }
    if (dv != nullptr)
    {
        *dv = d[1];
    }
    return std::isfinite(p.Length()) && std::isfinite(d[1].Length()) && std::isfinite(d[2].Length());
}

double edge_distance(SmEdge *e, const SmVector3d &p)
{
    SmPoint3d q;
    double t = 0.0;
    double d = HUGE_VAL;
    if (SmApiEdgeClosestPoint(e, SmPoint3d(p), q, t, d) != SM_SUCCESS || !std::isfinite(d))
    {
        return HUGE_VAL;
    }
    return d;
}

double nearest_edge(const SmTArray<SmEdge *> &edges, const SmVector3d &p, SmEdge **which)
{
    double best = HUGE_VAL;
    for (ULONG i = 0; i < edges.GetSize(); ++i)
    {
        double d = edge_distance(edges[i], p);
        if (d < best)
        {
            best = d;
            if (which != nullptr)
            {
                *which = edges[i];
            }
        }
    }
    return best;
}

constexpr int SAMPLES = 9;

// Whether side `s` of the face's box is a side of the face: its samples lie on the face's
// edges (false for a side that collapses to a point, a pole).
bool is_boundary(const SmSurface *surf, const SmExtent2d &dom, const Side &s, const SmTArray<SmEdge *> &edges, double tol, bool &collapsed)
{
    SmVector3d a;
    SmVector3d b;
    collapsed = false;
    if (!eval(surf, on_side(dom, s, 0.0), a, nullptr, nullptr) || !eval(surf, on_side(dom, s, 1.0), b, nullptr, nullptr))
    {
        return false;
    }
    SmVector3d m;
    if (!eval(surf, on_side(dom, s, 0.5), m, nullptr, nullptr))
    {
        return false;
    }
    if ((a - b).Length() <= tol && (a - m).Length() <= tol)
    {
        collapsed = true;
        return true;
    }
    for (int i = 0; i <= SAMPLES; ++i)
    {
        SmVector3d p;
        if (!eval(surf, on_side(dom, s, double(i) / SAMPLES), p, nullptr, nullptr) || nearest_edge(edges, p, nullptr) > tol)
        {
            return false;
        }
    }
    return true;
}

// The edge of `edges` that side `s` of `dom` runs along (nearest the side's middle), or null for
// a side collapsed to a point.
SmEdge *side_edge(const SmSurface *surf, const SmExtent2d &dom, const Side &s, const SmTArray<SmEdge *> &edges, double on)
{
    SmVector3d a;
    SmVector3d b;
    SmVector3d m;
    if (!eval(surf, on_side(dom, s, 0.0), a, nullptr, nullptr) || !eval(surf, on_side(dom, s, 1.0), b, nullptr, nullptr) ||
        !eval(surf, on_side(dom, s, 0.5), m, nullptr, nullptr) || ((a - b).Length() <= on && (a - m).Length() <= on))
    {
        return nullptr;
    }
    SmEdge *e = nullptr;
    return nearest_edge(edges, m, &e) <= on ? e : nullptr;
}

// `to` takes a copy of each attribute of `from` (as SMLib's copies pass them on).
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

} // namespace

extern "C" int32_t smlib_extend_sheet(const smlib_brep *sheet, const double *point, double tol, double distance, int32_t mode, smlib_brep **out)
{
    if (!valid(sheet) || point == nullptr || out == nullptr || !finite(point, 3) || !std::isfinite(tol) || !(tol > 0.0) ||
        !std::isfinite(distance) || !(distance > 0.0) || mode < SMLIB_EXTEND_NATURAL || mode > SMLIB_EXTEND_CURVATURE)
    {
        return SMLIB_ERR_ARG;
    }
    Call c("smlib_extend_sheet", handle(sheet) + ", " + vec(point, 3) + ", " + num(tol) + ", " + num(distance) + ", " + std::to_string(mode));
    c.input("sheet", sheet);
    return c.run(
        [&]
        {
            SmBrep *b = sheet->p;
            SmTArray<SmFace *> faces;
            b->GetFaces(faces);
            if (faces.GetSize() != 1)
            {
                return SMLIB_ERR_UNSUPPORTED;
            }
            SmFace *f = faces[0];
            SmSurface *surf = f->GetSurface();
            SmTArray<SmLoop *> loops;
            f->GetLoops(loops);
            if (surf == nullptr || loops.GetSize() != 1)
            {
                return SMLIB_ERR_UNSUPPORTED;
            }
            SmTArray<SmEdge *> edges;
            f->GetEdges(edges);
            SmEdge *picked = nullptr;
            SmVector3d at(point[0], point[1], point[2]);
            if (nearest_edge(edges, at, &picked) > tol || picked == nullptr)
            {
                return SMLIB_ERR_PICK_NONE;
            }
            double model_tol = std::max(static_cast<double>(b->GetTolerance()), 1e-9);
            double on = 10.0 * model_tol;
            SmExtent2d dom = f->GetUVDomain();
            // Every side of the box is a side of the face (the face is its whole box), and the
            // picked edge is all of one side.
            const Side *grow = nullptr;
            for (const Side &s : SIDES)
            {
                bool collapsed = false;
                if (!is_boundary(surf, dom, s, edges, on, collapsed))
                {
                    return SMLIB_ERR_UNSUPPORTED;
                }
                if (collapsed)
                {
                    continue;
                }
                bool all = true;
                for (int i = 0; i <= SAMPLES && all; ++i)
                {
                    SmVector3d p;
                    all = eval(surf, on_side(dom, s, double(i) / SAMPLES), p, nullptr, nullptr) && edge_distance(picked, p) <= on;
                }
                if (all)
                {
                    grow = &s;
                }
            }
            if (grow == nullptr)
            {
                return SMLIB_ERR_UNSUPPORTED;
            }
            const SmContext &ctx = *SmApiGetOrCreateContext();
            // SMLib extends the whole surface; the face must be all of it.
            SmExtent2d natural = surf->GetNaturalUVDomain();
            double ptol = 1e-9 * std::max(1.0, std::max(natural.GetMax().x - natural.GetMin().x, natural.GetMax().y - natural.GetMin().y));
            if (std::fabs(natural.GetMin().x - dom.GetMin().x) > ptol || std::fabs(natural.GetMin().y - dom.GetMin().y) > ptol ||
                std::fabs(natural.GetMax().x - dom.GetMax().x) > ptol || std::fabs(natural.GetMax().y - dom.GetMax().y) > ptol)
            {
                return SMLIB_ERR_UNSUPPORTED;
            }
            SmContinuityType ct = mode == SMLIB_EXTEND_TANGENT ? SM_CT_G1 : mode == SMLIB_EXTEND_CURVATURE ? SM_CT_G1_G2 : SM_CT_CINFINITY;
            SmSurface *extended = nullptr;
            int32_t s = status(static_cast<SmBSplineSurface *>(surf)->CreateExtendedSurface(ctx, grow->which, distance, ct, extended));
            std::unique_ptr<SmSurface> clean(extended);
            if (s != SMLIB_OK || extended == nullptr)
            {
                return s != SMLIB_OK ? s : SMLIB_ERR_EMPTY;
            }
            SmSurface *grown = nullptr;
            s = status(extended->CopyAndAddAnalytics(ctx, grown));
            if (s != SMLIB_OK || grown == nullptr)
            {
                delete grown;
                return s != SMLIB_OK ? s : SMLIB_ERR_EMPTY;
            }
            SmExtent2d box = grown->GetNaturalUVDomain();
            std::unique_ptr<SmBrep> r(new (ctx) SmBrep());
            r->SetTolerance(b->GetTolerance());
            SmFace *nf = nullptr;
            s = status(r->CreateFaceFromSurface(grown, box, nf));
            if (s != SMLIB_OK || nf == nullptr)
            {
                return s != SMLIB_OK ? s : SMLIB_ERR_EMPTY;
            }
            // The face and each side stand for the face and side they continue.
            copy_attributes(f, nf, ctx);
            SmTArray<SmEdge *> made;
            nf->GetEdges(made);
            double new_on = 10.0 * std::max(static_cast<double>(r->GetTolerance()), 1e-9);
            for (const Side &side : SIDES)
            {
                SmEdge *from = side_edge(surf, dom, side, edges, on);
                SmEdge *to = side_edge(grown, box, side, made, new_on);
                if (from != nullptr && to != nullptr)
                {
                    copy_attributes(from, to, ctx);
                }
            }
            // The new face faces the way the old one did.
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
            return take(r.release(), out);
        }
    );
}
