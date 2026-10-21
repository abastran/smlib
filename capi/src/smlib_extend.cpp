// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// smlib_extend_sheet: a one-face sheet bounded by its surface's parameter lines, continued past
// one of those sides. A B-spline surface is extended by SMLib (SmBSplineSurface::
// CreateExtendedSurface: its end spans continued, or a G1 or G2 run-off); an analytic surface
// keeps its own geometry and the face's parameter box grows on that side, where the distance
// is the same parameter step all along it (a plane, a cylinder or cone along its rulings or
// round its axis, a sphere or torus across a circle of constant speed). A plane continued at
// G1 or G2 is the plane itself.

#include "smlib_shim.h"

#include <SmApiGeneral.h>
#include <SmApiQueries.h>
#include <SmBSplineSurface.h>
#include <SmBrep.h>
#include <SmEdge.h>
#include <SmFace.h>
#include <SmFaceuse.h>
#include <SmLoop.h>
#include <SmPlane.h>

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

// The parameter step that moves `distance` across side `s` (outward), the same all along it;
// false where the speed across it varies along the side.
bool analytic_step(const SmSurface *surf, const SmExtent2d &dom, const Side &s, double distance, double &step)
{
    double first = 0.0;
    for (int i = 0; i <= SAMPLES; ++i)
    {
        SmVector3d p;
        SmVector3d du;
        SmVector3d dv;
        if (!eval(surf, on_side(dom, s, double(i) / SAMPLES), p, &du, &dv))
        {
            return false;
        }
        double speed = (s.u ? du : dv).Length();
        if (!(speed > 0.0) || !std::isfinite(speed))
        {
            return false;
        }
        if (i == 0)
        {
            first = speed;
        }
        else if (std::fabs(speed - first) > 1e-9 * first)
        {
            return false;
        }
    }
    step = distance / first;
    return std::isfinite(step) && step > 0.0;
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
            SmSurface *grown = nullptr;
            SmExtent2d box = dom;
            bool plane = surf->IsKindOf(SmPlane_TYPE);
            if (surf->IsKindOf(SmBSplineSurface_TYPE))
            {
                // SMLib extends the whole surface; the face must be all of it.
                SmExtent2d natural = surf->GetNaturalUVDomain();
                double ptol = 1e-9 * std::max(1.0, std::max(natural.GetMax().x - natural.GetMin().x, natural.GetMax().y - natural.GetMin().y));
                if (std::fabs(natural.GetMin().x - dom.GetMin().x) > ptol || std::fabs(natural.GetMin().y - dom.GetMin().y) > ptol ||
                    std::fabs(natural.GetMax().x - dom.GetMax().x) > ptol || std::fabs(natural.GetMax().y - dom.GetMax().y) > ptol)
                {
                    return SMLIB_ERR_UNSUPPORTED;
                }
                SmContinuityType ct = mode == SMLIB_EXTEND_TANGENT ? SM_CT_G1 : mode == SMLIB_EXTEND_CURVATURE ? SM_CT_G1_G2 : SM_CT_CINFINITY;
                int32_t s = status(static_cast<SmBSplineSurface *>(surf)->CreateExtendedSurface(ctx, grow->which, distance, ct, grown));
                if (s != SMLIB_OK || grown == nullptr)
                {
                    delete grown;
                    return s != SMLIB_OK ? s : SMLIB_ERR_EMPTY;
                }
                box = grown->GetNaturalUVDomain();
            }
            else
            {
                if (mode != SMLIB_EXTEND_NATURAL && !plane)
                {
                    return SMLIB_ERR_UNSUPPORTED;
                }
                double step = 0.0;
                if (!analytic_step(surf, dom, *grow, distance, step))
                {
                    return SMLIB_ERR_NOT_EXACT;
                }
                SmPoint2d lo = dom.GetMin();
                SmPoint2d hi = dom.GetMax();
                double &edge = grow->u ? (grow->high ? hi.x : lo.x) : (grow->high ? hi.y : lo.y);
                edge += grow->high ? step : -step;
                box.SetMinMax(lo, hi);
                // A parameter box may not run past the surface's own (a sphere's poles) nor wrap
                // round it (more than a full turn of a periodic direction).
                SmExtent2d natural = surf->GetNaturalUVDomain();
                double span = grow->u ? hi.x - lo.x : hi.y - lo.y;
                double nlo = grow->u ? natural.GetMin().x : natural.GetMin().y;
                double nhi = grow->u ? natural.GetMax().x : natural.GetMax().y;
                if (std::isfinite(nlo) && std::isfinite(nhi) && nhi - nlo < 1e12)
                {
                    bool periodic = surf->IsPeriodic(natural, grow->u ? SM_SP_U : SM_SP_V) == TRUE;
                    double ptol = 1e-12 * std::max(1.0, nhi - nlo);
                    double a = grow->u ? lo.x : lo.y;
                    double z = grow->u ? hi.x : hi.y;
                    if (periodic ? span > nhi - nlo + ptol : (a < nlo - ptol || z > nhi + ptol))
                    {
                        return SMLIB_ERR_NOT_EXACT;
                    }
                }
                if (surf->Copy(ctx, grown) != SM_SUCCESS || grown == nullptr)
                {
                    delete grown;
                    return SMLIB_ERR_EMPTY;
                }
            }
            std::unique_ptr<SmBrep> r(new (ctx) SmBrep());
            r->SetTolerance(b->GetTolerance());
            SmFace *nf = nullptr;
            int32_t s = status(r->CreateFaceFromSurface(grown, box, nf));
            if (s != SMLIB_OK || nf == nullptr)
            {
                return s != SMLIB_OK ? s : SMLIB_ERR_EMPTY;
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
