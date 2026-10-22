// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// Split Face and Split Body.
//
// smlib_split_faces cuts chosen faces where tool surfaces meet them, the way SMLib's own plane
// cutter (SmPlaneCutter::MergeIntersection, SmCutter.cpp) cuts a subset of faces: the face's
// surface is intersected with the tool over the face's UV domain (CreatePlanarSectionCurves for a
// plane, GlobalSurfaceIntersect otherwise), and the untrimmed curves go to
// SmBrep::MergeCurvesOnSurface with the tool as the intersecting surface, which trims them to the
// face, cuts the edges they end on and splits the face (SmBrep.h: "merge 3DCurves known to be on
// a surface into the surface->Face possibly creating new vertices, edges and splitting Face into
// faces"). Unlike the plane cutter, a failed merge fails the call.
//
// smlib_split_body slices a solid with a one-face sheet on the tool surface (SM_BO_SLICE, which
// keeps the solid on the side the sheet's normal points away from and closes it with the sheet's
// piece inside it), once with the sheet and once with it flipped.

#include "smlib_shim.h"

#include <SmApiBrep.h>
#include <SmApiQueries.h>
#include <SmExtent3d.h>
#include <SmCurve.h>
#include <SmFace.h>
#include <SmSurface.h>
#include <SmTol.h>

#include <memory>
#include <set>

using namespace smlib_shim;

namespace
{

using Owned = std::unique_ptr<SmBrep>;

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

// Owns the curves an intersection hands back (a curve may be listed twice).
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

std::string list(const int64_t *v, int64_t n)
{
    std::string s = "[";
    for (int64_t i = 0; i < n && i < 16; ++i)
    {
        s += (i ? "," : "") + std::to_string(v[i]);
    }
    return s + (n > 16 ? ",...]" : "]");
}

// Cuts `face` where `tool` meets it; the faces it was split into besides itself go to `made`.
int32_t cut_face(SmBrep &brep, SmFace *face, const smlib_surface &tool, double tol, SmTArray<SmFace *> &made)
{
    const SmContext &ctx = *brep.GetContext();
    SmSurface *surface = face->GetSurface();
    if (surface == nullptr)
    {
        return static_cast<int32_t>(SM_ERR);
    }
    const SmExtent2d domain = face->GetUVDomain();
    SmApproxTol3d approx = SM_XSECT_TO_APPROXTOL3D(tol);
    Curves cs;
    SmPoint3d at;
    SmVector3d normal;
    if (tool.s->IsKindOf(SmPlane_TYPE) && tool.s->IsPlanar(SM_EFF_ZERO, &at, &normal) == TRUE)
    {
        // As SmPlaneCutter: a face far from the plane is not sectioned.
        SmExtent3d box;
        int32_t s = status(face->CalculateBoundingBox(box));
        if (s != SMLIB_OK)
        {
            return s;
        }
        if (std::fabs(box.DistanceToPlane(at, normal)) > tol)
        {
            return SMLIB_OK;
        }
        s = status(surface->CreatePlanarSectionCurves(ctx, domain, at, normal, &approx, NULL, &cs.c3, &cs.uv1));
        if (s != SMLIB_OK)
        {
            return s;
        }
    }
    else
    {
        SmBoolean edges[2] = {TRUE, TRUE};
        SmTArray<SmTsectCurveType> types;
        int32_t s = status(surface->GlobalSurfaceIntersect(ctx, domain, *tool.s, tool.domain, edges, &approx, NULL, &cs.c3, &cs.uv1, &cs.uv2,
                                                           &types, NULL));
        if (s != SMLIB_OK)
        {
            return s;
        }
        // A single point of contact cuts nothing.
        SmTArray<SmCurve *> c3, uv1;
        for (ULONG i = 0; i < cs.c3.GetSize(); ++i)
        {
            if (i < types.GetSize() && types[i] == SM_TC_TOUCHING)
            {
                continue;
            }
            c3.Add(cs.c3[i]);
            uv1.Add(i < cs.uv1.GetSize() ? cs.uv1[i] : nullptr);
        }
        if (c3.GetSize() == 0)
        {
            return SMLIB_OK;
        }
        SmTArray<SmEdge *> on;
        const bool all_uv = [&] {
            for (ULONG i = 0; i < uv1.GetSize(); ++i)
            {
                if (uv1[i] == nullptr)
                {
                    return false;
                }
            }
            return true;
        }();
        return status(brep.MergeCurvesOnSurface(*surface, tol, c3, all_uv ? &uv1 : nullptr, made, on, FALSE, tool.s));
    }
    if (cs.c3.GetSize() == 0)
    {
        return SMLIB_OK;
    }
    SmTArray<SmEdge *> on;
    const bool uv = cs.uv1.GetSize() == cs.c3.GetSize();
    return status(brep.MergeCurvesOnSurface(*surface, tol, cs.c3, uv ? &cs.uv1 : nullptr, made, on, FALSE, tool.s));
}

} // namespace

extern "C" int32_t smlib_split_faces(const smlib_brep *body, const int64_t *faces, int64_t n_faces, const smlib_surface *const *tools,
                                     int64_t n_tools, smlib_brep **out, int64_t *made)
{
    if (!valid(body) || faces == nullptr || made == nullptr || out == nullptr || tools == nullptr || n_faces < 1 ||
        n_faces > SMLIB_SPLIT_MAX_FACES || n_tools < 1 || n_tools > SMLIB_SPLIT_MAX_TOOLS)
    {
        return SMLIB_ERR_ARG;
    }
    for (int64_t t = 0; t < n_tools; ++t)
    {
        if (tools[t] == nullptr || tools[t]->s == nullptr)
        {
            return SMLIB_ERR_ARG;
        }
    }
    Call c("smlib_split_faces", handle(body) + ", faces=" + list(faces, n_faces) + ", n_tools=" + std::to_string(n_tools));
    c.input("body", body);
    return c.run(
        [&]
        {
            const ULONG count = body->p->GetNumFaces();
            std::set<int64_t> seen;
            for (int64_t i = 0; i < n_faces; ++i)
            {
                if (faces[i] < 0 || faces[i] >= static_cast<int64_t>(count) || !seen.insert(faces[i]).second)
                {
                    return SMLIB_ERR_ARG;
                }
            }
            int32_t s = SMLIB_OK;
            Owned a = copy_of(body, s);
            if (s != SMLIB_OK)
            {
                return s;
            }
            SmTArray<SmFace *> all;
            a->GetFaces(all);
            if (all.GetSize() != count)
            {
                return static_cast<int32_t>(SM_ERR);
            }
            // Each chosen face and, as tools cut it, its pieces.
            std::vector<std::vector<SmFace *>> pieces(static_cast<size_t>(n_faces));
            for (int64_t i = 0; i < n_faces; ++i)
            {
                pieces[i].push_back(all[static_cast<ULONG>(faces[i])]);
                made[i] = 0;
            }
            SmTemporaryChangeValue<SmBoolean> editing(a->m_bEditingEnabled, TRUE);
            SmTemporaryChangeValue<SmBoolean> doing(SM_CONST_CAST(SmContext *, a->GetContext())->GetDoingBooleanRef(), TRUE);
            const double tol = a->GetTolerance();
            int64_t total = 0;
            for (int64_t t = 0; t < n_tools; ++t)
            {
                for (int64_t i = 0; i < n_faces; ++i)
                {
                    std::vector<SmFace *> now = pieces[i];
                    for (SmFace *f : now)
                    {
                        SmTArray<SmFace *> split;
                        s = cut_face(*a, f, *tools[t], tol, split);
                        if (s != SMLIB_OK)
                        {
                            return s;
                        }
                        for (ULONG k = 0; k < split.GetSize(); ++k)
                        {
                            pieces[i].push_back(split[k]);
                        }
                        made[i] += static_cast<int64_t>(split.GetSize());
                        total += static_cast<int64_t>(split.GetSize());
                    }
                }
            }
            if (total == 0)
            {
                return SMLIB_ERR_NO_CUT;
            }
            return take(a.release(), out);
        }
    );
}

extern "C" int32_t smlib_split_body(const smlib_brep *body, const smlib_surface *tool, smlib_brep **back, smlib_brep **front)
{
    if (!valid(body) || tool == nullptr || tool->s == nullptr || back == nullptr || front == nullptr)
    {
        return SMLIB_ERR_ARG;
    }
    Call c("smlib_split_body", handle(body));
    c.input("body", body);
    return c.run(
        [&]
        {
            if (body->p->IsManifoldSolid() != TRUE)
            {
                return static_cast<int32_t>(SMLIB_ERR_NOT_SOLID);
            }
            const SmContext &ctx = *SmApiGetOrCreateContext();
            // One side: the body sliced by the tool's sheet, its normal turned away from it.
            auto side = [&](bool flip, Owned &part) -> int32_t
            {
                int32_t s = SMLIB_OK;
                Owned a = copy_of(body, s);
                if (s != SMLIB_OK)
                {
                    return s;
                }
                SmSurface *surface = nullptr;
                if (tool->s->Copy(ctx, surface) != SM_SUCCESS || surface == nullptr)
                {
                    delete surface;
                    return static_cast<int32_t>(SM_ERR);
                }
                // The other side's sheet is the surface turned over (its normal reversed);
                // SMLib's slice reads the side from the surface, not the face's use.
                if (flip)
                {
                    int32_t st = status(surface->SwapUV());
                    if (st != SMLIB_OK)
                    {
                        delete surface;
                        return st;
                    }
                }
                Owned sheet(new (ctx) SmBrep());
                sheet->SetTolerance(a->GetTolerance());
                SmFace *face = nullptr;
                s = status(sheet->CreateFaceFromSurface(surface, surface->GetNaturalUVDomain(), face));
                if (s != SMLIB_OK || face == nullptr)
                {
                    return s != SMLIB_OK ? s : static_cast<int32_t>(SMLIB_ERR_EMPTY);
                }
                SmBrep *r = nullptr;
                // The sheet is consumed by the slice.
                s = status(SmApiBooleanWithOptions(a.get(), sheet.get(), SM_BO_SLICE, FALSE, FALSE, FALSE, FALSE, FALSE, r));
                if (s != SMLIB_OK)
                {
                    return s;
                }
                sheet.release();
                if (r != nullptr && r != a.get())
                {
                    a.release();
                    a.reset(r);
                }
                if (a == nullptr || a->GetNumFaces() == 0)
                {
                    part.reset();
                    return SMLIB_OK;
                }
                part = std::move(a);
                return SMLIB_OK;
            };
            Owned b, f;
            int32_t s = side(true, b);
            if (s != SMLIB_OK)
            {
                return s;
            }
            s = side(false, f);
            if (s != SMLIB_OK)
            {
                return s;
            }
            if (b == nullptr || f == nullptr)
            {
                return SMLIB_ERR_NO_CUT;
            }
            *back = new smlib_brep{b.release()};
            *front = new smlib_brep{f.release()};
            return SMLIB_OK;
        }
    );
}
