// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// smlib_brep_mass_properties and smlib_brep_area: SMLib's exact (integrated, not tessellated)
// area, volume, centroid and inertia; smlib_brep_bounds: the exact box.

#include "smlib_shim.h"

#include <SmAttribute.h>
#include <SmExtent3d.h>
#include <SmFace.h>
#include <SmRegion.h>
#include <SmVertex.h>

using namespace smlib_shim;

namespace
{
// Every face counted once (SmApiBrepComputeArea, with the [1e-8, 1e-1] envelope instead of its
// [1e-4, 1e-1]): the total area and, when `first` is given, the area's first moments about `o`.
int32_t faces_area(const SmBrep *b, double accuracy, const SmPoint3d &o, double &area, SmVector3d *first)
{
    area = 0.0;
    SmVector3d sum(0.0, 0.0, 0.0);
    SmTArray<SmFace *> faces;
    b->GetFaces(faces);
    const ULONG flags = first != nullptr ? (SM_PPF_AREA | SM_PPF_AREA_FIRST_MOMENTS) : SM_PPF_AREA;
    for (ULONG i = 0; i < faces.GetSize(); ++i)
    {
        const SmFace *f = faces[i];
        SmExtent3d box;
        int32_t s = status(f->CalculateBoundingBox(box));
        if (s != SMLIB_OK)
        {
            return s;
        }
        const SmVector3d size = box.GetSize();
        const double edge = std::max(size.z, std::max(size.x, size.y));
        double a = 0.0, volume = 0.0;
        SmTArray<SmVector3d> moments;
        s = status(f->ComputePreciseProperties(SM_OT_SAME, accuracy, o, edge * edge / 100.0, 0.0, a, volume, moments, flags));
        if (s != SMLIB_OK)
        {
            return s;
        }
        // The area and its moments share one integrand, so a trim boundary integrated the other
        // way round flips both.
        const double sign = a < 0.0 ? -1.0 : 1.0;
        area += sign * a;
        if (first != nullptr && moments.GetSize() == 8)
        {
            sum = sum + moments[0] * sign;
        }
    }
    if (first != nullptr)
    {
        *first = sum;
    }
    return std::isfinite(area) ? SMLIB_OK : static_cast<int32_t>(SM_ERR_INVALID_INPUT);
}
} // namespace

extern "C" int32_t smlib_brep_mass_properties(const smlib_brep *brep, double relative_accuracy, double out[11])
{
    if (!valid(brep) || out == nullptr || !std::isfinite(relative_accuracy))
    {
        return SMLIB_ERR_ARG;
    }
    Call c("smlib_brep_mass_properties", handle(brep) + ", relative_accuracy=" + num(relative_accuracy));
    c.input("brep", brep);
    return c.run(
        [&]
        {
            const SmBrep *b = brep->p;
            // SmBrep::ComputePreciseProperties is called directly: SmApiBrepComputeMassProperties
            // clamps the accuracy to [1e-4, 1e-1], which leaves curved faces off by 1e-5 to 1e-3,
            // and refuses any body that is not IsManifoldSolid. The kernel routine itself takes
            // the [1e-8, 1e-1] SmApiBrepComputeVolume allows and is defined by regions: its volume
            // comes from "faces which bound non void regions" and its area from "faces which lie
            // between a solid region and a void region" (SmMassPropertyBehavior, SmTopoTypes.h).
            // So a body of several material regions, touching along edges or points
            // (non-manifold), measures as the sum over those regions; none is refused for it.
            // A body with no material region (a sheet, an open shell) has no volume: it is
            // measured face by face, as SmApiBrepComputeArea does, with its area centroid and no
            // inertia.
            bool material = false;
            {
                SmTArray<SmRegion *> regions;
                b->GetRegions(regions);
                for (ULONG i = 0; i < regions.GetSize() && !material; ++i)
                {
                    material = regions[i] != nullptr && !regions[i]->IsVoid();
                }
            }
            // A per-entity mass attribute would override the unit density.
            ULONG attr = SM_AI_MASS_PROPERTIES;
            SmTArray<SmFace *> faces;
            SmTArray<SmEdge *> edges;
            SmTArray<SmVertex *> vertices;
            b->GetFaces(faces, &attr);
            b->GetEdges(edges, &attr);
            b->GetVertices(vertices, &attr);
            if (faces.GetSize() != 0 || edges.GetSize() != 0 || vertices.GetSize() != 0)
            {
                return static_cast<int32_t>(SM_ERR_INVALID_INPUT);
            }
            const double accuracy = std::min(std::max(relative_accuracy, 1e-8), 1e-1);
            // Moments are taken about the middle of the body's box and moved to the centroid
            // below: about the world origin, a body far from it would lose its inertia to
            // cancellation.
            SmExtent3d box;
            int32_t s = status(b->CalculateBoundingBox(box));
            if (s != SMLIB_OK)
            {
                return s;
            }
            const SmPoint3d o = box.GetMid();
            if (!std::isfinite(o.x) || !std::isfinite(o.y) || !std::isfinite(o.z))
            {
                return static_cast<int32_t>(SM_ERR_INVALID_INPUT);
            }
            if (!material)
            {
                double area = 0.0;
                SmVector3d first(0.0, 0.0, 0.0);
                s = faces_area(b, accuracy, o, area, &first);
                if (s != SMLIB_OK)
                {
                    return s;
                }
                if (!(area > 0.0))
                {
                    return static_cast<int32_t>(SM_ERR_INVALID_INPUT);
                }
                const double r[11] = {area, 0.0, o.x + first.x / area, o.y + first.y / area, o.z + first.z / area, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
                if (!finite(r, 11))
                {
                    return static_cast<int32_t>(SM_ERR_INVALID_INPUT);
                }
                std::memcpy(out, r, sizeof r);
                return SMLIB_OK;
            }
            double area = 0.0, volume = 0.0, mass = 0.0;
            SmTArray<SmVector3d> m;
            m.SetSize(8);
            s = status(b->ComputePreciseProperties(accuracy, o, SmMassPropertyBehavior(1.0), area, volume, mass, m));
            if (s != SMLIB_OK)
            {
                return s;
            }
            if (m.GetSize() != 8 || !std::isfinite(mass) || !(mass > 0.0) || !(volume > 0.0))
            {
                return static_cast<int32_t>(SM_ERR_INVALID_INPUT);
            }
            // m[4]: first moments about o; m[7]: (y^2+z^2, z^2+x^2, x^2+y^2) about o; m[6]: (yz, zx, xy)
            // about o. Parallel axes from o to the centroid (unit density: mass = volume).
            const double d[3] = {m[4].x / mass, m[4].y / mass, m[4].z / mass};
            const double r[11] = {
                area,
                volume,
                o.x + d[0],
                o.y + d[1],
                o.z + d[2],
                m[7].x - mass * (d[1] * d[1] + d[2] * d[2]),
                m[7].y - mass * (d[0] * d[0] + d[2] * d[2]),
                m[7].z - mass * (d[0] * d[0] + d[1] * d[1]),
                m[6].x - mass * d[1] * d[2],
                m[6].y - mass * d[2] * d[0],
                m[6].z - mass * d[0] * d[1],
            };
            if (!finite(r, 11))
            {
                return static_cast<int32_t>(SM_ERR_INVALID_INPUT);
            }
            std::memcpy(out, r, sizeof r);
            return SMLIB_OK;
        }
    );
}

extern "C" int32_t smlib_brep_area(const smlib_brep *brep, double relative_accuracy, double *out)
{
    if (!valid(brep) || out == nullptr || !std::isfinite(relative_accuracy))
    {
        return SMLIB_ERR_ARG;
    }
    Call c("smlib_brep_area", handle(brep) + ", relative_accuracy=" + num(relative_accuracy));
    c.input("brep", brep);
    return c.run(
        [&]
        {
            const double accuracy = std::min(std::max(relative_accuracy, 1e-8), 1e-1);
            double total = 0.0;
            const int32_t s = faces_area(brep->p, accuracy, SmPoint3d(0.0, 0.0, 0.0), total, nullptr);
            if (s != SMLIB_OK)
            {
                return s;
            }
            *out = total;
            return SMLIB_OK;
        }
    );
}

extern "C" int32_t smlib_brep_bounds(const smlib_brep *brep, double lo[3], double hi[3])
{
    if (!valid(brep) || lo == nullptr || hi == nullptr)
    {
        return SMLIB_ERR_ARG;
    }
    Call c("smlib_brep_bounds", handle(brep));
    c.input("brep", brep);
    return c.run(
        [&]
        {
            // The tight box: the extrema of every vertex, edge and face on their exact geometry,
            // not the control hulls or tolerance-grown boxes CalculateBoundingBox unions.
            SmExtent3d box;
            int32_t s = status(brep->p->CalculateTightBoundingBox(box));
            if (s != SMLIB_OK)
            {
                return s;
            }
            const double r[6] = {box.GetUMin(), box.GetVMin(), box.GetWMin(), box.GetUMax(), box.GetVMax(), box.GetWMax()};
            if (!finite(r, 6) || r[0] > r[3] || r[1] > r[4] || r[2] > r[5])
            {
                return static_cast<int32_t>(SM_ERR_INVALID_INPUT);
            }
            std::memcpy(lo, r, 3 * sizeof(double));
            std::memcpy(hi, r + 3, 3 * sizeof(double));
            return SMLIB_OK;
        }
    );
}
