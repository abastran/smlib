// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// smlib_brep_mass_properties and smlib_brep_area: SMLib's exact (integrated, not tessellated)
// area, volume, centroid and inertia.

#include "smlib_shim.h"

#include <SmAttribute.h>
#include <SmExtent3d.h>
#include <SmFace.h>
#include <SmVertex.h>

using namespace smlib_shim;

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
            // What SmApiBrepComputeMassProperties checks. It clamps the accuracy to [1e-4, 1e-1],
            // which leaves curved faces off by 1e-5 to 1e-3; the kernel call itself takes the
            // [1e-8, 1e-1] that SmApiBrepComputeVolume allows, so this calls it directly.
            if (!b->IsManifoldSolid())
            {
                return static_cast<int32_t>(SM_ERR_INVALID_INPUT);
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
            // SmApiBrepComputeArea face by face, with the [1e-8, 1e-1] envelope instead of its
            // [1e-4, 1e-1]. Every face counts once, so sheets and open shells are measured too.
            const double accuracy = std::min(std::max(relative_accuracy, 1e-8), 1e-1);
            SmTArray<SmFace *> faces;
            brep->p->GetFaces(faces);
            double total = 0.0;
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
                double area = 0.0, volume = 0.0;
                SmTArray<SmVector3d> moments;
                s = status(f->ComputePreciseProperties(SM_OT_SAME, accuracy, box.GetMid(), edge * edge / 100.0, 0.0, area,
                                                      volume, moments, SM_PPF_AREA));
                if (s != SMLIB_OK)
                {
                    return s;
                }
                total += std::fabs(area);
            }
            if (!std::isfinite(total))
            {
                return static_cast<int32_t>(SM_ERR_INVALID_INPUT);
            }
            *out = total;
            return SMLIB_OK;
        }
    );
}
