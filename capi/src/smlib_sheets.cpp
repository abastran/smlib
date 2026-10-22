// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// Sheets (include/smlib.h, "Sheets"): sewing a body's faces along the edges that coincide, and a
// sheet's offset along its faces' normals.
//
// smlib_stitch glues coincident vertices and edges with SmBrep::StitchFaces at the caller's
// tolerance, once, and then orients the faces with SmBrep::OrientTrimmedSurfaces: outward where
// the faces close into a solid, else only made to agree with each other across the edges they
// share (bJustUnifyNormals, which needs no ray firing). SmBrep::StitchAndOrient is not used: it
// caps the tolerance at 100 zone tolerances and stitches a second time at a larger one, so the
// caller's tolerance would not be the one that decides.
//
// smlib_offset_sheet is SM_API's offset (SmApiOffsetBrepFull, through
// SmOffsetExecutive::OffsetBrep) of a sheet: every face onto its surface's offset along its
// normal, as smlib_thicken moves a sheet back by half before a symmetric thicken.

#include "smlib_shim.h"

#include <SmApiBrep.h>
#include <SmApiGeneral.h>
#include <SmApiQueries.h>

#include <memory>

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

} // namespace

extern "C" int32_t smlib_stitch(const smlib_brep *brep, double tolerance, int32_t *out_closed, smlib_brep **out)
{
    if (!valid(brep) || out == nullptr || out_closed == nullptr || !std::isfinite(tolerance) || !(tolerance > 0.0) || tolerance > 1e3)
    {
        return SMLIB_ERR_ARG;
    }
    Call c("smlib_stitch", handle(brep) + ", tolerance=" + num(tolerance));
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
            ULONG stitched = 0;
            ULONG lamina = 0;
            double vertex_gap = 0.0;
            double edge_gap = 0.0;
            s = status(w->StitchFaces(tolerance, stitched, lamina, vertex_gap, edge_gap));
            if (s != SMLIB_OK)
            {
                return s;
            }
            SmBoolean maybe_open = FALSE;
            s = status(w->OrientTrimmedSurfaces(TRUE, maybe_open, lamina != 0 ? TRUE : FALSE));
            if (s != SMLIB_OK)
            {
                return s;
            }
            *out_closed = lamina == 0 && w->IsManifoldSolid() ? 1 : 0;
            settle_tolerances(w.get());
            return take(w.release(), out);
        }
    );
}

extern "C" int32_t smlib_offset_sheet(const smlib_brep *sheet, double distance, smlib_brep **out)
{
    if (!valid(sheet) || out == nullptr || !std::isfinite(distance) || distance == 0.0 || std::fabs(distance) > 1e6)
    {
        return SMLIB_ERR_ARG;
    }
    Call c("smlib_offset_sheet", handle(sheet) + ", distance=" + num(distance));
    c.input("sheet", sheet);
    return c.run(
        [&]
        {
            int32_t s = SMLIB_OK;
            Owned w = copy_of(sheet, s);
            if (s != SMLIB_OK)
            {
                return s;
            }
            if (w->IsManifoldSolid())
            {
                return SMLIB_ERR_ARG;
            }
            SmBrep *moved = nullptr;
            s = status(SmApiOffsetBrepFull(w.get(), distance, TRUE, TRUE, moved));
            if (s != SMLIB_OK || moved == nullptr)
            {
                if (moved != w.get())
                {
                    delete moved;
                }
                return s != SMLIB_OK ? s : SMLIB_ERR_EMPTY;
            }
            // OffsetBrep leaves its input allocated (possibly changed) beside the result.
            if (moved != w.get())
            {
                w.reset(moved);
            }
            settle_tolerances(w.get());
            return take(w.release(), out);
        }
    );
}
