// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// smlib_extend_sheet: a square plane sheet continued past one side in each mode, a bilinear
// B-spline sheet continued naturally, and the shapes and arguments it refuses.

#include "../src/smlib_shim.h"

#include <SmApiGeneral.h>
#include <SmApiPrimitives.h>
#include <SmApiSurfaces.h>
#include <SmBrep.h>
#include <SmSurface.h>

#include <cmath>
#include <cstdio>

static int failures = 0;

#define CHECK(cond)                                                                                \
    do                                                                                             \
    {                                                                                              \
        if (!(cond))                                                                               \
        {                                                                                          \
            std::fprintf(stderr, "%s:%d: CHECK(%s) failed\n", __FILE__, __LINE__, #cond);          \
            failures++;                                                                            \
        }                                                                                          \
    } while (0)

static smlib_brep *wrap(SmBrep *p)
{
    return p == nullptr ? nullptr : new smlib_brep{p};
}

static double area_of(const smlib_brep *b)
{
    double a = -1.0;
    return smlib_brep_area(b, 1e-7, &a) == SMLIB_OK ? a : -1.0;
}

int main()
{
    CHECK(smlib_init() == SMLIB_OK);
    const SmContext &ctx = *SmApiGetOrCreateContext();

    // The square [0, 10] x [0, 10] at z = 0, continued 5 past its side x = 10: 15 x 10 in each
    // mode (a plane's G1 and G2 run-offs are the plane).
    SmBrep *p = nullptr;
    CHECK(SmApiCreatePlane(SmVector3d(0.0, 0.0, 0.0), 10.0, 10.0, p) == SM_SUCCESS);
    smlib_brep *square = wrap(p);
    const double side[3] = {10.0, 5.0, 0.0};
    for (int32_t mode : {SMLIB_EXTEND_NATURAL, SMLIB_EXTEND_TANGENT, SMLIB_EXTEND_CURVATURE})
    {
        smlib_brep *r = nullptr;
        CHECK(smlib_extend_sheet(square, side, 1e-3, 5.0, mode, &r) == SMLIB_OK);
        CHECK(std::fabs(area_of(r) - 150.0) < 1e-6);
        int64_t n[3] = {0, 0, 0};
        CHECK(smlib_brep_counts(r, n) == SMLIB_OK && n[0] == 1 && n[1] == 4);
        smlib_brep_free(r);
    }
    // The input is unchanged.
    CHECK(std::fabs(area_of(square) - 100.0) < 1e-6);

    // A bilinear B-spline sheet (a flat square 10 x 10) continued 4 past its side y = 10.
    SmSurface *s = nullptr;
    CHECK(SmApiCreateSurfaceFromCornerPoints(SmPoint3d(0, 0, 0), SmPoint3d(10, 0, 0), SmPoint3d(0, 10, 0), SmPoint3d(10, 10, 0), s) ==
          SM_SUCCESS);
    SmBrep *bs = new (ctx) SmBrep();
    SmFace *f = nullptr;
    CHECK(s != nullptr && bs->CreateFaceFromSurface(s, s->GetNaturalUVDomain(), f) == SM_SUCCESS);
    smlib_brep *spline = wrap(bs);
    const double top[3] = {5.0, 10.0, 0.0};
    smlib_brep *grown = nullptr;
    CHECK(smlib_extend_sheet(spline, top, 1e-3, 4.0, SMLIB_EXTEND_NATURAL, &grown) == SMLIB_OK);
    CHECK(std::fabs(area_of(grown) - 140.0) < 1e-4);

    // Refused: a box (not one face), a pick on no edge, bad arguments.
    const double o[3] = {0, 0, 0}, sz[3] = {10, 10, 10};
    smlib_brep *box = nullptr;
    CHECK(smlib_make_box(o, sz, &box) == SMLIB_OK);
    smlib_brep *none = nullptr;
    CHECK(smlib_extend_sheet(box, side, 1e-3, 5.0, SMLIB_EXTEND_NATURAL, &none) == SMLIB_ERR_UNSUPPORTED && none == nullptr);
    const double mid[3] = {5.0, 5.0, 0.0};
    CHECK(smlib_extend_sheet(square, mid, 1e-3, 5.0, SMLIB_EXTEND_NATURAL, &none) == SMLIB_ERR_PICK_NONE);
    CHECK(smlib_extend_sheet(square, side, 1e-3, -1.0, SMLIB_EXTEND_NATURAL, &none) == SMLIB_ERR_ARG);
    CHECK(smlib_extend_sheet(square, side, 1e-3, 5.0, 7, &none) == SMLIB_ERR_ARG);
    CHECK(smlib_extend_sheet(nullptr, side, 1e-3, 5.0, SMLIB_EXTEND_NATURAL, &none) == SMLIB_ERR_ARG);

    for (smlib_brep *b : {square, spline, grown, box})
    {
        smlib_brep_free(b);
    }
    if (failures == 0)
    {
        std::printf("extend: ok\n");
    }
    return failures == 0 ? 0 : 1;
}
