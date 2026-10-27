// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// smlib_imprint: a box cut by a plane sheet through its middle, a sheet cut by a solid that
// crosses it, a tool that misses, and rejected arguments.

#include "../src/smlib_shim.h"

#include <SmApiGeneral.h>
#include <SmApiPrimitives.h>

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

static smlib_brep *plane_at_z(double z, double x0, double size)
{
    SmApiGetOrCreateContext();
    SmBrep *p = nullptr;
    if (SmApiCreatePlane(SmVector3d(x0, x0, 0.0), size, size, p) != SM_SUCCESS)
    {
        return nullptr;
    }
    SmApiTranslate(p, SmVector3d(0.0, 0.0, z));
    return wrap(p);
}

int main()
{
    CHECK(smlib_init() == SMLIB_OK);
    const double o[3] = {0, 0, 0}, s[3] = {10, 10, 10};
    smlib_brep *box = nullptr;
    CHECK(smlib_make_box(o, s, &box) == SMLIB_OK);

    // A plane through the middle cuts the four side faces: 6 + 4 faces, the volume unchanged.
    smlib_brep *mid = plane_at_z(5.0, -5.0, 20.0);
    CHECK(mid != nullptr);
    smlib_brep *cut = nullptr;
    CHECK(smlib_imprint(box, mid, &cut) == SMLIB_OK);
    int64_t n[3] = {0, 0, 0};
    CHECK(smlib_brep_counts(cut, n) == SMLIB_OK && n[0] == 10);
    double v = 0.0;
    CHECK(smlib_brep_volume(cut, 1e-7, &v) == SMLIB_OK && std::fabs(v - 1000.0) < 1e-6);
    // The inputs are unchanged.
    CHECK(smlib_brep_counts(box, n) == SMLIB_OK && n[0] == 6);
    CHECK(smlib_brep_counts(mid, n) == SMLIB_OK && n[0] == 1);

    // The sheet cut by the box: its one face splits where the box's walls cross it.
    smlib_brep *sheet_cut = nullptr;
    CHECK(smlib_imprint(mid, box, &sheet_cut) == SMLIB_OK);
    CHECK(smlib_brep_counts(sheet_cut, n) == SMLIB_OK && n[0] >= 2);
    double area = 0.0;
    CHECK(smlib_brep_area(sheet_cut, 1e-7, &area) == SMLIB_OK && std::fabs(area - 400.0) < 1e-6);

    // A plane that misses the box cuts nothing.
    smlib_brep *far = plane_at_z(50.0, -5.0, 20.0);
    smlib_brep *none = nullptr;
    CHECK(smlib_imprint(box, far, &none) == SMLIB_ERR_NO_CUT && none == nullptr);

    CHECK(smlib_imprint(box, box, &none) == SMLIB_ERR_ARG);
    CHECK(smlib_imprint(nullptr, mid, &none) == SMLIB_ERR_ARG);
    CHECK(smlib_imprint(box, mid, nullptr) == SMLIB_ERR_ARG);

    for (smlib_brep *b : {box, mid, cut, sheet_cut, far})
    {
        smlib_brep_free(b);
    }
    if (failures == 0)
    {
        std::printf("imprint: ok\n");
    }
    return failures == 0 ? 0 : 1;
}
