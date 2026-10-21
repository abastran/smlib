/* SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors */
/* SPDX-License-Identifier: Apache-2.0 */

/* The operations from plain C: booleans, fillet, chamfer and shell on boxes against analytic
 * volumes, the selection errors, and hostile arguments (NULL, NaN, out of range). Every call
 * must return a status; under ASan, none may touch freed or foreign memory. */

#include "smlib.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int failures = 0;

#define CHECK(cond)                                                                                \
    do                                                                                             \
    {                                                                                              \
        if (!(cond))                                                                               \
        {                                                                                          \
            fprintf(stderr, "%s:%d: CHECK(%s) failed\n", __FILE__, __LINE__, #cond);               \
            failures++;                                                                            \
        }                                                                                          \
    } while (0)

static smlib_brep *box(double x, double y, double z, double dx, double dy, double dz)
{
    const double o[3] = {x, y, z};
    const double s[3] = {dx, dy, dz};
    smlib_brep *b = NULL;
    CHECK(smlib_make_box(o, s, &b) == SMLIB_OK);
    return b;
}

static double volume(const smlib_brep *b)
{
    double v = -1.0;
    CHECK(smlib_brep_volume(b, 1e-7, &v) == SMLIB_OK);
    return v;
}

static int near(double got, double want)
{
    return fabs(got - want) <= 1e-6 * fabs(want);
}

int main(void)
{
    smlib_brep *a = box(0, 0, 0, 20, 20, 20);
    smlib_brep *b = box(10, 5, 5, 20, 10, 10);
    smlib_brep *r = NULL;
    int64_t n = 0;

    CHECK(smlib_boolean(a, b, SMLIB_BOOLEAN_UNION, &r) == SMLIB_OK);
    CHECK(near(volume(r), 9000.0));
    CHECK(smlib_brep_solid_count(r, &n) == SMLIB_OK && n == 1);
    smlib_brep_free(r);
    r = NULL;
    CHECK(smlib_boolean(a, b, SMLIB_BOOLEAN_SUBTRACT, &r) == SMLIB_OK);
    CHECK(near(volume(r), 7000.0));
    smlib_brep_free(r);
    r = NULL;
    CHECK(smlib_boolean(a, b, SMLIB_BOOLEAN_INTERSECT, &r) == SMLIB_OK);
    CHECK(near(volume(r), 1000.0));
    smlib_brep_free(r);
    r = NULL;
    CHECK(near(volume(a), 8000.0) && near(volume(b), 2000.0));

    smlib_brep *far = box(50, 0, 0, 10, 10, 10);
    CHECK(smlib_boolean(a, far, SMLIB_BOOLEAN_INTERSECT, &r) == SMLIB_ERR_EMPTY && r == NULL);
    /* The same body twice: its copies coincide everywhere. Any status, but a success is right. */
    if (smlib_boolean(a, a, SMLIB_BOOLEAN_UNION, &r) == SMLIB_OK)
    {
        CHECK(near(volume(r), 8000.0));
        smlib_brep_free(r);
        r = NULL;
    }
    /* Two cubes sharing an edge: their union is closed but not manifold (four faces on that
     * edge); smlib_boolean hands it out, two solids of the two volumes, and refuses it only
     * without SMLIB_RESULT_NONMANIFOLD. Cubes sharing only a vertex are manifold edge by edge. */
    {
        smlib_brep *p = box(0, 0, 0, 10, 10, 10);
        smlib_brep *q = box(10, 10, 0, 10, 10, 10);
        smlib_brep *v = box(10, 10, 10, 10, 10, 10);
        CHECK(smlib_boolean(p, q, SMLIB_BOOLEAN_UNION, &r) == SMLIB_OK && r != NULL);
        if (r != NULL)
        {
            CHECK(near(volume(r), 2000.0));
            CHECK(smlib_brep_solid_count(r, &n) == SMLIB_OK && n == 2);
            smlib_brep_free(r);
            r = NULL;
        }
        CHECK(smlib_boolean_with(p, q, SMLIB_BOOLEAN_UNION, 0, &r) == SMLIB_ERR_NOT_SOLID && r == NULL);
        CHECK(smlib_boolean_with(p, q, SMLIB_BOOLEAN_UNION, SMLIB_RESULT_NONMANIFOLD, &r) == SMLIB_OK);
        smlib_brep_free(r);
        r = NULL;
        CHECK(smlib_boolean(p, v, SMLIB_BOOLEAN_UNION, &r) == SMLIB_OK && r != NULL);
        if (r != NULL)
        {
            CHECK(near(volume(r), 2000.0));
            smlib_brep_free(r);
            r = NULL;
        }
        CHECK(smlib_boolean_with(p, q, SMLIB_BOOLEAN_UNION, 2, &r) == SMLIB_ERR_ARG && r == NULL);
        smlib_brep_free(p);
        smlib_brep_free(q);
        smlib_brep_free(v);
    }
    CHECK(smlib_boolean(a, b, 7, &r) == SMLIB_ERR_ARG);
    CHECK(smlib_boolean(NULL, b, SMLIB_BOOLEAN_UNION, &r) == SMLIB_ERR_ARG);
    CHECK(smlib_boolean(a, NULL, SMLIB_BOOLEAN_UNION, &r) == SMLIB_ERR_ARG);
    CHECK(smlib_boolean(a, b, SMLIB_BOOLEAN_UNION, NULL) == SMLIB_ERR_ARG);
    CHECK(smlib_brep_solid_count(NULL, &n) == SMLIB_ERR_ARG);
    CHECK(smlib_brep_solid_count(a, NULL) == SMLIB_ERR_ARG);

    const double edge[3] = {10, 20, 20};
    int64_t failed = -1;
    CHECK(smlib_chamfer(a, edge, 1, 0.1, 3.0, &r, &failed) == SMLIB_OK);
    CHECK(near(volume(r), 7910.0));
    smlib_brep_free(r);
    r = NULL;
    CHECK(smlib_fillet(a, edge, 1, 0.1, 2.0, &r, &failed) == SMLIB_OK);
    smlib_brep_free(r);
    r = NULL;

    /* smlib_blend: each edge at its point's size. Two parallel top edges rounded at 2 and 3
     * lose (1 - pi/4) r^2 L each. */
    const double tops[6] = {10, 20, 20, 10, 0, 20};
    const double radii[2] = {2.0, 3.0};
    const double quarter = 1.0 - 3.14159265358979323846 / 4.0;
    CHECK(smlib_blend(a, SMLIB_BLEND_ROUND, tops, radii, 2, 0.1, &r, &failed) == SMLIB_OK);
    CHECK(fabs(volume(r) - (8000.0 - quarter * (4.0 + 9.0) * 20.0)) < 1e-3);
    smlib_brep_free(r);
    r = NULL;
    CHECK(smlib_blend(a, SMLIB_BLEND_CHAMFER, tops, radii, 2, 0.1, &r, &failed) == SMLIB_OK);
    CHECK(fabs(volume(r) - (8000.0 - 0.5 * (4.0 + 9.0) * 20.0)) < 1e-3);
    smlib_brep_free(r);
    r = NULL;
    /* A chord of 2 sqrt(2) on a right-angled edge is the round of radius 2. */
    const double chord = 2.0 * sqrt(2.0);
    CHECK(smlib_blend(a, SMLIB_BLEND_CHORD, edge, &chord, 1, 0.1, &r, &failed) == SMLIB_OK);
    CHECK(fabs(volume(r) - (8000.0 - quarter * 4.0 * 20.0)) < 1e-2);
    smlib_brep_free(r);
    r = NULL;
    /* A curvature-continuous blend between the round's contact lines starts flat at each one, so
     * it runs outside the round and takes less, but something. */
    CHECK(smlib_blend(a, SMLIB_BLEND_CURVATURE, edge, radii, 1, 0.1, &r, &failed) == SMLIB_OK);
    double g2 = volume(r);
    CHECK(g2 > 8000.0 - quarter * 4.0 * 20.0 && g2 < 8000.0 - 1.0);
    smlib_brep_free(r);
    r = NULL;
    CHECK(smlib_blend(a, 9, edge, radii, 1, 0.1, &r, NULL) == SMLIB_ERR_ARG);
    CHECK(smlib_blend(a, SMLIB_BLEND_ROUND, edge, NULL, 1, 0.1, &r, NULL) == SMLIB_ERR_ARG);
    const double bad[2] = {1.0, -1.0};
    CHECK(smlib_blend(a, SMLIB_BLEND_ROUND, tops, bad, 2, 0.1, &r, NULL) == SMLIB_ERR_ARG);

    /* A radius running linearly from 3 at x = 0 to 1 at x = 20 along the edge at y = 20, z = 20:
     * (1 - pi/4) L (r0^2 + r0 r1 + r1^2) / 3 taken away, whichever way the edge runs. */
    const double knots[4] = {0, 0, 1, 1};
    const double law[2] = {3.0, 1.0};
    const double from0[3] = {0, 20, 20};
    const double from1[3] = {20, 20, 20};
    double linear = 8000.0 - quarter * 20.0 * (9.0 + 3.0 + 1.0) / 3.0;
    CHECK(smlib_fillet_variable(a, edge, 0.1, from0, 1, knots, 4, law, 2, &r) == SMLIB_OK);
    double v0 = volume(r);
    CHECK(fabs(v0 - linear) < 0.05);
    smlib_brep_free(r);
    r = NULL;
    CHECK(smlib_fillet_variable(a, edge, 0.1, from1, 1, knots, 4, law, 2, &r) == SMLIB_OK);
    CHECK(fabs(volume(r) - v0) < 1e-6 * 8000.0);
    smlib_brep_free(r);
    r = NULL;
    /* A constant law (2 -> 2, and the same as a quadratic) is the round of radius 2; a quadratic
     * with control radii 3, 2, 1 is the linear law above. */
    const double flat[2] = {2.0, 2.0};
    CHECK(smlib_fillet_variable(a, edge, 0.1, from0, 1, knots, 4, flat, 2, &r) == SMLIB_OK);
    CHECK(fabs(volume(r) - (8000.0 - quarter * 4.0 * 20.0)) < 1e-6 * 8000.0);
    smlib_brep_free(r);
    r = NULL;
    const double knots2[6] = {0, 0, 0, 1, 1, 1};
    const double flat3[3] = {2.0, 2.0, 2.0};
    CHECK(smlib_fillet_variable(a, edge, 0.1, from0, 2, knots2, 6, flat3, 3, &r) == SMLIB_OK);
    CHECK(fabs(volume(r) - (8000.0 - quarter * 4.0 * 20.0)) < 1e-6 * 8000.0);
    smlib_brep_free(r);
    r = NULL;
    const double falling[3] = {3.0, 2.0, 1.0};
    CHECK(smlib_fillet_variable(a, edge, 0.1, from0, 2, knots2, 6, falling, 3, &r) == SMLIB_OK);
    CHECK(fabs(volume(r) - v0) < 1e-6 * 8000.0);
    smlib_brep_free(r);
    r = NULL;
    const double unclamped[4] = {0, 0.5, 1, 1};
    CHECK(smlib_fillet_variable(a, edge, 0.1, from0, 1, unclamped, 4, law, 2, &r) == SMLIB_ERR_ARG);
    CHECK(smlib_fillet_variable(a, edge, 0.1, from0, 4, knots, 4, law, 2, &r) == SMLIB_ERR_ARG);
    CHECK(smlib_fillet_variable(a, edge, 0.1, NULL, 1, knots, 4, law, 2, &r) == SMLIB_ERR_ARG);
    CHECK(smlib_fillet_variable(a, edge, 0.1, from0, 1, knots, 3, law, 2, &r) == SMLIB_ERR_ARG);

    const double top[3] = {10, 10, 20};
    CHECK(smlib_shell(a, top, 1, 0.1, -2.0, &r, &failed) == SMLIB_OK);
    CHECK(near(volume(r), 3392.0));
    smlib_brep_free(r);
    r = NULL;
    CHECK(smlib_shell(a, NULL, 0, 0.1, -2.0, &r, NULL) == SMLIB_OK);
    smlib_brep_free(r);
    r = NULL;

    /* Selection errors name the point. */
    const double two[6] = {10, 20, 20, 100, 100, 100};
    CHECK(smlib_fillet(a, two, 2, 0.1, 1.0, &r, &failed) == SMLIB_ERR_PICK_NONE && failed == 1 && r == NULL);
    const double corner[3] = {20, 20, 20};
    failed = -1;
    CHECK(smlib_fillet(a, corner, 1, 0.1, 1.0, &r, &failed) == SMLIB_ERR_PICK_AMBIGUOUS && failed == 0);
    CHECK(smlib_shell(a, edge, 1, 0.1, -1.0, &r, NULL) == SMLIB_ERR_PICK_AMBIGUOUS);

    /* Hostile arguments. */
    const double nan3[3] = {NAN, 0, 0};
    CHECK(smlib_fillet(NULL, edge, 1, 0.1, 1.0, &r, NULL) == SMLIB_ERR_ARG);
    CHECK(smlib_fillet(a, NULL, 1, 0.1, 1.0, &r, NULL) == SMLIB_ERR_ARG);
    CHECK(smlib_fillet(a, edge, 0, 0.1, 1.0, &r, NULL) == SMLIB_ERR_ARG);
    CHECK(smlib_fillet(a, edge, -1, 0.1, 1.0, &r, NULL) == SMLIB_ERR_ARG);
    CHECK(smlib_fillet(a, edge, 1 << 20, 0.1, 1.0, &r, NULL) == SMLIB_ERR_ARG);
    CHECK(smlib_fillet(a, nan3, 1, 0.1, 1.0, &r, NULL) == SMLIB_ERR_ARG);
    CHECK(smlib_fillet(a, edge, 1, NAN, 1.0, &r, NULL) == SMLIB_ERR_ARG);
    CHECK(smlib_fillet(a, edge, 1, 0.1, NAN, &r, NULL) == SMLIB_ERR_ARG);
    CHECK(smlib_fillet(a, edge, 1, 0.1, -1.0, &r, NULL) == SMLIB_ERR_ARG);
    CHECK(smlib_fillet(a, edge, 1, 0.1, 1.0, NULL, NULL) == SMLIB_ERR_ARG);
    CHECK(smlib_chamfer(a, edge, 1, 0.1, INFINITY, &r, NULL) == SMLIB_ERR_ARG);
    CHECK(smlib_shell(a, top, 1, 0.1, 0.0, &r, NULL) == SMLIB_ERR_ARG);
    CHECK(smlib_shell(a, NULL, 1, 0.1, -1.0, &r, NULL) == SMLIB_ERR_ARG);
    CHECK(r == NULL);

    /* Sizes SMLib cannot build: any status but a crash; a success must be a real body. */
    if (smlib_fillet(a, edge, 1, 0.1, 50.0, &r, NULL) == SMLIB_OK)
    {
        smlib_brep_free(r);
        r = NULL;
    }
    if (smlib_chamfer(a, edge, 1, 0.1, 30.0, &r, NULL) == SMLIB_OK)
    {
        smlib_brep_free(r);
        r = NULL;
    }
    CHECK(smlib_shell(a, top, 1, 0.1, -15.0, &r, NULL) == SMLIB_ERR_NO_ROOM && r == NULL);
    CHECK(near(volume(a), 8000.0));

    smlib_brep_free(far);
    smlib_brep_free(b);
    smlib_brep_free(a);
    if (failures != 0)
    {
        fprintf(stderr, "%d check(s) failed\n", failures);
        return 1;
    }
    printf("ops: all checks passed\n");
    return 0;
}
