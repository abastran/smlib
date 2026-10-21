// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// Booleans that SMLib got wrong (tangency-labs/solidcore issues), each against its closed form
// and SMLib's own check. Run it under ASan/UBSan too (README, "Fuzzing"). With arguments, only
// the named cases run (`smlib_test_booleans tangent_cylinders`).

#include "../src/smlib_shim.h"

#include <SmApiPrimitives.h>
#include <SmVector3d.h>

#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>

namespace
{

int failures = 0;

#define CHECK(cond)                                                                                \
    do                                                                                             \
    {                                                                                              \
        if (!(cond))                                                                               \
        {                                                                                          \
            std::fprintf(stderr, "%s:%d: CHECK(%s) failed\n", __FILE__, __LINE__, #cond);          \
            failures++;                                                                            \
        }                                                                                          \
    } while (0)

const double PI = 3.14159265358979323846;

smlib_brep *box(double x0, double y0, double z0, double x1, double y1, double z1)
{
    const double o[3] = {x0, y0, z0};
    const double s[3] = {x1 - x0, y1 - y0, z1 - z0};
    smlib_brep *b = nullptr;
    CHECK(smlib_make_box(o, s, &b) == SMLIB_OK);
    return b;
}

double volume(const smlib_brep *b)
{
    double v = -1.0;
    CHECK(smlib_brep_volume(b, 1e-7, &v) == SMLIB_OK);
    return v;
}

// SMLib's own validity check, every test (level 2), must find nothing wrong with `b`.
bool sound(const smlib_brep *b)
{
    static char buf[1 << 16];
    int64_t len = 0;
    int64_t bad = -1;
    if (b == nullptr || smlib_brep_check(b, 2, buf, sizeof buf, &len, &bad) != SMLIB_OK)
    {
        std::fprintf(stderr, "  smlib_brep_check did not run\n");
        return false;
    }
    if (bad != 0)
    {
        std::fprintf(stderr, "  %lld failed validity tests:\n%s\n", static_cast<long long>(bad), buf);
    }
    return bad == 0;
}

bool near(double got, double want, double rel)
{
    bool ok = std::fabs(got - want) <= rel * std::fabs(want);
    if (!ok)
    {
        std::fprintf(stderr, "  got %.9g, want %.9g\n", got, want);
    }
    return ok;
}

// `a op b` must be one valid body of volume `want` (relative `rel`) and `solids` solids.
void expect(const char *what, const smlib_brep *a, const smlib_brep *b, int32_t op, double want, double rel, int64_t solids)
{
    smlib_brep *r = nullptr;
    int32_t s = smlib_boolean(a, b, op, &r);
    if (s != SMLIB_OK)
    {
        std::fprintf(stderr, "%s: status %d\n", what, static_cast<int>(s));
        failures++;
        return;
    }
    int f0 = failures;
    CHECK(sound(r));
    CHECK(near(volume(r), want, rel));
    int64_t n = -1;
    CHECK(smlib_brep_solid_count(r, &n) == SMLIB_OK && n == solids);
    if (failures != f0)
    {
        std::fprintf(stderr, "%s: wrong result above\n", what);
    }
    smlib_brep_free(r);
}

// A cylinder of SMLib's own along `axis` (0 x, 1 y, 2 z) from its base centre.
smlib_brep *cylinder(int axis, double x, double y, double z, double r, double h)
{
    const double o[3] = {x, y, z};
    smlib_brep *b = nullptr;
    CHECK(smlib_make_cylinder(o, r, h, &b) == SMLIB_OK);
    if (b != nullptr && axis != 2)
    {
        // +z turned onto +x (about +y by 90 degrees) or +y (about +x by -90 degrees).
        SmVector3d about = axis == 0 ? SmVector3d(0, 1, 0) : SmVector3d(1, 0, 0);
        CHECK(SmApiRotate(b->p, SmVector3d(x, y, z), about, axis == 0 ? 90.0 : -90.0) == SM_SUCCESS);
    }
    return b;
}

// The same body through the flat arrays, as SolidCraft hands its cylinders over: analytic
// surfaces with no UV trims on them, SMLib rebuilding its own.
smlib_brep *through_flat(const smlib_brep *b)
{
    smlib_flat *f = nullptr;
    smlib_brep *out = nullptr;
    CHECK(smlib_brep_export(b, &f) == SMLIB_OK);
    CHECK(smlib_brep_import(f, &out) == SMLIB_OK);
    smlib_flat_free(f);
    return out;
}

smlib_brep *torus(double x, double y, double z, double R, double r)
{
    const double c[3] = {x, y, z};
    smlib_brep *b = nullptr;
    CHECK(smlib_make_torus(c, R, r, &b) == SMLIB_OK);
    return b;
}

// The volume of a torus about z (centre c, radii R, r) inside the prism over [x0,x1] x [y0,y1]
// that holds the torus's whole height there: the integral of the tube's height 2 sqrt(r^2 -
// (rho - R)^2), midpoint rule.
double torus_in_prism(double cx, double cy, double R, double r, double x0, double x1, double y0, double y1)
{
    const int n = 3000;
    double dx = (x1 - x0) / n;
    double dy = (y1 - y0) / n;
    double sum = 0.0;
    for (int i = 0; i < n; ++i)
    {
        double x = x0 + (i + 0.5) * dx - cx;
        for (int j = 0; j < n; ++j)
        {
            double y = y0 + (j + 0.5) * dy - cy;
            double t = std::sqrt(x * x + y * y) - R;
            if (t * t < r * r)
            {
                sum += 2.0 * std::sqrt(r * r - t * t);
            }
        }
    }
    return sum * dx * dy;
}

// The integral of sqrt(r^2 - u^2) from a to b (|a|, |b| <= r).
double half_disc(double r, double a, double b)
{
    auto F = [r](double u) { return 0.5 * (u * std::sqrt(r * r - u * u) + r * r * std::asin(u / r)); };
    return F(b) - F(a);
}

// tangency-labs/solidcore#383: a cylinder or a torus tangent to a face of a box along a line or
// a circle (the seeds of SolidCraft's random boolean runner). The tangency stays an edge where
// the faces touch, shared by four faces: a closed body, non-manifold there, which SMLib's check
// refused (SmEdge test 9: the tangent faces' binormals tie, so it could not see their order).
// Each cylinder goes in twice: as SMLib makes it and through the flat arrays.
void tangent_cylinders()
{
    const double pi = PI;
    struct Cyl
    {
        int seed;
        int axis;
        double base[3], r, h;
        double box[6];
        int32_t op;
        bool cyl_first;
        double want;
        int64_t solids;
    };
    const double a74 = half_disc(10, -2.5, 7.5) * 15.0;
    const Cyl cyls[] = {
        {1, 2, {15, 25, 22.5}, 7.5, 22.5, {22.5, 22.5, 30, 27.5, 50, 45}, SMLIB_BOOLEAN_UNION, true,
         pi * 7.5 * 7.5 * 22.5 + 5.0 * 27.5 * 15.0, 2},
        {32, 2, {12.5, 20, 10}, 5, 30, {7.5, 15, 0, 25, 32.5, 27.5}, SMLIB_BOOLEAN_SUBTRACT, false,
         17.5 * 17.5 * 27.5 - pi * 25.0 * 17.5, 1},
        {36, 1, {5, 30, 30}, 2.5, 7.5, {2.5, 17.5, 20, 27.5, 42.5, 27.5}, SMLIB_BOOLEAN_UNION, true,
         25.0 * 25.0 * 7.5 + pi * 6.25 * 7.5, 2},
        {39, 2, {25, 12.5, 22.5}, 5, 25, {20, 5, 5, 35, 20, 32.5}, SMLIB_BOOLEAN_SUBTRACT, false,
         15.0 * 15.0 * 27.5 - pi * 25.0 * 10.0, 1},
        {74, 1, {10, 0, 22.5}, 10, 27.5, {7.5, 12.5, 12.5, 17.5, 32.5, 22.5}, SMLIB_BOOLEAN_INTERSECT, false, a74, 1},
        {74, 1, {10, 0, 22.5}, 10, 27.5, {7.5, 12.5, 12.5, 17.5, 32.5, 22.5}, SMLIB_BOOLEAN_SUBTRACT, false, 2000.0 - a74, 1},
        {85, 0, {15, 10, 30}, 5, 25, {17.5, 5, 7.5, 32.5, 25, 37.5}, SMLIB_BOOLEAN_INTERSECT, false, 15.0 * 25.0 * pi, 1},
        {85, 0, {15, 10, 30}, 5, 25, {17.5, 5, 7.5, 32.5, 25, 37.5}, SMLIB_BOOLEAN_SUBTRACT, false,
         9000.0 - 15.0 * 25.0 * pi, 1},
    };
    for (const Cyl &c : cyls)
    {
        smlib_brep *made = cylinder(c.axis, c.base[0], c.base[1], c.base[2], c.r, c.h);
        smlib_brep *flat = through_flat(made);
        smlib_brep *b = box(c.box[0], c.box[1], c.box[2], c.box[3], c.box[4], c.box[5]);
        int k = 0;
        for (smlib_brep *cyl : {made, flat})
        {
            char what[64];
            std::snprintf(what, sizeof what, "seed %d op %d %s", c.seed, static_cast<int>(c.op), k++ ? "flat" : "made");
            if (c.cyl_first)
            {
                expect(what, cyl, b, c.op, c.want, 1e-6, c.solids);
            }
            else
            {
                expect(what, b, cyl, c.op, c.want, 1e-6, c.solids);
            }
        }
        for (smlib_brep *x : {made, flat, b})
        {
            smlib_brep_free(x);
        }
    }
    // Tori about z, tangent to a box's face along a circle. Seed 86's tube runs across the box's
    // corner just under its top: the sliver between the tube, the top and the corner touches the
    // rest along the circle only, a second solid.
    {
        smlib_brep *t = torus(20, 22.5, 25, 10, 2.5);
        smlib_brep *b = box(22.5, 25, 22.5, 30, 32.5, 50);
        expect("seed 4 cut", b, t, SMLIB_BOOLEAN_SUBTRACT, 7.5 * 7.5 * 27.5 - torus_in_prism(20, 22.5, 10, 2.5, 22.5, 30, 25, 32.5),
               1e-4, 1);
        smlib_brep_free(t);
        smlib_brep_free(b);
    }
    {
        smlib_brep *t = torus(25, 7.5, 25, 10, 5);
        smlib_brep *b = box(12.5, 2.5, 30, 37.5, 32.5, 55);
        expect("seed 65 union", b, t, SMLIB_BOOLEAN_UNION, 25.0 * 30.0 * 25.0 + 2.0 * pi * pi * 10.0 * 25.0, 1e-5, 2);
        smlib_brep_free(t);
        smlib_brep_free(b);
    }
    {
        smlib_brep *t = torus(20, 17.5, 30, 12.5, 2.5);
        smlib_brep *b = box(30, 20, 15, 35, 30, 32.5);
        expect("seed 86 cut", b, t, SMLIB_BOOLEAN_SUBTRACT, 5.0 * 10.0 * 17.5 - torus_in_prism(20, 17.5, 12.5, 2.5, 30, 35, 20, 30),
               1e-4, 2);
        smlib_brep_free(t);
        smlib_brep_free(b);
    }
}

// tangency-labs/solidcore#387's cylinder: base (0,50,-60), axis z, r 40, h 120, as SolidCraft's
// record hands it over: two half-cylinder faces, each a rational quadratic in four spans
// (weights 1, 0.75, 0.5, 0.5, 0.75, 1; the halves meet at 35 and 215 degrees), and bilinear caps;
// every edgeuse on a curved face with its UV trim. The flat arrays, by id.
smlib_brep *rational_cylinder()
{
    struct I
    {
        int32_t id;
        std::vector<int64_t> v;
    };
    struct F
    {
        int32_t id;
        std::vector<double> v;
    };
    const I ints[] = {
        {1, {1, 1}},
        {2, {1, 0}},
        {3, {4, 4}},
        {4, {0, 0}},
        {5, {0, 0}},
        {6, {0, 1, 2, 3, 0, 1, 2, 3}},
        {7, {0, 1, 1, 1, 1, 0, 0, 0}},
        {8, {1, 1, 1, 1}},
        {9, {0, 0, 0, 0}},
        {10, {0, 0, 0, 0}},
        {11, {2, 4, 4, 2}},
        {12, {-1, -1, -1, -1}},
        {13, {1, 0, 1, 2, 3, 4, 0, 4, 5, 2, 3, 5}},
        {14, {1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1}},
        {15, {2, 6, 0, 9, 10, 7, 1, 5, 11, 3, 4, 8}},
        {16, {0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0}},
        {17, {0, 0, 1, 0, 1, 0}},
        {18, {0, 1, 1, 0, 0, 2, 3, 2, 1, 3, 2, 3}},
        {21, {1, 1, 1, 1}},
        {22, {6, 6, 6, 6}},
        {23, {3, 3, 3, 3}},
        {26, {0, 0, 10, 7, 10, 7, 10, 7, 10, 7, 0, 0}},
        {27, {0, 0, 4, 4, 4, 4, 4, 4, 4, 4, 0, 0}},
        {28, {2, 6, 6, 2}},
        {29, {2, 2, 2, 2}},
        {30, {2, 3, 3, 2}},
        {31, {2, 2, 2, 2}},
        {32, {1}},
    };
    const F reals[] = {
        {101, {0.00014422205101855955}},
        {102, {-32.71472772631909, 26.983775509593045, -60.0, 32.71472772631909, 73.01622449040696, 60.0}},
        {103, {0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0}},
        {104, {0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0}},
        {107, {-32.71472772631909, 26.983775509593045, -60.0, 32.71472772631909, 73.01622449040696, -60.0, -32.71472772631909, 26.983775509593045, 60.0, 32.71472772631909, 73.01622449040696, 60.0}},
        {108, {-32.71472772631909, 26.983775509593045, -60.0, -25.04265289618343, 16.078866267486678, -60.0, 6.658860627247413, 5.777160028477432, -60.0, 39.373588353566504, 28.79338451888438, -60.0, 40.38680255645475, 62.11131524830059, -60.0, 32.71472772631909, 73.01622449040696, -60.0, 32.71472772631909, 73.01622449040696, -60.0, 25.042652896183437, 83.92113373251331, -60.0, -6.65886062724741, 94.22283997152256, -60.0, -39.37358835356651, 71.20661548111562, -60.0, -40.38680255645474, 37.888684751699394, -60.0, -32.71472772631908, 26.983775509593034, -60.0, 32.71472772631909, 73.01622449040696, 60.0, 25.042652896183437, 83.92113373251331, 60.0, -6.65886062724741, 94.22283997152256, 60.0, -39.37358835356651, 71.20661548111562, 60.0, -40.38680255645474, 37.888684751699394, 60.0, -32.71472772631908, 26.983775509593034, 60.0, -32.71472772631909, 26.983775509593045, 60.0, -25.04265289618343, 16.078866267486678, 60.0, 6.658860627247413, 5.777160028477432, 60.0, 39.373588353566504, 28.79338451888438, 60.0, 40.38680255645475, 62.11131524830059, 60.0, 32.71472772631909, 73.01622449040696, 60.0}},
        {109, {0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0}},
        {110, {1.0, 0.75, 0.5, 0.5, 0.75, 1.0, 1.0, 0.7500000000000001, 0.49999999999999994, 0.49999999999999983, 0.7499999999999999, 1.0, 1.0, 0.7500000000000001, 0.49999999999999994, 0.49999999999999983, 0.7499999999999999, 1.0, 1.0, 0.75, 0.5, 0.5, 0.75, 1.0}},
        {115, {-32.71472772631909, 26.983775509593045, -60.0, 32.71472772631909, 73.01622449040696, -60.0}},
        {116, {0.0, 0.0, 120.0, 0.0, 0.0, 120.0}},
        {136, {0.0, 0.0, 0.08333333333333358, -2.898371010036392e-17, 0.1666666666666664, 7.760815199088004e-17, 0.33333333333333354, 2.0570757031889698e-17, 0.4166666666666665, -1.9480461156518785e-17, 0.5833333333333336, 4.508066525742703e-17, 0.6666666666666664, 1.4457586211611358e-17, 0.8333333333333336, 1.259182713794139e-17, 0.9166666666666664, 4.801456310251341e-17, 1.0, 0.0, 0.9999999999999999, 0.0, 0.9999999999999994, 0.08333333333333302, 1.0000000000000007, 0.2500000000000006, 0.9999999999999986, 0.4999999999999993, 1.000000000000002, 0.7500000000000014, 0.9999999999999986, 0.9166666666666659, 0.9999999999999999, 1.0, 1.0, 1.0, 0.9166666666666664, 0.9999999999999999, 0.8333333333333336, 1.0000000000000002, 0.6666666666666664, 0.9999999999999999, 0.5833333333333336, 0.9999999999999999, 0.4166666666666665, 1.0000000000000002, 0.33333333333333337, 1.0, 0.1666666666666667, 0.9999999999999997, 0.08333333333333344, 1.0000000000000007, 0.0, 1.0, 1.3782156035472658e-22, 1.0, -2.2755236685602588e-18, 0.9166666666666659, 5.114434962684454e-18, 0.7500000000000014, -8.890953646062393e-18, 0.49999999999999906, 3.010141179406148e-17, 0.2500000000000012, -1.5346166037687467e-17, 0.08333333333333263, 0.0, 0.0, 0.0, 0.0, 0.08333333333333287, 4.6012186791769217e-17, 0.16666666666666716, 5.70547663916104e-17, 0.33333333333333276, 1.4291952170213445e-16, 0.41666666666666746, 2.3669254995652233e-17, 0.5833333333333327, 4.99098811276788e-17, 0.6666666666666671, -3.30318755373146e-17, 0.833333333333333, 1.1433458095462838e-16, 0.9166666666666671, -7.9120602017972e-17, 1.0, 0.0, 1.0, 0.0, 0.9999999999999987, 0.08333333333333302, 1.000000000000002, 0.2500000000000006, 0.9999999999999981, 0.4999999999999993, 1.0000000000000024, 0.7500000000000014, 0.9999999999999986, 0.9166666666666659, 1.0, 1.0, 1.0, 1.0, 0.9166666666666664, 0.9999999999999999, 0.8333333333333336, 1.0000000000000002, 0.6666666666666664, 0.9999999999999999, 0.5833333333333336, 0.9999999999999999, 0.4166666666666665, 1.0000000000000002, 0.33333333333333337, 1.0, 0.1666666666666667, 0.9999999999999997, 0.0833333333333333, 1.0000000000000007, 0.0, 1.0, 3.728351866584152e-29, 1.0, 8.480203362103058e-18, 0.9166666666666659, -1.631636387510398e-17, 0.7500000000000014, 1.9531053419202524e-17, 0.4999999999999993, 2.3165968984404475e-17, 0.2500000000000006, -1.489064146855586e-17, 0.08333333333333302, 0.0, 0.0}},
        {137, {0.0, 0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, -0.75, -0.75, -0.5, -0.5, -0.25, -0.25, -0.0, -0.0, -0.0, -0.0, -1.0, -1.0, -1.0, -1.0, -0.75, -0.5, -0.25, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, -0.75, -0.75, -0.5, -0.5, -0.25, -0.25, -0.0, -0.0, -0.0, -0.0, -1.0, -1.0, -1.0, -1.0, -0.75, -0.5, -0.25, -0.0, -0.0, -0.0, -0.0}},
        {138, {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}},
        {139, {40.38680255645475, 94.22283997152256, -60.0, 40.38680255645475, 5.777160028477439, -60.0, -40.38680255645473, 94.22283997152256, -60.0, -40.38680255645473, 5.777160028477439, -60.0, 32.71472772631909, 73.01622449040696, -60.0, 32.71472772631909, 73.01622449040696, 60.0, 25.042652896183437, 83.92113373251331, -60.0, 25.042652896183437, 83.92113373251331, 60.0, -6.65886062724741, 94.22283997152256, -60.0, -6.65886062724741, 94.22283997152256, 60.0, -39.37358835356651, 71.20661548111562, -60.0, -39.37358835356651, 71.20661548111562, 60.0, -40.38680255645474, 37.888684751699394, -60.0, -40.38680255645474, 37.888684751699394, 60.0, -32.71472772631908, 26.983775509593034, -60.0, -32.71472772631908, 26.983775509593034, 60.0, -32.71472772631909, 26.983775509593045, -60.0, -32.71472772631909, 26.983775509593045, 60.0, -25.04265289618343, 16.078866267486678, -60.0, -25.04265289618343, 16.078866267486678, 60.0, 6.658860627247413, 5.777160028477432, -60.0, 6.658860627247413, 5.777160028477432, 60.0, 39.373588353566504, 28.79338451888438, -60.0, 39.373588353566504, 28.79338451888438, 60.0, 40.38680255645475, 62.11131524830059, -60.0, 40.38680255645475, 62.11131524830059, 60.0, 32.71472772631909, 73.01622449040696, -60.0, 32.71472772631909, 73.01622449040696, 60.0, 40.38680255645475, 94.22283997152256, 60.0, 40.38680255645475, 5.777160028477439, 60.0, -40.38680255645473, 94.22283997152256, 60.0, -40.38680255645473, 5.777160028477439, 60.0}},
        {140, {0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0}},
        {141, {0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0}},
        {142, {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.7500000000000001, 0.7500000000000001, 0.49999999999999994, 0.49999999999999994, 0.49999999999999983, 0.49999999999999983, 0.7499999999999999, 0.7499999999999999, 1.0, 1.0, 1.0, 1.0, 0.75, 0.75, 0.5, 0.5, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}},
        {164, {0.00014422205101855955, 0.00014422205101855955, 0.00014422205101855955, 0.00014422205101855955, 0.00014422205101855955, 0.00014422205101855955}},
        {165, {0.00014422205101855955, 0.00014422205101855955, 0.00014422205101855955, 0.00014422205101855955}},
    };
    smlib_flat *f = smlib_flat_new();
    for (const I &a : ints)
    {
        CHECK(smlib_flat_set_i64(f, a.id, a.v.data(), static_cast<int64_t>(a.v.size())) == SMLIB_OK);
    }
    for (const F &a : reals)
    {
        CHECK(smlib_flat_set_f64(f, a.id, a.v.data(), static_cast<int64_t>(a.v.size())) == SMLIB_OK);
    }
    smlib_brep *b = nullptr;
    CHECK(smlib_brep_import(f, &b) == SMLIB_OK);
    smlib_flat_free(f);
    return b;
}

// tangency-labs/solidcore#387: a box minus that cylinder, which cuts a concave arc into its
// y = 0..20 slab from the +y side. Splitting the half cylinder's face, SMLib copied its surface as
// a surface of revolution, whose NURBS form is SMLib's own two-span arcs: same shape, other
// parameterization, so the trims it kept pointed up to 5 mm away (edge tolerances 0.78 mm). The
// volume: 28 * integral over |x| <= 14 of (sqrt(40^2 - x^2) - 30) comes off.
void rational_cylinder_cut()
{
    smlib_brep *b = box(-14, 0, -14, 14, 20, 14);
    smlib_brep *c = rational_cylinder();
    if (c != nullptr)
    {
        CHECK(sound(c));
        CHECK(near(volume(c), PI * 1600.0 * 120.0, 1e-6));
        const double cut = 28.0 * (half_disc(40.0, -14.0, 14.0) - 30.0 * 28.0);
        expect("box minus the rational cylinder", b, c, SMLIB_BOOLEAN_SUBTRACT, 28.0 * 20.0 * 28.0 - cut, 1e-6, 1);
    }
    smlib_brep_free(b);
    smlib_brep_free(c);
}

} // namespace

int main(int argc, char **argv)
{
    struct Case
    {
        const char *name;
        void (*run)();
    };
    const Case cases[] = {
        {"tangent_cylinders", tangent_cylinders},
        {"rational_cylinder_cut", rational_cylinder_cut},
    };
    for (const Case &c : cases)
    {
        bool chosen = argc < 2;
        for (int i = 1; i < argc; ++i)
        {
            chosen = chosen || std::strcmp(argv[i], c.name) == 0;
        }
        if (chosen)
        {
            int f0 = failures;
            c.run();
            std::printf("%s: %s\n", c.name, failures == f0 ? "ok" : "FAILED");
        }
    }
    if (failures != 0)
    {
        std::fprintf(stderr, "%d check(s) failed\n", failures);
        return 1;
    }
    std::printf("booleans: all checks passed\n");
    return 0;
}
