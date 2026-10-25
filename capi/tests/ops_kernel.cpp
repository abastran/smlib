// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// The operations on bodies the C API cannot make yet (a cylinder through SM_API), and the
// SMLib defects found under ASan, each as a regression case. Run it under ASan/UBSan
// (README, "Fuzzing"): every case must stay clean there, not only return the right status.

#include "../src/smlib_shim.h"

#include <SmApiBrep.h>
#include <SmApiFillets.h>
#include <SmApiPrimitives.h>
#include <SmApiQueries.h>
#include <SmBSplineCurve.h>
#include <SmCylinder.h>
#include <SmExtent3d.h>
#include <SmFace.h>
#include <SmFilletExecutive.h>
#include <SmOffsetSurface.h>
#include <SmVector3d.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <utility>
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

smlib_brep *box(double x, double y, double z, double dx, double dy, double dz)
{
    const double o[3] = {x, y, z};
    const double s[3] = {dx, dy, dz};
    smlib_brep *b = nullptr;
    CHECK(smlib_make_box(o, s, &b) == SMLIB_OK);
    return b;
}

// A Z-axis cylinder through SM_API (the C API has no primitives besides the box yet).
smlib_brep *cylinder(double x, double y, double z, double r, double h)
{
    SmApiGetOrCreateContext();
    SmBrep *p = nullptr;
    CHECK(SmApiCreateCylinder(SmVector3d(x, y, z), r, h, p) == SM_SUCCESS && p != nullptr);
    return new smlib_brep{p};
}

double volume_of(const smlib_brep *b)
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
    int64_t failures = -1;
    if (smlib_brep_check(b, 2, buf, sizeof buf, &len, &failures) != SMLIB_OK)
    {
        std::fprintf(stderr, "  smlib_brep_check did not run\n");
        return false;
    }
    if (failures != 0)
    {
        std::fprintf(stderr, "  %lld failed validity tests:\n%s\n", static_cast<long long>(failures), buf);
    }
    return failures == 0;
}

// Checks a result and frees it.
void done(smlib_brep *&r)
{
    CHECK(r != nullptr && sound(r));
    smlib_brep_free(r);
    r = nullptr;
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

// A sphere of radius `r` about the origin in six patches, each on its own sphere surface turned
// to it, as an importer hands over a cube projected onto the sphere: the topology of a box, its
// faces sphere patches (reference direction to the patch's middle, longitude and latitude in
// [-45, 45] degrees), its edges the great-circle arcs between the corners.
smlib_brep *cube_sphere(double r)
{
    using V = std::array<double, 3>;
    auto unit = [](V a) {
        double l = std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
        return V{a[0] / l, a[1] / l, a[2] / l};
    };
    smlib_brep *b = box(-1, -1, -1, 2, 2, 2);
    smlib_flat *f = nullptr;
    CHECK(smlib_brep_export(b, &f) == SMLIB_OK);
    smlib_brep_free(b);
    const int64_t *ev = nullptr;
    const int64_t *loops = nullptr;
    const int64_t *uses = nullptr;
    const int64_t *kinds = nullptr;
    const double *vp = nullptr;
    int64_t n_ev = 0, n_loops = 0, n_uses = 0, n_faces = 0, n_vp = 0;
    CHECK(smlib_flat_get_i64(f, SMLIB_FLAT_EDGE_VERTICES, &ev, &n_ev) == SMLIB_OK);
    CHECK(smlib_flat_get_i64(f, SMLIB_FLAT_LOOP_EDGEUSE_COUNT, &loops, &n_loops) == SMLIB_OK);
    CHECK(smlib_flat_get_i64(f, SMLIB_FLAT_EDGEUSE_EDGE, &uses, &n_uses) == SMLIB_OK);
    CHECK(smlib_flat_get_i64(f, SMLIB_FLAT_FACE_SURFACE_TYPE, &kinds, &n_faces) == SMLIB_OK);
    CHECK(smlib_flat_get_f64(f, SMLIB_FLAT_VERTEX_POSITION, &vp, &n_vp) == SMLIB_OK);
    CHECK(n_faces == 6 && n_loops == 6 && n_ev == 24 && n_vp == 24);
    std::vector<double> position, centre, axis, ref, radius, range;
    for (int64_t i = 0; i < n_vp / 3; ++i)
    {
        V u = unit({vp[3 * i], vp[3 * i + 1], vp[3 * i + 2]});
        for (double x : u)
        {
            position.push_back(r * x);
        }
    }
    for (int64_t e = 0; e < n_ev / 2; ++e)
    {
        V a = unit({vp[3 * ev[2 * e]], vp[3 * ev[2 * e] + 1], vp[3 * ev[2 * e] + 2]});
        V z = unit({vp[3 * ev[2 * e + 1]], vp[3 * ev[2 * e + 1] + 1], vp[3 * ev[2 * e + 1] + 2]});
        V n = unit({a[1] * z[2] - a[2] * z[1], a[2] * z[0] - a[0] * z[2], a[0] * z[1] - a[1] * z[0]});
        for (int k = 0; k < 3; ++k)
        {
            centre.push_back(0.0);
            axis.push_back(n[k]);
            ref.push_back(a[k]);
        }
        radius.push_back(r);
        range.push_back(0.0);
        range.push_back(std::acos(a[0] * z[0] + a[1] * z[1] + a[2] * z[2]));
    }
    std::vector<double> s_centre, s_axis, s_ref, s_radius, s_range;
    int64_t at = 0;
    for (int64_t i = 0; i < n_faces && at < n_uses; ++i)
    {
        // The patch's direction: the middle of its loop's corners (one loop per face, in order).
        V mid{0, 0, 0};
        for (int64_t k = 0; k < loops[i]; ++k)
        {
            int64_t e = uses[at + k];
            for (int j = 0; j < 2; ++j)
            {
                for (int q = 0; q < 3; ++q)
                {
                    mid[q] += vp[3 * ev[2 * e + j] + q];
                }
            }
        }
        at += loops[i];
        V d = unit(mid);
        V up = std::fabs(d[2]) > 0.5 ? V{1, 0, 0} : V{0, 0, 1};
        for (int q = 0; q < 3; ++q)
        {
            s_centre.push_back(0.0);
            s_axis.push_back(up[q]);
            s_ref.push_back(d[q]);
        }
        s_radius.push_back(r);
        for (double x : {-PI / 4.0, -PI / 4.0, PI / 4.0, PI / 4.0})
        {
            s_range.push_back(x);
        }
    }
    std::vector<int64_t> circles(static_cast<size_t>(n_ev / 2), SMLIB_FLAT_CURVE_CIRCLE);
    std::vector<int64_t> spheres(static_cast<size_t>(n_faces), SMLIB_FLAT_SURFACE_SPHERE);
    CHECK(smlib_flat_set_i64(f, SMLIB_FLAT_EDGE_CURVE_TYPE, circles.data(), static_cast<int64_t>(circles.size())) == SMLIB_OK);
    CHECK(smlib_flat_set_i64(f, SMLIB_FLAT_FACE_SURFACE_TYPE, spheres.data(), static_cast<int64_t>(spheres.size())) == SMLIB_OK);
    for (int id : {SMLIB_FLAT_EDGE_LINE_ORIGIN, SMLIB_FLAT_EDGE_LINE_DIRECTION, SMLIB_FLAT_FACE_PLANE_ORIGIN, SMLIB_FLAT_FACE_PLANE_AXIS,
                   SMLIB_FLAT_FACE_PLANE_REF_DIRECTION})
    {
        CHECK(smlib_flat_set_f64(f, id, nullptr, 0) == SMLIB_OK);
    }
    const std::pair<int, const std::vector<double> *> arrays[] = {
        {SMLIB_FLAT_VERTEX_POSITION, &position},       {SMLIB_FLAT_EDGE_CIRCLE_CENTER, &centre},
        {SMLIB_FLAT_EDGE_CIRCLE_AXIS, &axis},          {SMLIB_FLAT_EDGE_CIRCLE_REF_DIRECTION, &ref},
        {SMLIB_FLAT_EDGE_CIRCLE_RADIUS, &radius},      {SMLIB_FLAT_EDGE_RANGE, &range},
        {SMLIB_FLAT_FACE_SPHERE_CENTER, &s_centre},    {SMLIB_FLAT_FACE_SPHERE_AXIS, &s_axis},
        {SMLIB_FLAT_FACE_SPHERE_REF_DIRECTION, &s_ref}, {SMLIB_FLAT_FACE_SPHERE_RADIUS, &s_radius},
        {SMLIB_FLAT_FACE_RANGE, &s_range}};
    for (const auto &a : arrays)
    {
        CHECK(smlib_flat_set_f64(f, a.first, a.second->data(), static_cast<int64_t>(a.second->size())) == SMLIB_OK);
    }
    smlib_brep *out = nullptr;
    CHECK(smlib_brep_import(f, &out) == SMLIB_OK);
    smlib_flat_free(f);
    return out;
}

} // namespace

int main()
{
    smlib_brep *cube = box(0, 0, 0, 20, 20, 20);
    smlib_brep *r = nullptr;

    // SMLib defect 1 (SmProtoEdge's destructor left a pointer to itself in a ProtoVertex that
    // bounds both its ends; ~SmProtoTopologyManager then wrote through it): a box minus a
    // through cylinder.
    smlib_brep *rod = cylinder(10, 10, -1, 4, 22);
    CHECK(smlib_boolean(cube, rod, SMLIB_BOOLEAN_SUBTRACT, &r) == SMLIB_OK);
    CHECK(near(volume_of(r), 8000.0 - PI * 16.0 * 20.0, 1e-6));
    CHECK(sound(r));
    smlib_brep *holed = r;
    r = nullptr;
    CHECK(smlib_boolean(cube, rod, SMLIB_BOOLEAN_UNION, &r) == SMLIB_OK);
    CHECK(near(volume_of(r), 8000.0 + PI * 16.0 * 2.0, 1e-6));
    done(r);
    // A blind hole, ending inside.
    smlib_brep *pin = cylinder(10, 10, 10, 4, 15);
    CHECK(smlib_boolean(cube, pin, SMLIB_BOOLEAN_SUBTRACT, &r) == SMLIB_OK);
    CHECK(near(volume_of(r), 8000.0 - PI * 16.0 * 10.0, 1e-6));
    done(r);

    // The rim of the hole is one closed circular edge.
    // At 45 degrees: the cylinder's seam meets the rim at 0 degrees, and a point there is
    // as near to the seam as to the rim.
    const double rim[3] = {10 + 4 * std::sqrt(0.5), 10 + 4 * std::sqrt(0.5), 20};
    int32_t s = smlib_fillet(holed, rim, 1, 0.1, 1.0, &r, nullptr);
    if (s == SMLIB_OK)
    {
        // A round of r = 1 on a hole of radius R = 4: the cross-section (1 − π/4)·r² goes, its
        // centroid at R + r·(10 − 3π)/(12 − 3π) from the axis (Pappus).
        double area = (1.0 - PI / 4.0);
        double c = 4.0 + (10.0 - 3.0 * PI) / (12.0 - 3.0 * PI);
        CHECK(near(volume_of(r), volume_of(holed) - 2.0 * PI * c * area, 1e-5));
        done(r);
    }
    else
    {
        std::fprintf(stderr, "fillet of the hole rim failed: %d\n", static_cast<int>(s));
        failures++;
    }
    s = smlib_chamfer(holed, rim, 1, 0.1, 1.0, &r, nullptr);
    if (s == SMLIB_OK)
    {
        // A 1 × 1 triangle, centroid at R + 1/3.
        CHECK(near(volume_of(r), volume_of(holed) - 2.0 * PI * (4.0 + 1.0 / 3.0) * 0.5, 1e-5));
        done(r);
    }
    else
    {
        std::fprintf(stderr, "chamfer of the hole rim failed: %d\n", static_cast<int>(s));
        failures++;
    }
    // tangency-labs/solidcore#379: SmApiChamferFillet (SMLib's constant-distance solver, a
    // chord of 1 across each chamfer) on every edge of the holed cube, the hole's rims too,
    // gives a closed solid, the one the shim's chamfer (constant radius, linear section) gives
    // at the setback a chord of 1 has on these square edges, 1 / sqrt(2).
    {
        SmBrep *all = nullptr;
        CHECK(SmApiBrepCopy(holed->p, all) == SM_SUCCESS && all != nullptr);
        if (all != nullptr)
        {
            CHECK(SmApiChamferFillet(all, 1.0) == SM_SUCCESS);
            smlib_brep *chamfered = new smlib_brep{all};
            CHECK(all->IsManifoldSolid() && sound(chamfered));
            double edges[42];
            int n = 0;
            for (double a : {0.0, 20.0})
            {
                for (double b : {0.0, 20.0})
                {
                    const double rows[3][3] = {{10, a, b}, {a, 10, b}, {a, b, 10}};
                    for (const auto &row : rows)
                    {
                        for (double x : row)
                        {
                            edges[n++] = x;
                        }
                    }
                }
            }
            for (double z : {0.0, 20.0})
            {
                edges[n++] = rim[0];
                edges[n++] = rim[1];
                edges[n++] = z;
            }
            CHECK(smlib_chamfer(holed, edges, 14, 0.1, 1.0 / std::sqrt(2.0), &r, nullptr) == SMLIB_OK);
            if (r != nullptr)
            {
                CHECK(near(volume_of(chamfered), volume_of(r), 1e-8));
                done(r);
            }
            smlib_brep_free(chamfered);
        }
    }
    const double top[3] = {2, 2, 20};
    if (smlib_shell(holed, top, 1, 0.1, -1.0, &r, nullptr) == SMLIB_OK)
    {
        done(r);
    }

    // Every edge of the cube, then a tangent chain.
    double mids[36];
    int k = 0;
    const double ends[4][2] = {{0, 0}, {20, 0}, {0, 20}, {20, 20}};
    for (auto &e : ends)
    {
        double rows[3][3] = {{10, e[0], e[1]}, {e[0], 10, e[1]}, {e[0], e[1], 10}};
        for (auto &row : rows)
        {
            for (double x : row)
            {
                mids[k++] = x;
            }
        }
    }
    CHECK(smlib_fillet(cube, mids, 12, 0.1, 2.0, &r, nullptr) == SMLIB_OK);
    done(r);
    const double vertical[3] = {20, 20, 10};
    smlib_brep *round = nullptr;
    CHECK(smlib_fillet(cube, vertical, 1, 0.1, 5.0, &round, nullptr) == SMLIB_OK && sound(round));
    const double straight[3] = {7.5, 20, 20};
    CHECK(smlib_fillet(round, straight, 1, 0.1, 1.0, &r, nullptr) == SMLIB_OK);
    done(r);
    CHECK(smlib_chamfer(round, straight, 1, 0.1, 1.0, &r, nullptr) == SMLIB_OK);
    done(r);

    // tangency-labs/solidcore#420: a chamfer as wide as the corner rounds it runs round. A box
    // 40 x 40 x 20 with its vertical edges rounded R1, its bottom loop chamfered 1: on each
    // corner the round's arc gets a cone closing to a point on the round's axis (SMLib traced
    // it as an offset intersection, which collapses to that point, and failed). Each straight
    // side loses a 1 x 1 triangle along 38, each corner the triangle turned a quarter about the
    // axis: (pi / 2) * (1 / 3).
    {
        smlib_brep *plate = box(-20, -20, 0, 40, 40, 20);
        const double corners[12] = {20, 20, 10, -20, 20, 10, 20, -20, 10, -20, -20, 10};
        smlib_brep *rounded = nullptr;
        CHECK(smlib_fillet(plate, corners, 4, 0.1, 1.0, &rounded, nullptr) == SMLIB_OK);
        // The whole bottom loop: its four sides and the four arcs of the corner rounds.
        const double h = 19.0 + std::sqrt(0.5);
        const double bottom[24] = {0, -20, 0, 20, 0, 0, 0, 20, 0, -20, 0, 0, h, h, 0, -h, h, 0, h, -h, 0, -h, -h, 0};
        if (rounded != nullptr && smlib_chamfer(rounded, bottom, 8, 0.1, 1.0, &r, nullptr) == SMLIB_OK)
        {
            CHECK(near(volume_of(r), volume_of(rounded) - 4.0 * 38.0 * 0.5 - 4.0 * (PI / 2.0) / 3.0, 1e-7));
            done(r);
        }
        else
        {
            CHECK(!"the bottom loop of the rounded plate chamfered 1");
        }
        smlib_brep_free(rounded);
        smlib_brep_free(plate);
    }

    // tangency-labs/solidcore#417: a box 297 x 100 x 60 with its vertical edges rounded R5, then
    // its top rim rounded R5: on each corner the round of the corner round's arc is a ball
    // centered on that round's axis, which SMLib's torus turned the long way round (the
    // section's rail points, taken where the ball's offset faces collapse onto the axis, carry
    // no angle; on a smaller box it happened to turn the right way). Each straight edge loses
    // (1 - pi/4) 25 along 287 or 90, each corner the square 5 x 5 less the quarter disc, turned
    // a quarter about the axis: (pi/2) (62.5 - 125/3).
    {
        smlib_brep *block = box(0, 0, 0, 297, 100, 60);
        const double corners[12] = {0, 0, 30, 297, 0, 30, 0, 100, 30, 297, 100, 30};
        smlib_brep *rounded = nullptr;
        CHECK(smlib_fillet(block, corners, 4, 0.1, 5.0, &rounded, nullptr) == SMLIB_OK);
        // The whole top rim: its four sides and the four arcs of the corner rounds.
        const double a = 5.0 - 5.0 * std::sqrt(0.5);
        const double rim[24] = {148.5, 0, 60, 297, 50, 60, 148.5, 100, 60, 0, 50, 60, a, a, 60, 297 - a, a, 60, a, 100 - a, 60, 297 - a, 100 - a, 60};
        if (rounded != nullptr && smlib_fillet(rounded, rim, 8, 0.1, 5.0, &r, nullptr) == SMLIB_OK)
        {
            double removed = 754.0 * 25.0 * (1.0 - PI / 4.0) + 4.0 * (PI / 2.0) * (62.5 - 125.0 / 3.0);
            CHECK(near(volume_of(r), volume_of(rounded) - removed, 1e-7));
            done(r);
        }
        else
        {
            CHECK(!"the top rim of the block with R5 corners rounded R5");
        }
        smlib_brep_free(rounded);
        smlib_brep_free(block);
    }

    // tangency-labs/solidcore#416: a box 25 x 25 x 10 with its bottom edges rounded R2, then its
    // vertical edges R2. Each vertical edge's tangent chain runs on down the seam where two
    // bottom rounds meet, between two cylinders as wide as the round: the ball touching both
    // sits where their axes meet, the corner's ball, and SMLib traced the offset cylinders'
    // intersection (their axes, crossing at that one point) and failed. Each vertical edge
    // loses (1 - pi/4) 4 along 8, each corner the octant of the two bottom rounds' crossing
    // (16/3, a Steinmetz solid's eighth) less the corner ball's (4 pi / 3).
    {
        smlib_brep *tray = box(0, 0, 0, 25, 25, 10);
        const double bottom[12] = {12.5, 0, 0, 25, 12.5, 0, 12.5, 25, 0, 0, 12.5, 0};
        smlib_brep *rounded = nullptr;
        CHECK(smlib_fillet(tray, bottom, 4, 0.1, 2.0, &rounded, nullptr) == SMLIB_OK);
        // Each vertical edge and the seam where the two bottom rounds below it meet, which
        // continues it tangentially (a point of it: x = y = z = 2 - sqrt 2 at the origin corner).
        const double s2 = 2.0 - std::sqrt(2.0);
        const double vertical[24] = {0, 0, 6, 25, 0, 6, 25, 25, 6, 0, 25, 6, s2, s2, s2, 25 - s2, s2, s2, 25 - s2, 25 - s2, s2, s2, 25 - s2, s2};
        if (rounded != nullptr && smlib_fillet(rounded, vertical, 8, 0.1, 2.0, &r, nullptr) == SMLIB_OK)
        {
            double removed = 4.0 * (32.0 * (1.0 - PI / 4.0) + 16.0 / 3.0 - 4.0 * PI / 3.0);
            CHECK(near(volume_of(r), volume_of(rounded) - removed, 1e-7));
            done(r);
        }
        else
        {
            CHECK(!"the vertical edges of the tray with R2 bottom rounds rounded R2");
        }
        smlib_brep_free(rounded);
        smlib_brep_free(tray);
    }

    // tangency-labs/solidcore#453: a 20 cube with a hole r5 along y whose wall touches the
    // bottom face along one line, its back rim chamfered 3. With the rim's tangent chain the
    // bottom edge came too, four filleted edges met at the touching point among five, a corner
    // no corner type handles, and the smooth-sequence walk called it as a corner (a crash). No
    // chamfer of that size exists there (the rim has no face below it): refused, every kind.
    {
        SmBrep *bore = nullptr;
        CHECK(SmApiCreateCylinder(SmVector3d(0, 0, 0), 5.0, 20.0, bore) == SM_SUCCESS && bore != nullptr);
        CHECK(SmApiRotate(bore, SmVector3d(0, 0, 0), SmVector3d(1, 0, 0), 90.0) == SM_SUCCESS);
        CHECK(SmApiTranslate(bore, SmVector3d(10, 20, 5)) == SM_SUCCESS);
        smlib_brep tool{bore};
        smlib_brep *holed_low = nullptr;
        CHECK(smlib_boolean(cube, &tool, SMLIB_BOOLEAN_SUBTRACT, &holed_low) == SMLIB_OK);
        CHECK(near(volume_of(holed_low), 8000.0 - PI * 25.0 * 20.0, 1e-6));
        const double back_rim[3] = {10, 20, 10};
        CHECK(smlib_chamfer(holed_low, back_rim, 1, 0.1, 3.0, &r, nullptr) != SMLIB_OK && r == nullptr);
        for (int32_t kind : {SMLIB_BLEND_ROUND, SMLIB_BLEND_CURVATURE, SMLIB_BLEND_CHORD, SMLIB_BLEND_CHAMFER})
        {
            const double size = 3.0;
            int64_t failed = -1;
            CHECK(smlib_blend(holed_low, kind, back_rim, &size, 1, 0.1, &r, &failed) != SMLIB_OK && r == nullptr);
        }
        smlib_brep_free(holed_low);
        delete bore;
    }

    // tangency-labs/solidcore#418: a bar of stadium section (25 x 19 and two half discs r9.5)
    // 14 long along x, a slot x 3.5..9.5 cutting away the y < 0 half of one round end, and the
    // edge between the flat top and the slot floor (6 long) rounded r. At each of its ends the
    // top runs on tangentially onto the round end; the corner past that tangent edge is the
    // slot wall, whose extension the round ends on. SMLib took the tangent edge for the side
    // edge and the round for the end face (1001). Each round takes (1 - pi/4) r^2 along 6, a
    // chamfer of setback 2 the triangle of 2 x 2.
    {
        const int64_t kinds[4] = {SMLIB_CURVE_LINE, SMLIB_CURVE_ARC, SMLIB_CURVE_LINE, SMLIB_CURVE_ARC};
        const double reals[36] = {0, 0,  -9.5, 0, 25, -9.5, 0, 25, 0,   1, 0, 0, 0, 0, -1, 9.5, 0, PI,
                                  0, 25, 9.5,  0, 0,  9.5,  0, 0,  0,   1, 0, 0, 0, 0, 1,  9.5, 0, PI};
        const double along[3] = {1, 0, 0};
        smlib_brep *bar = nullptr;
        CHECK(smlib_extrude(kinds, 4, reals, 36, along, 14.0, 0.0, SMLIB_DRAFT_CORNER_TANGENT, 1, &bar, nullptr) == SMLIB_OK);
        smlib_brep *slot = box(3.5, -12, -10, 6, 12, 20);
        smlib_brep *slotted = nullptr;
        CHECK(smlib_boolean(bar, slot, SMLIB_BOOLEAN_SUBTRACT, &slotted) == SMLIB_OK && sound(slotted));
        const double floor_edge[3] = {6.5, 0, 9.5};
        for (double radius : {0.5, 2.0, 4.0})
        {
            if (slotted != nullptr && smlib_fillet(slotted, floor_edge, 1, 0.1, radius, &r, nullptr) == SMLIB_OK)
            {
                CHECK(near(volume_of(r), volume_of(slotted) - 6.0 * radius * radius * (1.0 - PI / 4.0), 1e-8));
                done(r);
            }
            else
            {
                CHECK(!"the slot floor's top edge of the stadium bar rounded");
            }
        }
        if (slotted != nullptr && smlib_chamfer(slotted, floor_edge, 1, 0.1, 2.0, &r, nullptr) == SMLIB_OK)
        {
            CHECK(near(volume_of(r), volume_of(slotted) - 6.0 * 2.0, 1e-8));
            done(r);
        }
        else
        {
            CHECK(!"the slot floor's top edge of the stadium bar chamfered");
        }
        for (smlib_brep *b : {bar, slot, slotted})
        {
            smlib_brep_free(b);
        }
    }

    // A constant radius law (r -> r) on the top back edge of a 60 x 40 x 20 box whose vertical
    // edge at (60, 40) is rounded R5. The round runs on past the edge's end until it meets the R5
    // cylinder, so the law is evaluated past its end; built with no slope (the shim's law curve
    // had radius in x and 0 in y), it could not be extended there, and the round failed (1001).
    // Removed: (1 - pi/4) r^2 along 55, and past x = 55 the square r x r less the quarter disc,
    // as far as the cylinder reaches: the integral of (r - sqrt(r^2 - u^2)) sqrt(25 - (5 - r +
    // u)^2) over 0 < u < r.
    {
        smlib_brep *block = box(0, 0, 0, 60, 40, 20);
        const double corner_edge[3] = {60, 40, 10};
        smlib_brep *rounded = nullptr;
        CHECK(smlib_fillet(block, corner_edge, 1, 0.1, 5.0, &rounded, nullptr) == SMLIB_OK);
        const double back_top[3] = {27.5, 40, 20};
        const double from[3] = {0, 40, 20};
        const double knots[4] = {0, 0, 1, 1};
        for (double radius : {1.0, 2.0, 3.0})
        {
            const double law[2] = {radius, radius};
            if (rounded != nullptr && smlib_fillet_variable(rounded, back_top, 0.1, from, 1, knots, 4, law, 2, &r) == SMLIB_OK)
            {
                double past = 0.0;
                const int steps = 200000;
                for (int i = 0; i < steps; ++i)
                {
                    double u = (i + 0.5) * radius / steps;
                    double w = 5.0 - radius + u;
                    past += (radius - std::sqrt(radius * radius - u * u)) * std::sqrt(std::max(0.0, 25.0 - w * w)) * radius / steps;
                }
                CHECK(near(volume_of(r), volume_of(rounded) - 55.0 * radius * radius * (1.0 - PI / 4.0) - past, 1e-6));
                done(r);
            }
            else
            {
                CHECK(!"a constant law into the R5 corner");
            }
        }
        smlib_brep_free(rounded);
        smlib_brep_free(block);
    }

    // tangency-labs/solidcore#458: a variable round whose law SMLib cannot evaluate past the edge's
    // ends (a constant law curve laid out with no slope, radius 2 as x and 0 as y, as the shim
    // once built it; extending a curve of zero speed fails). SmVariableRadiusFS::SetupOffsetValues
    // and SmTemporaryRadiusChange::SetUp ignored the failed evaluation and offset the surfaces by
    // an uninitialised radius (it went on to fail, so only a memory checker sees the read). Now
    // the evaluation's failure stops it: the round fails, or is the round of radius 2.
    {
        SmBrep *b = nullptr;
        CHECK(SmApiBrepCopy(cube->p, b) == SM_SUCCESS && b != nullptr);
        SmTArray<SmEdge *> all;
        b->GetEdges(all);
        SmTArray<SmEdge *> top_edge;
        for (ULONG i = 0; i < all.GetSize(); ++i)
        {
            SmVertex *v0 = nullptr;
            SmVertex *v1 = nullptr;
            all[i]->GetVertices(v0, v1);
            SmPoint3d m = (v0->GetPoint() + v1->GetPoint()) * 0.5;
            if (std::fabs(m.x - 10) < 1e-9 && std::fabs(m.y - 20) < 1e-9 && std::fabs(m.z - 20) < 1e-9)
            {
                top_edge.Add(all[i]);
            }
        }
        CHECK(top_edge.GetSize() == 1);
        const SmContext &ctx = *SmApiGetOrCreateContext();
        SmTArray<SmPoint3d> cps;
        cps.Add(SmPoint3d(2, 0, 0));
        cps.Add(SmPoint3d(2, 0, 0));
        SmTArray<ULONG> mults;
        mults.Add(2);
        mults.Add(2);
        SmTArray<double> unique;
        unique.Add(0.0);
        unique.Add(1.0);
        SmBSplineCurve *flat = nullptr;
        CHECK(SmBSplineCurve::CreateCanonical(ctx, 3, 1, cps, SM_CF_UNSPECIFIED, mults, unique, SM_KT_UNSPECIFIED, nullptr, nullptr,
                                              flat) == SM_SUCCESS && flat != nullptr);
        if (top_edge.GetSize() == 1 && flat != nullptr)
        {
            SmBSplineFilletLaw law(flat, SmExtent1d(0.0, 1.0));
            SmTArray<SmFilletLaw *> laws;
            laws.Add(&law);
            SmTArray<double> values;
            values.Add(2.0);
            SmTArray<SmFilletSurfaceGeneratorType> sections;
            sections.Add(SM_FSG_CIRCULAR);
            SmTArray<SmFilletSolverType> solvers;
            solvers.Add(SM_FS_VARIABLE_RADIUS);
            SmFilletExecutive exec(ctx, b);
            SmStatus st = exec.SetFilletParameters(b, top_edge, sections, solvers, values, &laws, -1.0, -1.0, -1.0, -1.0, 1, 1.0);
            if (st == SM_SUCCESS)
            {
                st = exec.DoFillet();
            }
            if (st == SM_SUCCESS && b->StitchAndOrient() == SM_SUCCESS)
            {
                smlib_brep got{b};
                CHECK(near(volume_of(&got), 8000.0 - 4.0 * 20.0 * (1.0 - PI / 4.0), 1e-6));
            }
        }
        delete b;
    }

    // tangency-labs/solidcore#430: a variable round whose radius law has a corner. SMLib traced
    // the round smoothly across the corner, averaging the tangents of its two sides, and skinned
    // one smooth surface over it: it missed the law past the corner and failed its own check.
    // The trace now ends a step at the corner with the tangents of each side, the surface is
    // skinned piece by piece, and the face is split there. On the cube's edge y = z = 20, the
    // round takes (1 - pi/4) r^2 per unit length; for 4 -> 2 -> 4 that is the same as for a
    // straight 4 -> 2 along the whole edge.
    {
        const double PI4 = 1.0 - PI / 4.0;
        const double edge[3] = {10, 20, 20};
        const double from0[3] = {0, 20, 20};
        const double from1[3] = {20, 20, 20};
        const double v_knots[5] = {0, 0, 0.5, 1, 1};
        const double v_radii[3] = {4, 2, 4};
        const double removed = PI4 * 20.0 * (16.0 + 8.0 + 4.0) / 3.0;
        for (const double *from : {from0, from1})
        {
            CHECK(smlib_fillet_variable(cube, edge, 0.1, from, 1, v_knots, 5, v_radii, 3, &r) == SMLIB_OK);
            if (r != nullptr)
            {
                CHECK(near(8000.0 - volume_of(r), removed, 1e-4));
                done(r);
            }
        }
        // Two corners (4, 2, 4, 2 at 0, 0.25, 0.75, 1): the straight 4 -> 2 in three pieces.
        const double w_knots[6] = {0, 0, 0.25, 0.75, 1, 1};
        const double w_radii[4] = {4, 2, 4, 2};
        CHECK(smlib_fillet_variable(cube, edge, 0.1, from0, 1, w_knots, 6, w_radii, 4, &r) == SMLIB_OK);
        if (r != nullptr)
        {
            CHECK(near(8000.0 - volume_of(r), removed, 1e-4));
            done(r);
        }
        // The same law as a quadratic with a corner at 0.5 (the knot repeated twice): 4, 3, 2
        // then 2, 3, 4, evenly spaced, are the straight halves again.
        const double q_knots[8] = {0, 0, 0, 0.5, 0.5, 1, 1, 1};
        const double q_radii[5] = {4, 3, 2, 3, 4};
        CHECK(smlib_fillet_variable(cube, edge, 0.1, from0, 2, q_knots, 8, q_radii, 5, &r) == SMLIB_OK);
        if (r != nullptr)
        {
            CHECK(near(8000.0 - volume_of(r), removed, 1e-4));
            done(r);
        }
        // A radius of zero inside the edge would pinch the round to a point: refused.
        const double p_radii[3] = {4, 0, 4};
        CHECK(smlib_fillet_variable(cube, edge, 0.1, from0, 1, v_knots, 5, p_radii, 3, &r) == SMLIB_ERR_UNSUPPORTED && r == nullptr);
    }

    // A variable round 2 -> r1 on the top edge x = 10, z = 30 of a 10 x 20 x 30 box. The top is
    // 10 wide: up to r1 = 10 the round fits, and takes (1 - pi/4) L (a^2 + ab + b^2) / 3. Past
    // it the contact line on the top leaves it across its far edge x = 0 before the round's end:
    // the round is too large for the top and is refused, from either end, as is a constant R12.
    // (The corner without an end rib, for a short side edge, had built r1 = 20 and 30 across the
    // top's far edge as valid bodies.)
    {
        smlib_brep *slab = box(0, 0, 0, 10, 20, 30);
        const double edge[3] = {10, 10, 30};
        const double from0[3] = {10, 0, 30};
        const double from1[3] = {10, 20, 30};
        const double knots[4] = {0, 0, 1, 1};
        const double fits[2] = {2, 10};
        CHECK(smlib_fillet_variable(slab, edge, 0.1, from0, 1, knots, 4, fits, 2, &r) == SMLIB_OK);
        if (r != nullptr)
        {
            CHECK(near(6000.0 - volume_of(r), (1.0 - PI / 4.0) * 20.0 * (4.0 + 20.0 + 100.0) / 3.0, 1e-4));
            done(r);
        }
        for (double r1 : {11.0, 20.0, 30.0, 50.0})
        {
            const double law[2] = {2, r1};
            const double back[2] = {r1, 2};
            CHECK(smlib_fillet_variable(slab, edge, 0.1, from0, 1, knots, 4, law, 2, &r) == SMLIB_ERR_NO_ROOM && r == nullptr);
            CHECK(smlib_fillet_variable(slab, edge, 0.1, from1, 1, knots, 4, back, 2, &r) == SMLIB_ERR_NO_ROOM && r == nullptr);
        }
        CHECK(smlib_fillet(slab, edge, 1, 0.1, 12.0, &r, nullptr) == SMLIB_ERR_NO_ROOM && r == nullptr);
        smlib_brep_free(slab);
    }

    // A chamfer of setback 2 on the top edge y = 0 of a prism 20 high on (0,0) (40,0)
    // (40 - 20s, 20) (20s, 20): its end faces slope by s, so the chamfer's contact line on the
    // top leaves it across an end face's top edge (the top's side edge at that end) before the
    // cross section at the edge's end. That is past the end, where the end face cuts the
    // chamfer off, not a blend too large for the top: the triangle 2 x 2 along 40 - 4s/3 comes
    // off.
    for (double s : {-0.5, 0.5})
    {
        const double xy[8] = {0, 0, 40, 0, 40 - 20 * s, 20, 20 * s, 20};
        int64_t kinds[4];
        double reals[24];
        for (int i = 0; i < 4; ++i)
        {
            const int j = (i + 1) % 4;
            const double seg[6] = {xy[2 * i], xy[2 * i + 1], 0, xy[2 * j], xy[2 * j + 1], 0};
            kinds[i] = SMLIB_CURVE_LINE;
            std::copy(seg, seg + 6, reals + 6 * i);
        }
        const double up[3] = {0, 0, 1};
        smlib_brep *prism = nullptr;
        CHECK(smlib_extrude(kinds, 4, reals, 24, up, 20.0, 0.0, SMLIB_DRAFT_CORNER_CURVATURE, 1, &prism, nullptr) == SMLIB_OK);
        const double edge[3] = {20, 0, 20};
        if (prism != nullptr && smlib_chamfer(prism, edge, 1, 0.1, 2.0, &r, nullptr) == SMLIB_OK)
        {
            CHECK(near(volume_of(prism) - volume_of(r), 2.0 * (40.0 - 4.0 * s / 3.0), 1e-8));
            done(r);
        }
        else
        {
            CHECK(!"a chamfer between sloped end faces");
        }
        smlib_brep_free(prism);
    }

    // A 100 x 60 x 30 slab less a basin lofted from 20..80 x 4..50 at the top to 30..70 x 20..40
    // at z = 10, its rim edge y = 4 chamfered 2. Its contact line on the top runs on over the top
    // past the rim's concave corners and leaves it at the slab's sides x = 0 and 100, beyond the
    // cross sections at the edge's ends; the one on the sloped wall leaves it across the wall's
    // sloped corner edges, before them. Both are past the ends. The section's drop below the top
    // is 2 / sqrt(1 + 0.8^2) and it comes off along 60 - drop / 3.
    {
        smlib_brep *slab = box(0, 0, 0, 100, 60, 30);
        const double rims[2][5] = {{20, 4, 80, 50, 30}, {30, 20, 70, 40, 10}};
        int64_t kinds[8];
        double reals[48];
        for (int s = 0; s < 2; ++s)
        {
            const double *q = rims[s];
            const double xy[4][2] = {{q[0], q[1]}, {q[2], q[1]}, {q[2], q[3]}, {q[0], q[3]}};
            for (int i = 0; i < 4; ++i)
            {
                const int j = (i + 1) % 4;
                const double seg[6] = {xy[i][0], xy[i][1], q[4], xy[j][0], xy[j][1], q[4]};
                kinds[4 * s + i] = SMLIB_CURVE_LINE;
                std::copy(seg, seg + 6, reals + 6 * (4 * s + i));
            }
        }
        const int64_t sizes[2] = {4, 4};
        smlib_brep *pocket = nullptr;
        smlib_brep *basin = nullptr;
        CHECK(smlib_loft(kinds, 8, reals, 48, sizes, 2, 1, 1, &pocket) == SMLIB_OK);
        CHECK(pocket != nullptr && smlib_boolean(slab, pocket, SMLIB_BOOLEAN_SUBTRACT, &basin) == SMLIB_OK);
        const double rim[3] = {50, 4, 30};
        const double drop = 2.0 / std::sqrt(1.0 + 0.64);
        if (basin != nullptr && smlib_chamfer(basin, rim, 1, 0.1, 2.0, &r, nullptr) == SMLIB_OK)
        {
            CHECK(near(volume_of(r), volume_of(basin) - drop * (60.0 - drop / 3.0), 1e-8));
            done(r);
        }
        else
        {
            CHECK(!"a chamfer on the sloped basin's rim");
        }
        for (smlib_brep *b : {slab, pocket, basin})
        {
            smlib_brep_free(b);
        }
    }

    // SolidCraft ladder 07-fillets/016: a 60 x 40 x 20 box, its top back edge rounded R5, its top
    // right edge R3, then the vertical edge x = 60, y = 40 (z 0..15). The R3 round leaves the R5
    // round's end arc on the x = 60 face 0.46 long; the vertical round's contact line there rolls
    // over the R3 round's tangent edge past the round's end and never meets the arc. The corner
    // has no end rib and the round is cut off where it meets the body: removed, outside the
    // vertical round, the integral of min(R5, R3 top) over the corner column.
    {
        smlib_brep *block = box(0, 0, 0, 60, 40, 20);
        const double back_top[3] = {30, 40, 20};
        const double right_top[3] = {60, 20, 20};
        const double vertical[3] = {60, 40, 7.5};
        smlib_brep *r5 = nullptr;
        smlib_brep *r3 = nullptr;
        CHECK(smlib_fillet(block, back_top, 1, 0.1, 5.0, &r5, nullptr) == SMLIB_OK);
        CHECK(r5 != nullptr && smlib_fillet(r5, right_top, 1, 0.1, 3.0, &r3, nullptr) == SMLIB_OK);
        const double radii[3] = {7.0, 4.0, 2.0};
        const double volumes[3] = {47418.926881, 47546.786538, 47591.811067};
        for (int i = 0; i < 3 && r3 != nullptr; ++i)
        {
            if (smlib_fillet(r3, vertical, 1, 0.1, radii[i], &r, nullptr) == SMLIB_OK)
            {
                CHECK(near(volume_of(r), volumes[i], 1e-8));
                done(r);
            }
            else
            {
                CHECK(!"ladder 016's vertical round under the short R5 end arc");
            }
        }
        for (smlib_brep *b : {block, r5, r3})
        {
            smlib_brep_free(b);
        }
    }

    // tangency-labs/solidcore#423: a sphere r2 in six patches (each on its own sphere surface),
    // cut flat at z = 1, its rim rounded R1. The clean-up after filleting removed every
    // topological edge in the body, the seams between the sphere's own patches too, and merged
    // the five patches into one face that no patch's parameterization could carry: the result
    // failed the check (loop orientation, rails half a millimetre off). It now leaves the
    // target's own faces apart. The round's section is the square 1 x 1 less the quarter disc,
    // centred on the circle of radius 1 at z = 0: 2 pi - pi^2 / 2 comes off the 9 pi.
    {
        smlib_brep *ball = cube_sphere(2.0);
        CHECK(ball != nullptr && sound(ball) && near(volume_of(ball), 32.0 * PI / 3.0, 1e-6));
        smlib_brep *lid = box(-3, -3, 1, 6, 6, 3);
        smlib_brep *cut = nullptr;
        CHECK(smlib_boolean(ball, lid, SMLIB_BOOLEAN_SUBTRACT, &cut) == SMLIB_OK && sound(cut));
        // The rim's four arcs, one on each side patch.
        const double q = std::sqrt(3.0);
        const double rim[12] = {q, 0, 1, 0, q, 1, -q, 0, 1, 0, -q, 1};
        if (cut != nullptr && smlib_fillet(cut, rim, 4, 0.5, 1.0, &r, nullptr) == SMLIB_OK)
        {
            CHECK(near(volume_of(r), 7.0 * PI + PI * PI / 2.0, 5e-5));
            done(r);
        }
        else
        {
            CHECK(!"the rim of the cut sphere in patches rounded R1");
        }
        for (smlib_brep *b : {ball, lid, cut})
        {
            smlib_brep_free(b);
        }
    }

    // tangency-labs/solidcore#422: a plate 40 x 20 x 10 with a through hole r3 at (20, 14.5),
    // its top edge y = 20 rounded R4. The round's contact line on the top (y = 16) runs across
    // the hole: the round is larger than the top face allows there, and is refused
    // (SMLIB_ERR_NO_ROOM) rather than built across the hole.
    {
        smlib_brep *plate = box(0, 0, 0, 40, 20, 10);
        smlib_brep *bore = cylinder(20, 14.5, -1, 3, 12);
        smlib_brep *holed_plate = nullptr;
        CHECK(smlib_boolean(plate, bore, SMLIB_BOOLEAN_SUBTRACT, &holed_plate) == SMLIB_OK);
        const double edge[3] = {5, 20, 10};
        CHECK(holed_plate != nullptr && smlib_fillet(holed_plate, edge, 1, 0.1, 4.0, &r, nullptr) == SMLIB_ERR_NO_ROOM && r == nullptr);
        for (smlib_brep *b : {plate, bore, holed_plate})
        {
            smlib_brep_free(b);
        }
    }

    // A cylinder r10 h20, its top rim chamfered or rounded as wide as the top: each blend's rail
    // on the top is the top's centre, so the blend would take the whole top and end in a point
    // (the chamfer a cone, the round a sphere), and is refused (SMLIB_ERR_NO_ROOM) in every form:
    // the equal chamfer, by smlib_blend too, two distances with either face first, a distance
    // and an angle, the round. SMLib built them, the top gone (SolidCraft's
    // chamfer_sides_curved_tests::too_large_or_unsupported_is_refused). 9.99 fits and takes the
    // triangle (the round: the square less the quarter disc) of its section revolved about the
    // axis (Pappus).
    {
        smlib_brep *can = cylinder(0, 0, 0, 10, 20);
        const double whole = PI * 100.0 * 20.0;
        const double at = 100.0 * PI / 180.0;
        const double top_rim[3] = {10.0 * std::cos(at), 10.0 * std::sin(at), 20.0};
        const double up[3] = {0, 0, 1};
        const double down[3] = {0, 0, -1};
        // The ring the triangle with corners (r, z) sweeps about the axis.
        auto ring = [](double r0, double z0, double r1, double z1, double r2, double z2)
        {
            double area = 0.5 * std::fabs((r1 - r0) * (z2 - z0) - (r2 - r0) * (z1 - z0));
            return 2.0 * PI * (r0 + r1 + r2) / 3.0 * area;
        };
        for (double d : {10.0, 9.99})
        {
            auto expect = [&](const char *what, int32_t s, double removed)
            {
                if (d < 10.0)
                {
                    if (s != SMLIB_OK || !near(volume_of(r), whole - removed, 1e-7))
                    {
                        std::fprintf(stderr, "  %s %g: status %d\n", what, d, static_cast<int>(s));
                        failures++;
                    }
                    done(r);
                }
                else if (s != SMLIB_ERR_NO_ROOM || r != nullptr)
                {
                    std::fprintf(stderr, "  %s %g: status %d, not refused for want of room\n", what, d, static_cast<int>(s));
                    failures++;
                    smlib_brep_free(r);
                    r = nullptr;
                }
            };
            const double equal = ring(10, 20, 10 - d, 20, 10, 20 - d);
            const double sides = ring(10, 20, 10 - d, 20, 10, 17);
            const double slope = std::tan(PI / 6.0);
            expect("chamfer", smlib_chamfer(can, top_rim, 1, 0.1, d, &r, nullptr), equal);
            expect("blend chamfer", smlib_blend(can, SMLIB_BLEND_CHAMFER, top_rim, &d, 1, 0.1, &r, nullptr), equal);
            expect("chamfer across the top, 3 down", smlib_chamfer_sides(can, top_rim, 1, 0.1, d, 3.0, 0, up, &r, nullptr), sides);
            expect("chamfer 3 down, across the top", smlib_chamfer_sides(can, top_rim, 1, 0.1, 3.0, d, 0, down, &r, nullptr), sides);
            expect("chamfer across the top at 30 degrees", smlib_chamfer_sides(can, top_rim, 1, 0.1, d, PI / 6.0, 1, up, &r, nullptr),
                   ring(10, 20, 10 - d, 20, 10, 20 - d * slope));
            expect("chamfer 3 down at 60 degrees", smlib_chamfer_sides(can, top_rim, 1, 0.1, d / std::tan(PI / 3.0), PI / 3.0, 1, down, &r, nullptr),
                   ring(10, 20, 10 - d, 20, 10, 20 - d / std::tan(PI / 3.0)));
            // The round's section: the square d x d at the rim less the quarter disc, its
            // centroid d (10 - 3 pi) / (12 - 3 pi) in from the rim.
            expect("fillet", smlib_fillet(can, top_rim, 1, 0.1, d, &r, nullptr),
                   2.0 * PI * (10.0 - d * (10.0 - 3.0 * PI) / (12.0 - 3.0 * PI)) * (1.0 - PI / 4.0) * d * d);
        }
        smlib_brep_free(can);
    }

    // tangency-labs/solidcore#419: an offset surface evaluated beyond its
    // base surface's domain extends the base, and the extension's shape depends on its length.
    // It was made to reach just past the first point evaluated there, and extended again from
    // itself for the next, so two copies of one offset surface evaluated in another order (the
    // fillet solver's, which finds a round's end points, and the trace's) gave points 0.03 mm
    // apart at the same parameters, and the trace refused the end points. Now the extension's
    // reach depends only on how far out the point lies. A half cylinder r5 offset by 4, beyond
    // both ends of its arc: once straight at two points, once after two nearer points first.
    {
        const SmContext &ctx = *SmApiGetOrCreateContext();
        SmCylinder *half = new (ctx) SmCylinder(SmPoint3d(0, 0, 0), SmVector3d(1, 0, 0), SmVector3d(0, 1, 0), 5.0, 0.0, 180.0, 10.0,
                                                FALSE, FALSE, FALSE, &ctx);
        SmOffsetSurface *first = new (ctx) SmOffsetSurface(4.0, *half);
        SmOffsetSurface *second = new (ctx) SmOffsetSurface(4.0, *half);
        first->SetAllowExtension(TRUE);
        second->SetAllowExtension(TRUE);
        SmExtent2d d = half->GetNaturalUVDomain();
        double length = d.GetMax().x - d.GetMin().x;
        double v = 0.5 * (d.GetMin().y + d.GetMax().y);
        SmPoint3d p;
        const double nearer[2] = {d.GetMax().x + 0.067 * length, d.GetMin().x - 0.067 * length};
        for (double u : nearer)
        {
            CHECK(first->EvaluatePoint(SmPoint2d(u, v), p) == SM_SUCCESS);
        }
        const double farther[2] = {d.GetMax().x + 0.116 * length, d.GetMin().x - 0.116 * length};
        for (double u : farther)
        {
            SmPoint3d a;
            SmPoint3d b;
            CHECK(second->EvaluatePoint(SmPoint2d(u, v), b) == SM_SUCCESS);
            CHECK(first->EvaluatePoint(SmPoint2d(u, v), a) == SM_SUCCESS);
            CHECK(a.DistanceBetween(b) < 1e-9);
        }
        delete first;
        delete second;
        delete half;
    }

    // Booleans that split, vanish or touch.
    smlib_brep *cut = box(8, -1, -1, 4, 22, 22);
    CHECK(smlib_boolean(cube, cut, SMLIB_BOOLEAN_SUBTRACT, &r) == SMLIB_OK);
    done(r);
    smlib_brep *big = box(-5, -5, -5, 30, 30, 30);
    CHECK(smlib_boolean(cube, big, SMLIB_BOOLEAN_SUBTRACT, &r) == SMLIB_ERR_EMPTY);
    smlib_brep *far = box(50, 0, 0, 10, 10, 10);
    CHECK(smlib_boolean(cube, far, SMLIB_BOOLEAN_UNION, &r) == SMLIB_OK);
    done(r);
    smlib_brep *touching = box(20, 0, 0, 10, 20, 20);
    if (smlib_boolean(cube, touching, SMLIB_BOOLEAN_UNION, &r) == SMLIB_OK)
    {
        CHECK(near(volume_of(r), 12000.0, 1e-9));
        done(r);
    }

    // History through the boolean shim: tagged operands in, the result's history out. A slot
    // across the top of the cube: every result face comes from a face of one operand, the
    // slot's walls and floor from the tool (operand 1), the rest from the cube (operand 0).
    {
        smlib_brep *slot = box(5, -1, 10, 10, 22, 15);
        smlib_brep *ta = nullptr;
        smlib_brep *tb = nullptr;
        CHECK(smlib_history_tag(cube, 0, &ta) == SMLIB_OK);
        CHECK(smlib_history_tag(slot, 1, &tb) == SMLIB_OK);
        CHECK(smlib_boolean(ta, tb, SMLIB_BOOLEAN_SUBTRACT, &r) == SMLIB_OK);
        int64_t counts[3] = {0, 0, 0};
        CHECK(smlib_brep_counts(r, counts) == SMLIB_OK && counts[0] == 10);
        int64_t len = 0;
        CHECK(smlib_history_read(r, nullptr, 0, &len) == SMLIB_ERR_BUFFER && len > 3);
        std::vector<int64_t> h(static_cast<size_t>(len));
        CHECK(smlib_history_read(r, h.data(), len, &len) == SMLIB_OK);
        CHECK(h[0] == SMLIB_HISTORY_VERSION && h[1] == counts[0] && h[2] == counts[1]);
        int from[2] = {0, 0};
        size_t at = 3;
        for (int64_t f = 0; f < h[1] && at + 1 < h.size(); ++f)
        {
            int64_t n = h[at + 1];
            CHECK(n >= 1 && h[at] != SMLIB_CHANGE_GENERATED);
            for (int64_t i = 0; i < n && at + 3 + 2 * i < h.size(); ++i)
            {
                int64_t operand = h[at + 2 + 2 * i];
                CHECK(operand == 0 || operand == 1);
                from[operand == 1]++;
            }
            at += 2 + 2 * static_cast<size_t>(n);
        }
        // The cube's 6 faces give 7 (the top is cut in two); the tool gives the 3 slot faces.
        CHECK(from[0] == 7 && from[1] == 3);
        done(r);
        for (smlib_brep *b : {slot, ta, tb})
        {
            smlib_brep_free(b);
        }
    }

    // SMLib defects 2-5 (curves of a failed imprint freed under their Edge, a failed imprint
    // carried on, a squeezed vertex reused, a combined-away ProtoEdge left in its ProtoFaces):
    // inward walls as thick as the cube's inradius or thicker, refused before SMLib sees them
    // (no cavity can remain).
    const double lid[3] = {10, 10, 20};
    CHECK(smlib_shell(cube, lid, 1, 0.1, -10.0, &r, nullptr) == SMLIB_ERR_NO_ROOM && r == nullptr);
    CHECK(smlib_shell(cube, lid, 1, 0.1, -15.0, &r, nullptr) == SMLIB_ERR_NO_ROOM);
    CHECK(smlib_shell(cube, nullptr, 0, 0.1, -10.0, &r, nullptr) == SMLIB_ERR_NO_ROOM);
    if (smlib_shell(cube, lid, 1, 0.1, -9.5, &r, nullptr) == SMLIB_OK)
    {
        CHECK(near(volume_of(r), 8000.0 - 1.0 * 1.0 * 10.5, 1e-6));
        done(r);
    }
    // tangency-labs/solidcore#380: SMLib itself, past the shim's guard. A closed inward shell
    // whose offset leaves nothing handed the cube back unhollowed as a success; it fails now,
    // open or closed, and leaves the input as it was.
    {
        SmTArray<SmFace *> faces;
        cube->p->GetFaces(faces);
        SmTArray<const SmFace *> lidface;
        for (ULONG i = 0; i < faces.GetSize(); ++i)
        {
            SmPoint3d lo;
            SmPoint3d hi;
            if (SmApiFaceBoundingBox(faces[i], TRUE, lo, hi) == SM_SUCCESS && lo[2] > 19.9)
            {
                lidface.Add(faces[i]);
            }
        }
        CHECK(lidface.GetSize() == 1);
        SmTArray<const SmFace *> none;
        for (double wall : {10.0, 12.0, 15.0})
        {
            for (const SmTArray<const SmFace *> *open : {&none, &lidface})
            {
                SmBrep *res = nullptr;
                CHECK(SmApiShellBrepFull(cube->p, -wall, TRUE, TRUE, TRUE, *open, res) != SM_SUCCESS && res == nullptr);
                delete res;
            }
        }
        CHECK(near(volume_of(cube), 8000.0, 1e-9) && sound(cube));
    }

    // SMLib defect 6 (SmOffsetExecutive::RemoveExcessFaces deleted faces its face maps still
    // held, which PiecewiseMerge then read), found by fuzz_ops: a closed inward shell of the
    // cube with two of its edges chamfered 0.1 mm. Any status, never a crash.
    {
        const double two[6] = {0, 0, 10, 0, 20, 10};
        smlib_brep *nicked = nullptr;
        CHECK(smlib_chamfer(cube, two, 2, 0.01, 0.1, &nicked, nullptr) == SMLIB_OK && sound(nicked));
        if (smlib_shell(nicked, nullptr, 0, 0.01, -2.624267578125, &r, nullptr) == SMLIB_OK)
        {
            done(r);
        }
        smlib_brep_free(nicked);
    }

    // tangency-labs/solidcore#380: the guard against walls with no room inside is the body's
    // own inradius, not its box. An L of two 6 mm legs, 40 long and 20 tall: its box allows 10,
    // its inradius is 6 * sqrt(2) / (1 + sqrt(2)) = 3.51 (at the inner corner).
    {
        smlib_brep *legs[2] = {box(0, 0, 0, 40, 6, 20), box(0, 0, 0, 6, 40, 20)};
        smlib_brep *ell = nullptr;
        CHECK(smlib_boolean(legs[0], legs[1], SMLIB_BOOLEAN_UNION, &ell) == SMLIB_OK && sound(ell));
        for (double wall : {5.0, 4.0, 3.6})
        {
            CHECK(smlib_shell(ell, nullptr, 0, 0.1, -wall, &r, nullptr) == SMLIB_ERR_NO_ROOM && r == nullptr);
        }
        // Open at the top, the wall only needs room from the sides and the floor: the same.
        const double top[3] = {20, 3, 20};
        CHECK(smlib_shell(ell, top, 1, 0.1, -3.6, &r, nullptr) == SMLIB_ERR_NO_ROOM);
        if (smlib_shell(ell, top, 1, 0.1, -1.0, &r, nullptr) == SMLIB_OK)
        {
            // Floor 1 thick under a cavity 19 tall, its plan the L inset by 1: 2 * 38 * 4 - 4 * 4.
            CHECK(std::fabs(volume_of(r) - (volume_of(ell) - 19.0 * (2.0 * 38.0 * 4.0 - 16.0))) < 1e-4);
            done(r);
        }
        else
        {
            CHECK(!"the L opened at the top, 1 mm walls");
        }
        smlib_brep_free(ell);
        smlib_brep_free(legs[0]);
        smlib_brep_free(legs[1]);
    }

    // tangency-labs/solidcore#403: a hemisphere r10 (a sphere cut by a box, its rim two arcs of
    // 270 and 90 degrees) shelled 2 open at the flat face. The rim's offset is a closed circle,
    // which both arcs were trimmed to as the same 270 degrees: the flat face got the long arc
    // twice and the cavity was no solid. Each arc now takes the part of the circle it runs along.
    {
        const double centre[3] = {0, 0, 0};
        smlib_brep *ball = nullptr;
        CHECK(smlib_make_sphere(centre, 10.0, &ball) == SMLIB_OK);
        smlib_brep *upper = box(-20, -20, 0, 40, 40, 20);
        smlib_brep *half = nullptr;
        CHECK(smlib_boolean(ball, upper, SMLIB_BOOLEAN_INTERSECT, &half) == SMLIB_OK && sound(half));
        const double flat[3] = {3, 2, 0};
        CHECK(smlib_shell(half, flat, 1, 0.1, -2.0, &r, nullptr) == SMLIB_OK);
        if (r != nullptr)
        {
            CHECK(near(volume_of(r), 2.0 / 3.0 * PI * (1000.0 - 512.0), 1e-5));
            done(r);
        }
        for (smlib_brep *b : {ball, upper, half})
        {
            smlib_brep_free(b);
        }
    }

    // tangency-labs/solidcore#403: the same shell gives the same answer every time, wherever the
    // heap puts SMLib's objects. SMLib's pointer maps (SmMapPtrToPtr, SmRelation) handed out
    // their keys in hash order, which is address order, so the offset processed edges and faces
    // in an order that changed from call to call: a 47.5 x 45 x 30 box joined with a cylinder
    // r7.5 along x whose seam lies on the box's front face, shelled 1 open at the bottom,
    // failed with 1001 or 1006 by turns. Heap blocks of varying size between the calls move
    // the objects about.
    {
        smlib_brep *block = box(0, 0, 0, 47.5, 45, 30);
        SmBrep *bar = nullptr;
        CHECK(SmApiCreateCylinder(SmVector3d(0, 0, 0), 7.5, 30.0, bar) == SM_SUCCESS && bar != nullptr);
        CHECK(SmApiRotate(bar, SmVector3d(0, 0, 0), SmVector3d(0, 1, 0), 90.0) == SM_SUCCESS);
        CHECK(SmApiTranslate(bar, SmVector3d(47.5, 0, 10)) == SM_SUCCESS);
        smlib_brep side{bar};
        smlib_brep *joined = nullptr;
        CHECK(smlib_boolean(block, &side, SMLIB_BOOLEAN_UNION, &joined) == SMLIB_OK && sound(joined));
        const double bottom[3] = {23.75, 22.5, 0};
        std::vector<std::vector<char>> spacers;
        int32_t first = 0;
        double first_volume = 0.0;
        for (int i = 0; i < 6; ++i)
        {
            spacers.emplace_back(static_cast<size_t>(16 + 97 * i * i));
            int32_t got = smlib_shell(joined, bottom, 1, 0.1, -1.0, &r, nullptr);
            double v = got == SMLIB_OK ? volume_of(r) : 0.0;
            if (i == 0)
            {
                first = got;
                first_volume = v;
            }
            CHECK(got == first && v == first_volume);
            smlib_brep_free(r);
            r = nullptr;
        }
        delete bar;
        for (smlib_brep *b : {block, joined})
        {
            smlib_brep_free(b);
        }
    }

    // tangency-labs/solidcore#397: an inward shell of a box 30 x 20 x 10 whose edges are rounded
    // thinner than the 2 mm wall. Each round's inward offset collapses (SMLib's offset surface is
    // NULL), which failed the shell (1001); the faces beside it are now mated across it, so the
    // inner walls meet at sharp edges: the cavity is the box shrunk by the wall on every side
    // (but the open top), whatever the radius. All twelve edges rounded (the corner balls
    // collapse too) and only the four vertical ones (their ends meet the top and bottom faces,
    // where the offset edges along them shrink to a corner). A round thicker than the wall (r 3)
    // still shells as before, to a cavity rounded r - t.
    {
        const double a = 30, b = 20, c = 10, t = 2;
        // a x b x c with all edges rounded (corners balls), as the core box grown by the radius
        auto rounded = [](double x, double y, double z, double rad)
        {
            double cx = x - 2 * rad, cy = y - 2 * rad, cz = z - 2 * rad;
            return cx * cy * cz + 2 * rad * (cx * cy + cy * cz + cx * cz) + PI * rad * rad * (cx + cy + cz) + 4.0 / 3.0 * PI * rad * rad * rad;
        };
        std::vector<double> all, vertical;
        for (double y : {0.0, b})
        {
            for (double z : {0.0, c})
            {
                all.insert(all.end(), {a / 2, y, z});
            }
        }
        for (double x : {0.0, a})
        {
            for (double z : {0.0, c})
            {
                all.insert(all.end(), {x, b / 2, z});
            }
        }
        for (double x : {0.0, a})
        {
            for (double y : {0.0, b})
            {
                all.insert(all.end(), {x, y, c / 2});
                vertical.insert(vertical.end(), {x, y, c / 2});
            }
        }
        const double top[3] = {a / 2, b / 2, c};
        smlib_brep *block = box(0, 0, 0, a, b, c);
        struct Case
        {
            const std::vector<double> *edges;
            double radius;
        };
        for (Case k : {Case{&all, 1.0}, Case{&all, 0.5}, Case{&vertical, 1.0}, Case{&vertical, 0.5}, Case{&all, 3.0}})
        {
            const double rad = k.radius;
            const bool every = k.edges == &all;
            smlib_brep *round = nullptr;
            CHECK(smlib_fillet(block, k.edges->data(), static_cast<int64_t>(k.edges->size() / 3), 0.1, rad, &round, nullptr) == SMLIB_OK);
            if (round == nullptr)
            {
                continue;
            }
            const double outer = every ? rounded(a, b, c, rad) : a * b * c - 4.0 * (1.0 - PI / 4.0) * rad * rad * c;
            CHECK(near(volume_of(round), outer, 1e-7));
            const double inner = rad > t ? rounded(a - 2 * t, b - 2 * t, c - 2 * t, rad - t) : (a - 2 * t) * (b - 2 * t) * (c - 2 * t);
            std::printf("box %s rounded r%g, inward shell %g: closed", every ? "all edges" : "vertical edges", rad, t);
            int32_t got = smlib_shell(round, nullptr, 0, 0.1, -t, &r, nullptr);
            std::printf(" %d", static_cast<int>(got));
            CHECK(got == SMLIB_OK);
            if (r != nullptr)
            {
                // (a cavity rounded r - t is SMLib's approximation of the offset, to its tolerance)
                CHECK(near(volume_of(r), outer - inner, rad > t ? 1e-6 : 1e-7));
                done(r);
            }
            if (rad <= t)
            {
                got = smlib_shell(round, top, 1, 0.1, -t, &r, nullptr);
                std::printf(", open at the top %d", static_cast<int>(got));
                CHECK(got == SMLIB_OK);
                if (r != nullptr)
                {
                    CHECK(near(volume_of(r), outer - (a - 2 * t) * (b - 2 * t) * (c - t), 1e-7));
                    done(r);
                }
            }
            std::printf("\n");
            smlib_brep_free(round);
        }
        smlib_brep_free(block);
    }

    // A box 20 x 10 x 30 open at the top with walls 4 thick: the sides need room 4 from each
    // other (10 > 8), the floor's is far, so the shell goes through even though the box's
    // smallest half extent is 5; walls of 5 leave none.
    {
        smlib_brep *slab = box(0, 0, 0, 20, 10, 30);
        const double top[3] = {10, 5, 30};
        CHECK(smlib_shell(slab, top, 1, 0.1, -4.0, &r, nullptr) == SMLIB_OK);
        if (r != nullptr)
        {
            CHECK(std::fabs(volume_of(r) - (6000.0 - 12.0 * 2.0 * 26.0)) < 1e-6);
            done(r);
        }
        CHECK(smlib_shell(slab, top, 1, 0.1, -5.0, &r, nullptr) == SMLIB_ERR_NO_ROOM);
        smlib_brep_free(slab);
    }

    // Offset faces: the cube's top up by 5, then down by 5; its side and top together; a
    // cylinder's wall out by 2. Down by the cube's height leaves nothing and is refused.
    {
        const double top[3] = {10, 10, 20};
        CHECK(smlib_offset_faces(cube, top, 1, 0.1, 5.0, &r, nullptr) == SMLIB_OK);
        if (r != nullptr)
        {
            CHECK(std::fabs(volume_of(r) - 10000.0) < 1e-6);
            done(r);
        }
        CHECK(smlib_offset_faces(cube, top, 1, 0.1, -5.0, &r, nullptr) == SMLIB_OK);
        if (r != nullptr)
        {
            CHECK(std::fabs(volume_of(r) - 6000.0) < 1e-6);
            done(r);
        }
        const double two[6] = {10, 10, 20, 20, 10, 10};
        CHECK(smlib_offset_faces(cube, two, 2, 0.1, 2.0, &r, nullptr) == SMLIB_OK);
        if (r != nullptr)
        {
            CHECK(std::fabs(volume_of(r) - 22.0 * 20.0 * 22.0) < 1e-6);
            done(r);
        }
        CHECK(smlib_offset_faces(cube, top, 1, 0.1, -20.0, &r, nullptr) == SMLIB_ERR_NO_ROOM);
        CHECK(smlib_offset_faces(cube, top, 1, 0.1, 0.0, &r, nullptr) == SMLIB_ERR_ARG);
        const double wall[3] = {14, 10, 5};
        if (smlib_offset_faces(rod, wall, 1, 0.1, 2.0, &r, nullptr) == SMLIB_OK)
        {
            done(r);
        }
        // tangency-labs/solidcore#401: a sphere r10 out by 2 is a sphere r12 whose edges lie
        // within their tolerances (the gap was 5.015e-6 over 5e-6 before the results' tolerances
        // were settled).
        const double centre[3] = {0, 0, 0};
        smlib_brep *ball = nullptr;
        CHECK(smlib_make_sphere(centre, 10.0, &ball) == SMLIB_OK);
        const double pole[3] = {10, 0, 0};
        CHECK(smlib_offset_faces(ball, pole, 1, 0.1, 2.0, &r, nullptr) == SMLIB_OK);
        if (r != nullptr)
        {
            CHECK(near(volume_of(r), 4.0 / 3.0 * PI * 1728.0, 1e-6));
            done(r);
        }
        smlib_brep_free(ball);
    }

    // tangency-labs/solidcore#402: a face moved away from a tangent neighbour that keeps its
    // surface no longer meets it. SMLib's local operation dropped the edge across the gap and
    // returned a body whose edge tolerance had grown to the gap. A face moved along its normal
    // stays parallel to its tangent neighbour: a step, closed by a wall between the two. A boss
    // wall moved off its base round (a r10 boss on a 60 x 40 x 10 plate, a r2 round at its
    // base, the wall out by 1) keeps the round and gains a ledge on top of it: the 8 of wall
    // above the round grows from r10 to r11, pi * (121 - 100) * 8. One piece of a split top
    // moved along its normal (a 60 x 40 x 20 block, the top split at x = 30,
    // the piece x > 30 up by 1 adds 30 x 40 x 1, down by 1 takes it off).
    {
        smlib_brep *plate = box(0, 0, 0, 60, 40, 10);
        const double at[3] = {30, 20, 10};
        smlib_brep *boss = nullptr;
        smlib_brep *joined = nullptr;
        smlib_brep *rounded = nullptr;
        CHECK(smlib_make_cylinder(at, 10.0, 10.0, &boss) == SMLIB_OK);
        CHECK(smlib_boolean(plate, boss, SMLIB_BOOLEAN_UNION, &joined) == SMLIB_OK);
        const double base[3] = {30, 30, 10};
        CHECK(smlib_fillet(joined, base, 1, 0.1, 2.0, &rounded, nullptr) == SMLIB_OK && sound(rounded));
        const double wall[3] = {30, 30, 16};
        CHECK(smlib_offset_faces(rounded, wall, 1, 0.1, 1.0, &r, nullptr) == SMLIB_OK && sound(r));
        if (r != nullptr)
        {
            CHECK(near(volume_of(r), volume_of(rounded) + PI * 21.0 * 8.0, 1e-6));
            done(r);
        }
        smlib_brep *block = box(0, 0, 0, 60, 40, 20);
        smlib_brep *cutter = box(30, -1, 10, 40, 42, 20);
        smlib_brep *split = nullptr;
        CHECK(smlib_imprint(block, cutter, &split) == SMLIB_OK);
        const double piece[3] = {45, 20, 20};
        CHECK(smlib_offset_faces(split, piece, 1, 0.1, 1.0, &r, nullptr) == SMLIB_OK && sound(r));
        CHECK(near(volume_of(r), 48000.0 + 1200.0, 1e-6));
        smlib_brep_free(r);
        r = nullptr;
        CHECK(smlib_offset_faces(split, piece, 1, 0.1, -1.0, &r, nullptr) == SMLIB_OK && sound(r));
        CHECK(near(volume_of(r), 48000.0 - 1200.0, 1e-6));
        smlib_brep_free(r);
        r = nullptr;
        for (smlib_brep *b : {plate, boss, joined, rounded, block, cutter, split})
        {
            smlib_brep_free(b);
        }
    }

    // Draft: the cube's side x = 20 tilted 5 degrees about its line on z = 0, leaning in going
    // up (+z): a wedge of 20 x 20 tan(5 deg) / 2 x 20 comes off; -5 degrees adds it. The top
    // isn't parallel to the pull and the rod's wall isn't flat: not supported.
    {
        const double side[3] = {20, 10, 10};
        const double origin[3] = {0, 0, 0};
        const double up[3] = {0, 0, 1};
        const double deg = 3.14159265358979323846 / 180.0;
        const double wedge = 0.5 * 20.0 * 20.0 * std::tan(5.0 * deg) * 20.0;
        CHECK(smlib_draft(cube, side, 1, 0.1, origin, up, up, 5.0 * deg, &r, nullptr) == SMLIB_OK);
        if (r != nullptr)
        {
            CHECK(std::fabs(volume_of(r) - (8000.0 - wedge)) < 1e-6);
            done(r);
        }
        CHECK(smlib_draft(cube, side, 1, 0.1, origin, up, up, -5.0 * deg, &r, nullptr) == SMLIB_OK);
        if (r != nullptr)
        {
            CHECK(std::fabs(volume_of(r) - (8000.0 + wedge)) < 1e-6);
            done(r);
        }
        const double top[3] = {10, 10, 20};
        CHECK(smlib_draft(cube, top, 1, 0.1, origin, up, up, 5.0 * deg, &r, nullptr) == SMLIB_ERR_UNSUPPORTED);
        const double wall[3] = {14, 10, 5};
        CHECK(smlib_draft(rod, wall, 1, 0.1, origin, up, up, 5.0 * deg, &r, nullptr) == SMLIB_ERR_UNSUPPORTED);
        const double flat[3] = {1, 0, 0};
        CHECK(smlib_draft(cube, side, 1, 0.1, origin, flat, up, 5.0 * deg, &r, nullptr) == SMLIB_ERR_ARG);
        CHECK(smlib_draft(cube, side, 1, 0.1, origin, up, up, 0.0, &r, nullptr) == SMLIB_ERR_ARG);
        CHECK(smlib_draft(cube, side, 1, 0.1, origin, up, up, 90.0 * deg, &r, nullptr) == SMLIB_ERR_ARG);
        CHECK(smlib_draft(cube, side, 1, 0.1, nullptr, up, up, 5.0 * deg, &r, nullptr) == SMLIB_ERR_ARG);
    }

    // A cone frustum whose side is two half-cone faces (a profile of two half circles r 20,
    // drafted inward), shelled 2 mm open at the top (rig case shell/3103). The inner rim meets
    // the top plane where the offsets of the side and the top were already mated along it, and
    // the offset's self-intersection merge intersected them again: the plane/offset-cone curve
    // ended 0.0036 mm short of the inner seam, which split each half of the rim there (three
    // extra vertices and edges, then the shell's difference left partially coincident edges).
    // The two halves of each wall and their rims share edges; the cavity is the frustum moved
    // in by the wall normal to the cone: r(z) = 20 - z tan(a) - t / cos(a), floor at z = t.
    for (double deg : {9.4623222080, 15.0})
    {
        const double a = deg * PI / 180.0, t = 2.0, h = 30.0, r0 = 20.0;
        const int64_t kinds[2] = {SMLIB_CURVE_ARC, SMLIB_CURVE_ARC};
        const double reals[24] = {0, 0, 0, 0, 0, 1, 1, 0, 0, r0, 0, PI, 0, 0, 0, 0, 0, 1, 1, 0, 0, r0, PI, 2 * PI};
        const double up[3] = {0, 0, 1};
        smlib_brep *frustum = nullptr;
        CHECK(smlib_extrude(kinds, 2, reals, 24, up, h, a, SMLIB_DRAFT_CORNER_TANGENT, 1, &frustum, nullptr) == SMLIB_OK);
        auto cone = [](double ra, double rb, double len) { return PI * len / 3.0 * (ra * ra + ra * rb + rb * rb); };
        auto inner = [&](double z) { return r0 - z * std::tan(a) - t / std::cos(a); };
        const double want = cone(r0, r0 - h * std::tan(a), h) - cone(inner(t), inner(h), h - t);
        const double lid[3] = {0, 0, h};
        CHECK(frustum != nullptr && smlib_shell(frustum, lid, 1, 0.1, -t, &r, nullptr) == SMLIB_OK);
        if (r != nullptr)
        {
            int64_t counts[3] = {0, 0, 0};
            CHECK(smlib_brep_counts(r, counts) == SMLIB_OK);
            CHECK(counts[0] == 7 && counts[1] == 12 && counts[2] == 8);
            CHECK(near(volume_of(r), want, 1e-5));
            done(r);
        }
        smlib_brep_free(frustum);
    }

    // A 40 x D x 10 plate with a 4 mm through hole at (x, y), shelled t inward open at the top:
    // where the hole lies nearer than t to a side, the offset of its wall (radius 2 + t) crosses
    // the side's offset, and the offset floor's hole loop crossed its outer loop (SmFace test
    // #38). The merge failed on that face (SM_ERR_NULL_POINTER, SmEdgeDefinition::Imprint); the
    // face is now rebuilt from the stretches of its loops that bound it, and comes apart in two
    // when the hole crosses both sides (D = 8). The cavity is the inner box less the part of the
    // hole's offset disk inside the floor, over the height 10 - t; a hole whose offset only
    // touches a side's offset is not resolved: SM_ERR_BAD_INTERSECTIONS. At x = 4.99 the hole's
    // offset crosses the side's 0.01 mm deep, at 4.3 degrees: the cut-off rims of that crossing
    // were welded by the offset's stitch, whose tolerance escalated to the gap between them, and
    // the shell failed (1001).
    {
        struct Plate
        {
            double x, y, d, t;
        };
        // Area of a disk of radius r between the lines at signed distances a < b from its centre.
        auto band = [](double rr, double a, double b)
        {
            auto F = [rr](double u)
            {
                u = std::max(-rr, std::min(rr, u));
                return u * std::sqrt(rr * rr - u * u) + rr * rr * std::asin(u / rr);
            };
            return F(b) - F(a);
        };
        for (Plate p : {Plate{4.5, 20, 40, 1.5}, Plate{4.5, 20, 40, 1.0}, Plate{20, 4, 8, 1.5}, Plate{5, 20, 40, 1.5},
                        Plate{4.99, 20, 40, 1.5}})
        {
            const double o[3] = {0, 0, 0}, size[3] = {40, p.d, 10}, axis[3] = {p.x, p.y, -1}, lid[3] = {35, p.d / 2, 10};
            smlib_brep *box = nullptr, *pin = nullptr, *plate = nullptr;
            CHECK(smlib_make_box(o, size, &box) == SMLIB_OK && smlib_make_cylinder(axis, 2.0, 12.0, &pin) == SMLIB_OK);
            CHECK(smlib_boolean(box, pin, SMLIB_BOOLEAN_SUBTRACT, &plate) == SMLIB_OK);
            const double rr = 2.0 + p.t;
            int32_t s = plate != nullptr ? smlib_shell(plate, lid, 1, 0.1, -p.t, &r, nullptr) : SMLIB_ERR_ARG;
            if (p.x - rr == p.t)
            {
                CHECK(s == SM_ERR_BAD_INTERSECTIONS && r == nullptr);
            }
            else
            {
                double inside = band(rr, std::max(p.t - p.y, -rr), std::min(p.d - p.t - p.y, rr));
                if (p.x - rr < p.t)
                {
                    inside -= band(rr, -rr, p.t - p.x);
                }
                const double cavity = ((40.0 - 2.0 * p.t) * (p.d - 2.0 * p.t) - inside) * (10.0 - p.t);
                CHECK(s == SMLIB_OK);
                if (r != nullptr)
                {
                    CHECK(near(volume_of(r), 40.0 * p.d * 10.0 - PI * 4.0 * 10.0 - cavity, 1e-7));
                    done(r);
                }
            }
            for (smlib_brep *b : {box, pin, plate})
            {
                smlib_brep_free(b);
            }
        }
    }

    // SolidCraft's ladder part 09/011 (tests/data/part-chamfers-shell.txt, as smlib_shell was
    // handed it): a block with a bezier-arched window 1 mm below a round, an M5 hole 0.35 mm from
    // that round, a chamfer, holes and a slot, shelled 1.5 inward open at the bottom. The round's
    // offset is tangent to the wall's offset along z = 32 and to the top's along y = 33, so the
    // window's and the hole's offsets meet them tangentially where they cross: those contacts were
    // found inside the curves instead of at their ends, and their sliver edges failed the merge.
    // The result is checked point by point: an inward offset is sharp at every edge, so each face
    // lies on the part's boundary or 1.5 or more from its faces but the open one, exactly 1.5
    // except near the mitres of the part's concave edges; sampled at one interior point per face
    // (where on a face that point falls decides whether it reads 1.5 or more).
    {
        smlib_flat *f = nullptr;
        smlib_brep *part = nullptr;
        CHECK(smlib_flat_read(SMLIB_TEST_DATA "/part-chamfers-shell.txt", &f) == SMLIB_OK && smlib_brep_import(f, &part) == SMLIB_OK);
        smlib_flat_free(f);
        const double bottom[3] = {40.745625596812403, 30.55585248496563, 0};
        CHECK(part != nullptr && smlib_shell(part, bottom, 1, 0.080157097709407002, -1.5, &r, nullptr) == SMLIB_OK);
        if (r != nullptr)
        {
            SmBrep *walls = nullptr;
            CHECK(SmApiBrepCopy(part->p, walls) == SM_SUCCESS && walls != nullptr);
            walls->m_bEditingEnabled = TRUE;
            SmTArray<SmFace *> wf;
            walls->GetFaces(wf);
            for (ULONG i = 0; i < wf.GetSize(); i++)
            {
                SmExtent3d e;
                CHECK(wf[i]->CalculateBoundingBox(e, TRUE) == SM_SUCCESS);
                if (std::fabs(e.GetMin(2)) < 1e-9 && std::fabs(e.GetMax(2)) < 1e-9)
                {
                    CHECK(walls->DeleteFace(wf[i], TRUE, FALSE) == SM_SUCCESS);
                }
            }
            SmTArray<SmFace *> faces;
            r->p->GetFaces(faces);
            int on_part = 0, at_wall = 0, past_wall = 0;
            for (ULONG i = 0; i < faces.GetSize(); i++)
            {
                SmPoint3d p, q;
                double d_part = -1.0, d_walls = -1.0;
                CHECK(faces[i]->CalculateAnInternalPoint(p) == SM_SUCCESS);
                CHECK(SmApiGetClosestPoint(part->p, p, q, d_part) == SM_SUCCESS);
                CHECK(SmApiGetClosestPoint(walls, p, q, d_walls) == SM_SUCCESS);
                if (d_part <= 1e-9)
                {
                    on_part++;
                }
                else
                {
                    CHECK(d_walls >= 1.5 - 1e-5);
                    (std::fabs(d_walls - 1.5) <= 1e-5 ? at_wall : past_wall)++;
                }
            }
            std::printf("part 09/011 shell: %d faces on the part, %d at the wall, %d past it (at a mitre)\n", on_part, at_wall, past_wall);
            CHECK(on_part == 29 && at_wall + past_wall == 26 && at_wall > past_wall);
            delete walls;
            done(r);
        }
        smlib_brep_free(part);
    }

    // SolidCraft's cross-hole cut (tests/data/cross-pin.txt less cross-hole.txt, as SolidCraft
    // handed them to smlib_boolean): a pin of radius 2 and height 6 about z, a hole of radius 1
    // through it along y, both seamed at 0.6131 rad. Each intersection loop is a closed curve that
    // touches the hole's seam line at its start, which is its end: moving that contact to the
    // loop's other end passed a one-way measurement of the whole turn against the line's one
    // point, and the cut failed (1001). The cavity is 24 pi less the integral of
    // 4 sqrt(1 - x^2) sqrt(4 - x^2) over [-1, 1] (x = sin a: 4 cos^2 a sqrt(4 - sin^2 a), smooth,
    // by Simpson's rule).
    {
        smlib_flat *fa = nullptr, *fb = nullptr;
        smlib_brep *pin = nullptr, *hole = nullptr;
        CHECK(smlib_flat_read(SMLIB_TEST_DATA "/cross-pin.txt", &fa) == SMLIB_OK && smlib_brep_import(fa, &pin) == SMLIB_OK);
        CHECK(smlib_flat_read(SMLIB_TEST_DATA "/cross-hole.txt", &fb) == SMLIB_OK && smlib_brep_import(fb, &hole) == SMLIB_OK);
        smlib_flat_free(fa);
        smlib_flat_free(fb);
        const int n = 2000;
        double lens = 0.0;
        for (int i = 0; i <= n; i++)
        {
            const double a = -PI / 2 + PI * i / n, c = std::cos(a), s = std::sin(a);
            lens += (i == 0 || i == n ? 1.0 : (i % 2 ? 4.0 : 2.0)) * 4.0 * c * c * std::sqrt(4.0 - s * s);
        }
        lens *= PI / n / 3.0;
        CHECK(pin != nullptr && hole != nullptr && smlib_boolean(pin, hole, SMLIB_BOOLEAN_SUBTRACT, &r) == SMLIB_OK);
        if (r != nullptr)
        {
            CHECK(near(volume_of(r), 24.0 * PI - lens, 1e-6));
            done(r);
        }
        smlib_brep_free(pin);
        smlib_brep_free(hole);
    }

    // Thicken: a sheet (the cube's top face alone, SMLib's own way of making one is not in the
    // C API; a solid is refused).
    CHECK(smlib_thicken(cube, 1.0, 0, &r) == SMLIB_ERR_ARG);
    CHECK(smlib_thicken(nullptr, 1.0, 0, &r) == SMLIB_ERR_ARG);
    CHECK(smlib_thicken(cube, NAN, 0, &r) == SMLIB_ERR_ARG);

    // Two holes 6 apart cut one after the other through a disc: the second hole's wall meets the
    // first's (trimmed to the disc) in two lines. Parallel cylinders are intersected by
    // SmCone::IntersectWithCone, which once took the lines' span from the whole cylinders, so
    // the lines started below the trimmed wall and were dropped (the boolean failed, 1001).
    {
        smlib_brep *disc = cylinder(0, 0, 0, 20, 10);
        smlib_brep *first = cylinder(6, 0, -1, 5, 12);
        smlib_brep *second = cylinder(3, 3 * std::sqrt(3.0), -1, 5, 12);
        smlib_brep *one = nullptr;
        smlib_brep *two = nullptr;
        CHECK(smlib_boolean(disc, first, SMLIB_BOOLEAN_SUBTRACT, &one) == SMLIB_OK && sound(one));
        CHECK(smlib_boolean(one, second, SMLIB_BOOLEAN_SUBTRACT, &two) == SMLIB_OK && sound(two));
        // The union of two r5 discs 6 apart: twice a disc less their lens.
        const double lens = 2.0 * 25.0 * std::acos(0.6) - 3.0 * 8.0;
        CHECK(near(volume_of(two), (PI * 400.0 - (2.0 * PI * 25.0 - lens)) * 10.0, 1e-6));
        for (smlib_brep *b : {disc, first, second, one, two})
        {
            smlib_brep_free(b);
        }
    }

    for (smlib_brep *b : {cube, rod, holed, pin, round, cut, big, far, touching})
    {
        smlib_brep_free(b);
    }
    if (failures != 0)
    {
        std::fprintf(stderr, "%d check(s) failed\n", failures);
        return 1;
    }
    std::printf("ops_kernel: all checks passed\n");
    return 0;
}
