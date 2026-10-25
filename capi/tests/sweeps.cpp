// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// The sweeps (src/smlib_sweeps.cpp) against closed-form volumes: extrusions straight and drafted,
// revolutions whole, partial and touching the axis, sweeps along lines and arcs, lofts, screw
// sweeps (Pappus) and sweeps through stations (twisted, scaled, round a closed path), each checked
// as a valid closed solid; then hostile curve lists and station maps.

#include "smlib.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

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

static const double PI = 3.14159265358979323846;

// A curve list being built.
struct List
{
    std::vector<int64_t> kinds;
    std::vector<double> reals;

    List &line(double ax, double ay, double az, double bx, double by, double bz)
    {
        kinds.push_back(SMLIB_CURVE_LINE);
        reals.insert(reals.end(), {ax, ay, az, bx, by, bz});
        return *this;
    }

    List &arc(double cx, double cy, double cz, double nx, double ny, double nz, double rx, double ry, double rz, double r, double a0, double a1)
    {
        kinds.push_back(SMLIB_CURVE_ARC);
        reals.insert(reals.end(), {cx, cy, cz, nx, ny, nz, rx, ry, rz, r, a0, a1});
        return *this;
    }

    // A closed polygon through `pts` (x, y, z each).
    List &polygon(const std::vector<double> &pts)
    {
        size_t n = pts.size() / 3;
        for (size_t i = 0; i < n; ++i)
        {
            size_t j = (i + 1) % n;
            line(pts[3 * i], pts[3 * i + 1], pts[3 * i + 2], pts[3 * j], pts[3 * j + 1], pts[3 * j + 2]);
        }
        return *this;
    }

    List &append(const List &o)
    {
        kinds.insert(kinds.end(), o.kinds.begin(), o.kinds.end());
        reals.insert(reals.end(), o.reals.begin(), o.reals.end());
        return *this;
    }

    int64_t nk() const
    {
        return static_cast<int64_t>(kinds.size());
    }
    int64_t nr() const
    {
        return static_cast<int64_t>(reals.size());
    }
};

// A full circle of radius r about +z through (cx, cy, cz)... in the plane square to `n`.
static List circle(double cx, double cy, double cz, double nx, double ny, double nz, double rx, double ry, double rz, double r)
{
    List l;
    l.arc(cx, cy, cz, nx, ny, nz, rx, ry, rz, r, 0.0, 2 * PI);
    return l;
}

// Checks `*out` (status `s`) is one valid solid of volume `v` within `rel`; frees it. `out`, not
// the body: the operation's call is an argument too, and arguments are evaluated in no set order.
static void expect(const char *name, int32_t s, smlib_brep *const *out, double v, double rel)
{
    smlib_brep *b = *out;
    if (s != SMLIB_OK || b == nullptr)
    {
        std::fprintf(stderr, "%s: status %d\n", name, static_cast<int>(s));
        failures++;
        return;
    }
    double got = 0.0;
    CHECK(smlib_brep_volume(b, 1e-6, &got) == SMLIB_OK);
    int32_t manifold = 0;
    CHECK(smlib_brep_is_manifold_solid(b, &manifold) == SMLIB_OK && manifold == 1);
    int64_t len = 0;
    int64_t bad = -1;
    char buf[4096] = {0};
    int32_t cs = smlib_brep_check(b, 2, buf, sizeof buf, &len, &bad);
    bool valid = cs == SMLIB_OK && bad == 0;
    double err = std::fabs(got - v) / v;
    std::printf("%-34s volume %.6f expected %.6f (rel %.2e) check %s\n", name, got, v, err, valid ? "valid" : buf);
    if (!(err <= rel))
    {
        std::fprintf(stderr, "%s: volume %.9g, expected %.9g\n", name, got, v);
        failures++;
    }
    CHECK(valid);
    smlib_brep_free(b);
}

static const double Z[3] = {0, 0, 1};
static const double O[3] = {0, 0, 0};

int main()
{
    CHECK(smlib_init() == SMLIB_OK);
    smlib_brep *b = nullptr;

    // Extrusions: a 10 x 6 rectangle 5 up; a disc with a hole; a drafted square (a frustum of a
    // pyramid: h/3 (A1 + A2 + sqrt(A1 A2))).
    List rect;
    rect.polygon({0, 0, 0, 10, 0, 0, 10, 6, 0, 0, 6, 0});
    b = nullptr;
    expect("extrude rectangle", smlib_extrude(rect.kinds.data(), rect.nk(), rect.reals.data(), rect.nr(), Z, 5.0, 0.0, SMLIB_DRAFT_CORNER_TANGENT, 1, &b, nullptr), &b, 300.0, 1e-9);
    const double down[3] = {0, 0, -1};
    b = nullptr;
    expect("extrude rectangle down", smlib_extrude(rect.kinds.data(), rect.nk(), rect.reals.data(), rect.nr(), down, 5.0, 0.0, SMLIB_DRAFT_CORNER_TANGENT, 1, &b, nullptr), &b, 300.0,
           1e-9);
    List ring = circle(0, 0, 0, 0, 0, 1, 1, 0, 0, 10);
    ring.append(circle(0, 0, 0, 0, 0, -1, 1, 0, 0, 4));
    b = nullptr;
    expect("extrude annulus", smlib_extrude(ring.kinds.data(), ring.nk(), ring.reals.data(), ring.nr(), Z, 3.0, 0.0, SMLIB_DRAFT_CORNER_TANGENT, 1, &b, nullptr), &b,
           PI * (100.0 - 16.0) * 3.0, 1e-6);
    List sq;
    sq.polygon({-5, -5, 0, 5, -5, 0, 5, 5, 0, -5, 5, 0});
    double t = std::tan(10.0 * PI / 180.0);
    double a2 = 10.0 - 2.0 * 4.0 * t;
    b = nullptr;
    expect("extrude square drafted 10 deg", smlib_extrude(sq.kinds.data(), sq.nk(), sq.reals.data(), sq.nr(), Z, 4.0, 10.0 * PI / 180.0, SMLIB_DRAFT_CORNER_TANGENT, 1, &b, nullptr), &b,
           4.0 / 3.0 * (100.0 + a2 * a2 + 10.0 * a2), 1e-6);
    double a3 = 10.0 + 2.0 * 4.0 * t;
    b = nullptr;
    expect("extrude square drafted -10 deg", smlib_extrude(sq.kinds.data(), sq.nk(), sq.reals.data(), sq.nr(), Z, 4.0, -10.0 * PI / 180.0, SMLIB_DRAFT_CORNER_TANGENT, 1, &b, nullptr),
           &b, 4.0 / 3.0 * (100.0 + a3 * a3 + 10.0 * a3), 1e-6);
    b = nullptr;
    expect("extrude square drafted down 10 deg", smlib_extrude(sq.kinds.data(), sq.nk(), sq.reals.data(), sq.nr(), down, 4.0, 10.0 * PI / 180.0, SMLIB_DRAFT_CORNER_TANGENT, 1, &b, nullptr),
           &b, 4.0 / 3.0 * (100.0 + a2 * a2 + 10.0 * a2), 1e-6);
    b = nullptr;
    expect("extrude square drafted down -10 deg", smlib_extrude(sq.kinds.data(), sq.nk(), sq.reals.data(), sq.nr(), down, 4.0, -10.0 * PI / 180.0, SMLIB_DRAFT_CORNER_TANGENT, 1, &b, nullptr),
           &b, 4.0 / 3.0 * (100.0 + a3 * a3 + 10.0 * a3), 1e-6);

    // Revolutions about z of a profile in the xz plane: a square off the axis (a ring), a
    // quarter turn of it, and a half disc on the axis (a sphere).
    const double Y[3] = {0, 1, 0};
    List off;
    off.polygon({6, 0, 0, 8, 0, 0, 8, 0, 2, 6, 0, 2});
    b = nullptr;
    expect("revolve square full", smlib_revolve(off.kinds.data(), off.nk(), off.reals.data(), off.nr(), O, Z, 2 * PI, 1, &b), &b,
           2 * PI * 7.0 * 4.0, 1e-5);
    b = nullptr;
    expect("revolve square 90", smlib_revolve(off.kinds.data(), off.nk(), off.reals.data(), off.nr(), O, Z, PI / 2, 1, &b), &b,
           PI / 2 * 7.0 * 4.0, 1e-6);
    List half;
    half.line(0, 0, -5, 0, 0, 5).arc(0, 0, 0, 0, -1, 0, 0, 0, 1, 5, 0, PI);
    b = nullptr;
    expect("revolve half disc", smlib_revolve(half.kinds.data(), half.nk(), half.reals.data(), half.nr(), O, Z, 2 * PI, 1, &b), &b,
           4.0 / 3.0 * PI * 125.0, 1e-5);
    // Partial turns past half a turn: a half disc on the axis three quarters round (three
    // quarters of a ball), and the square off the axis likewise.
    for (double turn : {0.5, 0.75, 0.9})
    {
        char name[64];
        std::snprintf(name, sizeof name, "revolve half disc %.2f turn", turn);
        b = nullptr;
        expect(name, smlib_revolve(half.kinds.data(), half.nk(), half.reals.data(), half.nr(), O, Z, 2 * PI * turn, 1, &b), &b,
               turn * 4.0 / 3.0 * PI * 125.0, 1e-5);
        std::snprintf(name, sizeof name, "revolve square %.2f turn", turn);
        b = nullptr;
        expect(name, smlib_revolve(off.kinds.data(), off.nk(), off.reals.data(), off.nr(), O, Z, 2 * PI * turn, 1, &b), &b,
               2 * PI * turn * 7.0 * 4.0, 1e-5);
    }
    (void)Y;

    // Sweeps: a disc along a line (a cylinder), along a quarter arc (a quarter torus), and off
    // the path's start (the profile keeps its offset: Pappus about the arc's centre).
    List disc = circle(0, 0, 0, 0, 0, 1, 1, 0, 0, 2);
    List up;
    up.line(0, 0, 0, 0, 0, 10);
    b = nullptr;
    expect("sweep disc along line",
           smlib_sweep(disc.kinds.data(), disc.nk(), disc.reals.data(), disc.nr(), up.kinds.data(), up.nk(), up.reals.data(), up.nr(),
                       SMLIB_SWEEP_FOLLOW, 1, &b),
           &b, PI * 4.0 * 10.0, 1e-4);
    // Along a straight path (here two lines one after the other) a sweep is an extrusion, its
    // walls exact: a rectangle's are six planes, a disc's wall a cylinder.
    List straight;
    straight.line(0, 0, 0, 0, 0, 4).line(0, 0, 4, 0, 0, 10);
    for (const List *profile : {&rect, &disc})
    {
        b = nullptr;
        CHECK(smlib_sweep(profile->kinds.data(), profile->nk(), profile->reals.data(), profile->nr(), straight.kinds.data(), straight.nk(),
                          straight.reals.data(), straight.nr(), SMLIB_SWEEP_FOLLOW, 1, &b) == SMLIB_OK);
        smlib_flat *f = nullptr;
        const int64_t *kinds = nullptr;
        int64_t n = 0;
        CHECK(b != nullptr && smlib_brep_export(b, &f) == SMLIB_OK);
        CHECK(f != nullptr && smlib_flat_get_i64(f, SMLIB_FLAT_FACE_SURFACE_TYPE, &kinds, &n) == SMLIB_OK);
        int64_t planes = 0, cylinders = 0;
        for (int64_t i = 0; i < n; ++i)
        {
            planes += kinds[i] == SMLIB_FLAT_SURFACE_PLANE;
            cylinders += kinds[i] == SMLIB_FLAT_SURFACE_CYLINDER;
        }
        std::printf("%-34s faces %lld, planes %lld, cylinders %lld\n", profile == &rect ? "sweep rectangle along lines" : "sweep disc along lines",
                    static_cast<long long>(n), static_cast<long long>(planes), static_cast<long long>(cylinders));
        CHECK(profile == &rect ? planes == n && n == 6 : planes == 2 && cylinders == n - 2);
        smlib_flat_free(f);
        expect(profile == &rect ? "sweep rectangle along lines" : "sweep disc along lines", SMLIB_OK, &b, profile == &rect ? 600.0 : PI * 40.0,
               1e-6);
    }
    // Round a closed path (a whole circle, R 30) the ends meet: the walls close on themselves,
    // with no caps (Pappus: a disc r 3 and a 5 x 5 square, each centred on the path, in the
    // plane square to it at the path's start).
    List ring_path;
    ring_path.arc(0, 0, 0, 0, 0, 1, 1, 0, 0, 30, 0, 2 * PI);
    List round_section = circle(30, 0, 0, 0, 1, 0, 1, 0, 0, 3);
    List square_section;
    square_section.polygon({27.5, 0, -2.5, 32.5, 0, -2.5, 32.5, 0, 2.5, 27.5, 0, 2.5});
    for (const List *section : {&round_section, &square_section})
    {
        const bool round = section == &round_section;
        b = nullptr;
        int32_t s = smlib_sweep(section->kinds.data(), section->nk(), section->reals.data(), section->nr(), ring_path.kinds.data(),
                                ring_path.nk(), ring_path.reals.data(), ring_path.nr(), SMLIB_SWEEP_FOLLOW, 1, &b);
        int64_t counts[3] = {0, 0, 0};
        CHECK(s == SMLIB_OK && b != nullptr && smlib_brep_counts(b, counts) == SMLIB_OK);
        // Round arcs about one axis the walls are SM_API's surfaces of revolution, analytic for
        // each profile curve: the disc's a torus; the square's two planes and two cylinders (its
        // lines parallel to the axis).
        smlib_flat *f = nullptr;
        const int64_t *types = nullptr;
        int64_t n = 0;
        CHECK(b != nullptr && smlib_brep_export(b, &f) == SMLIB_OK);
        CHECK(f != nullptr && smlib_flat_get_i64(f, SMLIB_FLAT_FACE_SURFACE_TYPE, &types, &n) == SMLIB_OK);
        int64_t nurbs = 0;
        for (int64_t i = 0; i < n; ++i)
        {
            nurbs += types[i] == SMLIB_FLAT_SURFACE_NURBS;
        }
        smlib_flat_free(f);
        std::printf("%-34s faces %lld, nurbs %lld\n", round ? "sweep disc round a circle" : "sweep square round a circle",
                    static_cast<long long>(counts[0]), static_cast<long long>(nurbs));
        CHECK(counts[0] == (round ? 1 : 4) && nurbs == 0);
        expect(round ? "sweep disc round a circle" : "sweep square round a circle", s, &b, (round ? PI * 9.0 : 25.0) * 2 * PI * 30.0, 1e-6);
    }
    // Round a closed path of lines and arcs (a stadium: 20 straight, a half turn of radius 10,
    // back, and round) the ends meet the same way: one wall per path piece, no caps, A L.
    List stadium;
    stadium.line(0, -10, 0, 20, -10, 0).arc(20, 0, 0, 0, 0, 1, 0, -1, 0, 10, 0, PI);
    stadium.line(20, 10, 0, 0, 10, 0).arc(0, 0, 0, 0, 0, 1, 0, 1, 0, 10, 0, PI);
    List stadium_section = circle(0, -10, 0, 1, 0, 0, 0, 1, 0, 2);
    b = nullptr;
    int32_t closed_status = smlib_sweep(stadium_section.kinds.data(), stadium_section.nk(), stadium_section.reals.data(), stadium_section.nr(),
                                        stadium.kinds.data(), stadium.nk(), stadium.reals.data(), stadium.nr(), SMLIB_SWEEP_FOLLOW, 1, &b);
    int64_t stadium_counts[3] = {0, 0, 0};
    CHECK(closed_status == SMLIB_OK && b != nullptr && smlib_brep_counts(b, stadium_counts) == SMLIB_OK);
    std::printf("%-34s faces %lld\n", "sweep disc round a stadium", static_cast<long long>(stadium_counts[0]));
    CHECK(stadium_counts[0] == 4);
    expect("sweep disc round a stadium", closed_status, &b, PI * 4.0 * (40.0 + 20.0 * PI), 1e-4);
    List bend;
    bend.arc(-20, 0, 0, 0, -1, 0, 1, 0, 0, 20, 0, PI / 2);
    b = nullptr;
    expect("sweep disc along quarter arc",
           smlib_sweep(disc.kinds.data(), disc.nk(), disc.reals.data(), disc.nr(), bend.kinds.data(), bend.nk(), bend.reals.data(), bend.nr(),
                       SMLIB_SWEEP_FOLLOW, 1, &b),
           &b, PI * 4.0 * (PI / 2 * 20.0), 1e-3);
    List offset = circle(5, 0, 0, 0, 0, 1, 1, 0, 0, 2);
    b = nullptr;
    expect("sweep offset disc along arc",
           smlib_sweep(offset.kinds.data(), offset.nk(), offset.reals.data(), offset.nr(), bend.kinds.data(), bend.nk(), bend.reals.data(),
                       bend.nr(), SMLIB_SWEEP_FOLLOW, 1, &b),
           &b, PI * 4.0 * (PI / 2 * 25.0), 1e-3);

    // Profiles with a hole along a path that bends (solidcore#411): each loop's wall comes in a
    // face per path piece, and the caps take the faces at the path's ends (Pappus: the
    // profile's centroid is on the path, so area x path length).
    List annulus = circle(0, 0, 0, 0, 0, 1, 1, 0, 0, 3);
    annulus.append(circle(0, 0, 0, 0, 0, -1, 1, 0, 0, 2));
    List frame;
    frame.polygon({-3, -3, 0, 3, -3, 0, 3, 3, 0, -3, 3, 0});
    frame.append(circle(0, 0, 0, 0, 0, -1, 1, 0, 0, 2));
    List hook;
    hook.line(0, 0, 0, 0, 0, 10).arc(20, 0, 10, 0, 1, 0, -1, 0, 0, 20, 0, PI / 2);
    List crank = hook;
    crank.line(20, 0, 30, 35, 0, 30);
    const double hook_length = 10.0 + PI / 2 * 20.0;
    struct Hollow
    {
        const char *name;
        const List *profile;
        double area;
        const List *path;
        double length;
    };
    const Hollow hollows[] = {
        {"sweep annulus along line+arc", &annulus, PI * (9.0 - 4.0), &hook, hook_length},
        {"sweep annulus along line+arc+line", &annulus, PI * (9.0 - 4.0), &crank, hook_length + 15.0},
        {"sweep holed square along line+arc", &frame, 36.0 - PI * 4.0, &hook, hook_length},
    };
    for (const Hollow &h : hollows)
    {
        b = nullptr;
        expect(h.name,
               smlib_sweep(h.profile->kinds.data(), h.profile->nk(), h.profile->reals.data(), h.profile->nr(), h.path->kinds.data(), h.path->nk(),
                           h.path->reals.data(), h.path->nr(), SMLIB_SWEEP_FOLLOW, 1, &b),
               &b, h.area * h.length, 1e-4);
    }

    // Lofts: square to square (a box), circle to circle (a frustum, ruled).
    List boxes;
    boxes.polygon({0, 0, 0, 10, 0, 0, 10, 10, 0, 0, 10, 0}).polygon({0, 0, 8, 10, 0, 8, 10, 10, 8, 0, 10, 8});
    const int64_t four[2] = {4, 4};
    b = nullptr;
    expect("loft square to square", smlib_loft(boxes.kinds.data(), boxes.nk(), boxes.reals.data(), boxes.nr(), four, 2, 1, 1, &b), &b, 800.0,
           1e-6);
    List cone = circle(0, 0, 0, 0, 0, 1, 1, 0, 0, 5);
    cone.append(circle(0, 0, 6, 0, 0, 1, 1, 0, 0, 3));
    const int64_t one[2] = {1, 1};
    b = nullptr;
    expect("loft circle to circle", smlib_loft(cone.kinds.data(), cone.nk(), cone.reals.data(), cone.nr(), one, 2, 1, 1, &b), &b,
           PI * 6.0 / 3.0 * (25.0 + 9.0 + 15.0), 1e-4);

    // Lofts through sections with holes: each section's loops (outline, then holes) are one
    // profile, lofted loop to loop. A tube (r 20 out, 12 in, 20 long), a square frame (20 x 20
    // with a 10 x 10 hole, 10 long, both rings ruled).
    List tube = circle(0, 0, 0, 0, 0, 1, 1, 0, 0, 20);
    tube.append(circle(0, 0, 0, 0, 0, -1, 1, 0, 0, 12));
    tube.append(circle(0, 0, 20, 0, 0, 1, 1, 0, 0, 20));
    tube.append(circle(0, 0, 20, 0, 0, -1, 1, 0, 0, 12));
    const int64_t two[2] = {2, 2};
    b = nullptr;
    expect("loft ring to ring", smlib_loft(tube.kinds.data(), tube.nk(), tube.reals.data(), tube.nr(), two, 2, 1, 1, &b), &b,
           PI * (400.0 - 144.0) * 20.0, 1e-4);
    List square_ring;
    square_ring.polygon({-10, -10, 0, 10, -10, 0, 10, 10, 0, -10, 10, 0}).polygon({-5, -5, 0, -5, 5, 0, 5, 5, 0, 5, -5, 0});
    square_ring.polygon({-10, -10, 10, 10, -10, 10, 10, 10, 10, -10, 10, 10}).polygon({-5, -5, 10, -5, 5, 10, 5, 5, 10, 5, -5, 10});
    const int64_t eight[2] = {8, 8};
    b = nullptr;
    expect("loft frame to frame",
           smlib_loft(square_ring.kinds.data(), square_ring.nk(), square_ring.reals.data(), square_ring.nr(), eight, 2, 1, 1, &b), &b,
           300.0 * 10.0, 1e-6);
    // Two holes whose sizes swap from one section to the other: each is lofted to the hole
    // given in its place, whatever the sizes. A 40 x 20 plate 10 long, less two square frusta
    // between 4 and 8 wide (h/3 (16 + 64 + 32) each).
    List two_holes;
    two_holes.polygon({-20, -10, 0, 20, -10, 0, 20, 10, 0, -20, 10, 0});
    two_holes.polygon({-12, -2, 0, -12, 2, 0, -8, 2, 0, -8, -2, 0}).polygon({6, -4, 0, 6, 4, 0, 14, 4, 0, 14, -4, 0});
    two_holes.polygon({-20, -10, 10, 20, -10, 10, 20, 10, 10, -20, 10, 10});
    two_holes.polygon({-14, -4, 10, -14, 4, 10, -6, 4, 10, -6, -4, 10}).polygon({8, -2, 10, 8, 2, 10, 12, 2, 10, 12, -2, 10});
    const int64_t twelve[2] = {12, 12};
    b = nullptr;
    expect("loft two holes swapping sizes",
           smlib_loft(two_holes.kinds.data(), two_holes.nk(), two_holes.reals.data(), two_holes.nr(), twelve, 2, 1, 1, &b), &b,
           800.0 * 10.0 - 2.0 * 10.0 / 3.0 * (16.0 + 64.0 + 32.0), 1e-6);

    // Lofts to a point: a square base 20 x 20 to an apex 15 above its centre (A h / 3), and an
    // off-centre apex below a triangle (the same, whatever the apex's sideways place).
    List base;
    base.polygon({-10, -10, 0, 10, -10, 0, 10, 10, 0, -10, 10, 0});
    const double top[3] = {0, 0, 15};
    b = nullptr;
    expect("loft square to a point", smlib_loft_to_point(base.kinds.data(), base.nk(), base.reals.data(), base.nr(), top, &b), &b,
           400.0 * 15.0 / 3.0, 1e-9);
    List tri;
    tri.polygon({0, 0, 5, 12, 0, 5, 0, 9, 5});
    const double below[3] = {30, -7, -1};
    b = nullptr;
    expect("loft triangle to a point below", smlib_loft_to_point(tri.kinds.data(), tri.nk(), tri.reals.data(), tri.nr(), below, &b), &b,
           54.0 * 6.0 / 3.0, 1e-9);
    // An apex in the base's plane, or a curved base, is refused.
    const double flat[3] = {30, 30, 0};
    b = nullptr;
    CHECK(smlib_loft_to_point(base.kinds.data(), base.nk(), base.reals.data(), base.nr(), flat, &b) == SMLIB_ERR_ARG && b == nullptr);
    b = nullptr;
    CHECK(smlib_loft_to_point(cone.kinds.data(), 1, cone.reals.data(), 12, top, &b) == SMLIB_ERR_ARG && b == nullptr);

    // Screw sweeps: a 0.4 x 0.4 square at radius 4..4.4 in the xz plane, 3 turns of pitch 1,
    // up and down (Pappus: area x 2 pi x centroid radius x turns).
    List tooth;
    tooth.polygon({4, 0, 0, 4.4, 0, 0, 4.4, 0, 0.4, 4, 0, 0.4});
    double screw = 0.16 * 2 * PI * 4.2 * 3.0;
    b = nullptr;
    expect("screw sweep up", smlib_screw_sweep(tooth.kinds.data(), tooth.nk(), tooth.reals.data(), tooth.nr(), O, Z, 1.0, 3.0, 1e-4, &b), &b,
           screw, 2e-3);
    b = nullptr;
    expect("screw sweep down", smlib_screw_sweep(tooth.kinds.data(), tooth.nk(), tooth.reals.data(), tooth.nr(), O, Z, -1.0, 3.0, 1e-4, &b), &b,
           screw, 2e-3);

    // tangency-labs/solidcore#434: two walls 12 thick and 60 tall swept along half circles,
    // wall A r 44..56 over 0..180 degrees, wall B r 52..64 over 90..270, overlapping over a
    // quarter turn and sharing the planes z = 0 and z = 60. Wall A's frame started with its z in
    // the path's plane (a start tangent along +y), and its flat top came out 2.9e-4 off plane;
    // joined, the walls failed (1006). Flat, the union merged the tops onto wall B's surface
    // continued past its ends, where its radial lines run back over the sector through its
    // centre: wall A's edges dropped onto the wrong sheet, and the result had UV gaps and the
    // wrong volume. V = V_A + V_B - overlap.
    {
        List path_a;
        path_a.arc(0, 0, 0, 0, 0, 1, 1, 0, 0, 50.0, 0.0, PI);
        List wall_a;
        wall_a.polygon({44, 0, 0, 56, 0, 0, 56, 0, 60, 44, 0, 60});
        List path_b;
        path_b.arc(0, 0, 0, 0, 0, 1, 0, 1, 0, 58.0, 0.0, PI);
        List wall_b;
        wall_b.polygon({0, 52, 0, 0, 64, 0, 0, 64, 60, 0, 52, 60});
        const double v_a = 60.0 * (PI / 2.0) * (56.0 * 56.0 - 44.0 * 44.0);
        const double v_b = 60.0 * (PI / 2.0) * (64.0 * 64.0 - 52.0 * 52.0);
        const double overlap = 60.0 * (PI / 4.0) * (56.0 * 56.0 - 52.0 * 52.0);
        smlib_brep *a = nullptr;
        smlib_brep *w = nullptr;
        CHECK(smlib_sweep(wall_a.kinds.data(), wall_a.nk(), wall_a.reals.data(), wall_a.nr(), path_a.kinds.data(), path_a.nk(), path_a.reals.data(),
                          path_a.nr(), SMLIB_SWEEP_FOLLOW, 1, &a) == SMLIB_OK);
        CHECK(smlib_sweep(wall_b.kinds.data(), wall_b.nk(), wall_b.reals.data(), wall_b.nr(), path_b.kinds.data(), path_b.nk(), path_b.reals.data(),
                          path_b.nr(), SMLIB_SWEEP_FOLLOW, 1, &w) == SMLIB_OK);
        b = nullptr;
        int32_t joined = a != nullptr && w != nullptr ? smlib_boolean(a, w, SMLIB_BOOLEAN_UNION, &b) : SMLIB_ERR_ARG;
        expect("two swept arc walls joined", joined, &b, v_a + v_b - overlap, 1e-5);
        smlib_brep_free(a);
        smlib_brep_free(w);
    }

    // Circles on faces whose (u, v) run at different speeds. An extrusion's side face is a plane
    // spanning its edge's length by the height on [0,1]^2, and a circle's exact trim there is an
    // ellipse as many times as long as it is wide. The bodies are SolidCraft's: a cylinder's seam
    // at 0.6131 rad from its profile plane's x axis, that axis Z x n for an axis n off Z.
    {
        const double seam = 0.6131;
        // SolidCraft #955: a 2 mm drill 10 deep into the 60 x 4 rim of a 60 x 40 x 4 sheet, from
        // the hole's bottom out and 0.72 past the rim, clipped to the sheet: its end on the rim.
        // The polygon of the trim's first quarter arc, square on the rim, folds back by 171
        // degrees in (u, v). SmEdge::ValidateGeometry's reversal test, run by the check, took
        // that for a reversal and moved the arc's middle control point onto its chord: test #7,
        // a kink at the quarter knot.
        List sheet;
        sheet.polygon({0, 40, 0, 0, 0, 0, 60, 0, 0, 60, 40, 0});
        List bore;
        bore.arc(30, 10, 2, 0, -1, 0, 1, 0, 0, 1.0, seam, seam + 2 * PI);
        const double Y[3] = {0, -1, 0};
        smlib_brep *slab = nullptr, *drill = nullptr;
        CHECK(smlib_extrude(sheet.kinds.data(), sheet.nk(), sheet.reals.data(), sheet.nr(), Z, 4.0, 0.0, SMLIB_DRAFT_CORNER_TANGENT, 1, &slab, nullptr) ==
              SMLIB_OK);
        CHECK(smlib_extrude(bore.kinds.data(), bore.nk(), bore.reals.data(), bore.nr(), Y, 10.722218803410712, 0.0, SMLIB_DRAFT_CORNER_TANGENT, 1, &drill,
                            nullptr) == SMLIB_OK);
        b = nullptr;
        int32_t s = slab != nullptr && drill != nullptr ? smlib_boolean(drill, slab, SMLIB_BOOLEAN_INTERSECT, &b) : SMLIB_ERR_ARG;
        expect("drill clipped to a sheet's rim", s, &b, 10.0 * PI, 1e-6);
        smlib_brep_free(slab);
        smlib_brep_free(drill);

        // A 300 x 200 x 3 plate joined to a r 4 rod along its edge y = 0, z = 3: the rod's end
        // discs lie on the plate's 200 x 3 ends, where the trim of the r 4 circle is an ellipse
        // whose radius at its ends, 0.0003 in (u, v), test #13 took for a cusp (and the reversal
        // test kinked the trim, as above). Inside the plate is the rod's quarter y >= 0,
        // z <= 3 less the half segment below z = 0, (r^2 acos(3/r) - 3 sqrt(r^2 - 9)) / 2.
        List plate;
        plate.polygon({0, 0, 0, 300, 0, 0, 300, 200, 0, 0, 200, 0});
        List rod;
        rod.arc(0, 0, 3, 1, 0, 0, 0, 1, 0, 4.0, seam, seam + 2 * PI);
        const double X[3] = {1, 0, 0};
        smlib_brep *p = nullptr, *r = nullptr;
        CHECK(smlib_extrude(plate.kinds.data(), plate.nk(), plate.reals.data(), plate.nr(), Z, 3.0, 0.0, SMLIB_DRAFT_CORNER_TANGENT, 1, &p, nullptr) ==
              SMLIB_OK);
        CHECK(smlib_extrude(rod.kinds.data(), rod.nk(), rod.reals.data(), rod.nr(), X, 300.0, 0.0, SMLIB_DRAFT_CORNER_TANGENT, 1, &r, nullptr) ==
              SMLIB_OK);
        b = nullptr;
        s = p != nullptr && r != nullptr ? smlib_boolean(p, r, SMLIB_BOOLEAN_UNION, &b) : SMLIB_ERR_ARG;
        const double below = 0.5 * (16.0 * std::acos(0.75) - 3.0 * std::sqrt(7.0));
        expect("plate with a rod along its edge", s, &b, 300.0 * 200.0 * 3.0 + 300.0 * (0.75 * PI * 16.0 + below), 1e-6);
        smlib_brep_free(p);
        smlib_brep_free(r);

        // SolidCraft #1038: a 1 mm hole through a 200 x 200 block whose top is one cubic wave
        // across x (u on [0,1] over its 200 mm, v its length in mm): the hole's rim on the wave
        // is an ellipse some 200 times as long as it is wide in (u, v), of radius 1.3e-5 at its
        // ends (test #13). The wave is z = 40 at the hole's axis x = 100 with slope 0.2256
        // there, so the hole removes 10 pi within 1e-5 of it.
        List wave;
        wave.line(0, 200, 0, 200, 200, 0).line(200, 200, 0, 200, 200, 40);
        wave.kinds.insert(wave.kinds.end(), {SMLIB_CURVE_NURBS, 3, 4});
        wave.reals.insert(wave.reals.end(), {200, 200, 40, 1, 133, 200, 70, 1, 67, 200, 10, 1, 0, 200, 40, 1, 0, 0, 0, 0, 1, 1, 1, 1});
        wave.line(0, 200, 40, 0, 200, 0);
        List pin;
        pin.arc(100, 100, -10, 0, 0, 1, 1, 0, 0, 0.5, seam, seam + 2 * PI);
        const double minus_y[3] = {0, -1, 0};
        smlib_brep *block = nullptr, *tool = nullptr;
        CHECK(smlib_extrude(wave.kinds.data(), wave.nk(), wave.reals.data(), wave.nr(), minus_y, 200.0, 0.0, SMLIB_DRAFT_CORNER_TANGENT, 1, &block,
                            nullptr) == SMLIB_OK);
        CHECK(smlib_extrude(pin.kinds.data(), pin.nk(), pin.reals.data(), pin.nr(), Z, 100.0, 0.0, SMLIB_DRAFT_CORNER_TANGENT, 1, &tool, nullptr) ==
              SMLIB_OK);
        double whole = 0.0;
        CHECK(block != nullptr && smlib_brep_volume(block, 1e-7, &whole) == SMLIB_OK);
        b = nullptr;
        s = block != nullptr && tool != nullptr ? smlib_boolean(block, tool, SMLIB_BOOLEAN_SUBTRACT, &b) : SMLIB_ERR_ARG;
        double holed = 0.0;
        CHECK(b != nullptr && smlib_brep_volume(b, 1e-7, &holed) == SMLIB_OK);
        std::printf("%-34s removed %.9f expected %.9f\n", "small hole through a wave", whole - holed, 10.0 * PI);
        CHECK(std::fabs(whole - holed - 10.0 * PI) <= 1e-5 * 10.0 * PI);
        expect("small hole through a wave", s, &b, whole - 10.0 * PI, 1e-5 * 10.0 * PI / whole);
        smlib_brep_free(block);
        smlib_brep_free(tool);
    }

    // Modelled ISO metric threads, the groove (the space between the teeth, ISO 68-1 basic
    // profile, 60 degree flanks) screwed a pitch past each end of the threaded length and cut
    // away, as SolidCraft cuts them. External: from the root (5H/8 under the major diameter, P/4
    // wide) out past the shaft's surface. Internal (a nut, a 20 mm disc with a bore): from inside
    // the bore wall out to the major diameter (P/8 wide). Each cut is a valid solid and removes
    // the groove's volume within the material: per turn 2 pi * integral of the groove's width
    // * r dr over the radii it cuts (a screw motion keeps volume as a revolution does).
    // M10x1.5 on 12 mm, M33x3.5 on 14 mm and M44x1 on 4 mm are SolidCraft's thread cases 1250,
    // 2193 and 2209; their flank's intersection with an end face meets the shaft's rim at a
    // grazing angle, and once ran past its vertex along the rim (SmFace test #34, partially
    // coincident edges, solidcore#435).
    struct Thread
    {
        double d, p, length;
        bool internal;
        double bore;
    };
    const Thread threads[] = {
        {10.0, 1.5, 20.0, false, 0.0}, {8.0, 1.25, 20.0, false, 0.0}, {10.0, 1.5, 12.0, false, 0.0}, {33.0, 3.5, 14.0, false, 0.0},
        {44.0, 1.0, 4.0, false, 0.0},  {6.0, 1.0, 24.0, false, 0.0},  {10.0, 1.5, 12.0, true, 8.5},
    };
    for (const Thread &t : threads)
    {
        const double d = t.d;
        const double p = t.p;
        const double h = p * std::sqrt(3.0) / 2.0;
        const double slope = std::tan(PI / 6.0);
        const double basic = 5.0 * h / 8.0;
        // The groove's trapezoid: half width w0 at radius r0, w1 at r1; the material from rm0 to rm1.
        double r0, w0, r1, w1, rm0, rm1;
        if (t.internal)
        {
            r0 = std::min(t.bore / 2.0, d / 2.0 - basic) - 0.1 * p;
            w0 = (p / 8.0 + 2.0 * (d / 2.0 - r0) * slope) / 2.0;
            r1 = d / 2.0;
            w1 = p / 16.0;
            rm0 = t.bore / 2.0;
            rm1 = r1;
        }
        else
        {
            r0 = d / 2.0 - basic;
            w0 = p / 8.0;
            r1 = d / 2.0 + 0.05 * p;
            w1 = (p / 4.0 + 2.0 * (r1 - r0) * slope) / 2.0;
            rm0 = r0;
            rm1 = d / 2.0;
        }
        List groove;
        groove.polygon({r0, 0, -p - w0, r1, 0, -p - w1, r1, 0, -p + w1, r0, 0, -p + w0});
        smlib_brep *tool = nullptr;
        smlib_brep *body = nullptr;
        CHECK(smlib_screw_sweep(groove.kinds.data(), groove.nk(), groove.reals.data(), groove.nr(), O, Z, p, (t.length + 2.0 * p) / p, p * 1e-4,
                                &tool) == SMLIB_OK);
        double material = 0.0;
        if (t.internal)
        {
            smlib_brep *disc = nullptr;
            smlib_brep *bore = nullptr;
            const double below[3] = {0, 0, -1};
            CHECK(smlib_make_cylinder(O, 10.0, t.length, &disc) == SMLIB_OK);
            CHECK(smlib_make_cylinder(below, t.bore / 2.0, t.length + 2.0, &bore) == SMLIB_OK);
            CHECK(disc != nullptr && bore != nullptr && smlib_boolean(disc, bore, SMLIB_BOOLEAN_SUBTRACT, &body) == SMLIB_OK);
            smlib_brep_free(disc);
            smlib_brep_free(bore);
            material = PI * (100.0 - t.bore * t.bore / 4.0) * t.length;
        }
        else
        {
            CHECK(smlib_make_cylinder(O, d / 2.0, t.length, &body) == SMLIB_OK);
            material = PI * d * d / 4.0 * t.length;
        }
        const double slope_w = (w1 - w0) / (r1 - r0);
        auto f = [&](double r) { return 2.0 * ((w0 - slope_w * r0) * r * r / 2.0 + slope_w * r * r * r / 3.0); };
        const double per_turn = 2.0 * PI * (f(rm1) - f(rm0));
        char name[64];
        std::snprintf(name, sizeof name, "%s M%gx%g on %g mm", t.internal ? "internal" : "external", d, p, t.length);
        b = nullptr;
        int32_t cut = tool != nullptr && body != nullptr ? smlib_boolean(body, tool, SMLIB_BOOLEAN_SUBTRACT, &b) : SMLIB_ERR_ARG;
        // Every thread's edges end at their vertices: no vertex needed a tolerance grown past the
        // body's 1e-5 (1e-5 mm, written as twice that) to hold a gap. In the internal M10 a
        // flank's curve on the bore stopped 3.1e-5 short of its vertex on the bore's rim
        // (solidcore#410); it is continued to the vertex now.
        smlib_flat *flat = nullptr;
        if (cut == SMLIB_OK && smlib_brep_export(b, &flat) == SMLIB_OK)
        {
            const double *tol = nullptr;
            int64_t n = 0;
            double widest = 0.0;
            CHECK(smlib_flat_get_f64(flat, SMLIB_FLAT_VERTEX_TOLERANCE, &tol, &n) == SMLIB_OK && n > 0);
            for (int64_t k = 0; k < n; ++k)
            {
                widest = std::max(widest, tol[k]);
            }
            if (!(widest <= 2e-5 * (1.0 + 1e-9)))
            {
                std::fprintf(stderr, "%s: a vertex tolerance grew to %g\n", name, widest / 2.0);
                failures++;
            }
        }
        smlib_flat_free(flat);
        expect(name, cut, &b, material - t.length / p * per_turn, 1e-4);
        smlib_brep_free(tool);
        smlib_brep_free(body);
    }

    // Sweeps through stations: the profile placed by affine maps along its path.
    {
        // A 10 x 10 square up 40 along z, turning 90 degrees about it: the section area stays, so
        // V = 100 * 40. The skin interpolates the stations with free ends, its error falling as
        // the square of the step: 2e-7 of the volume at these 1.4 degree steps. The bound is the
        // volume's own: measured at 1e-6 accuracy, even the exact frustum below comes out 1.1e-6
        // off (3.6e-8 at 1e-8).
        const int n = 65;
        List square;
        square.polygon({-5, -5, 0, 5, -5, 0, 5, 5, 0, -5, 5, 0});
        std::vector<double> twist;
        for (int i = 0; i < n; ++i)
        {
            double f = static_cast<double>(i) / (n - 1);
            double a = 0.5 * PI * f;
            twist.insert(twist.end(), {std::cos(a), std::sin(a), 0, -std::sin(a), std::cos(a), 0, 0, 0, 1, 0, 0, 40.0 * f});
        }
        b = nullptr;
        expect("stations: twisted square", smlib_sweep_stations(square.kinds.data(), square.nk(), square.reals.data(), square.nr(), twist.data(), n,
                                                                 SMLIB_ENDS_CAPPED, &b),
               &b, 4000.0, 1e-5);
        // A circle of radius 5 grown to 10 over 30: a cone frustum, h/3 (A1 + A2 + sqrt(A1 A2)).
        List disc = circle(0, 0, 0, 0, 0, 1, 1, 0, 0, 5.0);
        std::vector<double> grow;
        for (int i = 0; i < n; ++i)
        {
            double f = static_cast<double>(i) / (n - 1);
            double k = 1.0 + f;
            grow.insert(grow.end(), {k, 0, 0, 0, k, 0, 0, 0, 1, 0, 0, 30.0 * f});
        }
        b = nullptr;
        expect("stations: cone frustum", smlib_sweep_stations(disc.kinds.data(), disc.nk(), disc.reals.data(), disc.nr(), grow.data(), n, SMLIB_ENDS_CAPPED, &b),
               &b, 10.0 * PI * (25.0 + 100.0 + 50.0), 1e-5);
        // An 8 x 8 square (x 26..34, z -4..4) carried round the circle of radius 30 about z,
        // turning 90 degrees about the path on the way: it ends on itself, a quarter turn on, and
        // the ends are joined. Pappus: V = 64 * 2 pi 30.
        const int m = 97;
        List ring;
        ring.polygon({26, 0, -4, 34, 0, -4, 34, 0, 4, 26, 0, 4});
        std::vector<double> round;
        for (int i = 0; i < m; ++i)
        {
            double f = static_cast<double>(i) / (m - 1);
            double phi = 2.0 * PI * f, psi = 0.5 * PI * f;
            double cp = std::cos(phi), sp = std::sin(phi), cs = std::cos(psi), ss = std::sin(psi);
            // M = Rz(phi) Ry(psi); t = Rz(phi) (c - Ry(psi) c), c = (30, 0, 0).
            double col0[3] = {cp * cs, sp * cs, -ss};
            double col1[3] = {-sp, cp, 0};
            double col2[3] = {cp * ss, sp * ss, cs};
            double d[3] = {30.0 - 30.0 * cs, 0.0, 30.0 * ss};
            double t[3] = {cp * d[0] - sp * d[1], sp * d[0] + cp * d[1], d[2]};
            round.insert(round.end(), {col0[0], col0[1], col0[2], col1[0], col1[1], col1[2], col2[0], col2[1], col2[2], t[0], t[1], t[2]});
        }
        b = nullptr;
        expect("stations: twisted ring, joined",
               smlib_sweep_stations(ring.kinds.data(), ring.nk(), ring.reals.data(), ring.nr(), round.data(), m, SMLIB_ENDS_JOINED, &b), &b,
               64.0 * 2.0 * PI * 30.0, 1e-5);
        // Hostile maps: a singular matrix, a NaN, one station, unknown ends.
        std::vector<double> flat = twist;
        flat[12 * 5 + 0] = 0.0;
        flat[12 * 5 + 1] = 0.0;
        flat[12 * 5 + 2] = 0.0;
        b = nullptr;
        CHECK(smlib_sweep_stations(square.kinds.data(), square.nk(), square.reals.data(), square.nr(), flat.data(), n, SMLIB_ENDS_CAPPED, &b) ==
                  SMLIB_ERR_ARG &&
              b == nullptr);
        std::vector<double> nan_map = twist;
        nan_map[7] = NAN;
        CHECK(smlib_sweep_stations(square.kinds.data(), square.nk(), square.reals.data(), square.nr(), nan_map.data(), n, SMLIB_ENDS_CAPPED, &b) ==
                  SMLIB_ERR_ARG &&
              b == nullptr);
        CHECK(smlib_sweep_stations(square.kinds.data(), square.nk(), square.reals.data(), square.nr(), twist.data(), 1, SMLIB_ENDS_CAPPED, &b) ==
                  SMLIB_ERR_ARG &&
              b == nullptr);
        CHECK(smlib_sweep_stations(square.kinds.data(), square.nk(), square.reals.data(), square.nr(), twist.data(), n, 7, &b) == SMLIB_ERR_ARG &&
              b == nullptr);
        CHECK(smlib_sweep_stations(square.kinds.data(), square.nk(), square.reals.data(), square.nr(), nullptr, n, SMLIB_ENDS_CAPPED, &b) ==
                  SMLIB_ERR_ARG &&
              b == nullptr);
    }

    // Hostile lists: nothing is built, nothing crashes.
    const int64_t bad_kind[1] = {9};
    const double six[6] = {0, 0, 0, 1, 0, 0};
    b = nullptr;
    CHECK(smlib_extrude(bad_kind, 1, six, 6, Z, 1.0, 0.0, SMLIB_DRAFT_CORNER_TANGENT, 1, &b, nullptr) == SMLIB_ERR_ARG && b == nullptr);
    const int64_t line_kind[1] = {SMLIB_CURVE_LINE};
    CHECK(smlib_extrude(line_kind, 1, six, 5, Z, 1.0, 0.0, SMLIB_DRAFT_CORNER_TANGENT, 1, &b, nullptr) == SMLIB_ERR_ARG && b == nullptr);
    const double nan6[6] = {0, NAN, 0, 1, 0, 0};
    CHECK(smlib_extrude(line_kind, 1, nan6, 6, Z, 1.0, 0.0, SMLIB_DRAFT_CORNER_TANGENT, 1, &b, nullptr) == SMLIB_ERR_ARG && b == nullptr);
    const int64_t nurbs_kind[3] = {SMLIB_CURVE_NURBS, 3, 1000000000};
    CHECK(smlib_extrude(nurbs_kind, 3, six, 6, Z, 1.0, 0.0, SMLIB_DRAFT_CORNER_TANGENT, 1, &b, nullptr) == SMLIB_ERR_ARG && b == nullptr);
    CHECK(smlib_extrude(rect.kinds.data(), rect.nk(), rect.reals.data(), rect.nr(), Z, -1.0, 0.0, SMLIB_DRAFT_CORNER_TANGENT, 1, &b, nullptr) == SMLIB_ERR_ARG);
    CHECK(smlib_revolve(off.kinds.data(), off.nk(), off.reals.data(), off.nr(), O, Z, 7.0, 1, &b) == SMLIB_ERR_ARG);
    CHECK(smlib_loft(boxes.kinds.data(), boxes.nk(), boxes.reals.data(), boxes.nr(), four, 1, 1, 1, &b) == SMLIB_ERR_ARG);
    const int64_t wrong[2] = {4, 3};
    CHECK(smlib_loft(boxes.kinds.data(), boxes.nk(), boxes.reals.data(), boxes.nr(), wrong, 2, 1, 1, &b) == SMLIB_ERR_ARG);
    CHECK(smlib_screw_sweep(tooth.kinds.data(), tooth.nk(), tooth.reals.data(), tooth.nr(), O, Z, 0.0, 3.0, 1e-4, &b) == SMLIB_ERR_ARG);
    CHECK(b == nullptr);

    if (failures != 0)
    {
        std::fprintf(stderr, "%d failure(s)\n", failures);
        return 1;
    }
    std::printf("sweeps: ok\n");
    return 0;
}
