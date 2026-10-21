// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// The sweeps (src/smlib_sweeps.cpp) against closed-form volumes: extrusions straight and drafted,
// revolutions whole, partial and touching the axis, sweeps along lines and arcs, lofts and screw
// sweeps (Pappus), each checked as a valid closed solid; then hostile curve lists.

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
    expect("extrude rectangle", smlib_extrude(rect.kinds.data(), rect.nk(), rect.reals.data(), rect.nr(), Z, 5.0, 0.0, 1, &b), &b, 300.0, 1e-9);
    const double down[3] = {0, 0, -1};
    b = nullptr;
    expect("extrude rectangle down", smlib_extrude(rect.kinds.data(), rect.nk(), rect.reals.data(), rect.nr(), down, 5.0, 0.0, 1, &b), &b, 300.0,
           1e-9);
    List ring = circle(0, 0, 0, 0, 0, 1, 1, 0, 0, 10);
    ring.append(circle(0, 0, 0, 0, 0, -1, 1, 0, 0, 4));
    b = nullptr;
    expect("extrude annulus", smlib_extrude(ring.kinds.data(), ring.nk(), ring.reals.data(), ring.nr(), Z, 3.0, 0.0, 1, &b), &b,
           PI * (100.0 - 16.0) * 3.0, 1e-6);
    List sq;
    sq.polygon({-5, -5, 0, 5, -5, 0, 5, 5, 0, -5, 5, 0});
    double t = std::tan(10.0 * PI / 180.0);
    double a2 = 10.0 - 2.0 * 4.0 * t;
    b = nullptr;
    expect("extrude square drafted 10 deg", smlib_extrude(sq.kinds.data(), sq.nk(), sq.reals.data(), sq.nr(), Z, 4.0, 10.0 * PI / 180.0, 1, &b), &b,
           4.0 / 3.0 * (100.0 + a2 * a2 + 10.0 * a2), 1e-6);
    double a3 = 10.0 + 2.0 * 4.0 * t;
    b = nullptr;
    expect("extrude square drafted -10 deg", smlib_extrude(sq.kinds.data(), sq.nk(), sq.reals.data(), sq.nr(), Z, 4.0, -10.0 * PI / 180.0, 1, &b),
           &b, 4.0 / 3.0 * (100.0 + a3 * a3 + 10.0 * a3), 1e-6);
    b = nullptr;
    expect("extrude square drafted down 10 deg", smlib_extrude(sq.kinds.data(), sq.nk(), sq.reals.data(), sq.nr(), down, 4.0, 10.0 * PI / 180.0, 1, &b),
           &b, 4.0 / 3.0 * (100.0 + a2 * a2 + 10.0 * a2), 1e-6);
    b = nullptr;
    expect("extrude square drafted down -10 deg", smlib_extrude(sq.kinds.data(), sq.nk(), sq.reals.data(), sq.nr(), down, 4.0, -10.0 * PI / 180.0, 1, &b),
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

    // Hostile lists: nothing is built, nothing crashes.
    const int64_t bad_kind[1] = {9};
    const double six[6] = {0, 0, 0, 1, 0, 0};
    b = nullptr;
    CHECK(smlib_extrude(bad_kind, 1, six, 6, Z, 1.0, 0.0, 1, &b) == SMLIB_ERR_ARG && b == nullptr);
    const int64_t line_kind[1] = {SMLIB_CURVE_LINE};
    CHECK(smlib_extrude(line_kind, 1, six, 5, Z, 1.0, 0.0, 1, &b) == SMLIB_ERR_ARG && b == nullptr);
    const double nan6[6] = {0, NAN, 0, 1, 0, 0};
    CHECK(smlib_extrude(line_kind, 1, nan6, 6, Z, 1.0, 0.0, 1, &b) == SMLIB_ERR_ARG && b == nullptr);
    const int64_t nurbs_kind[3] = {SMLIB_CURVE_NURBS, 3, 1000000000};
    CHECK(smlib_extrude(nurbs_kind, 3, six, 6, Z, 1.0, 0.0, 1, &b) == SMLIB_ERR_ARG && b == nullptr);
    CHECK(smlib_extrude(rect.kinds.data(), rect.nk(), rect.reals.data(), rect.nr(), Z, -1.0, 0.0, 1, &b) == SMLIB_ERR_ARG);
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
