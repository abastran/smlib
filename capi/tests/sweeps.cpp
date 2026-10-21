// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// The sweeps (src/smlib_sweeps.cpp) against closed-form volumes: extrusions straight and drafted,
// revolutions whole, partial and touching the axis, sweeps along lines and arcs, lofts and screw
// sweeps (Pappus), each checked as a valid closed solid; then hostile curve lists.

#include "smlib.h"

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

// Checks `b` (status `s`) is one valid solid of volume `v` within `rel`; frees it.
static void expect(const char *name, int32_t s, smlib_brep *b, double v, double rel)
{
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
    expect("extrude rectangle", smlib_extrude(rect.kinds.data(), rect.nk(), rect.reals.data(), rect.nr(), Z, 5.0, 0.0, 1, &b), b, 300.0, 1e-9);
    const double down[3] = {0, 0, -1};
    b = nullptr;
    expect("extrude rectangle down", smlib_extrude(rect.kinds.data(), rect.nk(), rect.reals.data(), rect.nr(), down, 5.0, 0.0, 1, &b), b, 300.0,
           1e-9);
    List ring = circle(0, 0, 0, 0, 0, 1, 1, 0, 0, 10);
    ring.append(circle(0, 0, 0, 0, 0, -1, 1, 0, 0, 4));
    b = nullptr;
    expect("extrude annulus", smlib_extrude(ring.kinds.data(), ring.nk(), ring.reals.data(), ring.nr(), Z, 3.0, 0.0, 1, &b), b,
           PI * (100.0 - 16.0) * 3.0, 1e-6);
    List sq;
    sq.polygon({-5, -5, 0, 5, -5, 0, 5, 5, 0, -5, 5, 0});
    double t = std::tan(10.0 * PI / 180.0);
    double a2 = 10.0 - 2.0 * 4.0 * t;
    b = nullptr;
    expect("extrude square drafted 10 deg", smlib_extrude(sq.kinds.data(), sq.nk(), sq.reals.data(), sq.nr(), Z, 4.0, 10.0 * PI / 180.0, 1, &b), b,
           4.0 / 3.0 * (100.0 + a2 * a2 + 10.0 * a2), 1e-6);
    double a3 = 10.0 + 2.0 * 4.0 * t;
    b = nullptr;
    expect("extrude square drafted -10 deg", smlib_extrude(sq.kinds.data(), sq.nk(), sq.reals.data(), sq.nr(), Z, 4.0, -10.0 * PI / 180.0, 1, &b),
           b, 4.0 / 3.0 * (100.0 + a3 * a3 + 10.0 * a3), 1e-6);
    b = nullptr;
    expect("extrude square drafted down 10 deg", smlib_extrude(sq.kinds.data(), sq.nk(), sq.reals.data(), sq.nr(), down, 4.0, 10.0 * PI / 180.0, 1, &b),
           b, 4.0 / 3.0 * (100.0 + a2 * a2 + 10.0 * a2), 1e-6);
    b = nullptr;
    expect("extrude square drafted down -10 deg", smlib_extrude(sq.kinds.data(), sq.nk(), sq.reals.data(), sq.nr(), down, 4.0, -10.0 * PI / 180.0, 1, &b),
           b, 4.0 / 3.0 * (100.0 + a3 * a3 + 10.0 * a3), 1e-6);

    // Revolutions about z of a profile in the xz plane: a square off the axis (a ring), a
    // quarter turn of it, and a half disc on the axis (a sphere).
    const double Y[3] = {0, 1, 0};
    List off;
    off.polygon({6, 0, 0, 8, 0, 0, 8, 0, 2, 6, 0, 2});
    b = nullptr;
    expect("revolve square full", smlib_revolve(off.kinds.data(), off.nk(), off.reals.data(), off.nr(), O, Z, 2 * PI, 1, &b), b,
           2 * PI * 7.0 * 4.0, 1e-6);
    b = nullptr;
    expect("revolve square 90", smlib_revolve(off.kinds.data(), off.nk(), off.reals.data(), off.nr(), O, Z, PI / 2, 1, &b), b,
           PI / 2 * 7.0 * 4.0, 1e-6);
    List half;
    half.line(0, 0, -5, 0, 0, 5).arc(0, 0, 0, 0, -1, 0, 0, 0, 1, 5, 0, PI);
    b = nullptr;
    expect("revolve half disc", smlib_revolve(half.kinds.data(), half.nk(), half.reals.data(), half.nr(), O, Z, 2 * PI, 1, &b), b,
           4.0 / 3.0 * PI * 125.0, 1e-5);
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
           b, PI * 4.0 * 10.0, 1e-4);
    List bend;
    bend.arc(-20, 0, 0, 0, -1, 0, 1, 0, 0, 20, 0, PI / 2);
    b = nullptr;
    expect("sweep disc along quarter arc",
           smlib_sweep(disc.kinds.data(), disc.nk(), disc.reals.data(), disc.nr(), bend.kinds.data(), bend.nk(), bend.reals.data(), bend.nr(),
                       SMLIB_SWEEP_FOLLOW, 1, &b),
           b, PI * 4.0 * (PI / 2 * 20.0), 1e-3);
    List offset = circle(5, 0, 0, 0, 0, 1, 1, 0, 0, 2);
    b = nullptr;
    expect("sweep offset disc along arc",
           smlib_sweep(offset.kinds.data(), offset.nk(), offset.reals.data(), offset.nr(), bend.kinds.data(), bend.nk(), bend.reals.data(),
                       bend.nr(), SMLIB_SWEEP_FOLLOW, 1, &b),
           b, PI * 4.0 * (PI / 2 * 25.0), 1e-3);

    // Lofts: square to square (a box), circle to circle (a frustum, ruled).
    List boxes;
    boxes.polygon({0, 0, 0, 10, 0, 0, 10, 10, 0, 0, 10, 0}).polygon({0, 0, 8, 10, 0, 8, 10, 10, 8, 0, 10, 8});
    const int64_t four[2] = {4, 4};
    b = nullptr;
    expect("loft square to square", smlib_loft(boxes.kinds.data(), boxes.nk(), boxes.reals.data(), boxes.nr(), four, 2, 1, 1, &b), b, 800.0,
           1e-6);
    List cone = circle(0, 0, 0, 0, 0, 1, 1, 0, 0, 5);
    cone.append(circle(0, 0, 6, 0, 0, 1, 1, 0, 0, 3));
    const int64_t one[2] = {1, 1};
    b = nullptr;
    expect("loft circle to circle", smlib_loft(cone.kinds.data(), cone.nk(), cone.reals.data(), cone.nr(), one, 2, 1, 1, &b), b,
           PI * 6.0 / 3.0 * (25.0 + 9.0 + 15.0), 1e-4);

    // Screw sweeps: a 0.4 x 0.4 square at radius 4..4.4 in the xz plane, 3 turns of pitch 1,
    // up and down (Pappus: area x 2 pi x centroid radius x turns).
    List tooth;
    tooth.polygon({4, 0, 0, 4.4, 0, 0, 4.4, 0, 0.4, 4, 0, 0.4});
    double screw = 0.16 * 2 * PI * 4.2 * 3.0;
    b = nullptr;
    expect("screw sweep up", smlib_screw_sweep(tooth.kinds.data(), tooth.nk(), tooth.reals.data(), tooth.nr(), O, Z, 1.0, 3.0, 1e-4, &b), b,
           screw, 2e-3);
    b = nullptr;
    expect("screw sweep down", smlib_screw_sweep(tooth.kinds.data(), tooth.nk(), tooth.reals.data(), tooth.nr(), O, Z, -1.0, 3.0, 1e-4, &b), b,
           screw, 2e-3);

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
