// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// smlib_brep_mass_properties and smlib_brep_area against closed-form values: boxes, cylinders,
// spheres, cones, frustums and tori, at the origin and far from it. Prints the worst relative
// error of each quantity.

#include "../src/smlib_shim.h"

#include <SmApiPrimitives.h>

#include <cmath>
#include <cstdio>
#include <functional>

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

// Closed-form values: area, volume, centroid, Ixx Iyy Izz and Iyz Izx Ixy about the centroid.
struct Expected
{
    const char *name;
    std::function<SmApiStatus(SmBrep *&)> make;
    double v[11];
};

static Expected box(const char *name, SmVector3d o, double a, double b, double c)
{
    double V = a * b * c;
    return {name,
            [=](SmBrep *&p) { return SmApiCreateBox(o, a, b, c, p); },
            {2 * (a * b + b * c + c * a), V, o.x + a / 2, o.y + b / 2, o.z + c / 2, V * (b * b + c * c) / 12,
             V * (a * a + c * c) / 12, V * (a * a + b * b) / 12, 0, 0, 0}};
}

static Expected cylinder(const char *name, SmVector3d o, double r, double h)
{
    double V = PI * r * r * h;
    double Ixx = V * (3 * r * r + h * h) / 12;
    return {name,
            [=](SmBrep *&p) { return SmApiCreateCylinder(o, r, h, p); },
            {2 * PI * r * r + 2 * PI * r * h, V, o.x, o.y, o.z + h / 2, Ixx, Ixx, V * r * r / 2, 0, 0, 0}};
}

static Expected sphere(const char *name, SmVector3d o, double r)
{
    double V = 4.0 / 3.0 * PI * r * r * r;
    double I = 0.4 * V * r * r;
    return {name, [=](SmBrep *&p) { return SmApiCreateSphere(o, r, p); }, {4 * PI * r * r, V, o.x, o.y, o.z, I, I, I, 0, 0, 0}};
}

// A frustum (or a cone, rt = 0) on its base centre o, radius rb at z = 0 to rt at z = h. Its
// volume integrals are polynomials in z of degree at most 4, so three-point Gauss-Legendre gives
// them exactly.
static Expected cone(const char *name, SmVector3d o, double rb, double rt, double h)
{
    const double gx[3] = {-std::sqrt(0.6), 0.0, std::sqrt(0.6)}, gw[3] = {5.0 / 9, 8.0 / 9, 5.0 / 9};
    double V = 0, Mz = 0, Izz = 0, Zz = 0;
    for (int i = 0; i < 3; ++i)
    {
        double z = h / 2 * (1 + gx[i]), w = gw[i] * h / 2;
        double rho = rb + (rt - rb) * z / h, disc = PI * rho * rho;
        V += w * disc;
        Mz += w * disc * z;
        Zz += w * disc * z * z;
        Izz += w * PI * rho * rho * rho * rho / 2;
    }
    double zc = Mz / V;
    double Ixx = Izz / 2 + Zz - V * zc * zc;
    double s = std::sqrt((rb - rt) * (rb - rt) + h * h);
    return {name,
            [=](SmBrep *&p) { return SmApiCreateCone(o, rb, rt, h, p); },
            {PI * rb * rb + PI * rt * rt + PI * (rb + rt) * s, V, o.x, o.y, o.z + zc, Ixx, Ixx, Izz, 0, 0, 0}};
}

// A torus about the z axis through o (Pappus for area and volume).
static Expected torus(const char *name, SmVector3d o, double R, double r)
{
    double V = 2 * PI * PI * R * r * r;
    double Ixx = V * (R * R / 2 + 5 * r * r / 8);
    return {name,
            [=](SmBrep *&p) { return SmApiCreateTorus(o, R, r, p); },
            {4 * PI * PI * R * r, V, o.x, o.y, o.z, Ixx, Ixx, V * (R * R + 3 * r * r / 4), 0, 0, 0}};
}

int main()
{
    CHECK(smlib_init() == SMLIB_OK);
    const SmVector3d at0(0, 0, 0), far(1000.0, -2000.0, 500.0);
    Expected cases[] = {
        // tests/oracle (SolidCraft): 01-box, 02-cylinder, 38-coincident-cylinders.
        box("oracle 01-box 40x30x20", at0, 40, 30, 20),
        cylinder("oracle 02-cylinder r15 h50", at0, 15, 50),
        cylinder("oracle 38 cylinder r10 h50", at0, 10, 50),
        // The primitives of oracle 40-box-sphere-torus, whole.
        sphere("sphere r12", at0, 12),
        torus("torus R22 r4", at0, 22, 4),
        cone("frustum r10..4 h20", at0, 10, 4, 20),
        cone("cone r10 h20", at0, 10, 0, 20),
        cone("inverted frustum r3..9 h7", at0, 3, 9, 7),
        // The same far from the origin, and thin or flat ones.
        box("box far", far, 3, 50, 0.5),
        cylinder("cylinder far", far, 0.25, 80),
        sphere("sphere far", far, 7),
        torus("torus far", far, 30, 0.5),
        cone("frustum far", far, 6, 2, 9),
    };
    // Relative errors: area, volume and inertia against themselves; the centroid against the
    // body's box diagonal (as tests/oracle measures it), the products against its largest moment.
    const char *labels[] = {"area", "volume", "centroid", "inertia", "products"};
    double worst[5] = {0, 0, 0, 0, 0};
    std::printf("%-30s %10s %10s %10s %10s %10s\n", "body", labels[0], labels[1], labels[2], labels[3], labels[4]);
    for (const Expected &e : cases)
    {
        SmBrep *p = nullptr;
        CHECK(e.make(p) == SM_SUCCESS && p != nullptr);
        smlib_brep b{p};
        double m[11] = {0};
        int32_t s = smlib_brep_mass_properties(&b, 1e-8, m);
        CHECK(s == SMLIB_OK);
        double area = 0;
        CHECK(smlib_brep_area(&b, 1e-8, &area) == SMLIB_OK);
        SmExtent3d box;
        CHECK(p->CalculateBoundingBox(box) == SM_SUCCESS);
        double size = box.GetSize().Length();
        double imax = std::fmax(e.v[5], std::fmax(e.v[6], e.v[7]));
        double err[5] = {std::fabs(m[0] / e.v[0] - 1), std::fabs(m[1] / e.v[1] - 1), 0, 0, 0};
        err[0] = std::fmax(err[0], std::fabs(area / e.v[0] - 1));
        for (int i = 2; i < 5; ++i)
        {
            err[2] = std::fmax(err[2], std::fabs(m[i] - e.v[i]) / size);
        }
        for (int i = 5; i < 8; ++i)
        {
            err[3] = std::fmax(err[3], std::fabs(m[i] / e.v[i] - 1));
        }
        for (int i = 8; i < 11; ++i)
        {
            err[4] = std::fmax(err[4], std::fabs(m[i] - e.v[i]) / imax);
        }
        std::printf("%-30s", e.name);
        for (int i = 0; i < 5; ++i)
        {
            std::printf(" %10.1e", err[i]);
            worst[i] = std::fmax(worst[i], err[i]);
        }
        std::printf("\n");
        delete p;
    }
    // Asked for 1e-8; tests/oracle accepts 1e-3.
    for (int i = 0; i < 5; ++i)
    {
        if (!(worst[i] < 1e-6))
        {
            std::fprintf(stderr, "%s: worst relative error %.2e\n", labels[i], worst[i]);
            failures++;
        }
    }

    // A sheet has an area but no mass properties.
    SmBrep *sheet = nullptr;
    CHECK(SmApiCreatePlane(at0, 10.0, 20.0, sheet) == SM_SUCCESS && sheet != nullptr);
    smlib_brep sb{sheet};
    double m[11], a = 0;
    CHECK(smlib_brep_mass_properties(&sb, 1e-4, m) == SM_ERR_INVALID_INPUT);
    CHECK(smlib_brep_area(&sb, 1e-4, &a) == SMLIB_OK && std::fabs(a - 200.0) < 1e-9);
    delete sheet;

    // Rejected arguments.
    CHECK(smlib_brep_mass_properties(nullptr, 1e-4, m) == SMLIB_ERR_ARG);
    CHECK(smlib_brep_area(nullptr, 1e-4, &a) == SMLIB_ERR_ARG);
    smlib_brep *bx = nullptr;
    const double o3[3] = {0, 0, 0}, s3[3] = {1, 2, 3};
    CHECK(smlib_make_box(o3, s3, &bx) == SMLIB_OK);
    CHECK(smlib_brep_mass_properties(bx, NAN, m) == SMLIB_ERR_ARG);
    CHECK(smlib_brep_mass_properties(bx, 1e-4, nullptr) == SMLIB_ERR_ARG);
    CHECK(smlib_brep_area(bx, INFINITY, &a) == SMLIB_ERR_ARG);
    CHECK(smlib_brep_area(bx, 1e-4, nullptr) == SMLIB_ERR_ARG);
    smlib_brep_free(bx);

    if (failures == 0)
    {
        std::printf("smlib_brep_mass_properties test: ok\n");
    }
    return failures == 0 ? 0 : 1;
}
