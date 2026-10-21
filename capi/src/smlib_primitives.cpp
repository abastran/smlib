// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// The primitive solids of include/smlib.h beyond the box: cylinder, sphere, cone and torus, each
// SM_API's own constructor behind the same argument checks and Call as smlib_make_box.

#include "smlib_shim.h"

#include <SmApiPrimitives.h>
#include <SmVector3d.h>

using namespace smlib_shim;

namespace
{

// Runs one SM_API constructor that returns a new body through its last argument.
template <typename F>
int32_t make(F &&f, smlib_brep **out)
{
    SmBrep *r = nullptr;
    int32_t s = status(f(r));
    if (s != SMLIB_OK)
    {
        delete r;
        return s;
    }
    return take(r, out);
}

bool positive(double x)
{
    return std::isfinite(x) && x > 0.0;
}

} // namespace

extern "C" int32_t smlib_make_cylinder(const double origin[3], double radius, double height, smlib_brep **out)
{
    if (origin == nullptr || out == nullptr || !finite(origin, 3) || !positive(radius) || !positive(height))
    {
        return SMLIB_ERR_ARG;
    }
    return Call("smlib_make_cylinder", "origin=" + vec(origin, 3) + ", radius=" + num(radius) + ", height=" + num(height))
        .run([&] { return make([&](SmBrep *&r) { return SmApiCreateCylinder(SmVector3d(origin[0], origin[1], origin[2]), radius, height, r); }, out); });
}

extern "C" int32_t smlib_make_sphere(const double center[3], double radius, smlib_brep **out)
{
    if (center == nullptr || out == nullptr || !finite(center, 3) || !positive(radius))
    {
        return SMLIB_ERR_ARG;
    }
    return Call("smlib_make_sphere", "center=" + vec(center, 3) + ", radius=" + num(radius))
        .run([&] { return make([&](SmBrep *&r) { return SmApiCreateSphere(SmVector3d(center[0], center[1], center[2]), radius, r); }, out); });
}

extern "C" int32_t smlib_make_cone(const double origin[3], double base_radius, double top_radius, double height, smlib_brep **out)
{
    if (origin == nullptr || out == nullptr || !finite(origin, 3) || !positive(base_radius) || !std::isfinite(top_radius) ||
        !(top_radius >= 0.0) || top_radius == base_radius || !positive(height))
    {
        return SMLIB_ERR_ARG;
    }
    return Call(
               "smlib_make_cone",
               "origin=" + vec(origin, 3) + ", base_radius=" + num(base_radius) + ", top_radius=" + num(top_radius) + ", height=" + num(height)
    )
        .run(
            [&]
            {
                return make(
                    [&](SmBrep *&r) { return SmApiCreateCone(SmVector3d(origin[0], origin[1], origin[2]), base_radius, top_radius, height, r); }, out
                );
            }
        );
}

extern "C" int32_t smlib_make_torus(const double center[3], double major_radius, double minor_radius, smlib_brep **out)
{
    if (center == nullptr || out == nullptr || !finite(center, 3) || !positive(major_radius) || !positive(minor_radius) ||
        !(minor_radius < major_radius))
    {
        return SMLIB_ERR_ARG;
    }
    return Call("smlib_make_torus", "center=" + vec(center, 3) + ", major_radius=" + num(major_radius) + ", minor_radius=" + num(minor_radius))
        .run(
            [&]
            {
                return make(
                    [&](SmBrep *&r) { return SmApiCreateTorus(SmVector3d(center[0], center[1], center[2]), major_radius, minor_radius, r); }, out
                );
            }
        );
}
