// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// The operations on bodies the C API cannot make yet (a cylinder through SM_API), and the
// SMLib defects found under ASan, each as a regression case. Run it under ASan/UBSan
// (README, "Fuzzing"): every case must stay clean there, not only return the right status.

#include "../src/smlib_shim.h"

#include <SmApiBrep.h>
#include <SmApiPrimitives.h>
#include <SmApiQueries.h>
#include <SmVector3d.h>

#include <cmath>
#include <cstdio>
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
    // returned a body whose edge tolerance had grown to the gap; now it fails. A boss r10 on a
    // 60 x 40 x 10 plate with a r2 round at its base, the boss wall out by 1; and a block whose
    // top is split in two, one piece up by 1.
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
        CHECK(smlib_offset_faces(rounded, wall, 1, 0.1, 1.0, &r, nullptr) != SMLIB_OK && r == nullptr);
        smlib_brep *block = box(0, 0, 0, 60, 40, 20);
        smlib_brep *cutter = box(30, -1, 10, 40, 42, 20);
        smlib_brep *split = nullptr;
        CHECK(smlib_imprint(block, cutter, &split) == SMLIB_OK);
        const double piece[3] = {45, 20, 20};
        CHECK(smlib_offset_faces(split, piece, 1, 0.1, 1.0, &r, nullptr) != SMLIB_OK && r == nullptr);
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

    // Thicken: a sheet (the cube's top face alone, SMLib's own way of making one is not in the
    // C API; a solid is refused).
    CHECK(smlib_thicken(cube, 1.0, 0, &r) == SMLIB_ERR_ARG);
    CHECK(smlib_thicken(nullptr, 1.0, 0, &r) == SMLIB_ERR_ARG);
    CHECK(smlib_thicken(cube, NAN, 0, &r) == SMLIB_ERR_ARG);

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
