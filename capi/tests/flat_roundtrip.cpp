// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// Flat B-rep arrays through the C API: for a set of solids, export -> import -> compare volume,
// counts and manifoldness, then export again and compare the arrays with the first export.
// Then hostile arrays: truncated, out of range, NaN, huge counts, random mutations; each must
// return an error status, never crash.
//
// The solids are built with SM_API directly and handed to the C API as native .smb files.
//
//   smlib_flat_roundtrip <scratch dir>

#include "smlib.h"

#include <SmApiBrep.h>
#include <SmApiGeneral.h>
#include <SmApiImportExport.h>
#include <SmApiPrimitives.h>
#include <SmApiQueries.h>
#include <SmBrep.h>
#include <SmEdge.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <map>
#include <string>
#include <vector>

namespace
{

int failures = 0;
std::string dir = ".";

#define CHECK(cond)                                                                                \
    do                                                                                             \
    {                                                                                              \
        if (!(cond))                                                                               \
        {                                                                                          \
            std::fprintf(stderr, "%s:%d: CHECK(%s) failed\n", __FILE__, __LINE__, #cond);          \
            failures++;                                                                            \
        }                                                                                          \
    } while (0)

// Every id: int64 ids 1..32, double ids 101..163.
std::vector<int32_t> IntIds()
{
    std::vector<int32_t> v;
    for (int32_t id = SMLIB_FLAT_REGION_SHELL_COUNT; id <= SMLIB_FLAT_SOURCE; ++id)
    {
        v.push_back(id);
    }
    return v;
}

std::vector<int32_t> RealIds()
{
    std::vector<int32_t> v;
    for (int32_t id = SMLIB_FLAT_TOLERANCE; id <= SMLIB_FLAT_VERTEX_TOLERANCE; ++id)
    {
        v.push_back(id);
    }
    return v;
}

// A copy of every array of a flat set.
struct Arrays
{
    std::map<int32_t, std::vector<int64_t>> i;
    std::map<int32_t, std::vector<double>> f;
};

Arrays Snapshot(const smlib_flat *flat)
{
    Arrays a;
    for (int32_t id : IntIds())
    {
        const int64_t *d = nullptr;
        int64_t n = 0;
        CHECK(smlib_flat_get_i64(flat, id, &d, &n) == SMLIB_OK);
        a.i[id] = std::vector<int64_t>(d, d + n);
    }
    for (int32_t id : RealIds())
    {
        const double *d = nullptr;
        int64_t n = 0;
        CHECK(smlib_flat_get_f64(flat, id, &d, &n) == SMLIB_OK);
        a.f[id] = std::vector<double>(d, d + n);
    }
    return a;
}

smlib_flat *FromArrays(const Arrays &a)
{
    smlib_flat *flat = smlib_flat_new();
    for (auto &kv : a.i)
    {
        CHECK(smlib_flat_set_i64(flat, kv.first, kv.second.data(), static_cast<int64_t>(kv.second.size())) == SMLIB_OK);
    }
    for (auto &kv : a.f)
    {
        CHECK(smlib_flat_set_f64(flat, kv.first, kv.second.data(), static_cast<int64_t>(kv.second.size())) == SMLIB_OK);
    }
    return flat;
}

// Ints exactly, reals within 1e-9 relative (1e-9 absolute below magnitude 1). Returns the
// largest relative real difference, or infinity on a mismatch of sizes or ints.
double Compare(const Arrays &a, const Arrays &b, const char *what)
{
    double worst = 0.0;
    for (auto &kv : a.i)
    {
        if (kv.second != b.i.at(kv.first))
        {
            std::fprintf(stderr, "  %s: int array %d differs\n", what, kv.first);
            return INFINITY;
        }
    }
    for (auto &kv : a.f)
    {
        const std::vector<double> &x = kv.second, &y = b.f.at(kv.first);
        // EXTENT is SMLib's loose (non-tight) bounding box, informational and ignored by import;
        // it can differ between a B-rep and its re-import (blend1: the blend face's box comes
        // out 4.6e-5 wider), with an all-NURBS export too. Its length still has to match.
        if (kv.first == SMLIB_FLAT_EXTENT && x.size() == y.size())
        {
            continue;
        }
        // An edge's or vertex's own tolerance is what SMLib finds it needs: building the edge
        // again from the exported curves may need more (blend1: 2e-5 then 1.1e-4), never less.
        if ((kv.first == SMLIB_FLAT_EDGE_TOLERANCE || kv.first == SMLIB_FLAT_VERTEX_TOLERANCE) && x.size() == y.size())
        {
            for (size_t k = 0; k < x.size(); ++k)
            {
                if (y[k] < x[k] * (1.0 - 1e-9))
                {
                    std::fprintf(stderr, "  %s: real array %d [%zu]: %.17g then less, %.17g\n", what, kv.first, k, x[k], y[k]);
                    return INFINITY;
                }
            }
            continue;
        }
        if (x.size() != y.size())
        {
            std::fprintf(stderr, "  %s: real array %d has %zu then %zu values\n", what, kv.first, x.size(), y.size());
            return INFINITY;
        }
        for (size_t k = 0; k < x.size(); ++k)
        {
            double d = std::fabs(x[k] - y[k]) / std::max({1.0, std::fabs(x[k]), std::fabs(y[k])});
            if (d > worst)
            {
                worst = d;
            }
            if (d > 1e-9)
            {
                std::fprintf(stderr, "  %s: real array %d [%zu]: %.17g vs %.17g\n", what, kv.first, k, x[k], y[k]);
                return INFINITY;
            }
        }
    }
    return worst;
}

// Hands an SM_API B-rep to the C API through a native file.
smlib_brep *ToCApi(SmBrep *b, const char *name)
{
    if (b == nullptr)
    {
        return nullptr;
    }
    std::string path = dir + "/" + name + ".smb";
    CHECK(SmApiWriteBrepToFile(b, path.c_str(), FALSE) == SM_SUCCESS);
    delete b;
    smlib_brep *r = nullptr;
    CHECK(smlib_brep_read(path.c_str(), 0, &r) == SMLIB_OK);
    return r;
}

struct Facts
{
    int64_t counts[3] = {0, 0, 0};
    double volume = 0.0;
    int32_t manifold = -1;
};

Facts Measure(const smlib_brep *b)
{
    Facts f;
    CHECK(smlib_brep_counts(b, f.counts) == SMLIB_OK);
    CHECK(smlib_brep_volume(b, 1e-7, &f.volume) == SMLIB_OK);
    CHECK(smlib_brep_is_manifold_solid(b, &f.manifold) == SMLIB_OK);
    return f;
}

// Expected curve kinds per solid (nurbs, line, circle, ellipse). Export keeps SMLib's own curve
// classes: SmApiCreateCylinder, SmApiCreateCone and SmApiCreateTorus build their rims and seam
// circles as plain SmBSplineCurves (exactly circular, but not SmCircles), so those stay NURBS,
// as do a boolean's intersection curves.
const std::map<std::string, std::vector<int64_t>> g_kinds = {
    {"box", {0, 12, 0, 0}},      // twelve SmLines
    {"cylinder", {2, 1, 0, 0}},  // B-spline rims, a line seam
    {"cone", {2, 1, 0, 0}},      // B-spline rims, a line seam
    {"sphere", {0, 0, 1, 0}},    // the seam: an SmCircle half circle
    {"torus", {2, 0, 0, 0}},     // B-spline seam circles
    {"fillet", {0, 24, 24, 0}},  // the box's trimmed lines, the fillets' SmCircle arcs
    {"blend", {48, 0, 0, 0}},    // blend-curve fillets: all B-splines
    {"blend1", {4, 11, 0, 0}},   // the blend face's four B-spline edges, the box's eleven lines
    {"boolean", {2, 13, 0, 0}},  // the hole's intersection curves are B-splines
    {"incline", {1, 1, 0, 2}},   // B-spline bottom rim, line seam, two SmEllipse arcs
};

size_t g_analyticTrims = 0;

// The round trip; returns the first export (caller frees) for the hostile cases.
//
// `normalizes`: the import rewrites this solid into an equivalent canonical form, so the first
// re-export differs from the first export. SmuConvert's import does this to cones with a negative
// semi-angle (SmApiCreateCone's narrowing cone): it reverses the axis and both STEP parameters
// and rebuilds the cone from Vmin. Then the volume (integrated at 1e-7 over a different
// parameterisation) need only agree to 1e-7, and the arrays must reach a fixed point: the second
// and third exports must match.
//
// `volumeTolerance`: 1e-9 unless an edge curve is rebuilt with another parameterisation, which
// SMLib then integrates (and trims) differently, by up to its 1e-7 accuracy: the cone (above),
// the sphere, whose seam SmCircle is rebuilt over [-90, 90] degrees instead of [0, 1] (its
// volume moves by 1.9e-8 at accuracy 1e-7 and 8.5e-8 at 1e-8, where an all-NURBS export keeps it
// exactly), and the inclined cut, whose SmEllipse arcs are rebuilt likewise and whose cylinder
// trims SMLib rebuilds. The geometry is the same: export checks the rebuilt curve against the
// original pointwise.
smlib_flat *RoundTrip(const char *name, smlib_brep *b, bool normalizes = false, double volumeTolerance = 1e-9)
{
    if (b == nullptr)
    {
        std::fprintf(stderr, "%s: no solid\n", name);
        failures++;
        return nullptr;
    }
    Facts f0 = Measure(b);

    smlib_flat *flat1 = nullptr;
    int32_t s = smlib_brep_export(b, &flat1);
    CHECK(s == SMLIB_OK);
    if (s != SMLIB_OK)
    {
        std::fprintf(stderr, "%s: export -> %d\n", name, s);
        return nullptr;
    }
    smlib_brep *b2 = nullptr;
    s = smlib_brep_import(flat1, &b2);
    CHECK(s == SMLIB_OK);
    if (s != SMLIB_OK)
    {
        std::fprintf(stderr, "%s: import -> %d\n", name, s);
        return flat1;
    }
    Facts f1 = Measure(b2);
    double v8a = 0.0, v8b = 0.0;
    smlib_brep_volume(b, 1e-8, &v8a);
    smlib_brep_volume(b2, 1e-8, &v8b);
    std::printf("%-9s volume at 1e-8: %.12g -> %.12g (rel %.2e)\n", name, v8a, v8b, std::fabs(v8b - v8a) / std::max(1e-300, std::fabs(v8a)));
    double dv = std::fabs(f1.volume - f0.volume) / std::fabs(f0.volume);
    CHECK(dv <= volumeTolerance);
    CHECK(f1.counts[0] == f0.counts[0] && f1.counts[1] == f0.counts[1] && f1.counts[2] == f0.counts[2]);
    CHECK(f0.manifold == 1 && f1.manifold == 1);

    smlib_flat *flat2 = nullptr;
    CHECK(smlib_brep_export(b2, &flat2) == SMLIB_OK);
    Arrays a1 = Snapshot(flat1), a2 = Snapshot(flat2);
    double worst = 0.0;
    if (normalizes)
    {
        smlib_brep *b3 = nullptr;
        smlib_flat *flat4 = nullptr;
        CHECK(smlib_brep_import(flat2, &b3) == SMLIB_OK);
        CHECK(b3 != nullptr && smlib_brep_export(b3, &flat4) == SMLIB_OK);
        if (flat4 != nullptr)
        {
            Facts f2 = Measure(b3);
            CHECK(std::fabs(f2.volume - f1.volume) <= 1e-9 * std::fabs(f1.volume));
            worst = Compare(a2, Snapshot(flat4), name);
        }
        smlib_flat_free(flat4);
        smlib_brep_free(b3);
    }
    else
    {
        worst = Compare(a1, a2, name);
    }
    CHECK(worst <= 1e-9);

    // The text form reads back to the same arrays, exactly.
    std::string path = dir + "/" + name + ".flat.txt";
    CHECK(smlib_flat_write(flat1, path.c_str()) == SMLIB_OK);
    smlib_flat *flat3 = nullptr;
    CHECK(smlib_flat_read(path.c_str(), &flat3) == SMLIB_OK);
    if (flat3 != nullptr)
    {
        Arrays a3 = Snapshot(flat3);
        CHECK(a3.i == a1.i && a3.f == a1.f);
    }

    // Curve kinds of the first export: NURBS, line, circle, ellipse (edges and wire edges).
    int64_t kinds[4] = {0, 0, 0, 0};
    for (int32_t id : {SMLIB_FLAT_EDGE_CURVE_TYPE, SMLIB_FLAT_WIRE_EDGE_CURVE_TYPE})
    {
        for (int64_t t : a1.i[id])
        {
            CHECK(t >= 0 && t < 4);
            if (t >= 0 && t < 4)
            {
                kinds[t]++;
            }
        }
    }
    auto expected = g_kinds.find(name);
    if (expected != g_kinds.end())
    {
        const std::vector<int64_t> &k = expected->second;
        if (!(kinds[0] == k[0] && kinds[1] == k[1] && kinds[2] == k[2] && kinds[3] == k[3]))
        {
            std::fprintf(stderr, "%s: curve kinds nurbs/line/circle/ellipse %lld/%lld/%lld/%lld, expected %lld/%lld/%lld/%lld\n", name,
                         (long long)kinds[0], (long long)kinds[1], (long long)kinds[2], (long long)kinds[3], (long long)k[0], (long long)k[1],
                         (long long)k[2], (long long)k[3]);
            failures++;
        }
    }

    // UV trims carried on edges exported as analytic curves (their knots were re-mapped).
    size_t nAnalyticTrims = 0;
    if (!a1.i[SMLIB_FLAT_EDGEUSE_NURBS_VERTEX_COUNT].empty())
    {
        for (size_t k = 0; k < a1.i[SMLIB_FLAT_EDGEUSE_EDGE].size(); ++k)
        {
            nAnalyticTrims += a1.i[SMLIB_FLAT_EDGEUSE_NURBS_VERTEX_COUNT][k] != 0 &&
                              a1.i[SMLIB_FLAT_EDGE_CURVE_TYPE][static_cast<size_t>(a1.i[SMLIB_FLAT_EDGEUSE_EDGE][k])] != SMLIB_FLAT_CURVE_NURBS;
        }
    }
    g_analyticTrims += nAnalyticTrims;

    size_t nNurbsFaces = a1.i[SMLIB_FLAT_FACE_NURBS_U_ORDER].size();
    size_t nUV = 0;
    for (int64_t n : a1.i[SMLIB_FLAT_EDGEUSE_NURBS_VERTEX_COUNT])
    {
        nUV += n != 0;
    }
    std::printf(
        "%-9s faces/edges/vertices %lld/%lld/%lld -> %lld/%lld/%lld, manifold %d -> %d, volume %.12g -> %.12g (rel %.2e), "
        "NURBS faces %zu, UV trims %zu (%zu on analytic edges), curves nurbs/line/circle/ellipse %lld/%lld/%lld/%lld, re-export max rel diff %.2e\n",
        name, (long long)f0.counts[0], (long long)f0.counts[1], (long long)f0.counts[2], (long long)f1.counts[0], (long long)f1.counts[1],
        (long long)f1.counts[2], f0.manifold, f1.manifold, f0.volume, f1.volume, dv, nNurbsFaces, nUV, nAnalyticTrims, (long long)kinds[0], (long long)kinds[1],
        (long long)kinds[2], (long long)kinds[3], worst
    );

    smlib_flat_free(flat3);
    smlib_flat_free(flat2);
    smlib_brep_free(b2);
    smlib_brep_free(b);
    return flat1;
}

// A region's shells may come in any order: import puts the one that encloses the others first,
// as SMLib requires (solidcore#388). An island inside a void: a cube of 50 hollowed by a cube of
// 40, joined with a cube of 20 inside the cavity, which a cube of 10 hollows in turn. Each order
// of the regions (the infinite one first) and of the shells within each region must import as a
// valid body of the same volume.
void NestedShellOrder()
{
    auto box = [](double lo, double size) {
        const double o[3] = {lo, lo, lo}, s[3] = {size, size, size};
        smlib_brep *b = nullptr;
        CHECK(smlib_make_box(o, s, &b) == SMLIB_OK);
        return b;
    };
    smlib_brep *b50 = box(-25, 50), *b40 = box(-20, 40), *b20 = box(-10, 20), *b10 = box(-5, 10);
    smlib_brep *hollow = nullptr, *joined = nullptr, *nested = nullptr;
    CHECK(smlib_boolean(b50, b40, SMLIB_BOOLEAN_SUBTRACT, &hollow) == SMLIB_OK);
    CHECK(hollow && smlib_boolean(hollow, b20, SMLIB_BOOLEAN_UNION, &joined) == SMLIB_OK);
    CHECK(joined && smlib_boolean(joined, b10, SMLIB_BOOLEAN_SUBTRACT, &nested) == SMLIB_OK);
    smlib_flat *flat = nullptr;
    CHECK(nested && smlib_brep_export(nested, &flat) == SMLIB_OK);
    for (smlib_brep *b : {b50, b40, b20, b10, hollow, joined, nested})
    {
        smlib_brep_free(b);
    }
    if (flat == nullptr)
    {
        return;
    }
    const double want = 50.0 * 50 * 50 - 40.0 * 40 * 40 + 20.0 * 20 * 20 - 10.0 * 10 * 10;
    const Arrays a = Snapshot(flat);
    smlib_flat_free(flat);

    // Each shell's block: its counts and its faceuses. Faceuse shells only here.
    struct Shell
    {
        int64_t faceuses, wires, point;
        std::vector<int64_t> face, orientation;
    };
    const std::vector<int64_t> &regionShells = a.i.at(SMLIB_FLAT_REGION_SHELL_COUNT);
    std::vector<std::vector<Shell>> regions;
    size_t shell = 0, fu = 0;
    for (int64_t n : regionShells)
    {
        regions.emplace_back();
        for (int64_t k = 0; k < n; ++k, ++shell)
        {
            Shell s;
            s.faceuses = a.i.at(SMLIB_FLAT_SHELL_FACEUSE_COUNT)[shell];
            s.wires = a.i.at(SMLIB_FLAT_SHELL_WIRE_EDGE_COUNT)[shell];
            s.point = a.i.at(SMLIB_FLAT_SHELL_POINT_TYPE)[shell];
            CHECK(s.faceuses > 0 && s.wires == 0 && s.point == 0);
            const std::vector<int64_t> &face = a.i.at(SMLIB_FLAT_FACEUSE_FACE), &orient = a.i.at(SMLIB_FLAT_FACEUSE_ORIENTATION);
            s.face.assign(face.begin() + static_cast<ptrdiff_t>(fu), face.begin() + static_cast<ptrdiff_t>(fu + s.faceuses));
            s.orientation.assign(orient.begin() + static_cast<ptrdiff_t>(fu), orient.begin() + static_cast<ptrdiff_t>(fu + s.faceuses));
            fu += static_cast<size_t>(s.faceuses);
            regions.back().push_back(s);
        }
    }
    // The infinite region with the 50 cube's outside; the 50 hollow; the void between 40 and 20;
    // the 20 hollow; the void inside 10.
    CHECK(regions.size() == 5);
    size_t multi = 0;
    for (auto &r : regions)
    {
        multi += r.size() > 1;
    }
    CHECK(multi == 3);

    // Every order: regions 1.. permuted, then the shells of every region.
    std::vector<size_t> regionOrder;
    for (size_t r = 1; r < regions.size(); ++r)
    {
        regionOrder.push_back(r);
    }
    int orders = 0, good = 0;
    do
    {
        std::vector<std::vector<size_t>> shellOrder(regions.size());
        for (size_t r = 0; r < regions.size(); ++r)
        {
            for (size_t k = 0; k < regions[r].size(); ++k)
            {
                shellOrder[r].push_back(k);
            }
        }
        for (;;)
        {
            Arrays p = a;
            std::vector<int64_t> &rsc = p.i[SMLIB_FLAT_REGION_SHELL_COUNT], &rt = p.i[SMLIB_FLAT_REGION_TYPE];
            std::vector<int64_t> &sfc = p.i[SMLIB_FLAT_SHELL_FACEUSE_COUNT], &swc = p.i[SMLIB_FLAT_SHELL_WIRE_EDGE_COUNT],
                                 &spt = p.i[SMLIB_FLAT_SHELL_POINT_TYPE], &ff = p.i[SMLIB_FLAT_FACEUSE_FACE],
                                 &fo = p.i[SMLIB_FLAT_FACEUSE_ORIENTATION];
            rsc.clear();
            rt.clear();
            sfc.clear();
            swc.clear();
            spt.clear();
            ff.clear();
            fo.clear();
            std::vector<size_t> order = {0};
            order.insert(order.end(), regionOrder.begin(), regionOrder.end());
            for (size_t r : order)
            {
                rsc.push_back(static_cast<int64_t>(regions[r].size()));
                rt.push_back(a.i.at(SMLIB_FLAT_REGION_TYPE)[r]);
                for (size_t k : shellOrder[r])
                {
                    const Shell &s = regions[r][k];
                    sfc.push_back(s.faceuses);
                    swc.push_back(s.wires);
                    spt.push_back(s.point);
                    ff.insert(ff.end(), s.face.begin(), s.face.end());
                    fo.insert(fo.end(), s.orientation.begin(), s.orientation.end());
                }
            }
            smlib_flat *f = FromArrays(p);
            smlib_brep *b = nullptr;
            int32_t s = smlib_brep_import(f, &b);
            double v = 0.0;
            int64_t len = 0, bad = -1;
            char buf[4096] = {0};
            if (s == SMLIB_OK)
            {
                CHECK(smlib_brep_volume(b, 1e-9, &v) == SMLIB_OK);
                CHECK(smlib_brep_check(b, 2, buf, sizeof buf, &len, &bad) == SMLIB_OK);
            }
            orders++;
            if (s == SMLIB_OK && bad == 0 && std::fabs(v - want) <= 1e-9 * want)
            {
                good++;
            }
            else
            {
                std::fprintf(stderr, "nested shells, order %d: import %d, volume %.12g (want %.12g), failed checks %lld %s\n", orders,
                             static_cast<int>(s), v, want, static_cast<long long>(bad), buf);
                failures++;
            }
            smlib_brep_free(b);
            smlib_flat_free(f);

            // the next shell order, odometer-wise over the regions
            size_t r = 0;
            while (r < shellOrder.size() && !std::next_permutation(shellOrder[r].begin(), shellOrder[r].end()))
            {
                r++;
            }
            if (r == shellOrder.size())
            {
                break;
            }
        }
    } while (std::next_permutation(regionOrder.begin(), regionOrder.end()));
    std::printf("nested shells: %d of %d region and shell orders import as a valid body of volume %.12g\n", good, orders, want);
}

// Imports `a`; expects a failure status (and no crash). Returns the status.
int32_t ExpectReject(const Arrays &a, const char *what)
{
    smlib_flat *flat = FromArrays(a);
    smlib_brep *b = nullptr;
    int32_t s = smlib_brep_import(flat, &b);
    if (s == SMLIB_OK)
    {
        std::fprintf(stderr, "  hostile case accepted: %s\n", what);
        failures++;
    }
    CHECK(s == SMLIB_OK || b == nullptr);
    smlib_brep_free(b);
    smlib_flat_free(flat);
    return s;
}

// Truncations, bad indices, bad values and huge counts, all of which must be rejected.
int Hostile(const char *name, const smlib_flat *flat)
{
    const Arrays base = Snapshot(flat);
    int cases = 0;
    char what[128];

    // Each non-empty array one element short (SOURCE may be empty: it is optional).
    for (auto &kv : base.i)
    {
        if (kv.second.empty() || kv.first == SMLIB_FLAT_SOURCE)
        {
            continue;
        }
        Arrays a = base;
        a.i[kv.first].pop_back();
        std::snprintf(what, sizeof what, "%s: int array %d truncated", name, kv.first);
        ExpectReject(a, what);
        cases++;
    }
    // (EDGE_ and VERTEX_TOLERANCE may be empty: they are optional.)
    for (auto &kv : base.f)
    {
        bool optional = kv.first == SMLIB_FLAT_EDGE_TOLERANCE || kv.first == SMLIB_FLAT_VERTEX_TOLERANCE;
        if (kv.second.empty() || (optional && kv.second.size() == 1))
        {
            continue;
        }
        Arrays a = base;
        a.f[kv.first].pop_back();
        std::snprintf(what, sizeof what, "%s: real array %d truncated", name, kv.first);
        ExpectReject(a, what);
        cases++;
    }

    // Index, kind and count arrays: values out of range and huge.
    const int32_t indexIds[] = {SMLIB_FLAT_REGION_SHELL_COUNT, SMLIB_FLAT_REGION_TYPE,       SMLIB_FLAT_SHELL_FACEUSE_COUNT,
                                SMLIB_FLAT_FACEUSE_FACE,       SMLIB_FLAT_FACEUSE_ORIENTATION, SMLIB_FLAT_FACE_LOOP_COUNT,
                                SMLIB_FLAT_FACE_SURFACE_TYPE,  SMLIB_FLAT_LOOP_EDGEUSE_COUNT, SMLIB_FLAT_LOOP_VERTEX,
                                SMLIB_FLAT_EDGEUSE_EDGE,       SMLIB_FLAT_EDGEUSE_NEXT_RADIAL, SMLIB_FLAT_EDGE_CURVE_TYPE,
                                SMLIB_FLAT_EDGE_VERTICES,      SMLIB_FLAT_EDGE_NURBS_VERTEX_COUNT, SMLIB_FLAT_EDGE_NURBS_ORDER,
                                SMLIB_FLAT_FACE_NURBS_U_VERTEX_COUNT, SMLIB_FLAT_FACE_NURBS_V_ORDER, SMLIB_FLAT_EDGEUSE_NURBS_ORDER,
                                SMLIB_FLAT_VERTEX_POINT_TYPE};
    const int64_t badValues[] = {-1, -2, 7, 99, int64_t(1) << 40, std::numeric_limits<int64_t>::max(), std::numeric_limits<int64_t>::min()};
    for (int32_t id : indexIds)
    {
        if (base.i.at(id).empty())
        {
            continue;
        }
        for (int64_t v : badValues)
        {
            Arrays a = base;
            int64_t &x = a.i[id][a.i[id].size() / 2];
            if (x == v || (id == SMLIB_FLAT_LOOP_VERTEX && v == -1))
            {
                continue;
            }
            x = v;
            std::snprintf(what, sizeof what, "%s: int array %d set to %lld", name, id, (long long)v);
            ExpectReject(a, what);
            cases++;
        }
    }

    // Non-finite and absurd reals everywhere.
    const double badReals[] = {NAN, INFINITY, -INFINITY, 1e300};
    for (auto &kv : base.f)
    {
        if (kv.second.empty() || kv.first == SMLIB_FLAT_EXTENT)
        {
            continue;
        }
        for (double v : badReals)
        {
            Arrays a = base;
            a.f[kv.first][0] = v;
            std::snprintf(what, sizeof what, "%s: real array %d [0] set to %g", name, kv.first, v);
            ExpectReject(a, what);
            cases++;
        }
    }
    {
        // Tolerance zero, an empty set, and a mis-sized extent.
        Arrays a = base;
        a.f[SMLIB_FLAT_TOLERANCE][0] = 0.0;
        ExpectReject(a, "zero tolerance");
        Arrays e;
        ExpectReject(e, "empty");
        Arrays x = base;
        x.f[SMLIB_FLAT_EXTENT].push_back(1.0);
        ExpectReject(x, "extent of 7");
        // Export writes a tolerance per edge and per vertex; one too many is refused.
        CHECK(base.f.at(SMLIB_FLAT_EDGE_TOLERANCE).size() == base.i.at(SMLIB_FLAT_EDGE_CURVE_TYPE).size());
        CHECK(base.f.at(SMLIB_FLAT_VERTEX_TOLERANCE).size() == base.i.at(SMLIB_FLAT_VERTEX_POINT_TYPE).size());
        Arrays t = base;
        t.f[SMLIB_FLAT_EDGE_TOLERANCE].push_back(1e-5);
        ExpectReject(t, "an edge tolerance too many");
        Arrays z = base;
        if (!z.f[SMLIB_FLAT_VERTEX_TOLERANCE].empty())
        {
            z.f[SMLIB_FLAT_VERTEX_TOLERANCE][0] = 0.0;
            ExpectReject(z, "a zero vertex tolerance");
        }
        cases += 5;
    }
    return cases;
}

// Random single-value mutations: any status is fine, a crash is not.
void Mutations(const smlib_flat *flat, uint64_t seed, int rounds, int &accepted)
{
    const Arrays base = Snapshot(flat);
    std::vector<int32_t> ints, reals;
    for (auto &kv : base.i)
    {
        if (!kv.second.empty())
        {
            ints.push_back(kv.first);
        }
    }
    for (auto &kv : base.f)
    {
        if (!kv.second.empty())
        {
            reals.push_back(kv.first);
        }
    }
    uint64_t r = seed;
    auto next = [&r]
    {
        r ^= r << 13;
        r ^= r >> 7;
        r ^= r << 17;
        return r;
    };
    for (int k = 0; k < rounds; ++k)
    {
        Arrays a = base;
        if (next() % 2 == 0)
        {
            std::vector<int64_t> &v = a.i[ints[next() % ints.size()]];
            int64_t &x = v[next() % v.size()];
            switch (next() % 4)
            {
            case 0: x += 1; break;
            case 1: x -= 1; break;
            case 2: x = static_cast<int64_t>(next() % 64); break;
            default: x = static_cast<int64_t>(next()); break;
            }
        }
        else
        {
            std::vector<double> &v = a.f[reals[next() % reals.size()]];
            double &x = v[next() % v.size()];
            switch (next() % 4)
            {
            case 0: x = -x; break;
            case 1: x *= 1.0 + 1e-3; break;
            case 2: x += static_cast<double>(next() % 1000) / 100.0 - 5.0; break;
            default: x = 0.0; break;
            }
        }
        smlib_flat *m = FromArrays(a);
        smlib_brep *b = nullptr;
        if (smlib_brep_import(m, &b) == SMLIB_OK)
        {
            accepted++;
            double vol = 0.0;
            smlib_brep_volume(b, 1e-3, &vol);
        }
        smlib_brep_free(b);
        smlib_flat_free(m);
    }
}

} // namespace

int main(int argc, char **argv)
{
    if (argc > 1)
    {
        dir = argv[1];
    }
    CHECK(smlib_abi_version() == SMLIB_ABI_VERSION && SMLIB_ABI_VERSION == 3);
    CHECK(smlib_init() == SMLIB_OK);

    // The solids. With SMLIB_FLAT_SOLIDS set to a directory of <name>.smb files from an earlier
    // run, they are read from there instead of built (for sanitizer runs: SMLib's boolean and
    // fillet code has memory errors of its own that would stop the run before the converter).
    std::vector<std::pair<const char *, smlib_brep *>> solids;
    if (const char *from = std::getenv("SMLIB_FLAT_SOLIDS"))
    {
        for (const char *name : {"box", "cylinder", "cone", "sphere", "torus", "fillet", "blend", "blend1", "boolean", "incline"})
        {
            smlib_brep *b = nullptr;
            CHECK(smlib_brep_read((std::string(from) + "/" + name + ".smb").c_str(), 0, &b) == SMLIB_OK);
            solids.emplace_back(name, b);
        }
    }
    else
    {
    {
        SmBrep *b = nullptr;
        SmApiCreateBox(SmVector3d(1, 2, 3), 10, 20, 30, b);
        solids.emplace_back("box", ToCApi(b, "box"));
    }
    {
        SmBrep *b = nullptr;
        SmApiCreateCylinder(SmVector3d(0, 0, 0), 5, 12, b);
        solids.emplace_back("cylinder", ToCApi(b, "cylinder"));
    }
    {
        SmBrep *b = nullptr;
        SmApiCreateCone(SmVector3d(1, 1, 1), 6, 2, 9, b);
        solids.emplace_back("cone", ToCApi(b, "cone"));
    }
    {
        SmBrep *b = nullptr;
        SmApiCreateSphere(SmVector3d(0, 0, 0), 7, b);
        solids.emplace_back("sphere", ToCApi(b, "sphere"));
    }
    {
        SmBrep *b = nullptr;
        SmApiCreateTorus(SmVector3d(0, 0, 0), 10, 3, b);
        solids.emplace_back("torus", ToCApi(b, "torus"));
    }
    {
        SmBrep *b = nullptr;
        SmApiCreateBox(SmVector3d(0, 0, 0), 10, 10, 10, b);
        SmTArray<SmEdge *> edges;
        SmApiGetEdges(b, edges);
        SmApiStatus s = SmApiFilletEdges(b, edges, 1.0, 1, 1, 1.0);
        if (s != SM_SUCCESS)
        {
            std::fprintf(stderr, "fillet: SmApiFilletEdges -> %d\n", (int)s);
            failures++;
        }
        solids.emplace_back("fillet", ToCApi(b, "fillet"));
    }
    {
        // Blend-curve cross sections give NURBS fillet faces, so UV trims are carried too.
        SmBrep *b = nullptr;
        SmApiCreateBox(SmVector3d(0, 0, 0), 10, 10, 10, b);
        SmTArray<SmEdge *> edges;
        SmApiGetEdges(b, edges);
        SmApiStatus s = SmApiFilletEdges(b, edges, 1.0, 2, 2, 1.0);
        if (s != SM_SUCCESS)
        {
            std::fprintf(stderr, "blend: SmApiFilletEdges -> %d\n", (int)s);
            failures++;
        }
        solids.emplace_back("blend", ToCApi(b, "blend"));
    }
    {
        // One blended edge: one NURBS face with UV trims. (Its four boundary edges come out as
        // B-splines, so this still carries no trim on an analytic edge.)
        SmBrep *b = nullptr;
        SmApiCreateBox(SmVector3d(0, 0, 0), 10, 10, 10, b);
        SmTArray<SmEdge *> edges, one;
        SmApiGetEdges(b, edges);
        one.Add(edges[0]);
        SmApiStatus s = SmApiFilletEdges(b, one, 2.0, 2, 2, 1.0);
        if (s != SM_SUCCESS)
        {
            std::fprintf(stderr, "blend1: SmApiFilletEdges -> %d\n", (int)s);
            failures++;
        }
        solids.emplace_back("blend1", ToCApi(b, "blend1"));
    }
    {
        SmBrep *box = nullptr, *cyl = nullptr, *r = nullptr;
        SmApiCreateBox(SmVector3d(0, 0, 0), 20, 20, 20, box);
        SmApiCreateCylinder(SmVector3d(10, 10, -5), 4, 30, cyl);
        SmApiStatus s = SmApiBooleanDifference(box, cyl, r);
        if (s != SM_SUCCESS)
        {
            std::fprintf(stderr, "boolean: SmApiBooleanDifference -> %d\n", (int)s);
            failures++;
            delete cyl;
            delete box;
            r = nullptr;
        }
        solids.emplace_back("boolean", ToCApi(r, "boolean"));
    }
    {
        // A cylinder minus a tilted box: the inclined rim is made of ellipse arcs. (SmApiCut with
        // an inclined plane also makes an ellipse, but leaves the body open: it does not cap.)
        SmBrep *cyl = nullptr, *box = nullptr, *r = nullptr;
        SmApiCreateCylinder(SmVector3d(0, 0, 0), 5, 20, cyl);
        SmApiCreateBox(SmVector3d(-20, -20, 10), 40, 40, 40, box);
        SmApiRotate(box, SmVector3d(0, 0, 10), SmVector3d(0, 1, 0), 25.0);
        SmApiStatus s = SmApiBooleanDifference(cyl, box, r);
        if (s != SM_SUCCESS)
        {
            std::fprintf(stderr, "incline: SmApiBooleanDifference -> %d\n", (int)s);
            failures++;
            delete cyl;
            delete box;
            r = nullptr;
        }
        solids.emplace_back("incline", ToCApi(r, "incline"));
    }
    }

    std::vector<std::pair<const char *, smlib_flat *>> exports;
    for (auto &s : solids)
    {
        const bool cone = std::strcmp(s.first, "cone") == 0;
        const bool reparameterised = cone || std::strcmp(s.first, "sphere") == 0 || std::strcmp(s.first, "incline") == 0;
        smlib_flat *f = RoundTrip(s.first, s.second, cone, reparameterised ? 1e-7 : 1e-9);
        if (f != nullptr)
        {
            exports.emplace_back(s.first, f);
        }
    }

    // An analytic face's angular range may start a turn past the principal one (solidcore#428):
    // the torus with its major, minor or both ranges moved on by 2 pi is the same torus. Its
    // tube circle was once rebuilt over a clamped, empty parameter interval, which gave NaN
    // knots, an overflow in N_IgesClampSrf (#427) and then a refused import.
    for (auto &e : exports)
    {
        if (std::strcmp(e.first, "torus") != 0)
        {
            continue;
        }
        const double pi = 3.14159265358979323846;
        for (int shift = 1; shift <= 3; ++shift)
        {
            Arrays a = Snapshot(e.second);
            const std::vector<int64_t> &type = a.i[SMLIB_FLAT_FACE_SURFACE_TYPE];
            std::vector<double> &range = a.f[SMLIB_FLAT_FACE_RANGE];
            int moved = 0;
            for (size_t k = 0; k < type.size() && 4 * k + 3 < range.size(); ++k)
            {
                if (type[k] != SMLIB_FLAT_SURFACE_TORUS)
                {
                    continue;
                }
                moved++;
                for (int c = 0; c < 2; ++c)
                {
                    if (shift & (1 << c))
                    {
                        range[4 * k + c] += 2.0 * pi;
                        range[4 * k + 2 + c] += 2.0 * pi;
                    }
                }
            }
            CHECK(moved == 1);
            smlib_flat *f = FromArrays(a);
            smlib_brep *b = nullptr;
            int32_t s = smlib_brep_import(f, &b);
            const double want = 2.0 * pi * pi * 10.0 * 9.0;
            double v = 0.0;
            int64_t len = 0, bad = -1;
            char buf[4096] = {0};
            if (s == SMLIB_OK)
            {
                v = Measure(b).volume;
                CHECK(smlib_brep_check(b, 2, buf, sizeof buf, &len, &bad) == SMLIB_OK);
            }
            std::printf("torus, %s range a turn on: import %d, volume %.12g (want %.12g), failed checks %lld %s\n",
                        shift == 1 ? "major" : shift == 2 ? "minor" : "both", static_cast<int>(s), v, want, static_cast<long long>(bad), buf);
            CHECK(s == SMLIB_OK && std::fabs(v - want) <= 1e-7 * want && bad == 0);
            smlib_brep_free(b);
            smlib_flat_free(f);
        }
    }

    NestedShellOrder();

    // The accessors' argument checks.
    {
        smlib_flat *f = smlib_flat_new();
        const int64_t *pi = nullptr;
        const double *pf = nullptr;
        int64_t n = 0;
        const int64_t one = 1;
        const double half = 0.5;
        for (int32_t id : {0, 33, 100, 166, -1, 1000})
        {
            CHECK(smlib_flat_get_i64(f, id, &pi, &n) == SMLIB_ERR_ARG);
            CHECK(smlib_flat_get_f64(f, id, &pf, &n) == SMLIB_ERR_ARG);
            CHECK(smlib_flat_set_i64(f, id, &one, 1) == SMLIB_ERR_ARG);
            CHECK(smlib_flat_set_f64(f, id, &half, 1) == SMLIB_ERR_ARG);
        }
        CHECK(smlib_flat_get_f64(f, SMLIB_FLAT_REGION_TYPE, &pf, &n) == SMLIB_ERR_ARG);
        CHECK(smlib_flat_get_i64(f, SMLIB_FLAT_TOLERANCE, &pi, &n) == SMLIB_ERR_ARG);
        CHECK(smlib_flat_set_f64(f, SMLIB_FLAT_REGION_TYPE, &half, 1) == SMLIB_ERR_ARG);
        CHECK(smlib_flat_set_i64(f, SMLIB_FLAT_TOLERANCE, &one, 1) == SMLIB_ERR_ARG);
        CHECK(smlib_flat_set_i64(f, SMLIB_FLAT_REGION_TYPE, &one, -1) == SMLIB_ERR_ARG);
        CHECK(smlib_flat_set_i64(f, SMLIB_FLAT_REGION_TYPE, &one, SMLIB_FLAT_MAX_LEN + 1) == SMLIB_ERR_ARG);
        CHECK(smlib_flat_set_i64(f, SMLIB_FLAT_REGION_TYPE, nullptr, 3) == SMLIB_ERR_ARG);
        CHECK(smlib_flat_set_i64(f, SMLIB_FLAT_REGION_TYPE, nullptr, 0) == SMLIB_OK);
        CHECK(smlib_flat_set_i64(f, SMLIB_FLAT_REGION_TYPE, &one, 1) == SMLIB_OK);
        CHECK(smlib_flat_get_i64(f, SMLIB_FLAT_REGION_TYPE, &pi, &n) == SMLIB_OK && n == 1 && pi[0] == 1);
        CHECK(smlib_flat_get_i64(nullptr, SMLIB_FLAT_REGION_TYPE, &pi, &n) == SMLIB_ERR_ARG);
        smlib_brep *b = nullptr;
        CHECK(smlib_brep_import(f, &b) == SM_ERR_INVALID_INPUT && b == nullptr);
        CHECK(smlib_brep_import(nullptr, &b) == SMLIB_ERR_ARG);
        CHECK(smlib_brep_export(nullptr, &f) == SMLIB_ERR_ARG);
        CHECK(smlib_flat_read((dir + "/no-such-file").c_str(), &f) == SM_ERR_INVALID_INPUT);
        smlib_flat_free(f);
        smlib_flat_free(nullptr);
    }

    // A failing import dumps its arrays; the dump reads back as the same arrays.
    if (!exports.empty())
    {
        std::string dumps = dir + "/dumps";
        std::string mk = "mkdir -p '" + dumps + "' && rm -f '" + dumps + "'/*";
        CHECK(std::system(mk.c_str()) == 0);
        CHECK(smlib_set_dump_dir(dumps.c_str()) == SMLIB_OK);
        Arrays bad = Snapshot(exports[0].second);
        bad.i[SMLIB_FLAT_EDGEUSE_EDGE][0] = 1000;
        smlib_flat *f = FromArrays(bad);
        smlib_brep *b = nullptr;
        CHECK(smlib_brep_import(f, &b) == SM_ERR_INVALID_INPUT);
        CHECK(smlib_set_dump_dir(nullptr) == SMLIB_OK);
        smlib_flat *replay = nullptr;
        bool found = false;
        for (int n = 1; n < 100000 && !found; ++n)
        {
            std::string path = dumps + "/" + std::to_string(n) + "-smlib_brep_import-flat.txt";
            if (smlib_flat_read(path.c_str(), &replay) == SMLIB_OK)
            {
                found = true;
            }
        }
        CHECK(found);
        if (found)
        {
            Arrays r = Snapshot(replay);
            CHECK(r.i == bad.i && r.f == bad.f);
            CHECK(smlib_brep_import(replay, &b) == SM_ERR_INVALID_INPUT);
        }
        smlib_flat_free(replay);
        smlib_flat_free(f);
    }

    int hostile = 0, accepted = 0, mutated = 0;
    uint64_t seed = 0x9E3779B97F4A7C15ull;
    for (auto &e : exports)
    {
        hostile += Hostile(e.first, e.second);
        Mutations(e.second, seed++, 300, accepted);
        mutated += 300;
    }
    std::printf("hostile cases: %d rejected as expected; random mutations: %d, of which %d still imported\n", hostile, mutated, accepted);

    for (auto &e : exports)
    {
        smlib_flat_free(e.second);
    }
    if (failures == 0)
    {
        std::printf("smlib flat round trip: ok\n");
    }
    return failures == 0 ? 0 : 1;
}
