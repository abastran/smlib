// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// smlib_history_tag / smlib_history_read on booleans of primitives whose history is known by
// hand. The booleans themselves run through SM_API directly (as the operation shims do), on the
// tagged copies.

#include "smlib_shim.h"

#include <SmApiBrep.h>
#include <SmApiPrimitives.h>
#include <SmApiQueries.h>
#include <SmFace.h>
#include <SmVector3d.h>

#include <cmath>
#include <cstdio>
#include <map>
#include <set>
#include <vector>

namespace
{

int g_failed = 0;

#define CHECK(cond)                                                                                                    \
    do                                                                                                                 \
    {                                                                                                                  \
        if (!(cond))                                                                                                   \
        {                                                                                                              \
            std::fprintf(stderr, "%s:%d: CHECK failed: %s\n", __FILE__, __LINE__, #cond);                            \
            ++g_failed;                                                                                                \
        }                                                                                                              \
    } while (0)

using Origin = std::pair<int64_t, int64_t>;

struct Face
{
    int64_t change;
    std::vector<Origin> origins;
};

struct Edge
{
    std::vector<int64_t> faces;
    std::vector<Origin> origins;
};

struct History
{
    std::vector<Face> faces;
    std::vector<Edge> edges;

    size_t count(int64_t change) const
    {
        size_t n = 0;
        for (const Face &f : faces)
        {
            n += f.change == change;
        }
        return n;
    }

    size_t meets() const
    {
        size_t n = 0;
        for (const Edge &e : edges)
        {
            n += e.origins.empty();
        }
        return n;
    }
};

// Reads the history through the C API and decodes it, checking the layout as it goes.
History read(const smlib_brep *r)
{
    History h;
    int64_t len = 0;
    CHECK(smlib_history_read(r, nullptr, 0, &len) == SMLIB_ERR_BUFFER);
    std::vector<int64_t> v(static_cast<size_t>(len));
    int64_t again = 0;
    CHECK(smlib_history_read(r, v.data(), len, &again) == SMLIB_OK && again == len);
    size_t k = 0;
    auto next = [&]() -> int64_t
    {
        if (k >= v.size())
        {
            ++g_failed;
            std::fprintf(stderr, "history ran out at %zu\n", k);
            return 0;
        }
        return v[k++];
    };
    CHECK(next() == SMLIB_HISTORY_VERSION);
    int64_t nf = next(), ne = next();
    for (int64_t i = 0; i < nf; ++i)
    {
        Face f;
        f.change = next();
        for (int64_t n = next(); n > 0; --n)
        {
            int64_t op = next();
            f.origins.emplace_back(op, next());
        }
        h.faces.push_back(f);
    }
    for (int64_t i = 0; i < ne; ++i)
    {
        Edge e;
        for (int64_t m = next(); m > 0; --m)
        {
            e.faces.push_back(next());
        }
        for (int64_t n = next(); n > 0; --n)
        {
            int64_t op = next();
            e.origins.emplace_back(op, next());
        }
        h.edges.push_back(e);
    }
    CHECK(k == v.size());
    int64_t counts[3];
    CHECK(smlib_brep_counts(r, counts) == SMLIB_OK && counts[0] == nf && counts[1] == ne);
    return h;
}

smlib_brep *box(double x, double y, double z, double dx, double dy, double dz)
{
    SmBrep *p = nullptr;
    CHECK(SmApiCreateBox(SmVector3d(x, y, z), dx, dy, dz, p) == SM_SUCCESS && p != nullptr);
    return new smlib_brep{p};
}

smlib_brep *cylinder(double x, double y, double z, double r, double h)
{
    SmBrep *p = nullptr;
    CHECK(SmApiCreateCylinder(SmVector3d(x, y, z), r, h, p) == SM_SUCCESS && p != nullptr);
    return new smlib_brep{p};
}

smlib_brep *tagged(smlib_brep *b, int32_t operand)
{
    smlib_brep *t = nullptr;
    CHECK(smlib_history_tag(b, operand, &t) == SMLIB_OK && t != nullptr);
    smlib_brep_free(b);
    return t;
}

enum Op
{
    UNION,
    DIFFERENCE,
    INTERSECTION
};

// a op b on tagged copies (operands 0 and 1); the result handle owns the result.
smlib_brep *boolean(smlib_brep *a, smlib_brep *b, Op op)
{
    smlib_brep *ta = tagged(a, 0), *tb = tagged(b, 1);
    SmBrep *r = nullptr;
    SmApiStatus s = op == UNION        ? SmApiBooleanUnion(ta->p, tb->p, r)
                    : op == DIFFERENCE ? SmApiBooleanDifference(ta->p, tb->p, r)
                                       : SmApiBooleanIntersection(ta->p, tb->p, r);
    CHECK(s == SM_SUCCESS && r == ta->p);
    // SMLib kept the first operand as the result and deleted the second.
    tb->p = nullptr;
    smlib_brep_free(tb);
    ta->p = r;
    return ta;
}

double area(const SmBrep *b, int64_t face)
{
    SmTArray<SmFace *> faces;
    b->GetFaces(faces);
    double a = 0.0;
    CHECK(face >= 0 && face < static_cast<int64_t>(faces.GetSize()));
    CHECK(SmApiFaceComputeArea(faces[static_cast<ULONG>(face)], 1e-6, a) == SM_SUCCESS);
    return a;
}

// Every face's origins are in range and each Kept face has its origin's area; every edge bounds
// two faces (the results here are manifold solids).
void consistent(const History &h, const smlib_brep *r, const std::vector<const SmBrep *> &operands)
{
    for (size_t i = 0; i < h.faces.size(); ++i)
    {
        const Face &f = h.faces[i];
        CHECK((f.change == SMLIB_CHANGE_GENERATED) == f.origins.empty());
        CHECK((f.change == SMLIB_CHANGE_MERGED) == (f.origins.size() > 1));
        for (const Origin &o : f.origins)
        {
            CHECK(o.first >= 0 && o.first < static_cast<int64_t>(operands.size()));
            if (o.first >= 0 && o.first < static_cast<int64_t>(operands.size()))
            {
                CHECK(o.second >= 0 && o.second < static_cast<int64_t>(operands[o.first]->GetNumFaces()));
            }
        }
        if (f.change == SMLIB_CHANGE_KEPT)
        {
            double a = area(r->p, static_cast<int64_t>(i)), b = area(operands[f.origins[0].first], f.origins[0].second);
            CHECK(std::fabs(a - b) <= 1e-6 * b);
        }
        if (f.change == SMLIB_CHANGE_MODIFIED || f.change == SMLIB_CHANGE_SPLIT)
        {
            CHECK(area(r->p, static_cast<int64_t>(i)) < area(operands[f.origins[0].first], f.origins[0].second));
        }
    }
    for (const Edge &e : h.edges)
    {
        // Two faces, or one for a seam (the hole's wall meets itself there).
        CHECK(e.faces.size() == 2 || (e.faces.size() == 1 && e.origins.size() == 1));
        for (int64_t f : e.faces)
        {
            CHECK(f >= 0 && f < static_cast<int64_t>(h.faces.size()));
        }
    }
}

// Copies of the operands as tagged, kept for checking against.
std::vector<const SmBrep *> keep(std::initializer_list<smlib_brep *> bs)
{
    std::vector<const SmBrep *> v;
    for (smlib_brep *b : bs)
    {
        SmBrep *c = nullptr;
        CHECK(SmApiBrepCopy(b->p, c) == SM_SUCCESS);
        v.push_back(c);
    }
    return v;
}

void drop(std::vector<const SmBrep *> &v)
{
    for (const SmBrep *b : v)
    {
        delete b;
    }
    v.clear();
}

void overlapping_boxes_union()
{
    smlib_brep *a = box(0, 0, 0, 10, 10, 10), *b = box(5, 5, 5, 10, 10, 10);
    auto ops = keep({a, b});
    smlib_brep *r = boolean(a, b, UNION);
    History h = read(r);
    consistent(h, r, ops);
    // Each box keeps its three faces away from the other whole and its other three L-shaped.
    CHECK(h.faces.size() == 12);
    CHECK(h.count(SMLIB_CHANGE_KEPT) == 6);
    CHECK(h.count(SMLIB_CHANGE_MODIFIED) == 6);
    // The six edges where the two boxes' faces cross are new; every other edge is (part of) one
    // of theirs.
    CHECK(h.meets() == 6);
    for (const Edge &e : h.edges)
    {
        if (e.origins.empty())
        {
            const Face &x = h.faces[static_cast<size_t>(e.faces[0])], &y = h.faces[static_cast<size_t>(e.faces[1])];
            CHECK(x.change == SMLIB_CHANGE_MODIFIED && y.change == SMLIB_CHANGE_MODIFIED);
            CHECK(x.origins[0].first != y.origins[0].first);
        }
        CHECK(e.origins.size() <= 1);
    }
    smlib_brep_free(r);
    drop(ops);
}

void side_by_side_boxes_union()
{
    // Overlapping along x with the same y and z extents: one 15 x 10 x 10 box.
    smlib_brep *a = box(0, 0, 0, 10, 10, 10), *b = box(5, 0, 0, 10, 10, 10);
    auto ops = keep({a, b});
    smlib_brep *r = boolean(a, b, UNION);
    History h = read(r);
    consistent(h, r, ops);
    CHECK(h.faces.size() == 6);
    // The four long faces each join a face of each box; the two ends are each box's own.
    CHECK(h.count(SMLIB_CHANGE_MERGED) == 4);
    CHECK(h.count(SMLIB_CHANGE_KEPT) == 2);
    for (const Face &f : h.faces)
    {
        if (f.change == SMLIB_CHANGE_MERGED)
        {
            CHECK(f.origins.size() == 2 && f.origins[0].first == 0 && f.origins[1].first == 1);
        }
    }
    // No new edges: the long edges join one edge of each box, the end edges are each box's.
    CHECK(h.meets() == 0);
    size_t joined = 0;
    for (const Edge &e : h.edges)
    {
        joined += e.origins.size() == 2;
    }
    CHECK(joined == 4);
    smlib_brep_free(r);
    drop(ops);
}

void hole_through_box()
{
    smlib_brep *a = box(0, 0, 0, 10, 10, 10), *b = cylinder(5, 5, -5, 2, 20);
    auto ops = keep({a, b});
    smlib_brep *r = boolean(a, b, DIFFERENCE);
    History h = read(r);
    consistent(h, r, ops);
    // Top and bottom get a hole (Modified), the four sides are untouched, and the cylinder's
    // side is trimmed to the hole's wall.
    CHECK(h.faces.size() == 7);
    CHECK(h.count(SMLIB_CHANGE_KEPT) == 4);
    CHECK(h.count(SMLIB_CHANGE_MODIFIED) == 3);
    size_t from_tool = 0;
    for (const Face &f : h.faces)
    {
        from_tool += f.origins[0].first == 1;
    }
    CHECK(from_tool == 1);
    // The two circles where the wall meets top and bottom are new; the seam is part of the
    // cylinder's.
    CHECK(h.meets() == 2);
    smlib_brep_free(r);
    drop(ops);
}

void pocket()
{
    // A 4 x 4 pocket 7 deep from the top of a 10 mm cube: the tool's floor face lies inside the
    // cube, whole.
    smlib_brep *a = box(0, 0, 0, 10, 10, 10), *b = box(3, 3, 3, 4, 4, 20);
    auto ops = keep({a, b});
    smlib_brep *r = boolean(a, b, DIFFERENCE);
    History h = read(r);
    consistent(h, r, ops);
    CHECK(h.faces.size() == 11);
    // Five cube faces and the floor are whole; the cube's top and the four walls are trimmed.
    CHECK(h.count(SMLIB_CHANGE_KEPT) == 6);
    CHECK(h.count(SMLIB_CHANGE_MODIFIED) == 5);
    CHECK(h.meets() == 4);
    smlib_brep_free(r);
    drop(ops);
}

void slot_splits_box()
{
    // A slot right through the middle: two pieces, so the four long faces are each split in two.
    smlib_brep *a = box(0, 0, 0, 30, 10, 10), *b = box(10, -1, -1, 10, 12, 12);
    auto ops = keep({a, b});
    smlib_brep *r = boolean(a, b, DIFFERENCE);
    History h = read(r);
    consistent(h, r, ops);
    CHECK(h.faces.size() == 12);
    CHECK(h.count(SMLIB_CHANGE_SPLIT) == 8);
    CHECK(h.count(SMLIB_CHANGE_KEPT) == 2);
    CHECK(h.count(SMLIB_CHANGE_MODIFIED) == 2);
    std::map<Origin, int> uses;
    for (const Face &f : h.faces)
    {
        if (f.change == SMLIB_CHANGE_SPLIT)
        {
            ++uses[f.origins[0]];
        }
    }
    CHECK(uses.size() == 4);
    for (auto &u : uses)
    {
        CHECK(u.first.first == 0 && u.second == 2);
    }
    smlib_brep_free(r);
    drop(ops);
}

void apart_boxes_union()
{
    smlib_brep *a = box(0, 0, 0, 10, 10, 10), *b = box(20, 0, 0, 10, 10, 10);
    auto ops = keep({a, b});
    smlib_brep *r = boolean(a, b, UNION);
    History h = read(r);
    consistent(h, r, ops);
    CHECK(h.faces.size() == 12 && h.count(SMLIB_CHANGE_KEPT) == 12);
    CHECK(h.edges.size() == 24 && h.meets() == 0);
    smlib_brep_free(r);
    drop(ops);
}

void retag_and_hostile_arguments()
{
    smlib_brep *a = box(0, 0, 0, 1, 1, 1), *t = nullptr, *tt = nullptr;
    CHECK(smlib_history_tag(a, 3, &t) == SMLIB_OK);
    // Tagging a tagged body replaces its tags.
    CHECK(smlib_history_tag(t, 5, &tt) == SMLIB_OK);
    History h = read(tt);
    for (const Face &f : h.faces)
    {
        CHECK(f.change == SMLIB_CHANGE_KEPT && f.origins.size() == 1 && f.origins[0].first == 5);
    }
    // An untagged body: every face new, every edge where its faces meet.
    History u = read(a);
    CHECK(u.count(SMLIB_CHANGE_GENERATED) == 6 && u.meets() == 12);
    int64_t len = 0, small[4];
    CHECK(smlib_history_tag(nullptr, 0, &t) == SMLIB_ERR_ARG);
    CHECK(smlib_history_tag(a, -1, &t) == SMLIB_ERR_ARG);
    CHECK(smlib_history_tag(a, 1 << 30, &t) == SMLIB_ERR_ARG);
    CHECK(smlib_history_tag(a, 0, nullptr) == SMLIB_ERR_ARG);
    CHECK(smlib_history_read(nullptr, small, 4, &len) == SMLIB_ERR_ARG);
    CHECK(smlib_history_read(a, nullptr, 4, &len) == SMLIB_ERR_ARG);
    CHECK(smlib_history_read(a, small, -1, &len) == SMLIB_ERR_ARG);
    CHECK(smlib_history_read(a, small, 4, nullptr) == SMLIB_ERR_ARG);
    CHECK(smlib_history_read(a, small, 4, &len) == SMLIB_ERR_BUFFER && len > 4);
    smlib_brep_free(a);
    smlib_brep_free(t);
    smlib_brep_free(tt);
}

} // namespace

int main()
{
    CHECK(smlib_init() == SMLIB_OK);
    overlapping_boxes_union();
    side_by_side_boxes_union();
    hole_through_box();
    pocket();
    slot_splits_box();
    apart_boxes_union();
    retag_and_hostile_arguments();
    if (g_failed != 0)
    {
        std::fprintf(stderr, "%d checks failed\n", g_failed);
        return 1;
    }
    std::printf("history: all checks passed\n");
    return 0;
}
