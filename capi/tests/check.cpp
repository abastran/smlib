// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// smlib_brep_check: SMLib's primitives pass every test at every level; a vertex moved off its
// edges fails, with a report; the buffer protocol and the argument checks hold; the check of a
// body of many faces takes time about linear in them (tangency-labs/solidcore#431).

#include "../src/smlib_shim.h"

#include <SmApiPrimitives.h>
#include <SmEdge.h>
#include <SmVertex.h>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
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

// The test makes its own bodies with SM_API, as the shims do, and wraps them in handles.
template <typename F>
static smlib_brep *make(F &&f)
{
    SmBrep *p = nullptr;
    CHECK(f(p) == SM_SUCCESS && p != nullptr);
    return new smlib_brep{p};
}

struct Report
{
    int32_t status;
    int64_t failures;
    std::string text;
};

static Report check(const smlib_brep *b, int32_t level)
{
    static char buf[1 << 16];
    int64_t len = -1, n = -1;
    int32_t s = smlib_brep_check(b, level, buf, sizeof buf, &len, &n);
    return {s, n, s == SMLIB_OK ? std::string(buf, static_cast<size_t>(len)) : std::string()};
}

// A plate 3 mm thick with `teeth` teeth along one side, extruded from one profile of lines:
// 4 * teeth + 5 planar faces, the top and bottom ones bounded by 4 * teeth + 3 edges each.
static smlib_brep *comb(int teeth)
{
    std::vector<int64_t> kinds;
    std::vector<double> reals;
    auto line = [&](double x0, double y0, double x1, double y1)
    {
        kinds.push_back(SMLIB_CURVE_LINE);
        const double r[6] = {x0, y0, 0.0, x1, y1, 0.0};
        reals.insert(reals.end(), r, r + 6);
    };
    const double x = 2.0 * teeth;
    line(0, 0, x, 0);
    line(x, 0, x, 5);
    for (int i = teeth - 1; i >= 0; --i)
    {
        const double x0 = 2.0 * i;
        line(x0 + 2, 5, x0 + 1, 5);
        line(x0 + 1, 5, x0 + 1, 7);
        line(x0 + 1, 7, x0, 7);
        line(x0, 7, x0, 5);
    }
    line(0, 5, 0, 0);
    const double up[3] = {0, 0, 1};
    smlib_brep *b = nullptr;
    CHECK(smlib_extrude(kinds.data(), static_cast<int64_t>(kinds.size()), reals.data(), static_cast<int64_t>(reals.size()), up, 3.0,
                        0.0, 1, &b) == SMLIB_OK);
    return b;
}

// Seconds smlib_brep_check takes on `b` at level 2, which must find nothing wrong with it.
static double seconds_to_check(const smlib_brep *b)
{
    const auto t0 = std::chrono::steady_clock::now();
    Report r = check(b, 2);
    const double t = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    CHECK(r.status == SMLIB_OK && r.failures == 0);
    return t;
}

int main()
{
    CHECK(smlib_init() == SMLIB_OK);
    const SmVector3d o(1.0, 2.0, 3.0);
    auto box = [&](SmBrep *&p) { return SmApiCreateBox(o, 10.0, 20.0, 30.0, p); };
    smlib_brep *bodies[] = {
        make(box),
        make([&](SmBrep *&p) { return SmApiCreateCylinder(o, 15.0, 50.0, p); }),
        make([&](SmBrep *&p) { return SmApiCreateSphere(o, 12.0, p); }),
        make([&](SmBrep *&p) { return SmApiCreateCone(o, 10.0, 4.0, 20.0, p); }),
        make([&](SmBrep *&p) { return SmApiCreateCone(o, 10.0, 0.0, 20.0, p); }),
        make([&](SmBrep *&p) { return SmApiCreateTorus(o, 22.0, 4.0, p); }),
    };
    const char *names[] = {"box", "cylinder", "sphere", "frustum", "cone", "torus"};
    for (size_t i = 0; i < sizeof bodies / sizeof bodies[0]; ++i)
    {
        for (int32_t level = 0; level <= 2; ++level)
        {
            Report r = check(bodies[i], level);
            if (r.status != SMLIB_OK || r.failures != 0 || !r.text.empty())
            {
                std::fprintf(stderr, "%s, level %d: status %d, %lld failures\n%s", names[i], level, r.status,
                             static_cast<long long>(r.failures), r.text.c_str());
            }
            CHECK(r.status == SMLIB_OK && r.failures == 0 && r.text.empty());
        }
    }

    // A box with one corner pulled 1 mm off its three edges and faces.
    smlib_brep *broken = make(box);
    SmTArray<SmVertex *> vertices;
    broken->p->GetVertices(vertices);
    CHECK(vertices.GetSize() == 8);
    SmPoint3d moved = vertices[0]->GetPoint();
    moved.x += 1.0;
    vertices[0]->SetPoint(moved);
    Report bad = check(broken, 2);
    CHECK(bad.status == SMLIB_OK && bad.failures > 0);
    CHECK(!bad.text.empty() && bad.text.back() == '\n');
    std::printf("moved vertex, level 2: %lld failed tests\n%s", static_cast<long long>(bad.failures), bad.text.c_str());

    // Too small a buffer: SMLIB_ERR_BUFFER with the needed length, and nothing else written.
    int64_t len = -1, n = -7;
    CHECK(smlib_brep_check(broken, 2, nullptr, 0, &len, &n) == SMLIB_ERR_BUFFER);
    CHECK(len == static_cast<int64_t>(bad.text.size()) && n == -7);
    char small[4] = {'x', 'x', 'x', 'x'};
    CHECK(smlib_brep_check(broken, 2, small, sizeof small, &len, &n) == SMLIB_ERR_BUFFER && small[0] == 'x');
    std::string exact(static_cast<size_t>(len) + 1, '\0');
    CHECK(smlib_brep_check(broken, 2, &exact[0], len + 1, &len, &n) == SMLIB_OK && n == bad.failures);
    CHECK(std::strcmp(exact.c_str(), bad.text.c_str()) == 0);

    // tangency-labs/solidcore#431: the check took time quadratic in the faces of a body (more
    // where one face has many edges): each face's vertices listed again for each point
    // classified on it, each list of topology walked for each of its members, every pair of
    // vertices, edges and faces compared. Eight times the teeth, eight times the faces, edges and
    // vertices, must take nowhere near 64 times as long (it took about 40 times as long before,
    // 6 to 10 times now).
    {
        smlib_brep *small = comb(60), *large = comb(480);
        CHECK(small != nullptr && large != nullptr);
        if (small != nullptr && large != nullptr)
        {
            int64_t counts[3] = {0, 0, 0};
            CHECK(smlib_brep_counts(large, counts) == SMLIB_OK && counts[0] == 4 * 480 + 5);
            // The best of three, so that a moment of load on the machine does not count.
            double ts = 1e30, tl = 1e30;
            for (int i = 0; i < 3; ++i)
            {
                ts = std::fmin(ts, seconds_to_check(small));
                tl = std::fmin(tl, seconds_to_check(large));
            }
            std::printf("check of a comb of 60 teeth: %.3f s, of 480 teeth: %.3f s (x%.1f)\n", ts, tl, tl / ts);
            CHECK(tl < 24.0 * ts);
        }
        // Coincident topology is still found on such a body: a corner moved onto another, and
        // a corner moved onto the middle of an edge it does not bound.
        if (small != nullptr)
        {
            SmTArray<SmTopology *> first, second;
            // (With every pair asked for, it returns TRUE whatever it finds: the pairs tell.)
            small->p->ContainsCoincidentTopology(-1.0, &first, &second, nullptr, SM_LEVEL_2, nullptr, TRUE);
            CHECK(first.GetSize() == 0 && second.GetSize() == 0);
            SmTArray<SmVertex *> vs;
            small->p->GetVertices(vs);
            CHECK(vs.GetSize() > 100);
            const SmPoint3d p0 = vs[0]->GetPoint();
            vs[0]->SetPoint(vs[100]->GetPoint());
            small->p->ContainsCoincidentTopology(-1.0, &first, &second, nullptr, SM_LEVEL_2, nullptr, TRUE);
            bool pair = false;
            for (ULONG i = 0; i < first.GetSize(); ++i)
            {
                pair |= (first[i] == vs[0] && second[i] == vs[100]) || (first[i] == vs[100] && second[i] == vs[0]);
            }
            CHECK(pair);
            vs[0]->SetPoint(p0);
            SmTArray<SmEdge *> es;
            small->p->GetEdges(es);
            SmEdge *e = es[es.GetSize() / 2];
            SmVertex *v = vs[vs.GetSize() / 3];
            CHECK(e->GetVertex() != v && e->GetOtherVertex(e->GetVertex()) != v);
            SmPoint3d mid;
            CHECK(e->GetCurve()->EvaluatePoint(e->GetInterval().Evaluate(0.5), mid) == SM_SUCCESS);
            v->SetPoint(mid);
            small->p->ContainsCoincidentTopology(-1.0, &first, &second, nullptr, SM_LEVEL_2, nullptr, TRUE);
            pair = false;
            for (ULONG i = 0; i < first.GetSize(); ++i)
            {
                pair |= first[i] == v && second[i] == e;
            }
            CHECK(pair);
        }
        smlib_brep_free(small);
        smlib_brep_free(large);
    }

    // Rejected arguments.
    char buf[16];
    CHECK(smlib_brep_check(nullptr, 2, buf, sizeof buf, &len, &n) == SMLIB_ERR_ARG);
    CHECK(smlib_brep_check(bodies[0], 3, buf, sizeof buf, &len, &n) == SMLIB_ERR_ARG);
    CHECK(smlib_brep_check(bodies[0], -1, buf, sizeof buf, &len, &n) == SMLIB_ERR_ARG);
    CHECK(smlib_brep_check(bodies[0], 2, nullptr, 8, &len, &n) == SMLIB_ERR_ARG);
    CHECK(smlib_brep_check(bodies[0], 2, buf, -1, &len, &n) == SMLIB_ERR_ARG);
    CHECK(smlib_brep_check(bodies[0], 2, buf, sizeof buf, nullptr, &n) == SMLIB_ERR_ARG);
    CHECK(smlib_brep_check(bodies[0], 2, buf, sizeof buf, &len, nullptr) == SMLIB_ERR_ARG);

    for (smlib_brep *b : bodies)
    {
        smlib_brep_free(b);
    }
    smlib_brep_free(broken);
    if (failures == 0)
    {
        std::printf("smlib_brep_check test: ok\n");
    }
    return failures == 0 ? 0 : 1;
}
