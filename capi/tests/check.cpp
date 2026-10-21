// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// smlib_brep_check: SMLib's primitives pass every test at every level; a vertex moved off its
// edges fails, with a report; the buffer protocol and the argument checks hold.

#include "../src/smlib_shim.h"

#include <SmApiPrimitives.h>
#include <SmVertex.h>

#include <cstdio>
#include <cstring>
#include <string>

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
