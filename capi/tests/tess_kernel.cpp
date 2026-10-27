// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// SMLib's tessellator (SmApiTessellate) on SM_API-made solids: every mesh is closed and
// consistently wound, every vertex is on the surface and every triangle lies within the chord
// tolerance of it, measured at its centroid, at every tolerance (tangency-labs/solidcore#390:
// slivers along a sphere's seam by the poles once reached 13 times it, and the triangles of
// spheres and tori 1.2 times).

#include <SmApiBrep.h>
#include <SmApiPrimitives.h>
#include <SmApiTessellationParams.h>
#include <SmBrep.h>
#include <SmPoly.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <functional>
#include <map>
#include <utility>

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

using Distance = std::function<double(const SmPoint3d &)>;

// Tessellates `b` (and deletes it) at `chord` and 20 degrees, checks the mesh is closed and
// consistently wound, and bounds the distance of every vertex and every triangle centroid from
// the body's surface, which `dist` gives for a point near it.
void Check(const char *name, SmBrep *b, double chord, const Distance &dist)
{
    CHECK(b != nullptr);
    if (b == nullptr)
    {
        return;
    }
    SmTessellationParams params;
    params.dChordHeightTol = chord;
    params.dCurveAngleTolDeg = 20.0;
    params.dSurfaceAngleTolDeg = 20.0;
    params.dMaxEdgeLength = 0.0;
    params.dMaxAspectRatio = 0.0;
    SmPolyBrep *poly = nullptr;
    auto t0 = std::chrono::steady_clock::now();
    CHECK(SmApiTessellate(b, poly, params) == SM_SUCCESS && poly != nullptr);
    double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    double worstVertex = 0.0, worstCentroid = 0.0;
    int over = 0, open = 0;
    long triangles = 0;
    if (poly != nullptr)
    {
        SmTArray<SmPolyFace *> faces;
        SmTArray<SmPolyEdge *> edges;
        poly->GetPolyFaces(faces);
        triangles = static_cast<long>(faces.GetSize());
        std::map<std::pair<const SmPolyVertex *, const SmPolyVertex *>, int> uses;
        for (ULONG i = 0; i < faces.GetSize(); ++i)
        {
            const SmPolyLoop *loop = faces[i]->GetOuterPolyLoop();
            CHECK(loop != nullptr);
            if (loop == nullptr)
            {
                continue;
            }
            loop->GetPolyEdges(edges);
            CHECK(edges.GetSize() == 3);
            if (edges.GetSize() != 3)
            {
                continue;
            }
            SmPoint3d p[3];
            for (ULONG k = 0; k < 3; ++k)
            {
                p[k] = edges[k]->GetStartPolyVertex()->GetPoint();
                uses[{edges[k]->GetStartPolyVertex(), edges[(k + 1) % 3]->GetStartPolyVertex()}]++;
                worstVertex = std::max(worstVertex, std::fabs(dist(p[k])));
            }
            SmPoint3d c((p[0].x + p[1].x + p[2].x) / 3, (p[0].y + p[1].y + p[2].y) / 3, (p[0].z + p[1].z + p[2].z) / 3);
            double d = std::fabs(dist(c));
            worstCentroid = std::max(worstCentroid, d);
            over += d > chord;
        }
        for (auto &e : uses)
        {
            auto back = uses.find({e.first.second, e.first.first});
            open += e.second != 1 || back == uses.end() || back->second != 1;
        }
        poly->SetOKBackPtrs(FALSE);
        delete poly;
    }
    delete b;
    std::printf("%-8s chord %-6g %6ld triangles in %6.1f ms, vertex off %.1e, worst centroid %.4f (%.2fx), %d over, %d open edges\n", name,
                chord, triangles, ms, worstVertex, worstCentroid, worstCentroid / chord, over, open);
    CHECK(triangles > 0);
    CHECK(open == 0);
    CHECK(worstVertex < 1e-6);
    CHECK(worstCentroid <= chord);
}

} // namespace

int main()
{
    const double chords[] = {0.1, 0.05, 0.02, 0.01, 0.005};
    for (double chord : chords)
    {
        SmBrep *b = nullptr;
        SmApiCreateSphere(SmVector3d(0, 0, 0), 12, b);
        Check("sphere", b, chord, [](const SmPoint3d &p) { return std::sqrt(p.x * p.x + p.y * p.y + p.z * p.z) - 12.0; });
    }
    for (double chord : chords)
    {
        SmBrep *b = nullptr;
        SmApiCreateTorus(SmVector3d(0, 0, 0), 10, 3, b);
        Check("torus", b, chord, [](const SmPoint3d &p) { return std::hypot(std::hypot(p.x, p.y) - 10.0, p.z) - 3.0; });
    }
    for (double chord : chords)
    {
        // The side and both caps: the distance to the nearest of the three.
        SmBrep *b = nullptr;
        SmApiCreateCylinder(SmVector3d(0, 0, 0), 5, 12, b);
        Check("cylinder", b, chord,
              [](const SmPoint3d &p)
              { return std::min({std::fabs(std::hypot(p.x, p.y) - 5.0), std::fabs(p.z), std::fabs(p.z - 12.0)}); });
    }
    if (failures == 0)
    {
        std::printf("tess kernel: ok\n");
    }
    return failures == 0 ? 0 : 1;
}
