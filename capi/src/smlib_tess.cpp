// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// smlib_brep_tessellate: SMLib's own tessellator (SmApiTessellate) as plain arrays: vertex
// positions, triangles and, for each triangle, the index of the B-rep face it came from.

#include "smlib_shim.h"

#include <SmApiBrep.h>
#include <SmApiQueries.h>
#include <SmApiTessellationParams.h>
#include <SmFace.h>
#include <SmPoly.h>

#include <memory>
#include <unordered_map>

using namespace smlib_shim;

struct smlib_mesh
{
    std::vector<double> positions;
    std::vector<int64_t> triangles;
    std::vector<int64_t> faces;
};

namespace
{

// Upper bound on the triangles of one mesh; a request that would produce more fails instead of
// exhausting memory.
constexpr size_t MAX_TRIANGLES = size_t(1) << 26;

// Owns a SmPolyBrep made from a B-rep that outlives it.
struct PolyDeleter
{
    void operator()(SmPolyBrep *p) const
    {
        if (p != nullptr)
        {
            p->SetOKBackPtrs(FALSE);
            delete p;
        }
    }
};

// Reads `poly`, whose faces carry their original faces of `brep`, into `out`.
int32_t read_mesh(const SmBrep *brep, const SmPolyBrep *poly, smlib_mesh &out)
{
    SmTArray<SmFace *> faces;
    brep->GetFaces(faces);
    std::unordered_map<const SmFace *, int64_t> face_index;
    for (ULONG i = 0; i < faces.GetSize(); ++i)
    {
        face_index.emplace(faces[i], static_cast<int64_t>(i));
    }

    SmTArray<SmPolyVertex *> vertices;
    poly->GetPolyVertices(vertices);
    std::unordered_map<const SmPolyVertex *, int64_t> vertex_index;
    out.positions.reserve(3 * size_t(vertices.GetSize()));
    for (ULONG i = 0; i < vertices.GetSize(); ++i)
    {
        const SmPoint3d &p = vertices[i]->GetPoint();
        if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z))
        {
            return static_cast<int32_t>(SM_ERR_INVALID_INPUT);
        }
        vertex_index.emplace(vertices[i], static_cast<int64_t>(i));
        out.positions.push_back(p.x);
        out.positions.push_back(p.y);
        out.positions.push_back(p.z);
    }

    SmTArray<SmPolyFace *> polys;
    poly->GetPolyFaces(polys);
    SmTArray<SmPolyEdge *> edges;
    SmTArray<SmPolyLoop *> loops;
    for (ULONG i = 0; i < polys.GetSize(); ++i)
    {
        const SmPolyFace *pf = polys[i];
        pf->GetPolyLoops(loops);
        const SmPolyLoop *loop = pf->GetOuterPolyLoop();
        // SmApiTessellate triangulates: every polygon is one loop of three edges.
        if (loop == nullptr || loops.GetSize() != 1)
        {
            return static_cast<int32_t>(SM_ERR);
        }
        loop->GetPolyEdges(edges);
        if (edges.GetSize() != 3)
        {
            return static_cast<int32_t>(SM_ERR);
        }
        auto f = face_index.find(pf->GetOriginalFace());
        if (f == face_index.end())
        {
            return static_cast<int32_t>(SM_ERR);
        }
        for (ULONG j = 0; j < 3; ++j)
        {
            auto v = vertex_index.find(edges[j]->GetStartPolyVertex());
            if (v == vertex_index.end())
            {
                return static_cast<int32_t>(SM_ERR);
            }
            out.triangles.push_back(v->second);
        }
        out.faces.push_back(f->second);
        if (out.faces.size() > MAX_TRIANGLES)
        {
            return static_cast<int32_t>(SM_ERR_INVALID_INPUT);
        }
    }
    return SMLIB_OK;
}

} // namespace

extern "C" int32_t smlib_brep_tessellate(const smlib_brep *brep, double chord_tolerance, double angle_tolerance_deg,
                                         smlib_mesh **out)
{
    if (!valid(brep) || out == nullptr || !std::isfinite(chord_tolerance) || !(chord_tolerance >= 0.0) ||
        !std::isfinite(angle_tolerance_deg) || !(angle_tolerance_deg >= 0.0) || !(angle_tolerance_deg < 180.0) ||
        (chord_tolerance == 0.0 && angle_tolerance_deg == 0.0))
    {
        return SMLIB_ERR_ARG;
    }
    Call c("smlib_brep_tessellate",
           handle(brep) + ", chord_tolerance=" + num(chord_tolerance) + ", angle_tolerance_deg=" + num(angle_tolerance_deg));
    c.input("brep", brep);
    return c.run(
        [&]
        {
            // SmApiTessellate takes a mutable B-rep and its mesh points back into it: work on a
            // copy, which lives until the mesh is read and freed.
            SmBrep *raw = nullptr;
            int32_t s = status(SmApiBrepCopy(brep->p, raw));
            std::unique_ptr<SmBrep> copy(raw);
            if (s != SMLIB_OK)
            {
                return s;
            }
            if (copy == nullptr)
            {
                return static_cast<int32_t>(SMLIB_ERR_EMPTY);
            }
            SmTessellationParams params;
            params.dChordHeightTol = chord_tolerance;
            params.dCurveAngleTolDeg = angle_tolerance_deg;
            params.dSurfaceAngleTolDeg = angle_tolerance_deg;
            params.dMaxEdgeLength = 0.0;
            params.dMaxAspectRatio = 0.0;
            SmPolyBrep *poly_raw = nullptr;
            s = status(SmApiTessellate(copy.get(), poly_raw, params));
            // Declared after the copy, so freed before it.
            std::unique_ptr<SmPolyBrep, PolyDeleter> poly(poly_raw);
            if (s != SMLIB_OK)
            {
                return s;
            }
            if (poly == nullptr)
            {
                return static_cast<int32_t>(SMLIB_ERR_EMPTY);
            }
            auto mesh = std::make_unique<smlib_mesh>();
            s = read_mesh(copy.get(), poly.get(), *mesh);
            if (s != SMLIB_OK)
            {
                return s;
            }
            if (mesh->faces.empty())
            {
                return static_cast<int32_t>(SMLIB_ERR_EMPTY);
            }
            *out = mesh.release();
            return SMLIB_OK;
        }
    );
}

extern "C" int32_t smlib_mesh_arrays(const smlib_mesh *mesh, const double **positions, int64_t *n_vertices,
                                     const int64_t **triangles, const int64_t **faces, int64_t *n_triangles)
{
    if (mesh == nullptr || positions == nullptr || n_vertices == nullptr || triangles == nullptr || faces == nullptr ||
        n_triangles == nullptr)
    {
        return SMLIB_ERR_ARG;
    }
    *positions = mesh->positions.data();
    *n_vertices = static_cast<int64_t>(mesh->positions.size() / 3);
    *triangles = mesh->triangles.data();
    *faces = mesh->faces.data();
    *n_triangles = static_cast<int64_t>(mesh->faces.size());
    return SMLIB_OK;
}

extern "C" void smlib_mesh_free(smlib_mesh *mesh)
{
    delete mesh;
}
