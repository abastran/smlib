// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// Choosing edges, faces and vertices by points (include/smlib.h, "Selecting edges and faces"), shared by
// the shims that take a selection: each point picks the edge (face) nearest to it, which must lie
// within the caller's tolerance and be nearer than any other by more than the body's modelling
// tolerance. Internal to the library.

#ifndef SMLIB_PICK_H
#define SMLIB_PICK_H

#include "smlib_shim.h"

#include <SmApiQueries.h>
#include <SmEdge.h>
#include <SmFace.h>
#include <SmVertex.h>

#include <vector>

namespace smlib_pick
{

// The entity nearest each point under the picking rule, each once; with `first`, the index of
// the point that first picked each. On failure, the index of the point that matched none or
// several goes to `failed`.
template <typename T, typename Dist>
inline int32_t pick(const SmTArray<T *> &all, const double *points, int64_t n, double tol, double model_tol, Dist dist,
             std::vector<T *> &picked, int64_t *failed, std::vector<int64_t> *first = nullptr)
{
    for (int64_t i = 0; i < n; ++i)
    {
        SmPoint3d p(points[3 * i], points[3 * i + 1], points[3 * i + 2]);
        T *best = nullptr;
        double d0 = HUGE_VAL;
        double d1 = HUGE_VAL;
        for (ULONG k = 0; k < all.GetSize(); ++k)
        {
            double d = 0.0;
            if (!dist(all[k], p, d) || !std::isfinite(d))
            {
                continue;
            }
            if (d < d0)
            {
                d1 = d0;
                d0 = d;
                best = all[k];
            }
            else if (d < d1)
            {
                d1 = d;
            }
        }
        int32_t s = SMLIB_OK;
        if (best == nullptr || d0 > tol)
        {
            s = SMLIB_ERR_PICK_NONE;
        }
        else if (d1 <= tol && d1 - d0 <= model_tol)
        {
            s = SMLIB_ERR_PICK_AMBIGUOUS;
        }
        if (s != SMLIB_OK)
        {
            if (failed != nullptr)
            {
                *failed = i;
            }
            return s;
        }
        if (std::find(picked.begin(), picked.end(), best) == picked.end())
        {
            picked.push_back(best);
            if (first != nullptr)
            {
                first->push_back(i);
            }
        }
    }
    return SMLIB_OK;
}

inline bool edge_distance(SmEdge *e, const SmPoint3d &p, double &d)
{
    SmPoint3d q;
    double t = 0.0;
    return SmApiEdgeClosestPoint(e, p, q, t, d) == SM_SUCCESS;
}

inline bool vertex_distance(SmVertex *v, const SmPoint3d &p, double &d)
{
    d = (v->GetPoint() - p).Length();
    return true;
}

inline bool face_distance(SmFace *f, const SmPoint3d &p, double &d)
{
    SmPoint3d q;
    return SmApiFaceClosestPoint(f, p, q, d) == SM_SUCCESS;
}

inline int32_t pick_edges(SmBrep *b, const double *points, int64_t n, double tol, std::vector<SmEdge *> &out, int64_t *failed,
                   std::vector<int64_t> *first = nullptr)
{
    SmTArray<SmEdge *> all;
    b->GetEdges(all);
    return pick(all, points, n, tol, static_cast<double>(b->GetTolerance()), edge_distance, out, failed, first);
}

inline int32_t pick_faces(SmBrep *b, const double *points, int64_t n, double tol, std::vector<SmFace *> &out, int64_t *failed)
{
    SmTArray<SmFace *> all;
    b->GetFaces(all);
    return pick(all, points, n, tol, static_cast<double>(b->GetTolerance()), face_distance, out, failed);
}

inline int32_t pick_vertices(SmBrep *b, const double *points, int64_t n, double tol, std::vector<SmVertex *> &out, int64_t *failed,
                              std::vector<int64_t> *first = nullptr)
{
    SmTArray<SmVertex *> all;
    b->GetVertices(all);
    return pick(all, points, n, tol, static_cast<double>(b->GetTolerance()), vertex_distance, out, failed, first);
}

} // namespace smlib_pick

#endif
