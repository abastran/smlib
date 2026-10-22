// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// Adapted from source/BREP_SM_USD/src/SmuConvert.cpp (NVIDIA, Apache-2.0) without USD types;
// SolidCraft change.
//
// SmBrep <-> flat arrays for ONE B-rep. The arrays are the per-brep slice of
// UsdBrepData::UsdBrepArrayData: TfToken values became small integers, GfVec2d/GfVec3d/GfVec2i
// became flattened runs of doubles or int64s, and everything about materials, prim paths,
// colours, SdfPath attributes and multiple breps is gone. With one B-rep, every "global" index
// of SmuConvert equals its "local" index, so MapGlobalToLocal / MapLocalToGlobal disappear.
//
// Where SmuConvert reads its arrays on trust (it relies on BrepReadFromUsdStage and
// UsdBrepArraySpans for some checks and on nothing for others), BrepFromFlat validates the
// whole input first (Validate below), so the port itself can index freely.

#include "smlib_flat.h"

#include <array>
#include <SmApiGeneral.h>
#include <SmAttribute.h>
#include <SmBSplineCurve.h>
#include <SmBSplineSurface.h>
#include <SmBrep.h>
#include <SmTrimmingTools.h>
#include <SmBrepData.h>
#include <SmCircle.h>
#include <SmCone.h>
#include <SmCurve.h>
#include <SmCylinder.h>
#include <SmEdge.h>
#include <SmEllipse.h>
#include <SmFace.h>
#include <SmLine.h>
#include <SmLoop.h>
#include <SmPlane.h>
#include <SmRegion.h>
#include <SmShell.h>
#include <SmSphere.h>
#include <SmSurface.h>
#include <SmTol.h>
#include <SmTorus.h>
#include <SmVertex.h>
#include <SmNurbsSrf.h>
#include <nurbs.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <unordered_map>

namespace smlib_flat_impl
{

/* ---- The field table ------------------------------------------------------------------------ */

const std::vector<Field> &Fields()
{
    static const std::vector<Field> fields = []
    {
        std::vector<Field> v;
#define SMLIB_FLAT_ROW_I(id, m) v.push_back(Field{id, #id + 11, &FlatBrep::m, nullptr});
#define SMLIB_FLAT_ROW_F(id, m) v.push_back(Field{id, #id + 11, nullptr, &FlatBrep::m});
        SMLIB_FLAT_FIELDS(SMLIB_FLAT_ROW_I, SMLIB_FLAT_ROW_F)
#undef SMLIB_FLAT_ROW_I
#undef SMLIB_FLAT_ROW_F
        return v;
    }();
    return fields;
}

const Field *FindField(int32_t id)
{
    for (const Field &f : Fields())
    {
        if (f.id == id)
        {
            return &f;
        }
    }
    return nullptr;
}

namespace
{

// Token values of USD as integers (see include/smlib.h).
enum : int64_t
{
    kSurfaceNurbs = SMLIB_FLAT_SURFACE_NURBS,
    kSurfacePlane = SMLIB_FLAT_SURFACE_PLANE,
    kSurfaceCylinder = SMLIB_FLAT_SURFACE_CYLINDER,
    kSurfaceCone = SMLIB_FLAT_SURFACE_CONE,
    kSurfaceSphere = SMLIB_FLAT_SURFACE_SPHERE,
    kSurfaceTorus = SMLIB_FLAT_SURFACE_TORUS,
    kCurveNurbs = SMLIB_FLAT_CURVE_NURBS,
    kCurveLine = SMLIB_FLAT_CURVE_LINE,
    kCurveCircle = SMLIB_FLAT_CURVE_CIRCLE,
    kCurveEllipse = SMLIB_FLAT_CURVE_ELLIPSE,
};

// A face whose arrays carry a NURBS surface block: a NURBS face, or an analytic face given with
// its SMLib NURBS form (FACE_NURBS_FORM), whose range and UV trims are in that NURBS's parameters.
bool FaceHasNurbs(const FlatBrep &F, size_t iFace)
{
    return F.face_surface_type[iFace] == kSurfaceNurbs || (!F.face_nurbs_form.empty() && F.face_nurbs_form[iFace] == 1);
}

// USD writes USDBREP_NO_OBJECT_INDEX (9999999) for "no vertex"; the flat arrays use -1.
constexpr int64_t kNoIndex = -1;

// Caps on hostile input. Topology counts stay below SMLib's index sentinels (SM_NO_OBJECT, and
// MakeTopologyFromData's `m_lUVCurve < 999999` test), NURBS degrees within NLib's NL_MAXDEG.
constexpr int64_t kMaxTopology = 999000;
constexpr int64_t kMaxOrder = NL_MAXDEG + 1;
constexpr double kMaxValue = 1.0e9;

void Push3(std::vector<double> &v, double x, double y, double z)
{
    v.push_back(x);
    v.push_back(y);
    v.push_back(z);
}

SmPoint3d Get3(const std::vector<double> &v, size_t i)
{
    return SmPoint3d(v[3 * i], v[3 * i + 1], v[3 * i + 2]);
}

SmVector3d GetVec3(const std::vector<double> &v, size_t i)
{
    return SmVector3d(v[3 * i], v[3 * i + 1], v[3 * i + 2]);
}

/* ---- Helpers ported from SmuConvert's anonymous namespace ------------------------------------ */

// An edge's or vertex's zone tolerance, or where larger the B-rep's own or the largest of the
// faces around it: SMLib holds an approximated edge (a fillet's rail) to its faces' and its
// B-rep's tolerance, which may be more than the edge's own, and the arrays' edge tolerance bounds
// the edge's distance from its faces.
template <typename T>
double ZoneTolWithFaces(const T *pTopology, double dBrepZoneTol)
{
    double dTol = std::max(static_cast<double>(SmTol::GetZoneTol3d(pTopology)), dBrepZoneTol);
    SmTArray<SmFace *> sFaces;
    pTopology->GetFaces(sFaces);
    for (ULONG ii = 0; ii < sFaces.GetSize(); ++ii)
    {
        const double dFace = static_cast<double>(SmTol::GetZoneTol3d(sFaces[ii]));
        if (std::isfinite(dFace) && dFace > dTol)
        {
            dTol = dFace;
        }
    }
    return dTol;
}

template <typename T>
void AccumulateMaximumTopologyTolerance(const SmTArray<T *> &rTopologies, double &dMaximumTolerance)
{
    for (ULONG ii = 0; ii < rTopologies.GetSize(); ++ii)
    {
        const double dTolerance = static_cast<double>(SmTol::GetZoneTol3d(rTopologies[ii]));
        if (std::isfinite(dTolerance) && dTolerance >= 0.0 && dTolerance > dMaximumTolerance)
        {
            dMaximumTolerance = dTolerance;
        }
    }
}

double GetMaxZoneTol3d(const SmBrep *pBrep)
{
    if (pBrep == nullptr)
    {
        return 0.0;
    }

    double dMaximumTolerance = SmTol::GetZoneTol3d(pBrep);

    SmTArray<SmFace *> sFaces;
    pBrep->GetFaces(sFaces);
    AccumulateMaximumTopologyTolerance(sFaces, dMaximumTolerance);

    SmTArray<SmEdge *> sEdges;
    pBrep->GetEdges(sEdges);
    AccumulateMaximumTopologyTolerance(sEdges, dMaximumTolerance);

    SmTArray<SmVertex *> sVertices;
    pBrep->GetVertices(sVertices);
    AccumulateMaximumTopologyTolerance(sVertices, dMaximumTolerance);

    return dMaximumTolerance;
}

double SnapSMLibAngularBoundaryDeg(double dDegrees)
{
    static constexpr double dSnapTolDeg = 1.0e-9;

    for (int iQuarterTurn = -4; iQuarterTurn <= 4; ++iQuarterTurn)
    {
        const double dBoundary = 90.0 * static_cast<double>(iQuarterTurn);
        if (std::fabs(dDegrees - dBoundary) <= dSnapTolDeg)
        {
            return dBoundary;
        }
    }

    return dDegrees;
}

double UsdRadiansToSMLibDegrees(double dRadians)
{
    return SnapSMLibAngularBoundaryDeg(SM_RAD2DEG(dRadians));
}

/*******************************************************************/ /**
 PURPOSE: rtn TRUE when a flat (one entry per multiplicity) knot vector is clamped, i.e. the
          first and last knot values each repeat at least (degree + 1) times. Some authoring
          kernels emit periodic curves/surfaces with an unclamped knot vector; SMLib's Bezier
          decomposition assumes clamped knots and reads out of bounds on unclamped input.
 ***********************************************************************/
bool smu_IsFlatKnotVectorClamped(const SmTArray<double> &crFlatKnots, uint32_t iDegree)
{
    uint32_t iNum = crFlatKnots.GetSize();

    // Need at least (degree+1) knots at each end to be clamped, i.e. iNum >= 2 * (iDegree + 1).
    // Written as a subtraction/division so a malformed, oversized iDegree cannot overflow uint32
    // in (iDegree + 1) * 2 and slip past the guard; this also guarantees iDegree and
    // (iNum - 1 - iDegree) index in bounds.
    if (iNum < 2 || iDegree > (iNum - 2) / 2)
    {
        return false;
    }

    // clamped => first value repeats (degree+1) times and last value repeats (degree+1) times.
    // Check every entry in each end run rather than just the endpoint pair: comparing only
    // crFlatKnots[0]==crFlatKnots[iDegree] would misclassify a malformed (non-monotonic) flat
    // knot vector as clamped, bypassing CreateCanonical() on exactly the bad-data path we are
    // guarding.
    for (uint32_t ii = 1; ii <= iDegree; ++ii)
    {
        if (crFlatKnots[ii] != crFlatKnots[0] || crFlatKnots[iNum - 1 - ii] != crFlatKnots[iNum - 1])
        {
            return false;
        }
    }

    return true;
}

/*******************************************************************/ /**
 PURPOSE: collapse a flat knot vector (one entry per multiplicity) into the unique knot values
          plus their multiplicities, as expected by SmBSpline{Curve,Surface}::CreateCanonical().
 ***********************************************************************/
void smu_CompressFlatKnots(const SmTArray<double> &crFlatKnots, SmTArray<double> &rUniqueKnots, SmTArray<ULONG> &rMultiplicities)
{
    rUniqueKnots.SetSize(0);
    rMultiplicities.SetSize(0);

    for (uint32_t ii = 0; ii < crFlatKnots.GetSize(); ii++)
    {
        double dKnot = crFlatKnots[ii];
        if (rUniqueKnots.GetSize() != 0 && rUniqueKnots[rUniqueKnots.GetSize() - 1] == dKnot)
        {
            rMultiplicities[rMultiplicities.GetSize() - 1] += 1;
        }
        else
        {
            rUniqueKnots.Add(dKnot);
            rMultiplicities.Add(1);
        }
    }
}

/*******************************************************************/ /**
 PURPOSE: rtn SmBSplineCurve Obj built from BrepCurve array data
 ***********************************************************************/
SmBSplineCurve *smu_CreateSmBSplineCurve_FromArrays(
    const SmContext &crContext,
    const SmTArray<SmPoint3d> &crControlPoints,
    const SmTArray<double> &crWeights,
    const SmTArray<double> &crKnots,
    uint32_t iDegree,
    uint32_t iDim
)
{
    uint32_t lDim = iDim;

    // Allocate memory for curve given degree, number of knots and number of control points
    uint32_t lNumKnots = crKnots.GetSize();
    uint32_t lNumControlPoints = crControlPoints.GetSize();

    // Periodic curves (e.g. UV trim curves on periodic surfaces) arrive with an unclamped knot
    // vector. SMLib's Bezier decomposition assumes clamped knots. When unclamped, route through
    // CreateCanonical(), which clamps the curve to its valid parametric domain before
    // constructing it.
    if (!smu_IsFlatKnotVectorClamped(crKnots, iDegree))
    {
        SmTArray<double> sUniqueKnots;
        SmTArray<ULONG> sMults;
        smu_CompressFlatKnots(crKnots, sUniqueKnots, sMults);

        SmBSplineCurve *pClampedCurve = NULL;
        SmStatus eStatus = SmBSplineCurve::CreateCanonical(
            crContext, iDim, iDegree, crControlPoints, SM_CF_UNSPECIFIED, sMults, sUniqueKnots, SM_KT_UNSPECIFIED, &crWeights, NULL,
            pClampedCurve
        );

        if (eStatus != SM_SUCCESS)
        {
            delete pClampedCurve;
            return NULL;
        }

        return pClampedCurve;
    }

    gw_CURVE *cur = sm_AllocateNurbCurve(lNumControlPoints - 1, (gw_DEGREE)iDegree, lNumKnots - 1);

    if (cur == NULL)
    {
        return NULL;
    }

    gw_CPOINT *Pw;
    gw_REAL *U;
    N_CrvGetCPtsAndKnots(cur, &Pw, &U);

    if (iDim == 3)
    {
        for (uint32_t ii = 0; ii < lNumControlPoints; ii++)
        {
            SmPoint3d sCPt = crControlPoints[ii];
            if (crWeights[ii] == 1.0)
            {
                N_CPtFromWxWyWz(sCPt.x, sCPt.y, sCPt.z, crWeights[ii], &Pw[ii]);
            }
            else
            {
                N_CPtFromWxWyWz(crWeights[ii] * sCPt.x, crWeights[ii] * sCPt.y, crWeights[ii] * sCPt.z, crWeights[ii], &Pw[ii]);
            }
        }
    }
    else
    {
        for (uint32_t ii = 0; ii < lNumControlPoints; ii++)
        {
            SmPoint3d sCPt = crControlPoints[ii];
            if (crWeights[ii] == 1.0)
            {
                N_CPtFromWxWyWz(sCPt.x, sCPt.y, NL_NOZ, 1.0, &Pw[ii]);
            }
            else
            {
                N_CPtFromWxWyWz(crWeights[ii] * sCPt.x, crWeights[ii] * sCPt.y, NL_NOZ, crWeights[ii], &Pw[ii]);
            }
        }
    }

    for (uint32_t ii = 0; ii < lNumKnots; ii++)
    {
        U[ii] = crKnots[ii];
    }

    // SetFromGwNurb copies the temporary NURB, so release that temporary on both success and
    // failure instead of leaking one allocation per imported NURBS curve.
    SmBSplineCurve *pSmCurve = new (crContext) SmBSplineCurve(lDim);
    const SmStatus sSetStatus = pSmCurve->SetFromGwNurb(0, cur);
    smos_Free(cur);
    if (sSetStatus != SM_SUCCESS)
    {
        delete pSmCurve;
        return NULL;
    }

    return pSmCurve;
}

/*******************************************************************/ /**
 PURPOSE: rtn SmBSplineSurface Obj built from BrepSurface array data
 ***********************************************************************/
SmBSplineSurface *smu_CreateSmBSplineSurface_FromArrays(
    const SmContext &crSmContext,
    const SmTArray<SmPoint3d> &crCPts,
    const SmTArray<double> &crWeights,
    const SmTArray<double> &crUKnots,
    const SmTArray<double> &crVKnots,
    uint32_t iUDegree,
    uint32_t iVDegree
)
{
    // USD stores 'Deg+1' knots at each end, for clamped knots.
    uint32_t iNumUKnots = crUKnots.GetSize();
    uint32_t iNumVKnots = crVKnots.GetSize();

    // Control Point Counts:
    uint32_t iNumCPtsU = iNumUKnots - iUDegree - 1;
    uint32_t iNumCPtsV = iNumVKnots - iVDegree - 1;

    // Periodic surfaces arrive with an unclamped knot vector. SMLib's Bezier decomposition
    // (sm_DecomposeSrf) assumes clamped knots and reads out of bounds on unclamped input. When
    // unclamped, route through CreateCanonical(), which clamps the surface to its valid
    // parametric domain before constructing it.
    if (!smu_IsFlatKnotVectorClamped(crUKnots, iUDegree) || !smu_IsFlatKnotVectorClamped(crVKnots, iVDegree))
    {
        SmTArray<double> sUUniqueKnots, sVUniqueKnots;
        SmTArray<ULONG> sUMults, sVMults;
        smu_CompressFlatKnots(crUKnots, sUUniqueKnots, sUMults);
        smu_CompressFlatKnots(crVKnots, sVUniqueKnots, sVMults);

        SmBSplineSurface *pClampedSurf = NULL;
        SmStatus eStatus = SmBSplineSurface::CreateCanonical(
            crSmContext, iUDegree, iVDegree, crCPts, SM_SF_UNSPECIFIED, sUMults, sVMults, sUUniqueKnots, sVUniqueKnots, SM_KT_UNSPECIFIED,
            &crWeights, NULL, pClampedSurf
        );

        if (eStatus != SM_SUCCESS)
        {
            delete pClampedSurf;
            return NULL;
        }

        return pClampedSurf;
    }

    gw_SURFACE *sur = sm_AllocateNurbSurface(iNumCPtsU - 1, iNumCPtsV - 1, iUDegree, iVDegree, iNumUKnots - 1, iNumVKnots - 1);
    if (sur == NULL)
    {
        return NULL;
    }

    gw_CPOINT **Pw;
    gw_REAL *U, *V;
    N_SrfGetCPtsAndKnots(sur, &Pw, &U, &V);

    // Gather all control points move from vector to array
    //  Vec crCPts = [Puov0, Puov1, ..., PuovM, ..., PuMv0, PuMv1, ..., PuMvN]
    // maps to Array Pw = [ P[u0][v0], P[u0][v1], ..., P[u0][vN], ]
    //                    [ . . .,                                ]
    //                    [ P[uM][v0], P[uM][V1], ..., P[uM][vN]  ]
    for (uint32_t ii = 0, usdIdx = 0; ii < iNumCPtsU; ii++)
    {
        for (uint32_t jj = 0; jj < iNumCPtsV; jj++, usdIdx++)
        {
            NL_REAL x(crCPts[usdIdx].x);
            NL_REAL y(crCPts[usdIdx].y);
            NL_REAL z(crCPts[usdIdx].z);
            NL_REAL w(crWeights[usdIdx]);

            N_CPtFromWxWyWz(w * x, w * y, w * z, w, &Pw[ii][jj]);
        }
    }

    for (uint32_t ii = 0; ii < iNumUKnots; ii++)
    {
        U[ii] = crUKnots[ii];
    }

    for (uint32_t ii = 0; ii < iNumVKnots; ii++)
    {
        V[ii] = crVKnots[ii];
    }

    SmBSplineSurface *pNewBSSurf = new (crSmContext) SmBSplineSurface(sur, FALSE, &crSmContext);

    return pNewBSSurf;
}

// The running read positions of SmuConvert's UsdBrepArraySpans geometry start indices.
struct Cursor
{
    size_t nurbs = 0; // index into the per-curve (or per-surface) count/order arrays
    size_t cv = 0;    // control point (and weight) index
    size_t knot = 0;  // knot index (U knots for surfaces)
    size_t vknot = 0; // V knot index (surfaces)
    size_t circle = 0, line = 0, ellipse = 0;
};

// The arrays of one NURBS edge-curve family: edges or wire edges.
struct CurveArrays
{
    const std::vector<int64_t> &count;
    const std::vector<int64_t> &order;
    const std::vector<double> &cvs;
    const std::vector<double> &knots;
    const std::vector<double> &weights;
};

/*******************************************************************/ /**
 PURPOSE: rtn 3d SmBSplineCurve built from the arrays for the next NURBS edge (or wire edge)
          or NULL for errors. Port of smu_CreateSmBSplineCurve3d_ForEdgeIndex and
          smu_CreateSmBSplineCurve3d_ForWireEdgeIndex, which are identical apart from their
          arrays.
 ***********************************************************************/
SmBSplineCurve *smu_CreateSmBSplineCurve3d(const SmContext &crContext, const CurveArrays &a, Cursor &c)
{
    size_t iNurbIndex = c.nurbs;
    uint32_t iDimension = 3;
    uint32_t iDeg = static_cast<uint32_t>(a.order[iNurbIndex] - 1);
    uint32_t iVertexCount = static_cast<uint32_t>(a.count[iNurbIndex]);
    uint32_t iKnotCount = iVertexCount + iDeg + 1;

    SmTArray<SmPoint3d> sCPts(iVertexCount);
    SmTArray<double> sWeights(iVertexCount);
    SmTArray<double> sKnots(iKnotCount);

    for (uint32_t jj = 0; jj < iVertexCount; ++jj)
    {
        sCPts.Add(Get3(a.cvs, c.cv + jj));
    }

    for (uint32_t jj = 0; jj < iKnotCount; ++jj)
    {
        sKnots.Add(a.knots[c.knot + jj]);
    }

    // Set weights, default to 1.0
    if (a.weights.size() * 3 != a.cvs.size())
    {
        sWeights.InsertAt(0, 1.0, iVertexCount);
    }
    else
    {
        for (uint32_t jj = 0; jj < iVertexCount; ++jj)
        {
            sWeights.Add(a.weights[c.cv + jj]);
        }
    }

    SmBSplineCurve *pCurve = smu_CreateSmBSplineCurve_FromArrays(crContext, sCPts, sWeights, sKnots, iDeg, iDimension);

    // increment start indices
    c.nurbs++;
    c.cv += iVertexCount;
    c.knot += iKnotCount;

    return pCurve;
}

/*******************************************************************/ /**
 PURPOSE: rtn 2d SmBSplineCurve built from the arrays for edgeuse iEdgeuseIndex, or NULL when
          it has none, or when its face is analytic.

 NOTES: port of smu_CreateSmBSplineCurve2d_ForEdgeuseIndex. Validate() has already rejected
        out-of-range packed arrays and invalid intervals (which SmuConvert skips with an
        assert), so only the analytic-face case remains here.
 ***********************************************************************/
SmBSplineCurve *smu_CreateSmBSplineCurve2d_ForEdgeuseIndex(
    const SmContext &crContext,
    const FlatBrep &F,
    Cursor &c,
    size_t iEdgeuseIndex,
    bool bFaceHasNurbs,             // in : the owner face carries a NURBS block; only its UV curves are imported
    SmBoolean bEdgeuseSameOrientation // in : TRUE when the UV curve follows the edge's curve direction
)
{
    // no work - no UVTrimCurve
    if (F.edgeuse_nurbs_order[iEdgeuseIndex] == 0 || F.edgeuse_nurbs_vertex_count[iEdgeuseIndex] == 0)
    {
        return NULL;
    }

    uint32_t iDimension = 2;
    uint32_t iDeg = static_cast<uint32_t>(F.edgeuse_nurbs_order[iEdgeuseIndex] - 1);
    uint32_t iVertexCount = static_cast<uint32_t>(F.edgeuse_nurbs_vertex_count[iEdgeuseIndex]);
    uint32_t iKnotCount = iVertexCount + iDeg + 1;
    const size_t lVertexIndex = c.cv;
    const size_t lKnotIndex = c.knot;

    const auto fAdvancePackedCurve = [&]()
    {
        c.cv += iVertexCount;
        c.knot += iKnotCount;
    };

    // An analytic face given without its NURBS form has no parameters its UV curves could be in:
    // discard them and let SMLib construct trims from the 3D edges when they are needed.
    if (!bFaceHasNurbs)
    {
        fAdvancePackedCurve();
        return NULL;
    }

    SmTArray<SmPoint3d> sCPts(iVertexCount);
    SmTArray<double> sWeights(iVertexCount);
    SmTArray<double> sKnots(iKnotCount);

    for (uint32_t jj = 0; jj < iVertexCount; ++jj)
    {
        const size_t lSourceVertexIndex = lVertexIndex + (bEdgeuseSameOrientation ? jj : iVertexCount - 1 - jj);
        sCPts.Add(SmPoint3d(F.edgeuse_nurbs_cvs[2 * lSourceVertexIndex], F.edgeuse_nurbs_cvs[2 * lSourceVertexIndex + 1], 0));
    }

    // SMLib requires every UV trim curve to use its edge's parameter interval. USD permits an
    // independent NURBS knot interval, and analytic edge ranges are radians in USD but degrees
    // in SMLib, so affinely map the UV knot domain to the converted edge range while preserving
    // shape.
    const size_t iEdge = static_cast<size_t>(F.edgeuse_edge[iEdgeuseIndex]);
    double dTargetMin = F.edge_range[2 * iEdge];
    double dTargetMax = F.edge_range[2 * iEdge + 1];
    if (F.edge_curve_type[iEdge] == kCurveCircle || F.edge_curve_type[iEdge] == kCurveEllipse)
    {
        dTargetMin = UsdRadiansToSMLibDegrees(dTargetMin);
        dTargetMax = UsdRadiansToSMLibDegrees(dTargetMax);
    }
    const double dSourceMin = F.edgeuse_nurbs_knots[lKnotIndex + iDeg];
    const double dSourceMax = F.edgeuse_nurbs_knots[lKnotIndex + iVertexCount];
    const double dParameterScale = (dTargetMax - dTargetMin) / (dSourceMax - dSourceMin);
    for (uint32_t jj = 0; jj < iKnotCount; ++jj)
    {
        const double dSourceKnot = bEdgeuseSameOrientation ? F.edgeuse_nurbs_knots[lKnotIndex + jj]
                                                           : dSourceMin + dSourceMax - F.edgeuse_nurbs_knots[lKnotIndex + iKnotCount - 1 - jj];
        sKnots.Add(dTargetMin + (dSourceKnot - dSourceMin) * dParameterScale);
    }

    // Set weights, default to 1.0;
    if (F.edgeuse_nurbs_weights.size() * 2 != F.edgeuse_nurbs_cvs.size())
    {
        sWeights.InsertAt(0, 1.0, iVertexCount);
    }
    else
    {
        for (uint32_t jj = 0; jj < iVertexCount; ++jj)
        {
            const size_t lSourceWeightIndex = lVertexIndex + (bEdgeuseSameOrientation ? jj : iVertexCount - 1 - jj);
            sWeights.Add(F.edgeuse_nurbs_weights[lSourceWeightIndex]);
        }
    }

    SmBSplineCurve *pBSplineCurve = smu_CreateSmBSplineCurve_FromArrays(crContext, sCPts, sWeights, sKnots, iDeg, iDimension);

    fAdvancePackedCurve();

    return pBSplineCurve;
}

/*******************************************************************/ /**
 PURPOSE: rtn 3d SmBSplineSurface built from the arrays for the next NURBS face, or NULL.
          Port of smu_CreateSmBSplineSurface3d_ForFaceIndex.
 ***********************************************************************/
SmBSplineSurface *smu_CreateSmBSplineSurface3d(const SmContext &crContext, const FlatBrep &F, Cursor &c)
{
    const size_t lSurfaceIndex = c.nurbs;
    uint32_t iUDeg = static_cast<uint32_t>(F.face_nurbs_u_order[lSurfaceIndex] - 1);
    uint32_t iVDeg = static_cast<uint32_t>(F.face_nurbs_v_order[lSurfaceIndex] - 1);
    uint32_t iUVertexCount = static_cast<uint32_t>(F.face_nurbs_u_vertex_count[lSurfaceIndex]);
    uint32_t iVVertexCount = static_cast<uint32_t>(F.face_nurbs_v_vertex_count[lSurfaceIndex]);
    uint32_t iUKnotCount = iUVertexCount + iUDeg + 1;
    uint32_t iVKnotCount = iVVertexCount + iVDeg + 1;

    uint32_t iNumCPts = iUVertexCount * iVVertexCount;
    SmTArray<SmPoint3d> sCPts(iNumCPts);
    SmTArray<double> sWeights(iNumCPts);
    SmTArray<double> sUKnots(iUKnotCount);
    SmTArray<double> sVKnots(iVKnotCount);

    for (uint32_t jj = 0; jj < iNumCPts; ++jj)
    {
        sCPts.Add(Get3(F.face_nurbs_cvs, c.cv + jj));
    }

    for (uint32_t jj = 0; jj < iUKnotCount; ++jj)
    {
        sUKnots.Add(F.face_nurbs_u_knots[c.knot + jj]);
    }

    for (uint32_t jj = 0; jj < iVKnotCount; ++jj)
    {
        sVKnots.Add(F.face_nurbs_v_knots[c.vknot + jj]);
    }

    // Set weights, default to 1.0
    if (F.face_nurbs_weights.size() * 3 != F.face_nurbs_cvs.size())
    {
        sWeights.InsertAt(0, 1.0, iNumCPts);
    }
    else
    {
        for (uint32_t jj = 0; jj < iNumCPts; ++jj)
        {
            sWeights.Add(F.face_nurbs_weights[c.cv + jj]);
        }
    }

    SmBSplineSurface *pBSplineSurface = smu_CreateSmBSplineSurface_FromArrays(crContext, sCPts, sWeights, sUKnots, sVKnots, iUDeg, iVDeg);

    // increment start indices
    c.nurbs++;
    c.cv += iNumCPts;
    c.knot += iUKnotCount;
    c.vknot += iVKnotCount;

    return pBSplineSurface;
}

/*******************************************************************/ /**
 PURPOSE: append SmBSplineCurve2d data (Edgeuse->UVTrimCurves); NULL appends an empty curve.
          Port of smu_AddSmBSplineCurve2d_ToBrepArrayData.
 ***********************************************************************/
SmStatus smu_AddSmBSplineCurve2d(SmBSplineCurve *pBSplineCurve2d, SmBoolean bSameOrientation, FlatBrep &F, const double *pRemap = nullptr, const double *pAffine = nullptr)
{
    // For NULL UVTrimCurves
    if (pBSplineCurve2d == NULL)
    {
        // add VertexCount = 0 curve to preserve the edgeuse/UVTrimCurve one-to-one relationship
        F.edgeuse_nurbs_vertex_count.push_back(0);
        F.edgeuse_nurbs_order.push_back(0);
        return SM_SUCCESS;
    }

    gw_CURVE *nurb = pBSplineCurve2d->GetOrCreateGwNurbPointer();
    if (nurb == NULL)
    {
        return SM_ERR;
    }
    const bool bReverse = bSameOrientation != TRUE; // reverse the curve if the edgeuse runs in the opposite direction as its edge
    NL_CURVE sReversedNurb;
    NL_STACKS SC;
    if (bReverse)
    {
        N_InitNurbs(&SC);
        N_CrvInitArrays(&sReversedNurb);
        if (NL_YES == N_CrvReverse(nurb, &sReversedNurb, &SC))
        {
            N_EndNurbs(&SC);
            return SM_ERR;
        }
        nurb = &sReversedNurb;
    }

    NL_REAL *knts = nurb->knt->U;
    uint32_t lNumKnots = nurb->knt->m + 1;
    NL_DEGREE lDegree = nurb->p;
    NL_CPOINT *cpts = nurb->pol->Pw;
    uint32_t lNumCpts = nurb->pol->n + 1;

    F.edgeuse_nurbs_vertex_count.push_back(lNumCpts);
    F.edgeuse_nurbs_order.push_back(lDegree + 1);

    for (uint32_t i = 0; i < lNumKnots; i++)
    {
        double k = static_cast<double>(knts[i]);
        if (pRemap != nullptr)
        {
            // SolidCraft: the edge is exported as an analytic curve whose imported parameter
            // runs over [pRemap[2], pRemap[3]] instead of the edge's [pRemap[0], pRemap[1]];
            // map the knots the same way so the trim's parameter keeps matching its edge's.
            k = k == pRemap[1] ? pRemap[3] : pRemap[2] + (k - pRemap[0]) / (pRemap[1] - pRemap[0]) * (pRemap[3] - pRemap[2]);
        }
        F.edgeuse_nurbs_knots.push_back(k);
    }

    for (uint32_t i = 0; i < lNumCpts; ++i)
    {
        double x = cpts[i].x, y = cpts[i].y;
        double sWeight = 1.0;
        if (cpts[i].w != NL_NOW)
        {
            x /= cpts[i].w;
            y /= cpts[i].w;
            sWeight = cpts[i].w;
        }
        if (pAffine != nullptr)
        {
            // An affine map of the face's parameters (a plane's NURBS (u, v) onto the arrays'):
            // exact for a B-spline, its weights unchanged.
            const double dX = pAffine[0] + pAffine[1] * x + pAffine[2] * y;
            y = pAffine[3] + pAffine[4] * x + pAffine[5] * y;
            x = dX;
        }
        F.edgeuse_nurbs_cvs.push_back(x);
        F.edgeuse_nurbs_cvs.push_back(y);
        F.edgeuse_nurbs_weights.push_back(sWeight);
    }

    if (bReverse)
    {
        N_EndNurbs(&SC);
    }

    return SM_SUCCESS;
}

/*******************************************************************/ /**
 PURPOSE: append SmBSplineCurve3d data (Edge->Curves and WireEdge->Curves).
          Port of smu_AddSmBSplineCurve3d_ToBrepArrayData.
 ***********************************************************************/
SmStatus smu_AddSmBSplineCurve3d(
    SmBSplineCurve &rBSplineCurve3d,
    std::vector<int64_t> &count,
    std::vector<int64_t> &order,
    std::vector<double> &cvs,
    std::vector<double> &knots,
    std::vector<double> &weights
)
{
    gw_CURVE *nurb = rBSplineCurve3d.GetOrCreateGwNurbPointer();
    if (nurb == NULL)
    {
        return SM_ERR;
    }
    NL_REAL *knts = nurb->knt->U;
    uint32_t lNumKnots = nurb->knt->m + 1;
    NL_DEGREE lDegree = nurb->p;
    NL_CPOINT *cpts = nurb->pol->Pw;
    uint32_t lNumCpts = nurb->pol->n + 1;

    count.push_back(lNumCpts);
    order.push_back(lDegree + 1);

    for (uint32_t i = 0; i < lNumKnots; i++)
    {
        knots.push_back(static_cast<double>(knts[i]));
    }

    for (uint32_t i = 0; i < lNumCpts; ++i)
    {
        double x = cpts[i].x, y = cpts[i].y, z = cpts[i].z;
        double sWeight = 1.0;
        if (cpts[i].w != NL_NOW)
        {
            x /= cpts[i].w;
            y /= cpts[i].w;
            z /= cpts[i].w;
            sWeight = cpts[i].w;
        }
        Push3(cvs, x, y, z);
        weights.push_back(sWeight);
    }

    return SM_SUCCESS;
}

/*******************************************************************/ /**
 PURPOSE: append an SMLib surface. Port of smu_AddSmBSplineSurface_ToBrepArrayData.
 ***********************************************************************/
SmStatus smu_AddSmBSplineSurface(const SmBSplineSurface &crBSplineSurface, FlatBrep &F)
{
    NL_SURFACE *surf = ((SmBSplineSurface &)crBSplineSurface).GetOrCreateGwNurbPointer();
    if (surf == NULL)
    {
        return SM_ERR;
    }

    NL_CNET *cnet = surf->net;
    NL_REAL *uKnts = surf->knu->U;
    NL_REAL *vKnts = surf->knv->U;
    uint32_t lNumUKnots = surf->knu->m + 1;
    uint32_t lNumVKnots = surf->knv->m + 1;
    uint32_t degreeU = surf->p;
    uint32_t degreeV = surf->q;
    uint32_t numControlPointsInU = cnet->n + 1;
    uint32_t numControlPointsInV = cnet->m + 1;

    // Check state: surfaces have nonDegenerate knot vectors - else abort
    if (!(uKnts[0] < uKnts[lNumUKnots - 1]) || !(vKnts[0] < vKnts[lNumVKnots - 1]))
    {
        return SM_ERR; // "NURBS surface has degenerate knots and can't be converted"
    }

    F.face_nurbs_u_vertex_count.push_back(numControlPointsInU);
    F.face_nurbs_v_vertex_count.push_back(numControlPointsInV);

    // order = degree + 1 values
    F.face_nurbs_u_order.push_back(degreeU + 1);
    F.face_nurbs_v_order.push_back(degreeV + 1);

    for (uint32_t i = 0; i < lNumUKnots; i++)
    {
        F.face_nurbs_u_knots.push_back(static_cast<double>(uKnts[i]));
    }

    for (uint32_t i = 0; i < lNumVKnots; i++)
    {
        F.face_nurbs_v_knots.push_back(static_cast<double>(vKnts[i]));
    }

    // Control points U-major: P[u0][v0], P[u0][v1], ..., P[uM][vN]
    for (uint32_t u = 0; u < numControlPointsInU; ++u)
    {
        for (uint32_t v = 0; v < numControlPointsInV; ++v)
        {
            NL_CPOINT control_point = cnet->Pw[u][v];
            double x = control_point.x, y = control_point.y, z = control_point.z;
            double dWeight = 1.0;

            if (control_point.w != NL_NOW)
            {
                dWeight = control_point.w;
                x /= dWeight;
                y /= dWeight;
                z /= dWeight;
            }

            Push3(F.face_nurbs_cvs, x, y, z);
            F.face_nurbs_weights.push_back(dWeight);
        }
    }

    return SM_SUCCESS;
}

/*******************************************************************/ /**
 PURPOSE: Return the compacted index for a SmBrepData index, given the sorted SmBrepData
          indices removed from the list. Port of smu_CompactLocalIndex (which loops over every
          removed index; this binary-searches the sorted list, same result).
 ***********************************************************************/
int64_t smu_CompactLocalIndex(int64_t lBrepData_LocalIndex, const std::vector<int64_t> &rSortedRemoved)
{
    return lBrepData_LocalIndex - (std::lower_bound(rSortedRemoved.begin(), rSortedRemoved.end(), lBrepData_LocalIndex) - rSortedRemoved.begin());
}

// An edge or wire-edge curve as NURBS: the curve itself when it is a B-spline (lines,
// circles and ellipses are B-spline subclasses, so this is their exact NURBS form), else an
// approximation within sApproxTol3d over the edge interval, owned by `clean`.
SmBSplineCurve *EdgeCurveAsBSpline(SmCurve *pCurve, const SmExtent1d &crInterval, SmApproxTol3d sApproxTol3d, SmObjDelete &clean)
{
    if (pCurve == NULL)
    {
        return NULL;
    }
    // for now - force every curve to be a BSplineCurve - later may add support for other curve types
    if (pCurve->IsKindOf(SmBSplineCurve_TYPE))
    {
        return SM_REINTERPRET_CAST(SmBSplineCurve *, pCurve);
    }
    SmTArray<double> vBreaks;
    vBreaks.Add(crInterval.GetMin());
    vBreaks.Add(crInterval.GetMax());
    double dAchievedTol = 0.0;
    SmBSplineCurve *pBSplineCurve = NULL;
    SmStatus stat = pCurve->ApproximateCurve(
        *(pCurve->GetContext()), // in : memory context for new object
        vBreaks,                 // in : input curve params exactly interpolated by approx curve.
        sApproxTol3d,            // in : Max dist allowed between approx and original curves.
        dAchievedTol,            // out: Max dist between approx and original curves seen at sample points.
        pBSplineCurve            // out: the new BSpline (always Dim==3)
    );
    clean.SetObj(pBSplineCurve);
    // A curve SMLib could not approximate to the tolerance is none: the caller refuses NULL.
    return stat == SM_SUCCESS ? pBSplineCurve : NULL;
}

// The 3D curve SmuConvert's import builds for a LINE, CIRCLE or ELLIPSE edge: an unbounded line
// P + s * V, or a circle or ellipse over [dMin, dMax] degrees whose NURBS parameter is then set
// to that same interval (EditParameterization). dMin/dMax are the edge range in SMLib units.
// Export uses it too, to check what an import will rebuild.
SmStatus BuildAnalyticEdgeCurve(
    const SmContext &ctx,
    int64_t type,
    double dMin,
    double dMax,
    const SmPoint3d &center,
    const SmVector3d &axis, // the line's direction for a line
    const SmVector3d &ref,
    double dRadius1,
    double dRadius2,
    SmCurve *&out
)
{
    out = NULL;
    if (type == kCurveLine)
    {
        SmLine *pLine = NULL;
        SmStatus stat = SmLine::CreateCanonical(ctx, center, axis, pLine);
        if (stat != SM_SUCCESS || pLine == NULL)
        {
            delete pLine;
            return stat != SM_SUCCESS ? stat : SM_ERR;
        }
        out = pLine;
        return SM_SUCCESS;
    }
    // Build coordinate system: X = refDirection, Y = axis cross refDirection
    SmVector3d sYAxis = axis * ref;
    SmAxis2Placement sPlacement;
    SmStatus placed = sPlacement.SetCanonical(center, ref, sYAxis);
    if (placed != SM_SUCCESS)
    {
        return placed;
    }
    SmExtent1d sDegreeInterval(dMin, dMax);
    if (type == kCurveCircle)
    {
        SmCircle *pCircle = NULL;
        SmStatus stat = SmCircle::CreateCanonical(ctx, sPlacement, dRadius1, pCircle, &sDegreeInterval);
        if (stat != SM_SUCCESS || pCircle == NULL)
        {
            delete pCircle;
            return stat != SM_SUCCESS ? stat : SM_ERR;
        }
        stat = pCircle->EditParameterization(sDegreeInterval, FALSE);
        if (stat != SM_SUCCESS)
        {
            delete pCircle;
            return stat;
        }
        out = pCircle;
        return SM_SUCCESS;
    }
    SmEllipse *pEllipse = NULL;
    SmStatus stat = SmEllipse::CreateCanonical(ctx, sPlacement, dRadius1, dRadius2, pEllipse, &sDegreeInterval);
    if (stat != SM_SUCCESS || pEllipse == NULL)
    {
        delete pEllipse;
        return stat != SM_SUCCESS ? stat : SM_ERR;
    }
    stat = pEllipse->EditParameterization(sDegreeInterval, FALSE);
    if (stat != SM_SUCCESS)
    {
        delete pEllipse;
        return stat;
    }
    out = pEllipse;
    return SM_SUCCESS;
}

// How one edge's curve is exported (SolidCraft change: SmuConvert exports every edge curve as
// NURBS; this exports lines, circles and ellipses as their own kinds, in the arrays and units
// its import reads).
struct EdgeExport
{
    int64_t type = kCurveNurbs;
    double lo = 0.0, hi = 0.0;     // the exported range: radians for circles and ellipses
    double smLo = 0.0, smHi = 0.0; // the imported edge's interval (SMLib units: degrees for circles and ellipses)
    SmPoint3d center;              // circle/ellipse center, line origin
    SmVector3d axis;               // circle/ellipse axis, line direction
    SmVector3d ref;                // circle/ellipse reference direction
    double r1 = 0.0, r2 = 0.0;     // radius; x and y radius for ellipses
    // The imported curve's parameter is an affine image of this edge's, so UV trims can be
    // carried by mapping their knots; otherwise they are dropped and SMLib rebuilds them.
    bool affine = false;
};

// Decides how to export the curve of an edge with NURBS interval `ivl`. The most-derived class
// is tested first: SmCircle derives from SmEllipse, and every analytic curve from
// SmBSplineCurve. The edge's interval ends are converted with the curve's own
// ConvertTFromNURBSToSTEP; a curve whose STEP parameter runs against its NURBS parameter
// (InsideOut) is written with its axis (or direction) reversed and its range negated, which
// is the same point set traversed the same way. The analytic form is used only when the
// curve the import builds from it matches this one pointwise (at the edge's ends and along
// the STEP parameter); anything else stays NURBS.
EdgeExport ClassifyEdgeCurve(const SmContext &ctx, SmCurve *pCurve, const SmExtent1d &ivl)
{
    EdgeExport e;
    EdgeExport nurbs;
    if (pCurve == NULL)
    {
        return nurbs;
    }
    if (pCurve->IsKindOf(SmCircle_TYPE))
    {
        SmAxis2Placement sOrigin;
        double dRadius = 0.0;
        if (SM_REINTERPRET_CAST(SmCircle *, pCurve)->GetCanonical(sOrigin, dRadius) != SM_SUCCESS)
        {
            return nurbs;
        }
        e.type = kCurveCircle;
        e.center = sOrigin.GetOrigin();
        e.axis = sOrigin.GetZAxis();
        e.ref = sOrigin.GetXAxis();
        e.r1 = e.r2 = dRadius;
    }
    else if (pCurve->IsKindOf(SmEllipse_TYPE))
    {
        SmAxis2Placement sOrigin;
        double dRadius1 = 0.0, dRadius2 = 0.0;
        if (SM_REINTERPRET_CAST(SmEllipse *, pCurve)->GetCanonical(sOrigin, dRadius1, dRadius2) != SM_SUCCESS)
        {
            return nurbs;
        }
        e.type = kCurveEllipse;
        e.center = sOrigin.GetOrigin();
        e.axis = sOrigin.GetZAxis();
        e.ref = sOrigin.GetXAxis();
        e.r1 = dRadius1;
        e.r2 = dRadius2;
    }
    else if (pCurve->IsKindOf(SmLine_TYPE))
    {
        SmPoint3d sPoint;
        SmVector3d sVector; // Line(s) = Pt + s * Vec, Vec not unit length
        if (SM_REINTERPRET_CAST(SmLine *, pCurve)->GetCanonical(sPoint, sVector) != SM_SUCCESS)
        {
            return nurbs;
        }
        e.type = kCurveLine;
        e.center = sPoint;
        e.axis = sVector;
    }
    else
    {
        return nurbs;
    }

    double s0 = 0.0, s1 = 0.0;
    if (pCurve->ConvertTFromNURBSToSTEP(ivl.GetMin(), s0) != SM_SUCCESS || pCurve->ConvertTFromNURBSToSTEP(ivl.GetMax(), s1) != SM_SUCCESS ||
        !std::isfinite(s0) || !std::isfinite(s1))
    {
        return nurbs;
    }
    // sNew = sign * sOld + shift
    double sign = 1.0, shift = 0.0;
    if (s0 > s1)
    {
        sign = -1.0;
        s0 = -s0;
        s1 = -s1;
        e.axis = -e.axis;
    }
    if (!(s1 > s0))
    {
        return nurbs;
    }
    if (e.type == kCurveLine)
    {
        if (std::fabs(s0) > SM_INFINITE_PARAMETER || std::fabs(s1) > SM_INFINITE_PARAMETER)
        {
            return nurbs;
        }
        e.lo = e.smLo = s0;
        e.hi = e.smHi = s1;
    }
    else
    {
        // SmCircle/SmEllipse::CreateCanonical take intervals within [-360, 360] degrees.
        if (s1 - s0 > 360.0 + 1e-9)
        {
            return nurbs;
        }
        while (s1 > 360.0 + 1e-9)
        {
            s0 -= 360.0;
            s1 -= 360.0;
            shift -= 360.0;
        }
        while (s0 < -360.0 - 1e-9)
        {
            s0 += 360.0;
            s1 += 360.0;
            shift += 360.0;
        }
        if (s1 > 360.0 + 1e-9)
        {
            return nurbs;
        }
        e.lo = SM_DEG2RAD(s0);
        e.hi = SM_DEG2RAD(s1);
        e.smLo = UsdRadiansToSMLibDegrees(e.lo); // exactly what the import computes
        e.smHi = UsdRadiansToSMLibDegrees(e.hi);
    }

    SmCurve *pNew = NULL;
    if (BuildAnalyticEdgeCurve(ctx, e.type, e.smLo, e.smHi, e.center, e.axis, e.ref, e.r1, e.r2, pNew) != SM_SUCCESS || pNew == NULL)
    {
        delete pNew;
        return nurbs;
    }
    SmObjDelete sCleanNew(pNew);

    SmPoint3d a, b;
    if (pCurve->EvaluatePoint(ivl.GetMin(), a) != SM_SUCCESS || pCurve->EvaluatePoint(ivl.GetMax(), b) != SM_SUCCESS)
    {
        return nurbs;
    }
    const double dScale = 1.0 + e.r1 + e.r2 + a.GetMaxDimension() + b.GetMaxDimension();
    const double dTol = 1e-9 * dScale;
    bool bAffine = true;
    for (int k = 0; k <= 8; ++k)
    {
        const double f = k / 8.0;
        // the same point by the STEP parameter of each curve
        const double sNew = e.smLo + f * (e.smHi - e.smLo);
        const double sOld = (sNew - shift) * sign;
        double tOld = 0.0, tNew = 0.0;
        SmPoint3d pOld, pNew2, pAffOld, pAffNew;
        if (pCurve->ConvertTFromSTEPToNURBS(sOld, tOld) != SM_SUCCESS || pNew->ConvertTFromSTEPToNURBS(sNew, tNew) != SM_SUCCESS ||
            pCurve->EvaluatePoint(tOld, pOld) != SM_SUCCESS || pNew->EvaluatePoint(tNew, pNew2) != SM_SUCCESS || pOld.DistanceBetween(pNew2) > dTol)
        {
            return nurbs;
        }
        // the same point by each edge's NURBS parameter, affinely matched
        if (pCurve->EvaluatePoint(ivl.GetMin() + f * (ivl.GetMax() - ivl.GetMin()), pAffOld) != SM_SUCCESS ||
            pNew->EvaluatePoint(e.smLo + f * (e.smHi - e.smLo), pAffNew) != SM_SUCCESS)
        {
            return nurbs;
        }
        const double d = pAffOld.DistanceBetween(pAffNew);
        if ((k == 0 || k == 8) && d > dTol)
        {
            return nurbs; // the ends must match exactly as the edge's vertices do
        }
        bAffine = bAffine && d <= dTol;
    }
    e.affine = bAffine;
    return e;
}

// The arrays of one edge family (edges or wire edges) that export writes.
struct EdgeArraysOut
{
    std::vector<int64_t> &type;
    std::vector<double> &range;
    std::vector<int64_t> &nurbsCount, &nurbsOrder;
    std::vector<double> &nurbsCvs, &nurbsKnots, &nurbsWeights;
    std::vector<double> &circleCenter, &circleAxis, &circleRef, &circleRadius;
    std::vector<double> &lineOrigin, &lineDirection;
    std::vector<double> &ellipseCenter, &ellipseAxis, &ellipseRef, &ellipseX, &ellipseY;
};

// Writes one edge's curve kind, range and geometry.
SmStatus WriteEdgeCurve(EdgeArraysOut &o, const EdgeExport &e, SmCurve *pCurve, const SmExtent1d &ivl, SmApproxTol3d sApproxTol3d)
{
    o.type.push_back(e.type);
    if (e.type == kCurveNurbs)
    {
        o.range.push_back(ivl.GetMin());
        o.range.push_back(ivl.GetMax());
        SmObjDelete sCleanBSplineCurve;
        SmBSplineCurve *pBSplineCurve = EdgeCurveAsBSpline(pCurve, ivl, sApproxTol3d, sCleanBSplineCurve);
        if (pBSplineCurve == NULL)
        {
            return SM_ERR; // SmuConvert asserts and dereferences NULL
        }
        return smu_AddSmBSplineCurve3d(*pBSplineCurve, o.nurbsCount, o.nurbsOrder, o.nurbsCvs, o.nurbsKnots, o.nurbsWeights);
    }
    o.range.push_back(e.lo);
    o.range.push_back(e.hi);
    if (e.type == kCurveLine)
    {
        Push3(o.lineOrigin, e.center.x, e.center.y, e.center.z);
        Push3(o.lineDirection, e.axis.x, e.axis.y, e.axis.z);
    }
    else if (e.type == kCurveCircle)
    {
        Push3(o.circleCenter, e.center.x, e.center.y, e.center.z);
        Push3(o.circleAxis, e.axis.x, e.axis.y, e.axis.z);
        Push3(o.circleRef, e.ref.x, e.ref.y, e.ref.z);
        o.circleRadius.push_back(e.r1);
    }
    else
    {
        Push3(o.ellipseCenter, e.center.x, e.center.y, e.center.z);
        Push3(o.ellipseAxis, e.axis.x, e.axis.y, e.axis.z);
        Push3(o.ellipseRef, e.ref.x, e.ref.y, e.ref.z);
        o.ellipseX.push_back(e.r1);
        o.ellipseY.push_back(e.r2);
    }
    return SM_SUCCESS;
}

} // namespace

/*******************************************************************/ /**
 PURPOSE: SmBrep -> FlatBrep. Port of SMU_BrepConvert::BrepAppend_SMLibToUsdBrep for one B-rep
          with bExportUVCurves = true and bBoundUnboundedFaceRanges = false (SmuConvert's
          defaults; bounding face ranges may attach UV trims to the source B-rep, and the C API
          never changes its inputs).

 NOTES: SmBrepData vs. flat organization (from SmuConvert):
   BrepData            | flat arrays
   Region Objects      | Region arrays
   Shell Objects       | Shell arrays
   Faceuse Objects     | Faceuse arrays
   Face Objects        | Face arrays
   Loop Objects        | Loop arrays
   Edgeuse Objects     | Edgeuse arrays, culled of wireEdge-Edgeuses
     m_vEdgeuses       |   == [Edgeuses, 1 Edgeuse per WireEdge]
   Edge Objects        | divided into Edge and WireEdge arrays
     m_vEdges          |   == [FaceEdges, WireEdges]
   Vertex Objects      | divided into ShellPoint and Vertex arrays
     m_vVertices       |   == [ShellVertices, Vertices]
   so every stored vertex, edge and edgeuse index is compacted to account for the split.
 ***********************************************************************/
SmStatus FlatFromBrep(const SmBrep &crBrep, FlatBrep &F)
{
    F = FlatBrep();
    const SmContext *cpContext = crBrep.GetContext();
    if (cpContext == NULL)
    {
        return SM_ERR;
    }
    const SmBrep *pBrep_ii = &crBrep;
    SmTArray<SmAttribute *> sAttributes;    // FromBrep's attribute list; attributes are not exported
    std::vector<int64_t> sVertexShellIdx;   // BrepData vertex indices of VertexShells (removed from the vertex arrays)
    std::vector<int64_t> sWireEdgeIdx;      // BrepData edge indices of WireEdges (removed from the edge arrays)
    std::vector<int64_t> sNonFaceEdgeuseIdx; // BrepData edgeuse indices of WireEdge edgeuses (removed from the edgeuse arrays)

    // A curved face SMLib made without UV trims (MakeFaceWithCurves attaches them only when its
    // 3D curves were made from them) gets SMLib's own, made as SMLib makes them when it needs
    // them, so every exported trim is SMLib's.
    {
        SmTArray<SmFace *> sFaces;
        pBrep_ii->GetFaces(sFaces);
        for (ULONG ii = 0; ii < sFaces.GetSize(); ++ii)
        {
            const SmFace *pFace = sFaces[ii];
            if (pFace == NULL || pFace->GetSurface() == NULL || pFace->GetSurface()->IsKindOf(SmPlane_TYPE) ||
                SmTrimmingTools::AreUVCurvesPresent(pFace) == TRUE)
            {
                continue;
            }
            double dMeanGap = 0.0, dMaxGap = 0.0, dMeanVertexGap = 0.0, dMaxVertexGap = 0.0, dMeanCurveGap = 0.0, dMaxCurveGap = 0.0;
            SmStatus made = pFace->CreateUVTrimCurves(FALSE, NULL, NULL, dMeanGap, dMaxGap, dMeanVertexGap, dMaxVertexGap, dMeanCurveGap, dMaxCurveGap);
            if (made != SM_SUCCESS)
            {
                return made;
            }
        }
    }

    // temp BrepData_ii built from Brep_ii; each SmBrepData topology object points back to its
    // source SmBrep topology object (example: SmVertexData::m_pVertex).
    SmBrepData *pBrepData_ii = new (*cpContext) SmBrepData(FALSE, SM_DS_SMLIB);
    SmObjDelete sClean(pBrepData_ii);
    SmStatus sFrom = pBrepData_ii->FromBrep(*pBrep_ii, sAttributes);
    if (sFrom != SM_SUCCESS)
    {
        return sFrom;
    }

    // Get number of each object in this Brep
    uint32_t lNumRegions = pBrepData_ii->m_lNumRegions;
    uint32_t lNumShells = pBrepData_ii->m_vShells.GetSize();
    uint32_t lNumFaceuses = pBrepData_ii->m_vFaceuses.GetSize();
    uint32_t lNumFaces = pBrepData_ii->m_vFaces.GetSize();
    uint32_t lNumLoops = pBrepData_ii->m_vLoops.GetSize();
    uint32_t lNumEdgeuses = pBrepData_ii->m_vEdgeuses.GetSize();
    uint32_t lNumEdges = pBrepData_ii->m_vEdges.GetSize();
    uint32_t lNumVertices = pBrepData_ii->m_vVertices.GetSize();
    std::vector<SmBoolean> sLoopExportsUVCurves(lNumLoops, FALSE);
    // Per loop: the affine map of its plane face's NURBS (u, v) onto the arrays' plane (u, v), or none.
    std::vector<std::array<double, 6>> sLoopPlaneMap(lNumLoops);
    std::vector<SmBoolean> sLoopHasPlaneMap(lNumLoops, FALSE);

    // The arrays store the tolerance needed to preserve the source topology: the largest
    // effective face, edge, or vertex tolerance.
    double sZoneTol3d = GetMaxZoneTol3d(pBrep_ii);
    const double dBrepZoneTol = static_cast<double>(SmTol::GetZoneTol3d(pBrep_ii));
    SmXSectTol3d sXSectTol3d = SM_ZONE_TO_XSECTTOL3D(sZoneTol3d);
    SmApproxTol3d sApproxTol3d = SmTol::GetApproxTol3d(sXSectTol3d);

    // Brep_ii bounding box (== Usd Extent)
    SmExtent3d sLocalBBox;
    SmStatus boxed = pBrep_ii->CalculateBoundingBox(sLocalBBox, FALSE, FALSE, FALSE);
    if (boxed != SM_SUCCESS)
    {
        return boxed;
    }

    // Brep data
    F.tolerance.push_back(sXSectTol3d);
    F.extent = {sLocalBBox.GetUMin(), sLocalBBox.GetVMin(), sLocalBBox.GetWMin(), sLocalBBox.GetUMax(), sLocalBBox.GetVMax(), sLocalBBox.GetWMax()};
    F.source.push_back(1);

    // Edge pointer -> BrepData edge index (SmuConvert searches linearly: smu_FindEdgeDataObj_ForSmEdge)
    std::unordered_map<const SmEdge *, int64_t> sEdgeIndex;
    for (uint32_t jj = 0; jj < lNumEdges; ++jj)
    {
        sEdgeIndex.emplace(pBrepData_ii->m_vEdges[jj].m_pEdge, jj);
    }

    // How each edge's curve is exported (needed by the edgeuses' UV trims before the edges).
    std::vector<EdgeExport> sEdgeExports(lNumEdges);
    for (uint32_t jj = 0; jj < lNumEdges; ++jj)
    {
        SmEdgeData &rEdgeData = pBrepData_ii->m_vEdges[jj];
        if (rEdgeData.m_pEdge != NULL)
        {
            sEdgeExports[jj] = ClassifyEdgeCurve(*cpContext, rEdgeData.m_pEdge->GetCurve(), rEdgeData.m_vInterval);
        }
    }
    EdgeArraysOut sEdgeOut{F.edge_curve_type, F.edge_range, F.edge_nurbs_vertex_count, F.edge_nurbs_order, F.edge_nurbs_cvs,
                           F.edge_nurbs_knots, F.edge_nurbs_weights, F.edge_circle_center, F.edge_circle_axis, F.edge_circle_ref_direction,
                           F.edge_circle_radius, F.edge_line_origin, F.edge_line_direction, F.edge_ellipse_center, F.edge_ellipse_axis,
                           F.edge_ellipse_ref_direction, F.edge_ellipse_x_radius, F.edge_ellipse_y_radius};
    EdgeArraysOut sWireEdgeOut{F.wire_edge_curve_type, F.wire_edge_range, F.wire_edge_nurbs_vertex_count, F.wire_edge_nurbs_order,
                               F.wire_edge_nurbs_cvs, F.wire_edge_nurbs_knots, F.wire_edge_nurbs_weights, F.wire_edge_circle_center,
                               F.wire_edge_circle_axis, F.wire_edge_circle_ref_direction, F.wire_edge_circle_radius, F.wire_edge_line_origin,
                               F.wire_edge_line_direction, F.wire_edge_ellipse_center, F.wire_edge_ellipse_axis, F.wire_edge_ellipse_ref_direction,
                               F.wire_edge_ellipse_x_radius, F.wire_edge_ellipse_y_radius};

    // Regions
    for (uint32_t jj = 0; jj < lNumRegions; ++jj)
    {
        SmRegionData &rRegionData = pBrepData_ii->m_vRegions[jj];
        // rRegionData.m_bIsVoidFlag uses UNSURE for notSet; region:type has no notSet value.
        F.region_shell_count.push_back(rRegionData.m_lNumShells);
        F.region_type.push_back(rRegionData.m_bIsVoidFlag ? 1 : 0);
    }

    // Shells, WireEdges, WireEdge->Curves, and ShellVertex Points
    std::vector<int64_t> sWireEdgeRawVertices; // BrepData vertex indices, compacted below
    for (uint32_t jj = 0; jj < lNumShells; ++jj)
    {
        SmShellData &rShellData = pBrepData_ii->m_vShells[jj];
        SmShell *pShell = rShellData.m_pShell1;
        if (pShell == NULL)
        {
            return SM_ERR;
        }

        SmTArray<SmFaceuse *> sShellFaceuses;
        SmTArray<SmEdge *> sShellWireEdges;
        pShell->GetFaceuses(sShellFaceuses);
        SmStatus listed = pShell->GetWireEdges(sShellWireEdges);
        if (listed != SM_SUCCESS)
        {
            return listed;
        }

        uint32_t lFaceuseCnt = sShellFaceuses.GetSize();
        uint32_t lWireEdgeCnt = sShellWireEdges.GetSize();

        // shell point type (in the future may add option for MultiPoints)
        const bool bPointShell = lFaceuseCnt == 0 && lWireEdgeCnt == 0;

        F.shell_faceuse_count.push_back(lFaceuseCnt);
        F.shell_wire_edge_count.push_back(lWireEdgeCnt);
        F.shell_point_type.push_back(bPointShell ? 1 : 0);

        // for every shell WireEdge - append wireEdge and wireEdge->Curve3d data
        for (uint32_t kk = 0; kk < lWireEdgeCnt; kk++)
        {
            SmEdge *pWireEdge = sShellWireEdges[kk];
            auto found = sEdgeIndex.find(pWireEdge);
            if (found == sEdgeIndex.end())
            {
                return SM_ERR; // SmuConvert dereferences the missing SmEdgeData here
            }
            SmEdgeData *pWireEdgeData = &pBrepData_ii->m_vEdges[static_cast<ULONG>(found->second)];

            // Remember BrepData wireEdge indices, to map BrepData edge indices into compacted
            // edge indices. (SmuConvert records the wire edge's m_lCurve here and later compares
            // it with edge indices; that equals the edge index only while no curve is shared.)
            sWireEdgeIdx.push_back(found->second);

            // WireEdge VertexIndices are compacted after sVertexShellIdx is built.
            sWireEdgeRawVertices.push_back(pWireEdgeData->m_lStartVertex);
            sWireEdgeRawVertices.push_back(pWireEdgeData->m_lEndVertex);

            SmStatus s = WriteEdgeCurve(
                sWireEdgeOut, sEdgeExports[static_cast<size_t>(found->second)], pWireEdge->GetCurve(), pWireEdgeData->m_vInterval, sApproxTol3d
            );
            if (s != SM_SUCCESS)
            {
                return s;
            }
        }

        // when Shell is a Vertex shell - append ShellVertex position
        if (bPointShell)
        {
            if (pShell->GetVertex() == NULL)
            {
                return SM_ERR;
            }
            sVertexShellIdx.push_back(rShellData.m_lVertex);
            SmPoint3d sVertexShellPoint = pShell->GetVertex()->GetPoint();
            Push3(F.shell_point_position, sVertexShellPoint.x, sVertexShellPoint.y, sVertexShellPoint.z);
        }
    }
    std::sort(sVertexShellIdx.begin(), sVertexShellIdx.end());
    std::sort(sWireEdgeIdx.begin(), sWireEdgeIdx.end());

    // compact the WireEdge->VertexIndices now that sVertexShellIdx is built
    for (int64_t v : sWireEdgeRawVertices)
    {
        F.wire_edge_vertices.push_back(smu_CompactLocalIndex(v, sVertexShellIdx));
    }

    // Faceuses
    for (uint32_t jj = 0; jj < lNumFaceuses; ++jj)
    {
        SmFaceuseData &rFaceuseData = pBrepData_ii->m_vFaceuses[jj];
        F.faceuse_face.push_back(rFaceuseData.m_lFace);
        // m_bOrientation uses UNSURE for notSet; faceuse:orientationType has no notSet value, defaults to opposite
        F.faceuse_orientation.push_back(rFaceuseData.m_bOrientation == TRUE ? 1 : 0);
    }

    // Faces and Face->Surfaces
    for (uint32_t jj = 0; jj < lNumFaces; ++jj)
    {
        SmFaceData &rFaceData = pBrepData_ii->m_vFaces[jj];
        SmFace *pFace = rFaceData.m_pFace;
        if (pFace == NULL || pFace->GetSurface() == NULL || rFaceData.m_lSurface >= pBrepData_ii->m_vSurfaces.GetSize())
        {
            return SM_ERR; // SmuConvert writes surface type "none" here and has no surface to read
        }

        // m_bRectangularTrim uses UNSURE for notSet; face:trimType has no notSet value, defaults to general
        const int64_t tFaceTrimType = (rFaceData.m_bRectangularTrim == TRUE) ? 1 : 0;

        // face surfaceType - detect analytic types first, then fall back to NURB
        int64_t tFaceSurfaceType = kSurfaceNurbs;
        if (pFace->GetSurface()->IsKindOf(SmSphere_TYPE))
        {
            tFaceSurfaceType = kSurfaceSphere;
        }
        else if (pFace->GetSurface()->IsKindOf(SmPlane_TYPE))
        {
            tFaceSurfaceType = kSurfacePlane;
        }
        else if (pFace->GetSurface()->IsKindOf(SmCylinder_TYPE))
        {
            tFaceSurfaceType = kSurfaceCylinder;
        }
        else if (pFace->GetSurface()->IsKindOf(SmCone_TYPE))
        {
            tFaceSurfaceType = kSurfaceCone;
        }
        else if (pFace->GetSurface()->IsKindOf(SmTorus_TYPE))
        {
            tFaceSurfaceType = kSurfaceTorus;
        }

        // A curved analytic face is also given with its NURBS form (FACE_NURBS_FORM), the NURBS
        // SMLib works on: its range and UV trims are in that NURBS's parameters, as SMLib holds
        // them. A plane's are in its own (u, v), which carries no trims.
        const bool bNurbsForm = tFaceSurfaceType != kSurfaceNurbs && tFaceSurfaceType != kSurfacePlane;
        // A plane's trims are in its NURBS (u, v); the arrays' plane (u, v) are distances along
        // its unit axes (its STEP parameters scaled by GetUVScale), an affine image of them.
        std::array<double, 6> sPlaneMap{};
        SmBoolean bPlaneMap = FALSE;
        if (tFaceSurfaceType == kSurfacePlane)
        {
            const SmSurface *pPlaneSurface = pFace->GetSurface();
            const SmVector2d sScale = static_cast<const SmPlane *>(pPlaneSurface)->GetUVScale();
            SmPoint2d q[3];
            const SmPoint2d at[3] = {SmPoint2d(0.0, 0.0), SmPoint2d(1.0, 0.0), SmPoint2d(0.0, 1.0)};
            bPlaneMap = TRUE;
            for (int k = 0; k < 3 && bPlaneMap; ++k)
            {
                SmPoint2d st;
                bPlaneMap = pPlaneSurface->ConvertUVFromNURBSToSTEP(at[k], st) == SM_SUCCESS ? TRUE : FALSE;
                q[k] = SmPoint2d(st.x * sScale.x, st.y * sScale.y);
            }
            sPlaneMap = {q[0].x, q[1].x - q[0].x, q[2].x - q[0].x, q[0].y, q[1].y - q[0].y, q[2].y - q[0].y};
        }
        for (uint32_t kk = 0; kk < rFaceData.m_lNumLoops; ++kk)
        {
            const uint32_t lLoopIndex = rFaceData.m_lStartLoop + kk;
            if (lLoopIndex >= sLoopExportsUVCurves.size())
            {
                return SM_ERR; // "face loop index is out of bounds"
            }
            sLoopExportsUVCurves[lLoopIndex] = tFaceSurfaceType == kSurfaceNurbs || bNurbsForm || bPlaneMap;
            sLoopPlaneMap[lLoopIndex] = sPlaneMap;
            sLoopHasPlaneMap[lLoopIndex] = bPlaneMap;
        }

        // face range - must match the UV conventions declared by each surface API's schema:
        //   Plane    : (u, v) linear
        //   Sphere   : (u longitude, v latitude)    in radians
        //   Cylinder : (u angle,     v linear)      u in radians
        //   Cone     : (u angle,     v linear)      u in radians
        //   Torus    : (u major ang, v minor ang)   in radians
        //   Nurb     : native Nurb (u, v)
        // For analytic surfaces we source bounds from the surface's STEP UV domain (SmLib degrees
        // for angles) and convert angular parameters to radians. The Face range need not be
        // tight to the trimmed region. For NURB surfaces we keep the Face's native (Nurb) UV
        // domain.
        const SmExtent2d &sFaceRangeDomain = rFaceData.m_vUVDomain;
        double sFaceUVMin[2] = {sFaceRangeDomain.GetUMin(), sFaceRangeDomain.GetVMin()};
        double sFaceUVMax[2] = {sFaceRangeDomain.GetUMax(), sFaceRangeDomain.GetVMax()};
        if (tFaceSurfaceType != kSurfaceNurbs && !bNurbsForm)
        {
            SmExtent2d sStepDom = pFace->GetSurface()->GetSTEPUVDomain();
            double dUMin = sStepDom.GetUMin();
            double dVMin = sStepDom.GetVMin();
            double dUMax = sStepDom.GetUMax();
            double dVMax = sStepDom.GetVMax();

            if (tFaceSurfaceType == kSurfaceSphere || tFaceSurfaceType == kSurfaceTorus)
            {
                sFaceUVMin[0] = SM_DEG2RAD(dUMin);
                sFaceUVMin[1] = SM_DEG2RAD(dVMin);
                sFaceUVMax[0] = SM_DEG2RAD(dUMax);
                sFaceUVMax[1] = SM_DEG2RAD(dVMax);
            }
            else if (tFaceSurfaceType == kSurfaceCylinder || tFaceSurfaceType == kSurfaceCone)
            {
                sFaceUVMin[0] = SM_DEG2RAD(dUMin);
                sFaceUVMin[1] = dVMin;
                sFaceUVMax[0] = SM_DEG2RAD(dUMax);
                sFaceUVMax[1] = dVMax;
            }
            else // plane
            {
                // Plane parameters are distances along unit placement axes.
                SmVector2d sScale = static_cast<const SmPlane *>(pFace->GetSurface())->GetUVScale();
                sStepDom.Scale(sScale.x, &sScale.y);
                sFaceUVMin[0] = sStepDom.GetUMin();
                sFaceUVMin[1] = sStepDom.GetVMin();
                sFaceUVMax[0] = sStepDom.GetUMax();
                sFaceUVMax[1] = sStepDom.GetVMax();
            }
        }

        F.face_loop_count.push_back(rFaceData.m_lNumLoops);
        F.face_surface_type.push_back(tFaceSurfaceType);
        F.face_trim_type.push_back(tFaceTrimType);
        F.face_nurbs_form.push_back(bNurbsForm ? 1 : 0);
        F.face_range.insert(F.face_range.end(), {sFaceUVMin[0], sFaceUVMin[1], sFaceUVMax[0], sFaceUVMax[1]});

        // face surface
        SmSurface *pSurface = pBrepData_ii->m_vSurfaces[rFaceData.m_lSurface];

        if (tFaceSurfaceType == kSurfaceSphere)
        {
            SmSphere *pSphere = SM_REINTERPRET_CAST(SmSphere *, pSurface);
            SmAxis2Placement sOrigin;
            double dRadius = 0.0;
            SmStatus canonical = pSphere->GetCanonical(sOrigin, dRadius);
            if (canonical != SM_SUCCESS)
            {
                return canonical;
            }

            SmPoint3d sCenter = sOrigin.GetOrigin();
            SmVector3d sZAxis = sOrigin.GetZAxis();
            SmVector3d sXAxis = sOrigin.GetXAxis();

            Push3(F.face_sphere_center, sCenter.x, sCenter.y, sCenter.z);
            Push3(F.face_sphere_axis, sZAxis.x, sZAxis.y, sZAxis.z);
            Push3(F.face_sphere_ref_direction, sXAxis.x, sXAxis.y, sXAxis.z);
            F.face_sphere_radius.push_back(dRadius);
        }
        else if (tFaceSurfaceType == kSurfacePlane)
        {
            SmPlane *pPlane = SM_REINTERPRET_CAST(SmPlane *, pSurface);
            SmAxis2Placement sPlacement;
            SmStatus canonical = pPlane->GetCanonical(sPlacement);
            if (canonical != SM_SUCCESS)
            {
                return canonical;
            }

            SmPoint3d sOriginPt = sPlacement.GetOrigin();
            SmVector3d sZAxis = sPlacement.GetZAxis();
            SmVector3d sXAxis = sPlacement.GetXAxis();

            Push3(F.face_plane_origin, sOriginPt.x, sOriginPt.y, sOriginPt.z);
            Push3(F.face_plane_axis, sZAxis.x, sZAxis.y, sZAxis.z);
            Push3(F.face_plane_ref_direction, sXAxis.x, sXAxis.y, sXAxis.z);
        }
        else if (tFaceSurfaceType == kSurfaceCylinder)
        {
            SmCylinder *pCylinder = SM_REINTERPRET_CAST(SmCylinder *, pSurface);
            SmAxis2Placement sPlacement;
            double dRadius = 0.0;
            SmStatus canonical = pCylinder->GetCanonical(sPlacement, dRadius);
            if (canonical != SM_SUCCESS)
            {
                return canonical;
            }

            SmPoint3d sOriginPt = sPlacement.GetOrigin();
            SmVector3d sZAxis = sPlacement.GetZAxis();
            SmVector3d sXAxis = sPlacement.GetXAxis();

            Push3(F.face_cylinder_origin, sOriginPt.x, sOriginPt.y, sOriginPt.z);
            Push3(F.face_cylinder_axis, sZAxis.x, sZAxis.y, sZAxis.z);
            Push3(F.face_cylinder_ref_direction, sXAxis.x, sXAxis.y, sXAxis.z);
            F.face_cylinder_radius.push_back(dRadius);
        }
        else if (tFaceSurfaceType == kSurfaceCone)
        {
            SmCone *pCone = SM_REINTERPRET_CAST(SmCone *, pSurface);
            SmAxis2Placement sPlacement;
            double dRadius = 0.0;
            double dSemiAngleDeg = 0.0;
            SmStatus canonical = pCone->GetCanonical(sPlacement, dRadius, dSemiAngleDeg);
            if (canonical != SM_SUCCESS)
            {
                return canonical;
            }

            SmPoint3d sOriginPt = sPlacement.GetOrigin();
            SmVector3d sZAxis = sPlacement.GetZAxis();
            SmVector3d sXAxis = sPlacement.GetXAxis();

            Push3(F.face_cone_origin, sOriginPt.x, sOriginPt.y, sOriginPt.z);
            Push3(F.face_cone_axis, sZAxis.x, sZAxis.y, sZAxis.z);
            Push3(F.face_cone_ref_direction, sXAxis.x, sXAxis.y, sXAxis.z);
            F.face_cone_radius.push_back(dRadius);
            F.face_cone_semi_angle.push_back(SM_DEG2RAD(dSemiAngleDeg));
        }
        else if (tFaceSurfaceType == kSurfaceTorus)
        {
            SmTorus *pTorus = SM_REINTERPRET_CAST(SmTorus *, pSurface);
            SmAxis2Placement sPlacement;
            double dMajorRadius = 0.0;
            double dMinorRadius = 0.0;
            SmStatus canonical = pTorus->GetCanonical(sPlacement, dMajorRadius, dMinorRadius);
            if (canonical != SM_SUCCESS)
            {
                return canonical;
            }

            SmPoint3d sOriginPt = sPlacement.GetOrigin();
            SmVector3d sZAxis = sPlacement.GetZAxis();
            SmVector3d sXAxis = sPlacement.GetXAxis();

            Push3(F.face_torus_origin, sOriginPt.x, sOriginPt.y, sOriginPt.z);
            Push3(F.face_torus_axis, sZAxis.x, sZAxis.y, sZAxis.z);
            Push3(F.face_torus_ref_direction, sXAxis.x, sXAxis.y, sXAxis.z);
            F.face_torus_major_radius.push_back(dMajorRadius);
            F.face_torus_minor_radius.push_back(dMinorRadius);
        }
        if (bNurbsForm)
        {
            // Every analytic class is an SmBSplineSurface: its own NURBS, in face order with the
            // NURBS faces' surfaces.
            SmStatus s = smu_AddSmBSplineSurface(*SM_REINTERPRET_CAST(SmBSplineSurface *, pSurface), F);
            if (s != SM_SUCCESS)
            {
                return s;
            }
        }
        else if (tFaceSurfaceType == kSurfaceNurbs)
        {
            // NURB or approximation path
            SmBSplineSurface *pBSplineSurface = NULL;
            SmObjDelete sCleanBSplineSurface;

            if (pFace->GetSurface()->IsKindOf(SmBSplineSurface_TYPE))
            {
                pBSplineSurface = SM_REINTERPRET_CAST(SmBSplineSurface *, pSurface);
            }
            else
            {
                double dAchievedTol = 0.0;
                SmStatus approximated = pSurface->ApproximateSurface(
                    *pSurface->GetContext(), // in : context for new object construction
                    sApproxTol3d,            // in : max allowed approximation distance
                    dAchievedTol,            // out: Max tolerance achieved.
                    5,                       // in : How many times to split spans before giving up
                    pBSplineSurface,         // out: approx surf with the smallest number of subdivisions possible within tol
                    TRUE                     // in : TRUE = if this Surface is a BSpline, copy it FALSE= approximate the surface
                );
                sCleanBSplineSurface.SetObj(pBSplineSurface);
                if (approximated != SM_SUCCESS)
                {
                    return approximated;
                }
            }

            if (pBSplineSurface == NULL)
            {
                // SmuConvert logs "unable to convert face surface to BSpline" and continues,
                // leaving a NURBS face without a surface block; the arrays would be misaligned.
                return SM_ERR;
            }

            SmStatus s = smu_AddSmBSplineSurface(*pBSplineSurface, F);
            if (s != SM_SUCCESS)
            {
                return s; // SmuConvert ignores this status, with the same misalignment
            }
        }
    }

    // Loops
    for (uint32_t jj = 0; jj < lNumLoops; ++jj)
    {
        SmLoopData &rLoopData = pBrepData_ii->m_vLoops[jj];
        uint32_t lLoopType = rLoopData.m_lLoopType;
        // Loop must be either edgeuse_TYPE (lLoopType==0) or vertex_TYPE (lLoopType==1) - else abort
        if (lLoopType == 0 && rLoopData.m_lNumEU > 0)
        {
            F.loop_edgeuse_count.push_back(rLoopData.m_lNumEU);
            F.loop_vertex.push_back(kNoIndex);
        }
        else if (lLoopType == 1 && rLoopData.m_lVertex != SM_NO_OBJECT && rLoopData.m_lVertex != SM_UNDEF_ULONG)
        {
            F.loop_edgeuse_count.push_back(0);
            F.loop_vertex.push_back(smu_CompactLocalIndex(rLoopData.m_lVertex, sVertexShellIdx));
        }
        else
        {
            return SM_ERR; // "Incorrect LoopType"
        }
    }

    // Edgeuses: only edgeuses that connect edges to faces; wireEdge edgeuses are removed.
    for (int64_t jj = 0; jj < static_cast<int64_t>(lNumEdgeuses); ++jj)
    {
        if (pBrepData_ii->m_vEdgeuses[static_cast<ULONG>(jj)].m_lEUType != 0)
        {
            sNonFaceEdgeuseIdx.push_back(jj);
        }
    }
    for (int64_t jj = 0; jj < static_cast<int64_t>(lNumEdgeuses); ++jj)
    {
        SmEUData &rEUData = pBrepData_ii->m_vEdgeuses[static_cast<ULONG>(jj)];

        // ignore edgeuses used for wireEdges, ShellVertices, and VertexLoops
        if (rEUData.m_lEUType != 0) // 0 == Edgeuse connects an Edge to a Face
        {
            continue;
        }

        F.edgeuse_edge.push_back(smu_CompactLocalIndex(rEUData.m_lEdge, sWireEdgeIdx));
        // m_bOrientation uses UNSURE for notSet; edgeuse:orientationType has no notSet value, defaults to opposite
        F.edgeuse_orientation.push_back(rEUData.m_bOrientation == TRUE ? 1 : 0);

        // radial edgeuse list:
        //  SMLib: each SmEUData obj represents a mated top/bot pair of SmEdgeuse objects. The
        //         radial edgeuse order is stored as a pair of Next Edgeuse indices
        //           SmEUData::m_lNextEU     = the topEdgeuse's next edgeuse in the radial list.
        //           SmEUData::m_lMateNextEU = the botEdgeuse's next edgeuse in the radial list.
        //         one of the top/bot pair points to its mate and the other to a bottom or top
        //         edgeuse member of another SmEUData object.
        //           a nextIndex value:[positive ii] = the topEdgeuse of the (ii)th    SmEUData
        //           a nextIndex value:[negative ii] = the botEdgeuse of the (-ii-1)th SmEUData
        //  flat: each edgeuse record also represents a mated top/bot pair.
        //         EDGEUSE_NEXT_RADIAL  = the next edgeuse record on the radial list
        //         EDGEUSE_RADIAL_ENTRY = the radial order of this pair:
        //              1 (topEntry)    = { topEdgeuse, botEdgeuse }
        //              0 (bottomEntry) = { botEdgeuse, topEdgeuse }
        const int64_t lNextEU = static_cast<int64_t>(rEUData.m_lNextEU);
        const int64_t lMateNextEU = static_cast<int64_t>(rEUData.m_lMateNextEU);
        const int64_t lTopEUNext = (lNextEU >= 0) ? lNextEU : (-lNextEU - 1);
        const int64_t lBotEUNext = (lMateNextEU >= 0) ? lMateNextEU : (-lMateNextEU - 1);

        // when thisEdgeuse->topEdgeuseNext points to its own mated botEdgeuse, the entry is topEntry.
        const bool bTopEntry = lTopEUNext == jj;
        F.edgeuse_radial_entry.push_back(bTopEntry ? 1 : 0);

        // topEntry: next = SMLib's BotEUNext; bottomEntry: next = SMLib's TopEUNext. Compacted
        // as an edgeuse index once every removed edgeuse is known (it is, above).
        F.edgeuse_next_radial.push_back(smu_CompactLocalIndex(bTopEntry ? lBotEUNext : lTopEUNext, sNonFaceEdgeuseIdx));

        // Optional UV trim curve. Analytic-face curves are represented as missing because their
        // coordinates use SMLib's internal NURBS parameterization, not the analytic UV.
        if (rEUData.m_lLoop >= sLoopExportsUVCurves.size())
        {
            return SM_ERR; // "edgeuse loop index is out of bounds"
        }
        const SmBoolean bExportThisUVCurve = sLoopExportsUVCurves[rEUData.m_lLoop];
        SmCurve *pUVCurve = bExportThisUVCurve && rEUData.m_lUVCurve != SM_NO_OBJECT && rEUData.m_lUVCurve < pBrepData_ii->m_vUVCurves.GetSize()
                                ? pBrepData_ii->m_vUVCurves[rEUData.m_lUVCurve]
                                : NULL;
        SmBSplineCurve *pBSplineCurve2d = (pUVCurve != NULL && pUVCurve->IsKindOf(SmBSplineCurve_TYPE)) ? SM_REINTERPRET_CAST(SmBSplineCurve *, pUVCurve) : NULL;
        // SolidCraft: on an edge exported as an analytic curve, the trim's knots follow the
        // imported curve's parameter, or the trim is dropped when that is not an affine image of
        // the edge's (SMLib then rebuilds it from the 3D edge, as for analytic faces).
        double sRemap[4];
        const double *pRemap = nullptr;
        if (pBSplineCurve2d != NULL && rEUData.m_lEdge < sEdgeExports.size() && sEdgeExports[rEUData.m_lEdge].type != kCurveNurbs)
        {
            const EdgeExport &rE = sEdgeExports[rEUData.m_lEdge];
            const SmExtent1d &rIvl = pBrepData_ii->m_vEdges[rEUData.m_lEdge].m_vInterval;
            if (rE.affine)
            {
                sRemap[0] = rIvl.GetMin();
                sRemap[1] = rIvl.GetMax();
                sRemap[2] = rE.smLo;
                sRemap[3] = rE.smHi;
                pRemap = sRemap;
            }
            else
            {
                pBSplineCurve2d = NULL;
            }
        }
        const size_t lCountBefore = F.edgeuse_nurbs_vertex_count.size();
        const SmStatus sUvCurveStatus = smu_AddSmBSplineCurve2d(pBSplineCurve2d, rEUData.m_bOrientation, F, pRemap, sLoopHasPlaneMap[rEUData.m_lLoop] ? sLoopPlaneMap[rEUData.m_lLoop].data() : nullptr);
        if (sUvCurveStatus != SM_SUCCESS && F.edgeuse_nurbs_vertex_count.size() == lCountBefore)
        {
            // UV trims are optional. Preserve the edgeuse/UV-curve one-to-one record alignment
            // and let the reader reconstruct this trim from 3D topology.
            F.edgeuse_nurbs_vertex_count.push_back(0);
            F.edgeuse_nurbs_order.push_back(0);
        }
    }
    // Edges - only FaceEdges here, WireEdges are done in the Shell block
    for (uint32_t jj = 0; jj < lNumEdges; ++jj)
    {
        SmEdgeData &rEdgeData = pBrepData_ii->m_vEdges[jj];
        SmEdge *pEdge = rEdgeData.m_pEdge;

        if (std::binary_search(sWireEdgeIdx.begin(), sWireEdgeIdx.end(), static_cast<int64_t>(jj)))
        {
            continue;
        }
        if (pEdge == NULL)
        {
            return SM_ERR;
        }

        F.edge_vertices.push_back(smu_CompactLocalIndex(rEdgeData.m_lStartVertex, sVertexShellIdx));
        F.edge_vertices.push_back(smu_CompactLocalIndex(rEdgeData.m_lEndVertex, sVertexShellIdx));
        F.edge_tolerance.push_back(SM_ZONE_TO_XSECTTOL3D(ZoneTolWithFaces(pEdge, dBrepZoneTol)));

        SmStatus s = WriteEdgeCurve(sEdgeOut, sEdgeExports[jj], pEdge->GetCurve(), rEdgeData.m_vInterval, sApproxTol3d);
        if (s != SM_SUCCESS)
        {
            return s;
        }
    }

    // Vertices - only Vertices here, ShellVertices are done in the Shell block
    for (uint32_t jj = 0; jj < lNumVertices; ++jj)
    {
        if (std::binary_search(sVertexShellIdx.begin(), sVertexShellIdx.end(), static_cast<int64_t>(jj)))
        {
            continue;
        }
        const SmPoint3d &crPt = pBrepData_ii->m_vVertices[jj].m_vPoint;
        F.vertex_point_type.push_back(1); // BrepPointAPI
        Push3(F.vertex_position, crPt.x, crPt.y, crPt.z);
        const SmVertex *pVertex = pBrepData_ii->m_vVertices[jj].m_pVertex;
        double dVertexTol = sXSectTol3d;
        if (pVertex != NULL)
        {
            dVertexTol = SM_ZONE_TO_XSECTTOL3D(ZoneTolWithFaces(pVertex, dBrepZoneTol));
        }
        F.vertex_tolerance.push_back(dVertexTol);
    }

    return SM_SUCCESS;
}

/* ---- Import ---------------------------------------------------------------------------------- */

namespace
{

// What Validate learns about valid input.
struct Shape
{
    size_t regions = 0, shells = 0, faceuses = 0, faces = 0, loops = 0, edgeuses = 0, edges = 0, wire_edges = 0, vertices = 0;
    size_t point_shells = 0;
    bool uv_curves = false;
};

bool Bits(const std::vector<int64_t> &v)
{
    return std::all_of(v.begin(), v.end(), [](int64_t x) { return x == 0 || x == 1; });
}

bool InRange(const std::vector<int64_t> &v, int64_t lo, int64_t hi)
{
    return std::all_of(v.begin(), v.end(), [=](int64_t x) { return x >= lo && x <= hi; });
}

// Sum of non-negative counts, or -1 when one is negative or the sum exceeds kMaxTopology.
int64_t SumCounts(const std::vector<int64_t> &v)
{
    int64_t s = 0;
    for (int64_t x : v)
    {
        if (x < 0 || x > kMaxTopology)
        {
            return -1;
        }
        s += x;
        if (s > kMaxTopology)
        {
            return -1;
        }
    }
    return s;
}

double Length(const SmVector3d &v)
{
    return std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

// An axis and a reference direction that span a plane, as SmAxis2Placement::SetCanonical needs.
bool GoodFrame(const std::vector<double> &axis, const std::vector<double> &ref, size_t i)
{
    SmVector3d a = GetVec3(axis, i), r = GetVec3(ref, i);
    double la = Length(a), lr = Length(r);
    if (!(la > 1e-12) || !(lr > 1e-12))
    {
        return false;
    }
    SmVector3d c = a * r;
    return Length(c) > 1e-9 * la * lr;
}

// A flat knot vector of `count` points and order `order`: non-decreasing with a non-empty
// interval [k[order-1], k[count]].
bool GoodKnots(const double *k, int64_t count, int64_t order)
{
    for (int64_t i = 1; i < count + order; ++i)
    {
        if (!(k[i] >= k[i - 1]))
        {
            return false;
        }
    }
    return k[order - 1] < k[count];
}

bool GoodWeights(const std::vector<double> &w)
{
    return std::all_of(w.begin(), w.end(), [](double x) { return x > 0.0; });
}

// One family of NURBS curves (edges, wire edges or UV curves): per-curve counts and orders,
// packed points of `dim` coordinates, knots and weights. Skips curves with count == 0 when
// `optional` (UV curves), where count == 0 must come with order == 0.
bool GoodNurbsCurves(
    const std::vector<int64_t> &count,
    const std::vector<int64_t> &order,
    const std::vector<double> &cvs,
    const std::vector<double> &knots,
    const std::vector<double> &weights,
    int dim,
    bool optional
)
{
    if (count.size() != order.size())
    {
        return false;
    }
    int64_t points = 0, knotTotal = 0;
    for (size_t i = 0; i < count.size(); ++i)
    {
        if (optional && count[i] == 0 && order[i] == 0)
        {
            continue;
        }
        if (order[i] < 2 || order[i] > kMaxOrder || count[i] < order[i] || count[i] > SMLIB_FLAT_MAX_LEN)
        {
            return false;
        }
        points += count[i];
        knotTotal += count[i] + order[i];
        if (points > SMLIB_FLAT_MAX_LEN || knotTotal > SMLIB_FLAT_MAX_LEN)
        {
            return false;
        }
    }
    if (static_cast<int64_t>(cvs.size()) != dim * points || static_cast<int64_t>(knots.size()) != knotTotal)
    {
        return false;
    }
    if (!weights.empty() && (static_cast<int64_t>(weights.size()) != points || !GoodWeights(weights)))
    {
        return false;
    }
    int64_t k = 0;
    for (size_t i = 0; i < count.size(); ++i)
    {
        if (optional && count[i] == 0)
        {
            continue;
        }
        if (!GoodKnots(knots.data() + k, count[i], order[i]))
        {
            return false;
        }
        k += count[i] + order[i];
    }
    return true;
}

// Ranges of one edge family against its curves.
bool GoodEdges(
    const std::vector<int64_t> &type,
    const std::vector<int64_t> &vertices,
    const std::vector<double> &range,
    size_t n,
    size_t nVertices,
    const CurveArrays &nurbs,
    const std::vector<double> &circleCenter,
    const std::vector<double> &circleAxis,
    const std::vector<double> &circleRef,
    const std::vector<double> &circleRadius,
    const std::vector<double> &lineOrigin,
    const std::vector<double> &lineDirection,
    const std::vector<double> &ellipseCenter,
    const std::vector<double> &ellipseAxis,
    const std::vector<double> &ellipseRef,
    const std::vector<double> &ellipseX,
    const std::vector<double> &ellipseY
)
{
    if (type.size() != n || vertices.size() != 2 * n || range.size() != 2 * n || !InRange(type, kCurveNurbs, kCurveEllipse) ||
        !InRange(vertices, 0, static_cast<int64_t>(nVertices) - 1))
    {
        return false;
    }
    size_t nN = 0, nC = 0, nL = 0, nE = 0;
    for (int64_t t : type)
    {
        nN += t == kCurveNurbs;
        nC += t == kCurveCircle;
        nL += t == kCurveLine;
        nE += t == kCurveEllipse;
    }
    if (nurbs.count.size() != nN || !GoodNurbsCurves(nurbs.count, nurbs.order, nurbs.cvs, nurbs.knots, nurbs.weights, 3, false))
    {
        return false;
    }
    if (circleCenter.size() != 3 * nC || circleAxis.size() != 3 * nC || circleRef.size() != 3 * nC || circleRadius.size() != nC ||
        lineOrigin.size() != 3 * nL || lineDirection.size() != 3 * nL || ellipseCenter.size() != 3 * nE || ellipseAxis.size() != 3 * nE ||
        ellipseRef.size() != 3 * nE || ellipseX.size() != nE || ellipseY.size() != nE)
    {
        return false;
    }
    for (size_t i = 0; i < nC; ++i)
    {
        if (!(circleRadius[i] > 0.0) || !GoodFrame(circleAxis, circleRef, i))
        {
            return false;
        }
    }
    for (size_t i = 0; i < nL; ++i)
    {
        if (!(Length(GetVec3(lineDirection, i)) > 1e-12))
        {
            return false;
        }
    }
    for (size_t i = 0; i < nE; ++i)
    {
        if (!(ellipseX[i] > 0.0) || !(ellipseY[i] > 0.0) || !GoodFrame(ellipseAxis, ellipseRef, i))
        {
            return false;
        }
    }
    // Each edge's range: non-empty; inside its NURBS curve's knot interval; at most one turn
    // for circles and ellipses (radians).
    size_t iNurbs = 0;
    int64_t knot = 0;
    for (size_t i = 0; i < n; ++i)
    {
        double lo = range[2 * i], hi = range[2 * i + 1];
        if (!(lo < hi))
        {
            return false;
        }
        if (type[i] == kCurveNurbs)
        {
            int64_t c = nurbs.count[iNurbs], o = nurbs.order[iNurbs];
            double k0 = nurbs.knots[knot + o - 1], k1 = nurbs.knots[knot + c];
            double tol = 1e-9 * std::max(1.0, k1 - k0);
            if (lo < k0 - tol || hi > k1 + tol)
            {
                return false;
            }
            knot += c + o;
            iNurbs++;
        }
        else if ((type[i] == kCurveCircle || type[i] == kCurveEllipse) &&
                 (hi - lo > 2.0 * SM_PI + 1e-9 || std::fabs(lo) > 2.0 * SM_PI * (1.0 + 1e-9) || std::fabs(hi) > 2.0 * SM_PI * (1.0 + 1e-9)))
        {
            return false; // SmCircle/SmEllipse::CreateCanonical take [-360, 360] degrees
        }
        else if (type[i] == kCurveLine && (std::fabs(lo) > SM_INFINITE_PARAMETER || std::fabs(hi) > SM_INFINITE_PARAMETER))
        {
            return false; // an imported line is unbounded over +-SM_INFINITE_PARAMETER
        }
    }
    return true;
}

// Every check BrepFromFlat relies on. Returns false on anything malformed.
bool Validate(const FlatBrep &F, Shape &S)
{
    // Every double is finite and of sane magnitude.
    for (const Field &f : Fields())
    {
        if (f.f != nullptr)
        {
            for (double x : F.*(f.f))
            {
                if (!std::isfinite(x) || std::fabs(x) > kMaxValue)
                {
                    return false;
                }
            }
        }
    }

    // B-rep data
    if (F.tolerance.size() != 1 || !(F.tolerance[0] > 0.0) || F.tolerance[0] > 1.0)
    {
        return false;
    }
    if (!F.extent.empty() && F.extent.size() != 6)
    {
        return false;
    }
    if (F.source.size() > 1 || !Bits(F.source))
    {
        return false;
    }

    // Regions
    S.regions = F.region_shell_count.size();
    int64_t shells = SumCounts(F.region_shell_count);
    if (S.regions < 1 || S.regions > static_cast<size_t>(kMaxTopology) || F.region_type.size() != S.regions || !Bits(F.region_type) || shells < 0)
    {
        return false;
    }
    S.shells = static_cast<size_t>(shells);

    // Shells
    if (F.shell_faceuse_count.size() != S.shells || F.shell_wire_edge_count.size() != S.shells || F.shell_point_type.size() != S.shells ||
        !Bits(F.shell_wire_edge_count) || !Bits(F.shell_point_type))
    {
        return false;
    }
    int64_t faceuses = SumCounts(F.shell_faceuse_count);
    if (faceuses < 0 || faceuses % 2 != 0)
    {
        return false;
    }
    for (size_t i = 0; i < S.shells; ++i)
    {
        // a faceuse shell, a wire-edge shell (one edge), or a point shell
        if (F.shell_faceuse_count[i] > 0)
        {
            if (F.shell_wire_edge_count[i] != 0)
            {
                return false;
            }
        }
        else if (F.shell_wire_edge_count[i] == 0)
        {
            if (F.shell_point_type[i] != 1)
            {
                return false;
            }
            S.point_shells++;
        }
    }
    S.wire_edges = static_cast<size_t>(SumCounts(F.shell_wire_edge_count));
    if (F.shell_point_position.size() != 3 * S.point_shells)
    {
        return false;
    }

    // Faceuses: each face has exactly two, one on each side.
    S.faceuses = static_cast<size_t>(faceuses);
    S.faces = S.faceuses / 2;
    if (F.faceuse_face.size() != S.faceuses || F.faceuse_orientation.size() != S.faceuses || !Bits(F.faceuse_orientation) ||
        !InRange(F.faceuse_face, 0, static_cast<int64_t>(S.faces) - 1))
    {
        return false;
    }
    {
        std::vector<uint8_t> sides(S.faces, 0);
        for (size_t i = 0; i < S.faceuses; ++i)
        {
            uint8_t bit = F.faceuse_orientation[i] ? 1 : 2;
            uint8_t &s = sides[static_cast<size_t>(F.faceuse_face[i])];
            if (s & bit)
            {
                return false;
            }
            s |= bit;
        }
    }

    // Faces
    if (F.face_loop_count.size() != S.faces || F.face_surface_type.size() != S.faces || F.face_trim_type.size() != S.faces ||
        F.face_range.size() != 4 * S.faces || !Bits(F.face_trim_type) || !InRange(F.face_surface_type, kSurfaceNurbs, kSurfaceTorus) ||
        !InRange(F.face_loop_count, 1, kMaxTopology))
    {
        return false;
    }
    if (!F.face_nurbs_form.empty())
    {
        if (F.face_nurbs_form.size() != S.faces || !Bits(F.face_nurbs_form))
        {
            return false;
        }
        for (size_t i = 0; i < S.faces; ++i)
        {
            if (F.face_nurbs_form[i] == 1 && (F.face_surface_type[i] == kSurfaceNurbs || F.face_surface_type[i] == kSurfacePlane))
            {
                return false;
            }
        }
    }
    for (size_t i = 0; i < S.faces; ++i)
    {
        const double u0 = F.face_range[4 * i], v0 = F.face_range[4 * i + 1], u1 = F.face_range[4 * i + 2], v1 = F.face_range[4 * i + 3];
        if (!(u0 < u1) || !(v0 < v1))
        {
            return false;
        }
        // Angular parameters (radians): at most one turn, near the principal range; a sphere's
        // latitude within [-pi/2, pi/2]. A face given with its NURBS form has a NURBS range.
        const int64_t t = FaceHasNurbs(F, i) ? kSurfaceNurbs : F.face_surface_type[i];
        const bool uAngle = t == kSurfaceSphere || t == kSurfaceTorus || t == kSurfaceCylinder || t == kSurfaceCone;
        const bool vAngle = t == kSurfaceSphere || t == kSurfaceTorus;
        const double turn = 2.0 * SM_PI * (1.0 + 1e-9), far = 4.0 * SM_PI;
        if ((uAngle && (u1 - u0 > turn || std::fabs(u0) > far || std::fabs(u1) > far)) ||
            (vAngle && (v1 - v0 > turn || std::fabs(v0) > far || std::fabs(v1) > far)) ||
            (t == kSurfaceSphere && (v0 < -0.5 * SM_PI - 1e-9 || v1 > 0.5 * SM_PI + 1e-9)))
        {
            return false;
        }
    }
    int64_t loops = SumCounts(F.face_loop_count);
    if (loops < 0)
    {
        return false;
    }
    S.loops = static_cast<size_t>(loops);

    // Vertices
    S.vertices = F.vertex_point_type.size();
    if (S.vertices > static_cast<size_t>(kMaxTopology) || S.vertices + S.point_shells > static_cast<size_t>(kMaxTopology) ||
        !InRange(F.vertex_point_type, 1, 1) || F.vertex_position.size() != 3 * S.vertices)
    {
        return false;
    }
    // Per-edge and per-vertex tolerances: none, or one each, as TOLERANCE is.
    auto tolerances_ok = [](const std::vector<double> &t, size_t n)
    {
        return t.empty() || (t.size() == n && std::all_of(t.begin(), t.end(), [](double x) { return std::isfinite(x) && x > 0.0 && x <= 1.0; }));
    };
    if (!tolerances_ok(F.vertex_tolerance, S.vertices))
    {
        return false;
    }

    // Loops
    if (F.loop_edgeuse_count.size() != S.loops || F.loop_vertex.size() != S.loops)
    {
        return false;
    }
    int64_t edgeuses = SumCounts(F.loop_edgeuse_count);
    if (edgeuses < 0)
    {
        return false;
    }
    S.edgeuses = static_cast<size_t>(edgeuses);
    for (size_t i = 0; i < S.loops; ++i)
    {
        if (F.loop_edgeuse_count[i] == 0 ? !(F.loop_vertex[i] >= 0 && F.loop_vertex[i] < static_cast<int64_t>(S.vertices))
                                         : F.loop_vertex[i] != kNoIndex)
        {
            return false;
        }
    }

    // Edges and wire edges
    S.edges = F.edge_curve_type.size();
    if (S.edges > static_cast<size_t>(kMaxTopology) || S.edges + S.wire_edges > static_cast<size_t>(kMaxTopology) || S.wire_edges + S.edgeuses > static_cast<size_t>(kMaxTopology))
    {
        return false;
    }
    if (!tolerances_ok(F.edge_tolerance, S.edges))
    {
        return false;
    }
    if (!GoodEdges(
            F.edge_curve_type, F.edge_vertices, F.edge_range, S.edges, S.vertices,
            CurveArrays{F.edge_nurbs_vertex_count, F.edge_nurbs_order, F.edge_nurbs_cvs, F.edge_nurbs_knots, F.edge_nurbs_weights}, F.edge_circle_center,
            F.edge_circle_axis, F.edge_circle_ref_direction, F.edge_circle_radius, F.edge_line_origin, F.edge_line_direction, F.edge_ellipse_center,
            F.edge_ellipse_axis, F.edge_ellipse_ref_direction, F.edge_ellipse_x_radius, F.edge_ellipse_y_radius
        ) ||
        !GoodEdges(
            F.wire_edge_curve_type, F.wire_edge_vertices, F.wire_edge_range, S.wire_edges, S.vertices,
            CurveArrays{F.wire_edge_nurbs_vertex_count, F.wire_edge_nurbs_order, F.wire_edge_nurbs_cvs, F.wire_edge_nurbs_knots, F.wire_edge_nurbs_weights},
            F.wire_edge_circle_center, F.wire_edge_circle_axis, F.wire_edge_circle_ref_direction, F.wire_edge_circle_radius, F.wire_edge_line_origin,
            F.wire_edge_line_direction, F.wire_edge_ellipse_center, F.wire_edge_ellipse_axis, F.wire_edge_ellipse_ref_direction,
            F.wire_edge_ellipse_x_radius, F.wire_edge_ellipse_y_radius
        ))
    {
        return false;
    }

    // Edgeuses
    if (F.edgeuse_edge.size() != S.edgeuses || F.edgeuse_orientation.size() != S.edgeuses || F.edgeuse_next_radial.size() != S.edgeuses ||
        F.edgeuse_radial_entry.size() != S.edgeuses || !Bits(F.edgeuse_orientation) || !Bits(F.edgeuse_radial_entry) ||
        !InRange(F.edgeuse_edge, 0, static_cast<int64_t>(S.edges) - 1) || !InRange(F.edgeuse_next_radial, 0, static_cast<int64_t>(S.edgeuses) - 1))
    {
        return false;
    }
    {
        // The radial lists: NEXT_RADIAL is a permutation that stays on each edge and cycles
        // once through all of each edge's edgeuses; every edge has at least one edgeuse.
        std::vector<int64_t> perEdge(S.edges, 0);
        std::vector<uint8_t> hit(S.edgeuses, 0);
        std::vector<int64_t> first(S.edges, -1);
        for (size_t i = 0; i < S.edgeuses; ++i)
        {
            size_t e = static_cast<size_t>(F.edgeuse_edge[i]);
            size_t n = static_cast<size_t>(F.edgeuse_next_radial[i]);
            if (F.edgeuse_edge[n] != F.edgeuse_edge[i] || hit[n])
            {
                return false;
            }
            hit[n] = 1;
            perEdge[e]++;
            if (first[e] < 0)
            {
                first[e] = static_cast<int64_t>(i);
            }
        }
        for (size_t e = 0; e < S.edges; ++e)
        {
            if (first[e] < 0)
            {
                return false;
            }
            int64_t steps = 0, cur = first[e];
            do
            {
                cur = F.edgeuse_next_radial[static_cast<size_t>(cur)];
                steps++;
            } while (cur != first[e] && steps <= perEdge[e]);
            if (steps != perEdge[e])
            {
                return false;
            }
        }
    }

    // Every vertex is used; every edgeuse loop is a closed chain of its edges' vertices.
    {
        std::vector<uint8_t> used(S.vertices, 0);
        for (int64_t v : F.edge_vertices)
        {
            used[static_cast<size_t>(v)] = 1;
        }
        for (int64_t v : F.wire_edge_vertices)
        {
            used[static_cast<size_t>(v)] = 1;
        }
        size_t eu = 0;
        for (size_t l = 0; l < S.loops; ++l)
        {
            int64_t n = F.loop_edgeuse_count[l];
            if (n == 0)
            {
                used[static_cast<size_t>(F.loop_vertex[l])] = 1;
                continue;
            }
            for (int64_t k = 0; k < n; ++k)
            {
                size_t a = eu + static_cast<size_t>(k), b = eu + static_cast<size_t>((k + 1) % n);
                size_t ea = static_cast<size_t>(F.edgeuse_edge[a]), eb = static_cast<size_t>(F.edgeuse_edge[b]);
                int64_t endA = F.edge_vertices[2 * ea + (F.edgeuse_orientation[a] ? 1 : 0)];
                int64_t startB = F.edge_vertices[2 * eb + (F.edgeuse_orientation[b] ? 0 : 1)];
                if (endA != startB)
                {
                    return false;
                }
            }
            eu += static_cast<size_t>(n);
        }
        if (!std::all_of(used.begin(), used.end(), [](uint8_t u) { return u != 0; }))
        {
            return false;
        }
    }

    // UV trim curves: none, or one record per edgeuse.
    if (F.edgeuse_nurbs_vertex_count.size() != F.edgeuse_nurbs_order.size() ||
        (!F.edgeuse_nurbs_vertex_count.empty() && F.edgeuse_nurbs_vertex_count.size() != S.edgeuses))
    {
        return false;
    }
    S.uv_curves = !F.edgeuse_nurbs_vertex_count.empty();
    for (size_t i = 0; i < F.edgeuse_nurbs_vertex_count.size(); ++i)
    {
        if ((F.edgeuse_nurbs_vertex_count[i] == 0) != (F.edgeuse_nurbs_order[i] == 0))
        {
            return false;
        }
    }
    if (!GoodNurbsCurves(F.edgeuse_nurbs_vertex_count, F.edgeuse_nurbs_order, F.edgeuse_nurbs_cvs, F.edgeuse_nurbs_knots, F.edgeuse_nurbs_weights, 2, true))
    {
        return false;
    }

    // Surfaces
    size_t nT[6] = {0, 0, 0, 0, 0, 0};
    for (int64_t t : F.face_surface_type)
    {
        nT[static_cast<size_t>(t)]++;
    }
    size_t nN = 0;
    for (size_t i = 0; i < S.faces; ++i)
    {
        nN += FaceHasNurbs(F, i) ? 1 : 0;
    }
    const size_t nP = nT[kSurfacePlane], nCy = nT[kSurfaceCylinder], nCo = nT[kSurfaceCone], nS = nT[kSurfaceSphere],
                 nTo = nT[kSurfaceTorus];
    if (F.face_nurbs_u_vertex_count.size() != nN || F.face_nurbs_v_vertex_count.size() != nN || F.face_nurbs_u_order.size() != nN ||
        F.face_nurbs_v_order.size() != nN)
    {
        return false;
    }
    {
        int64_t points = 0, uk = 0, vk = 0;
        for (size_t i = 0; i < nN; ++i)
        {
            int64_t cu = F.face_nurbs_u_vertex_count[i], cv = F.face_nurbs_v_vertex_count[i];
            int64_t ou = F.face_nurbs_u_order[i], ov = F.face_nurbs_v_order[i];
            if (ou < 2 || ou > kMaxOrder || ov < 2 || ov > kMaxOrder || cu < ou || cv < ov || cu > SMLIB_FLAT_MAX_LEN || cv > SMLIB_FLAT_MAX_LEN ||
                cu * cv > SMLIB_FLAT_MAX_LEN)
            {
                return false;
            }
            points += cu * cv;
            uk += cu + ou;
            vk += cv + ov;
            if (points > SMLIB_FLAT_MAX_LEN || uk > SMLIB_FLAT_MAX_LEN || vk > SMLIB_FLAT_MAX_LEN)
            {
                return false;
            }
        }
        if (static_cast<int64_t>(F.face_nurbs_cvs.size()) != 3 * points || static_cast<int64_t>(F.face_nurbs_u_knots.size()) != uk ||
            static_cast<int64_t>(F.face_nurbs_v_knots.size()) != vk)
        {
            return false;
        }
        if (!F.face_nurbs_weights.empty() && (static_cast<int64_t>(F.face_nurbs_weights.size()) != points || !GoodWeights(F.face_nurbs_weights)))
        {
            return false;
        }
        uk = vk = 0;
        for (size_t i = 0; i < nN; ++i)
        {
            int64_t cu = F.face_nurbs_u_vertex_count[i], cv = F.face_nurbs_v_vertex_count[i];
            int64_t ou = F.face_nurbs_u_order[i], ov = F.face_nurbs_v_order[i];
            if (!GoodKnots(F.face_nurbs_u_knots.data() + uk, cu, ou) || !GoodKnots(F.face_nurbs_v_knots.data() + vk, cv, ov))
            {
                return false;
            }
            uk += cu + ou;
            vk += cv + ov;
        }
        // Each NURBS face's range overlaps its surface's knot domain (MakeTopologyFromData
        // intersects the two).
        size_t iNurbs = 0;
        uk = vk = 0;
        for (size_t i = 0; i < S.faces; ++i)
        {
            if (!FaceHasNurbs(F, i))
            {
                continue;
            }
            int64_t cu = F.face_nurbs_u_vertex_count[iNurbs], cv = F.face_nurbs_v_vertex_count[iNurbs];
            int64_t ou = F.face_nurbs_u_order[iNurbs], ov = F.face_nurbs_v_order[iNurbs];
            const double *ku = F.face_nurbs_u_knots.data() + uk, *kv = F.face_nurbs_v_knots.data() + vk;
            if (!(std::min(F.face_range[4 * i + 2], ku[cu]) > std::max(F.face_range[4 * i], ku[ou - 1])) ||
                !(std::min(F.face_range[4 * i + 3], kv[cv]) > std::max(F.face_range[4 * i + 1], kv[ov - 1])))
            {
                return false;
            }
            uk += cu + ou;
            vk += cv + ov;
            iNurbs++;
        }
    }
    if (F.face_plane_origin.size() != 3 * nP || F.face_plane_axis.size() != 3 * nP || F.face_plane_ref_direction.size() != 3 * nP ||
        F.face_cylinder_origin.size() != 3 * nCy || F.face_cylinder_axis.size() != 3 * nCy || F.face_cylinder_ref_direction.size() != 3 * nCy ||
        F.face_cylinder_radius.size() != nCy || F.face_cone_origin.size() != 3 * nCo || F.face_cone_axis.size() != 3 * nCo ||
        F.face_cone_ref_direction.size() != 3 * nCo || F.face_cone_radius.size() != nCo || F.face_cone_semi_angle.size() != nCo ||
        F.face_sphere_center.size() != 3 * nS || F.face_sphere_axis.size() != 3 * nS || F.face_sphere_ref_direction.size() != 3 * nS ||
        F.face_sphere_radius.size() != nS || F.face_torus_origin.size() != 3 * nTo || F.face_torus_axis.size() != 3 * nTo ||
        F.face_torus_ref_direction.size() != 3 * nTo || F.face_torus_major_radius.size() != nTo || F.face_torus_minor_radius.size() != nTo)
    {
        return false;
    }
    for (size_t i = 0; i < nP; ++i)
    {
        if (!GoodFrame(F.face_plane_axis, F.face_plane_ref_direction, i))
        {
            return false;
        }
    }
    for (size_t i = 0; i < nCy; ++i)
    {
        if (!(F.face_cylinder_radius[i] > 0.0) || !GoodFrame(F.face_cylinder_axis, F.face_cylinder_ref_direction, i))
        {
            return false;
        }
    }
    for (size_t i = 0; i < nCo; ++i)
    {
        if (!(F.face_cone_radius[i] >= 0.0) || !(std::fabs(F.face_cone_semi_angle[i]) < 0.5 * SM_PI * (1.0 - 1e-9)) ||
            !GoodFrame(F.face_cone_axis, F.face_cone_ref_direction, i))
        {
            return false;
        }
    }
    for (size_t i = 0; i < nS; ++i)
    {
        if (!(F.face_sphere_radius[i] > 0.0) || !GoodFrame(F.face_sphere_axis, F.face_sphere_ref_direction, i))
        {
            return false;
        }
    }
    for (size_t i = 0; i < nTo; ++i)
    {
        if (!(F.face_torus_major_radius[i] > 0.0) || !(F.face_torus_minor_radius[i] > 0.0) || !GoodFrame(F.face_torus_axis, F.face_torus_ref_direction, i))
        {
            return false;
        }
    }
    return true;
}

// Whether two B-spline surfaces are the same NURBS: equal sizes and orders, knots equal to within
// rounding, control points within dTol and weights equal to within rounding.
bool SameNurbs(const SmBSplineSurface &crA, const SmBSplineSurface &crB, double dTol)
{
    FlatBrep A, B;
    if (smu_AddSmBSplineSurface(crA, A) != SM_SUCCESS || smu_AddSmBSplineSurface(crB, B) != SM_SUCCESS)
    {
        return false;
    }
    if (A.face_nurbs_u_vertex_count != B.face_nurbs_u_vertex_count || A.face_nurbs_v_vertex_count != B.face_nurbs_v_vertex_count ||
        A.face_nurbs_u_order != B.face_nurbs_u_order || A.face_nurbs_v_order != B.face_nurbs_v_order)
    {
        return false;
    }
    const auto fClose = [](const std::vector<double> &x, const std::vector<double> &y, double t)
    {
        if (x.size() != y.size())
        {
            return false;
        }
        for (size_t i = 0; i < x.size(); ++i)
        {
            if (!(std::fabs(x[i] - y[i]) <= t * smos_Max(1.0, std::fabs(x[i]))))
            {
                return false;
            }
        }
        return true;
    };
    return fClose(A.face_nurbs_u_knots, B.face_nurbs_u_knots, 1e-12) && fClose(A.face_nurbs_v_knots, B.face_nurbs_v_knots, 1e-12) &&
           fClose(A.face_nurbs_weights, B.face_nurbs_weights, 1e-12) && fClose(A.face_nurbs_cvs, B.face_nurbs_cvs, dTol);
}

// The analytic surface of kind tKind (cylinder, cone, sphere or torus) that SMLib's own test finds
// in crNurbs, keeping crNurbs as its NURBS form; NULL when SMLib finds none, or makes one with
// another NURBS (SmCone::IsNurbSurfaceCone rebuilds its NURBS and reparameterises it, where the
// others adopt the one tested). The face's trims are in crNurbs's parameters, so only that NURBS
// will do.
SmSurface *AnalyticOfNurbs(const SmContext &ctx, const SmBSplineSurface &crNurbs, int64_t tKind, double dTol)
{
    SmBSplineSurface *pAnalytic = NULL;
    SmBoolean bFound = FALSE;
    if (tKind == kSurfaceCylinder)
    {
        SmCylinder *p = NULL;
        bFound = SmCylinder::IsNurbSurfaceCylinder(ctx, &crNurbs, p);
        pAnalytic = p;
    }
    else if (tKind == kSurfaceCone)
    {
        SmCone *p = NULL;
        bFound = SmCone::IsNurbSurfaceCone(ctx, &crNurbs, p);
        pAnalytic = p;
    }
    else if (tKind == kSurfaceSphere)
    {
        SmSphere *p = NULL;
        bFound = SmSphere::IsNurbSurfaceSphere(ctx, &crNurbs, p);
        pAnalytic = p;
    }
    else if (tKind == kSurfaceTorus)
    {
        SmTorus *p = NULL;
        bFound = SmTorus::IsNurbSurfaceTorus(ctx, &crNurbs, p);
        pAnalytic = p;
    }
    // IsNurbSurfaceCone can return an error status as its SmBoolean: only TRUE is a yes.
    if (bFound != TRUE || pAnalytic == NULL || !SameNurbs(crNurbs, *pAnalytic, dTol))
    {
        delete pAnalytic;
        return NULL;
    }
    return pAnalytic;
}

// An SmAxis2Placement (`out`) from an origin, an axis (z) and a reference direction (x):
// X = refDirection, Y = axis cross refDirection.
SmStatus Placement(const SmPoint3d &origin, const SmVector3d &zAxis, const SmVector3d &xAxis, SmAxis2Placement &out)
{
    SmVector3d sYAxis = zAxis * xAxis; // cross product
    return out.SetCanonical(origin, xAxis, sYAxis);
}

// The 3D curve of one edge (or wire edge), as SmuConvert builds it for its curve type.
// dMin/dMax are the edge range already converted to SMLib units (degrees for circles and
// ellipses).
SmStatus MakeEdgeCurve(
    const SmContext &ctx,
    int64_t type,
    double dMin,
    double dMax,
    const CurveArrays &nurbs,
    Cursor &c,
    const std::vector<double> &circleCenter,
    const std::vector<double> &circleAxis,
    const std::vector<double> &circleRef,
    const std::vector<double> &circleRadius,
    const std::vector<double> &lineOrigin,
    const std::vector<double> &lineDirection,
    const std::vector<double> &ellipseCenter,
    const std::vector<double> &ellipseAxis,
    const std::vector<double> &ellipseRef,
    const std::vector<double> &ellipseX,
    const std::vector<double> &ellipseY,
    SmCurve *&out
)
{
    out = NULL;
    if (type == kCurveNurbs)
    {
        out = smu_CreateSmBSplineCurve3d(ctx, nurbs, c);
        return out != NULL ? SM_SUCCESS : SM_ERR;
    }
    if (type == kCurveCircle)
    {
        size_t idx = c.circle++;
        return BuildAnalyticEdgeCurve(
            ctx, type, dMin, dMax, Get3(circleCenter, idx), GetVec3(circleAxis, idx), GetVec3(circleRef, idx), circleRadius[idx], circleRadius[idx], out
        );
    }
    if (type == kCurveLine)
    {
        size_t idx = c.line++;
        return BuildAnalyticEdgeCurve(ctx, type, dMin, dMax, Get3(lineOrigin, idx), GetVec3(lineDirection, idx), SmVector3d(), 0.0, 0.0, out);
    }
    size_t idx = c.ellipse++;
    return BuildAnalyticEdgeCurve(
        ctx, type, dMin, dMax, Get3(ellipseCenter, idx), GetVec3(ellipseAxis, idx), GetVec3(ellipseRef, idx), ellipseX[idx], ellipseY[idx], out
    );
}

// The SMLib form of a USD edge range: radians -> degrees for circles and ellipses.
void EdgeRangeToSMLib(int64_t type, double &dMin, double &dMax)
{
    if (type == kCurveCircle || type == kCurveEllipse)
    {
        dMin = UsdRadiansToSMLibDegrees(dMin);
        dMax = UsdRadiansToSMLibDegrees(dMax);
    }
}

/*******************************************************************/ /**
 PURPOSE: FlatBrep -> SmBrepData. Port of SMU_BrepConvert::BrepMove_OneUsdBrepToSMLibData.

 NOTES: packing of SmBrepData arrays (from SmuConvert):
      m_vEdges    = [FaceEdges, WireEdges]
      m_vVertices = [ShellVertices, Vertices]
      m_vEdgeuses = [new WireEdge-Edgeuses, FaceEdge-Edgeuses] (no SmEUData for LoopVertices
                    or ShellVertices), with WireEdge_ii-Edgeuse_EdgeIndex = EdgeCount + ii
 ***********************************************************************/
SmStatus FlatToBrepData(const FlatBrep &F, const Shape &S, const SmContext &crSmContext, SmBrepData *&prSmBrepData_ii)
{
    prSmBrepData_ii = new (crSmContext) SmBrepData(FALSE, SM_DS_USD);
    SmObjDelete sClean(prSmBrepData_ii);
    SmBrepData *pData = prSmBrepData_ii;

    // Restore the packed maximum as the topology zone tolerance. Reconstruct the corresponding
    // SMLib intersection tolerance only for approximation.
    SmXSectTol3d sXSectTol3d = F.tolerance[0];
    SmZoneTol3d sZoneTol3d = SmTol::GetZoneTol3d(sXSectTol3d);
    SmApproxTol3d sApproxTol3d = SmTol::GetApproxTol3d(sXSectTol3d);

    pData->m_bManifold = UNSURE;                                     // determined during brep construction
    pData->m_sZoneTol3d = sZoneTol3d;
    pData->m_dThisModelSizeEstimate = SM_MODEL_SIZE_ESTIMATE;        // not in the arrays - default value
    pData->m_dThisLargeSmallSizeRatio = SM_LARGE_SMALL_SIZE_RATIO;   // not in the arrays - default value
    pData->m_sApproxTol3d = sApproxTol3d;
    pData->m_lStartRegion = 0; // First region is always the infinite region
    pData->m_lNumRegions = static_cast<ULONG>(S.regions);
    pData->m_vColor.Set(0.0, 0.0, 0.0);
    pData->m_bIsGeomBorrowed = FALSE; // curves & surfaces in arrays are internally allocated.

    const ULONG nEdges = static_cast<ULONG>(S.edges);
    const ULONG nWire = static_cast<ULONG>(S.wire_edges);
    const ULONG nPointShells = static_cast<ULONG>(S.point_shells);

    // Regions
    pData->m_vRegions.SetSize(static_cast<ULONG>(S.regions));
    ULONG lShellStartIndex = 0;
    for (size_t i = 0; i < S.regions; ++i)
    {
        SmRegionData &rRegionData = pData->m_vRegions[static_cast<ULONG>(i)];
        rRegionData.ReSet();
        rRegionData.m_lStartShell = lShellStartIndex;
        rRegionData.m_lNumShells = static_cast<ULONG>(F.region_shell_count[i]);
        rRegionData.m_bIsVoidFlag = F.region_type[i] == 1 ? TRUE : FALSE;
        // make unused values look like SMLib unused values (simplifies debugging)
        rRegionData.m_lFlags = 0;
        rRegionData.m_lUserIndex1 = 0;
        rRegionData.m_lUserIndex2 = 0;
        lShellStartIndex += rRegionData.m_lNumShells;
    }

    // Shells - set m_vShells and the wire-edge entries of m_vEdgeuses
    pData->m_vShells.SetSize(static_cast<ULONG>(S.shells));
    pData->m_vEdgeuses.SetSize(nWire + static_cast<ULONG>(S.edgeuses));
    pData->m_vVertices.SetSize(nPointShells + static_cast<ULONG>(S.vertices));

    ULONG lFaceuseStartIndex = 0;
    ULONG lWireEdgeIndex = 0;
    ULONG lShellVertexPointIndex = 0;
    for (ULONG iLocal = 0; iLocal < S.shells; iLocal++)
    {
        SmShellData &rShellData = pData->m_vShells[iLocal];
        rShellData.ReSet();
        const int64_t nFaceuses = F.shell_faceuse_count[iLocal];
        const int64_t nWireEdges = F.shell_wire_edge_count[iLocal];
        rShellData.m_lShellType = nFaceuses > 0 ? 0 : nWireEdges > 0 ? 1 : 2; // faceuse, wireEdge or Vertex shell
        rShellData.m_lNumFaceuses = 0;
        rShellData.m_lFlags = 0;
        rShellData.m_lUserIndex1 = 0;
        rShellData.m_lUserIndex2 = 0;

        if (nFaceuses > 0)
        {
            rShellData.m_lFaceuseStart = lFaceuseStartIndex;
            rShellData.m_lNumFaceuses = static_cast<ULONG>(nFaceuses);
            lFaceuseStartIndex += static_cast<ULONG>(nFaceuses);
        }

        // Shells with WireEdges: m_vEdges = [FaceEdges, WireEdges], m_vEdgeuses = [WireEdge-Edgeuses, FaceEdge-Edgeuses]
        if (nWireEdges > 0)
        {
            rShellData.m_lEdge = lWireEdgeIndex + nEdges;
            for (int64_t ii = 0; ii < nWireEdges; ii++)
            {
                // 1 new WireEdge-Edgeuse per WireEdge, placed first on m_vEdgeuses to match SMLib use
                SmEUData &rEdgeuseData = pData->m_vEdgeuses[lWireEdgeIndex];
                rEdgeuseData.ReSet();
                rEdgeuseData.m_lEUType = 2; // 2 = WireEdge
                rEdgeuseData.m_lEdge = nEdges + lWireEdgeIndex;
                rEdgeuseData.m_lShell = iLocal;
                rEdgeuseData.m_bOrientation = TRUE;
                rEdgeuseData.m_lFlags = 0;
                rEdgeuseData.m_lUserIndex1 = 0;
                rEdgeuseData.m_lUserIndex2 = 0;
                rEdgeuseData.m_lMateFlags = 0;
                rEdgeuseData.m_lMateUserIndex1 = 0;
                rEdgeuseData.m_lMateUserIndex2 = 0;
                rEdgeuseData.m_lUVCurve = SM_NO_OBJECT; // not like SMLib but consistent with all other edgeuse values
                lWireEdgeIndex++;
            }
        }

        // VertexShells: their vertices come before the other vertices
        if (rShellData.m_lShellType == 2)
        {
            rShellData.m_lVertex = lShellVertexPointIndex;
            SmVertexData &rShell_VertexData = pData->m_vVertices[rShellData.m_lVertex];
            rShell_VertexData.ReSet();
            rShell_VertexData.m_sZoneTol3d = sZoneTol3d;
            rShell_VertexData.m_lFlags = 0;
            rShell_VertexData.m_lUserIndex1 = 0;
            rShell_VertexData.m_lUserIndex2 = 0;
            rShell_VertexData.m_vPoint = Get3(F.shell_point_position, lShellVertexPointIndex);
            lShellVertexPointIndex++;
        }
    }

    // Faceuses
    pData->m_vFaceuses.SetSize(static_cast<ULONG>(S.faceuses));
    for (ULONG iLocal = 0; iLocal < S.faceuses; iLocal++)
    {
        SmFaceuseData &rFaceuseData = pData->m_vFaceuses[iLocal];
        rFaceuseData.ReSet();
        rFaceuseData.m_lFace = static_cast<ULONG>(F.faceuse_face[iLocal]);
        rFaceuseData.m_bOrientation = F.faceuse_orientation[iLocal] == 1 ? TRUE : FALSE;
        rFaceuseData.m_lFlags = 0;
        rFaceuseData.m_lUserIndex1 = 0;
        rFaceuseData.m_lUserIndex2 = 0;
    }

    // Faces
    pData->m_vFaces.SetSize(static_cast<ULONG>(S.faces));
    ULONG lLoopIndex = 0;
    Cursor sNurbsSurfaces;
    size_t iSphere = 0, iPlane = 0, iCylinder = 0, iCone = 0, iTorus = 0;
    for (ULONG iLocal = 0; iLocal < S.faces; iLocal++)
    {
        SmFaceData &rFaceData = pData->m_vFaces[iLocal];
        rFaceData.ReSet();
        const int64_t tSurfaceType = F.face_surface_type[iLocal];
        const bool bNurbsForm = tSurfaceType != kSurfaceNurbs && FaceHasNurbs(F, iLocal);

        // Mapping face:trimType -> SMLib m_bRectangularTrim (tri-state):
        //   "rectangular" -> TRUE   (definitively naturally/rectangularly trimmed)
        //   "general"     -> UNSURE (there is no way to distinguish FALSE from UNSURE; treat as
        //                            "not determined" so round-trips of freshly-built SMLib
        //                            geometry (all UNSURE) are lossless, and the tessellator
        //                            re-computes as needed)
        SmBoolean bRectangularTrim = F.face_trim_type[iLocal] == 1 ? TRUE : UNSURE;
        SmExtent2d sDomain(F.face_range[4 * iLocal], F.face_range[4 * iLocal + 1], F.face_range[4 * iLocal + 2], F.face_range[4 * iLocal + 3]);

        // Analytic surface face:range angular parameters are radians, but SMLib analytic
        // surfaces (SmSphere, SmCylinder, SmCone, SmTorus) use degrees for their STEP
        // parameterization. Convert angular parameters from radians to degrees.
        if (bNurbsForm)
        {
            // The range is in the NURBS form's parameters already.
        }
        else if (tSurfaceType == kSurfaceSphere || tSurfaceType == kSurfaceTorus)
        {
            // Sphere/Torus: both U (rotation) and V (latitude/tube angle) are angular
            sDomain = SmExtent2d(
                UsdRadiansToSMLibDegrees(sDomain.GetMin().x), UsdRadiansToSMLibDegrees(sDomain.GetMin().y), UsdRadiansToSMLibDegrees(sDomain.GetMax().x),
                UsdRadiansToSMLibDegrees(sDomain.GetMax().y)
            );
        }
        else if (tSurfaceType == kSurfaceCylinder || tSurfaceType == kSurfaceCone)
        {
            // Cylinder/Cone: only U (rotation) is angular, V is distance along axis
            sDomain = SmExtent2d(UsdRadiansToSMLibDegrees(sDomain.GetMin().x), sDomain.GetMin().y, UsdRadiansToSMLibDegrees(sDomain.GetMax().x), sDomain.GetMax().y);
        }

        rFaceData.m_sZoneTol3d = sZoneTol3d; // OldTol: face's tolerant neighborhood offset size
        rFaceData.m_lSurface = pData->m_vSurfaces.GetSize();
        rFaceData.m_lStartLoop = lLoopIndex;
        rFaceData.m_lNumLoops = static_cast<ULONG>(F.face_loop_count[iLocal]);
        if (rFaceData.m_vUVDomain.Set(sDomain) != SM_SUCCESS)
        {
            prSmBrepData_ii = NULL;
            return SM_ERR_INVALID_INPUT;
        }
        rFaceData.m_bRectangularTrim = bRectangularTrim;
        rFaceData.m_lFlags = 0;
        rFaceData.m_lUserIndex1 = 0;
        rFaceData.m_lUserIndex2 = 0;

        const SmContext &ctx = *pData->GetContext();
        SmSurface *pNewSurface = NULL;
        SmStatus stat = SM_SUCCESS;
        if (bNurbsForm)
        {
            // An analytic face given with its NURBS form: that NURBS, as the analytic class the
            // arrays name (its canonical entry is passed over; the NURBS defines the surface).
            iSphere += tSurfaceType == kSurfaceSphere ? 1 : 0;
            iCylinder += tSurfaceType == kSurfaceCylinder ? 1 : 0;
            iCone += tSurfaceType == kSurfaceCone ? 1 : 0;
            iTorus += tSurfaceType == kSurfaceTorus ? 1 : 0;
            SmBSplineSurface *pNurbs = smu_CreateSmBSplineSurface3d(ctx, F, sNurbsSurfaces);
            SmObjDelete sNurbsCleanup(pNurbs);
            pNewSurface = pNurbs != NULL ? AnalyticOfNurbs(ctx, *pNurbs, tSurfaceType, F.tolerance[0]) : NULL;
            stat = pNewSurface != NULL ? SM_SUCCESS : SM_ERR_INVALID_INPUT;
        }
        else if (tSurfaceType == kSurfaceSphere)
        {
            size_t i = iSphere++;
            SmAxis2Placement sOrigin;
            if (Placement(Get3(F.face_sphere_center, i), GetVec3(F.face_sphere_axis, i), GetVec3(F.face_sphere_ref_direction, i), sOrigin) != SM_SUCCESS)
            {
                prSmBrepData_ii = NULL;
                return SM_ERR_INVALID_INPUT;
            }
            SmSphere *pSphere = NULL;
            stat = SmSphere::CreateCanonical(ctx, sOrigin, F.face_sphere_radius[i], pSphere);
            pNewSurface = pSphere;
        }
        else if (tSurfaceType == kSurfacePlane)
        {
            size_t i = iPlane++;
            SmAxis2Placement sPlacement;
            if (Placement(Get3(F.face_plane_origin, i), GetVec3(F.face_plane_axis, i), GetVec3(F.face_plane_ref_direction, i), sPlacement) != SM_SUCCESS)
            {
                prSmBrepData_ii = NULL;
                return SM_ERR_INVALID_INPUT;
            }
            SmPlane *pPlane = NULL;
            stat = SmPlane::CreateCanonical(ctx, sPlacement, pPlane);
            // CreateCanonical sets the plane's UV domain to [-SM_INFINITE_PARAMETER,
            // SM_INFINITE_PARAMETER]. The face's trim curves use sDomain coordinates (e.g.
            // [0,1]), so the plane's analytical domain must match. AdjustSTEPUVDomain (below)
            // restricts the domain without changing the parameterization (UVScale stays 1,1),
            // unlike Reparameterize which would rescale the plane so that sDomain covers the
            // entire infinite 3D extent.
            pNewSurface = pPlane;
        }
        else if (tSurfaceType == kSurfaceCylinder)
        {
            size_t i = iCylinder++;
            SmAxis2Placement sPlacement;
            if (Placement(Get3(F.face_cylinder_origin, i), GetVec3(F.face_cylinder_axis, i), GetVec3(F.face_cylinder_ref_direction, i), sPlacement) != SM_SUCCESS)
            {
                prSmBrepData_ii = NULL;
                return SM_ERR_INVALID_INPUT;
            }
            SmCylinder *pCylinder = NULL;
            stat = SmCylinder::CreateCanonical(ctx, sPlacement, F.face_cylinder_radius[i], pCylinder);
            pNewSurface = pCylinder;
        }
        else if (tSurfaceType == kSurfaceCone)
        {
            size_t i = iCone++;
            SmPoint3d sConeOrigin = Get3(F.face_cone_origin, i);
            SmVector3d sZAxis = GetVec3(F.face_cone_axis, i);
            SmVector3d sXAxis = GetVec3(F.face_cone_ref_direction, i);
            double dRadius = F.face_cone_radius[i];
            double dSemiAngleDeg = SM_RAD2DEG(F.face_cone_semi_angle[i]);

            // SmCone::CreateCanonical requires a non-negative semi-angle. A cone with a negative
            // semi-angle has an exactly equivalent positive-angle form when its axis and both
            // STEP parameters are reversed. Reversing both parameters preserves the surface
            // normal as well as the evaluated points.
            if (dSemiAngleDeg < 0.0)
            {
                dSemiAngleDeg = -dSemiAngleDeg;
                sZAxis = -sZAxis;
                if (sDomain.SetMinMax(-sDomain.GetUMax(), -sDomain.GetVMax(), -sDomain.GetUMin(), -sDomain.GetVMin()) != SM_SUCCESS)
                {
                    prSmBrepData_ii = NULL;
                    return SM_ERR_INVALID_INPUT;
                }
            }

            SmCone *pCone = NULL;
            stat = SM_ERR;
            if (sDomain.IsBounded())
            {
                // AdjustSTEPUVDomain on an initially unbounded SmCone does not preserve the
                // surface when the requested axial interval starts away from zero. Construct the
                // equivalent finite segment directly at Vmin and normalize its SMLib axial
                // parameter to [0, height]. Analytic-face UV curves are regenerated from the 3D
                // edges, so no authored trim coordinates depend on the original V offset here.
                const double dVMin = sDomain.GetVMin();
                const double dVMax = sDomain.GetVMax();
                const double dHeight = dVMax - dVMin;
                const double dTanSemiAngle = std::tan(SM_DEG2RAD(dSemiAngleDeg));
                const double dBottomRadius = dRadius + dVMin * dTanSemiAngle;
                const double dTopRadius = dRadius + dVMax * dTanSemiAngle;
                if (dHeight > SM_EFF_ZERO && dBottomRadius >= -SM_EFF_ZERO && dTopRadius >= -SM_EFF_ZERO)
                {
                    sConeOrigin += dVMin * sZAxis;
                    SmAxis2Placement sPlacement;
                    if (Placement(sConeOrigin, sZAxis, sXAxis, sPlacement) != SM_SUCCESS)
                    {
                        prSmBrepData_ii = NULL;
                        return SM_ERR_INVALID_INPUT;
                    }
                    stat = SmCone::CreateCanonical(ctx, sPlacement, smos_Max(0.0, dBottomRadius), smos_Max(0.0, dTopRadius), dHeight, pCone);
                    if (sDomain.SetMinMax(sDomain.GetUMin(), 0.0, sDomain.GetUMax(), dHeight) != SM_SUCCESS)
                    {
                        prSmBrepData_ii = NULL;
                        return SM_ERR_INVALID_INPUT;
                    }
                }
            }
            else
            {
                SmAxis2Placement sPlacement;
                if (Placement(sConeOrigin, sZAxis, sXAxis, sPlacement) != SM_SUCCESS)
                {
                    prSmBrepData_ii = NULL;
                    return SM_ERR_INVALID_INPUT;
                }
                stat = SmCone::CreateCanonical(ctx, sPlacement, dRadius, dSemiAngleDeg, pCone);
            }
            pNewSurface = pCone;
        }
        else if (tSurfaceType == kSurfaceTorus)
        {
            size_t i = iTorus++;
            SmAxis2Placement sPlacement;
            if (Placement(Get3(F.face_torus_origin, i), GetVec3(F.face_torus_axis, i), GetVec3(F.face_torus_ref_direction, i), sPlacement) != SM_SUCCESS)
            {
                prSmBrepData_ii = NULL;
                return SM_ERR_INVALID_INPUT;
            }
            SmTorus *pTorus = NULL;
            stat = SmTorus::CreateCanonical(ctx, sPlacement, F.face_torus_major_radius[i], F.face_torus_minor_radius[i], pTorus);
            pNewSurface = pTorus;
        }
        else // NURBS
        {
            pNewSurface = smu_CreateSmBSplineSurface3d(ctx, F, sNurbsSurfaces);
            stat = pNewSurface != NULL ? SM_SUCCESS : SM_ERR;
        }

        SmObjDelete sSurfaceCleanup(pNewSurface);
        if (stat != SM_SUCCESS || pNewSurface == NULL)
        {
            prSmBrepData_ii = NULL;
            return stat != SM_SUCCESS ? stat : SM_ERR;
        }

        const bool bAnalyticSurface = tSurfaceType != kSurfaceNurbs && !bNurbsForm;
        if (bAnalyticSurface && sDomain.IsBounded())
        {
            SmStatus adjStat = pNewSurface->AdjustSTEPUVDomain(sDomain);
            if (adjStat != SM_SUCCESS)
            {
                prSmBrepData_ii = NULL;
                return adjStat;
            }
        }
        pData->m_vSurfaces.Add(pNewSurface);
        sSurfaceCleanup.Clear();

        // SmFaceData stores the NURBS parameter domain consumed by MakeTopologyFromData, while
        // face:range for analytic surfaces is a STEP domain (angles are converted to degrees
        // above). AdjustSTEPUVDomain rebuilds the analytic surface over that STEP range, but its
        // natural domain is the resulting NURBS interval. Convert the face domain as well so a
        // shifted angular range (for example 270..360 degrees) is not intersected with a NURBS
        // interval beginning at zero and cleared as empty.
        if (bAnalyticSurface && sDomain.IsBounded())
        {
            SmExtent2d sNurbDomain;
            SmStatus domainStat = pNewSurface->ConvertDomainFromSTEPToNURBS(sDomain, sNurbDomain);
            if (domainStat != SM_SUCCESS)
            {
                prSmBrepData_ii = NULL;
                return domainStat;
            }
            if (rFaceData.m_vUVDomain.Set(sNurbDomain) != SM_SUCCESS)
            {
                prSmBrepData_ii = NULL;
                return SM_ERR_INVALID_INPUT;
            }
        }

        lLoopIndex += static_cast<ULONG>(F.face_loop_count[iLocal]);
    }

    // Loops; m_vEdgeuses = [ WireEdge-Edgeuses, FaceEdge-Edgeuses ]. Each edgeuse's loop and
    // face follow from the loop counts (SmuConvert walks them with running stop indices; the
    // validated counts make that walk exactly this one).
    pData->m_vLoops.SetSize(static_cast<ULONG>(S.loops));
    std::vector<ULONG> sEdgeuseLoop(S.edgeuses), sEdgeuseFace(S.edgeuses);
    {
        ULONG lEdgeuseIndex = nWire;
        size_t eu = 0, face = 0;
        ULONG faceLoopEnd = S.faces > 0 ? static_cast<ULONG>(F.face_loop_count[0]) : 0;
        for (ULONG iLocal = 0; iLocal < S.loops; iLocal++)
        {
            while (iLocal >= faceLoopEnd && face + 1 < S.faces)
            {
                face++;
                faceLoopEnd += static_cast<ULONG>(F.face_loop_count[face]);
            }
            SmLoopData &rLoopData = pData->m_vLoops[iLocal];
            rLoopData.ReSet();
            const int64_t nEU = F.loop_edgeuse_count[iLocal];
            rLoopData.m_lLoopType = nEU > 0 ? 0 : 1; // edgeuse loop : vertex loop
            rLoopData.m_lNumEU = static_cast<ULONG>(nEU);
            rLoopData.m_lFlags = 0;
            rLoopData.m_lUserIndex1 = 0;
            rLoopData.m_lUserIndex2 = 0;
            if (nEU > 0)
            {
                rLoopData.m_lStartEU = lEdgeuseIndex;
                lEdgeuseIndex += static_cast<ULONG>(nEU);
                for (int64_t k = 0; k < nEU; ++k, ++eu)
                {
                    sEdgeuseLoop[eu] = iLocal;
                    sEdgeuseFace[eu] = static_cast<ULONG>(face);
                }
            }
            else
            {
                // SmBrepData->m_vEdgeuses does not contain a SmEUData for each LoopVertex
                rLoopData.m_lVertex = nPointShells + static_cast<ULONG>(F.loop_vertex[iLocal]);
            }
        }
    }

    // Edgeuses (face edges only; wire-edge edgeuses were made in the Shell block)
    std::vector<long> sEdgePrimEU(S.edges, SM_NO_OBJECT);
    Cursor sUVCurves;
    for (size_t k = 0; k < S.edgeuses; ++k)
    {
        const ULONG iLocal = nWire + static_cast<ULONG>(k);
        const size_t iNext = static_cast<size_t>(F.edgeuse_next_radial[k]);
        const long iLocalNextRadialEUIndex = static_cast<long>(nWire + iNext);
        const SmBoolean bThisTopEntry = F.edgeuse_radial_entry[k] == 1 ? TRUE : FALSE;
        const SmBoolean bNextTopEntry = F.edgeuse_radial_entry[iNext] == 1 ? TRUE : FALSE;
        const long iLocalNextRadialEUFlag = bNextTopEntry ? iLocalNextRadialEUIndex : -iLocalNextRadialEUIndex - 1;
        const SmBoolean bOrientation = F.edgeuse_orientation[k] == 1 ? TRUE : FALSE;

        SmEUData &rEdgeuseData = pData->m_vEdgeuses[iLocal];
        rEdgeuseData.ReSet();
        rEdgeuseData.m_lEUType = 0;
        rEdgeuseData.m_lEdge = static_cast<ULONG>(F.edgeuse_edge[k]);
        rEdgeuseData.m_lLoop = sEdgeuseLoop[k];
        rEdgeuseData.m_bOrientation = bOrientation;
        rEdgeuseData.m_bDownwardLU = false;
        rEdgeuseData.m_lFlags = 0;
        rEdgeuseData.m_lUserIndex1 = 0;
        rEdgeuseData.m_lUserIndex2 = 0;

        // map from {NextRadialIndex, RadialEntryType} to {NextRadialEUFlag, MateNextRadialEUFlag}
        //  topEntry:    From Top - to Bot - to NextRadialEUFlag
        //  bottomEntry: From Bot - to Top - to NextRadialEUFlag
        //  where NextRadialEUFlag = enter the next edgeuse's TopEdgeuse when Flag >= 0, its
        //  BotEdgeuse when Flag < 0
        if (bThisTopEntry)
        {
            rEdgeuseData.m_lNextEU = -static_cast<long>(iLocal) - 1; // TopEdgeuse
            rEdgeuseData.m_lMateNextEU = iLocalNextRadialEUFlag;     // BotEdgeuse
        }
        else
        {
            rEdgeuseData.m_lNextEU = iLocalNextRadialEUFlag;       // TopEdgeuse
            rEdgeuseData.m_lMateNextEU = static_cast<long>(iLocal); // BotEdgeuse
        }

        // Attach UV TrimCurves to edgeuses
        rEdgeuseData.m_lUVCurve = SM_NO_OBJECT;
        if (S.uv_curves)
        {
            SmBSplineCurve *pBSplineCurve = smu_CreateSmBSplineCurve2d_ForEdgeuseIndex(
                *pData->GetContext(), F, sUVCurves, k, FaceHasNurbs(F, sEdgeuseFace[k]), bOrientation
            );
            if (pBSplineCurve != NULL)
            {
                rEdgeuseData.m_lUVCurve = pData->m_vUVCurves.GetSize();
                pData->m_vUVCurves.Add(pBSplineCurve);
            }
        }

        // these will be future edge->PrimEU values. The PrimaryEdgeuse is the first edgeuse in
        // the radial list and the 2nd edgeuse must be that edgeuse's mate.
        //    for TopEntry EUs use the top Edgeuse index:[iLocal]
        //    for BotEntry EUs use the bot Edgeuse index:[-ilocal-1]
        if (sEdgePrimEU[rEdgeuseData.m_lEdge] == SM_NO_OBJECT)
        {
            sEdgePrimEU[rEdgeuseData.m_lEdge] = bThisTopEntry ? static_cast<long>(iLocal) : -static_cast<long>(iLocal) - 1;
        }
    }

    // Edges - m_vEdges = [FaceEdges, WireEdges]
    pData->m_vEdges.SetSize(nEdges + nWire);
    {
        Cursor c;
        const CurveArrays nurbs{F.edge_nurbs_vertex_count, F.edge_nurbs_order, F.edge_nurbs_cvs, F.edge_nurbs_knots, F.edge_nurbs_weights};
        for (ULONG iLocal = 0; iLocal < nEdges; iLocal++)
        {
            SmEdgeData &rEdgeData = pData->m_vEdges[iLocal];
            rEdgeData.ReSet();
            rEdgeData.m_sZoneTol3d = sZoneTol3d;
            if (!F.edge_tolerance.empty())
            {
                SmXSectTol3d sEdgeTol = F.edge_tolerance[iLocal];
                rEdgeData.m_sZoneTol3d = SmTol::GetZoneTol3d(sEdgeTol);
            }
            double dEdgeMin = F.edge_range[2 * iLocal];
            double dEdgeMax = F.edge_range[2 * iLocal + 1];
            // circle/ellipse edge range is radians; SMLib SmCircle/SmEllipse use degrees
            EdgeRangeToSMLib(F.edge_curve_type[iLocal], dEdgeMin, dEdgeMax);
            if (rEdgeData.m_vInterval.SetMinMax(dEdgeMin, dEdgeMax) != SM_SUCCESS)
            {
                prSmBrepData_ii = NULL;
                return SM_ERR_INVALID_INPUT;
            }
            rEdgeData.m_lCurve = pData->m_v3DCurves.GetSize();
            rEdgeData.m_lStartVertex = nPointShells + static_cast<ULONG>(F.edge_vertices[2 * iLocal]);
            rEdgeData.m_lEndVertex = nPointShells + static_cast<ULONG>(F.edge_vertices[2 * iLocal + 1]);
            rEdgeData.m_lPrimEU = sEdgePrimEU[iLocal];
            rEdgeData.m_lFlags = 0;
            rEdgeData.m_lUserIndex1 = 0;
            rEdgeData.m_lUserIndex2 = 0;

            SmCurve *pCurve = NULL;
            SmStatus stat = MakeEdgeCurve(
                *pData->GetContext(), F.edge_curve_type[iLocal], dEdgeMin, dEdgeMax, nurbs, c, F.edge_circle_center, F.edge_circle_axis,
                F.edge_circle_ref_direction, F.edge_circle_radius, F.edge_line_origin, F.edge_line_direction, F.edge_ellipse_center, F.edge_ellipse_axis,
                F.edge_ellipse_ref_direction, F.edge_ellipse_x_radius, F.edge_ellipse_y_radius, pCurve
            );
            if (stat != SM_SUCCESS)
            {
                prSmBrepData_ii = NULL;
                return stat;
            }
            pData->m_v3DCurves.Add(pCurve);
        }
    }

    // WireEdges
    {
        Cursor c;
        const CurveArrays nurbs{F.wire_edge_nurbs_vertex_count, F.wire_edge_nurbs_order, F.wire_edge_nurbs_cvs, F.wire_edge_nurbs_knots, F.wire_edge_nurbs_weights};
        for (ULONG w = 0; w < nWire; w++)
        {
            const ULONG iLocal = nEdges + w;
            SmEdgeData &rWireEdgeData = pData->m_vEdges[iLocal];
            rWireEdgeData.ReSet();
            rWireEdgeData.m_sZoneTol3d = sZoneTol3d;
            double dWireMin = F.wire_edge_range[2 * w];
            double dWireMax = F.wire_edge_range[2 * w + 1];
            EdgeRangeToSMLib(F.wire_edge_curve_type[w], dWireMin, dWireMax);
            if (rWireEdgeData.m_vInterval.SetMinMax(dWireMin, dWireMax) != SM_SUCCESS)
            {
                prSmBrepData_ii = NULL;
                return SM_ERR_INVALID_INPUT;
            }
            rWireEdgeData.m_lCurve = pData->m_v3DCurves.GetSize();
            rWireEdgeData.m_lStartVertex = nPointShells + static_cast<ULONG>(F.wire_edge_vertices[2 * w]);
            rWireEdgeData.m_lEndVertex = nPointShells + static_cast<ULONG>(F.wire_edge_vertices[2 * w + 1]);
            rWireEdgeData.m_lPrimEU = static_cast<long>(w); // WireEdge_Edgeuse index = WireEdge index - EdgeCount
            rWireEdgeData.m_lFlags = 0;
            rWireEdgeData.m_lUserIndex1 = 0;
            rWireEdgeData.m_lUserIndex2 = 0;

            SmCurve *pCurve = NULL;
            SmStatus stat = MakeEdgeCurve(
                *pData->GetContext(), F.wire_edge_curve_type[w], dWireMin, dWireMax, nurbs, c, F.wire_edge_circle_center, F.wire_edge_circle_axis,
                F.wire_edge_circle_ref_direction, F.wire_edge_circle_radius, F.wire_edge_line_origin, F.wire_edge_line_direction, F.wire_edge_ellipse_center,
                F.wire_edge_ellipse_axis, F.wire_edge_ellipse_ref_direction, F.wire_edge_ellipse_x_radius, F.wire_edge_ellipse_y_radius, pCurve
            );
            if (stat != SM_SUCCESS)
            {
                prSmBrepData_ii = NULL;
                return stat;
            }
            pData->m_v3DCurves.Add(pCurve);
        }
    }

    // Vertices - m_vVertices = [ShellVertices, Vertices]
    for (ULONG v = 0; v < S.vertices; v++)
    {
        SmVertexData &rVertexData = pData->m_vVertices[nPointShells + v];
        rVertexData.ReSet();
        rVertexData.m_sZoneTol3d = sZoneTol3d;
        if (!F.vertex_tolerance.empty())
        {
            SmXSectTol3d sVertexTol = F.vertex_tolerance[v];
            rVertexData.m_sZoneTol3d = SmTol::GetZoneTol3d(sVertexTol);
        }
        rVertexData.m_lFlags = 0;
        rVertexData.m_lUserIndex1 = 0;
        rVertexData.m_lUserIndex2 = 0;
        rVertexData.m_vPoint = Get3(F.vertex_position, v);
    }

    sClean.Clear();
    return SM_SUCCESS;
}

} // namespace

/*******************************************************************/ /**
 PURPOSE: FlatBrep -> new SmBrep. Port of SMU_BrepConvert::BrepMove_OneUsdBrepToSMLib, with
          the healer enabled (SmuConvert's default).
 ***********************************************************************/
// SMLib takes the first shell of every region but the infinite one as the shell that encloses
// the region, the others lying inside it (SmBrep::AssertValid test 12; BrepData has no other way
// to say which one it is). The arrays may list a region's shells in any order, so the one that
// encloses all the others is put first, decided by SMLib's own containment test
// (SmShell::IsInnerShellOf, the ray casts test 12 runs). Candidates go largest bounding box
// first, as an enclosing shell's box holds every other's, so the first one tried is nearly
// always it. A region where no shell can be shown to enclose the others keeps its order.
void PutEnclosingShellsFirst(SmBrep &rBrep)
{
    SmTArray<SmRegion *> sRegions;
    rBrep.GetRegions(sRegions);
    const SmRegion *pInfinite = rBrep.GetInfiniteRegion();
    for (ULONG ii = 0; ii < sRegions.GetSize(); ++ii)
    {
        SmRegion *pRegion = sRegions[ii];
        if (pRegion == NULL || pRegion == pInfinite)
        {
            continue;
        }
        SmTArray<SmShell *> sShells;
        pRegion->GetShells(sShells);
        if (sShells.GetSize() < 2)
        {
            continue;
        }
        std::vector<std::pair<double, SmShell *>> vCandidates;
        for (ULONG jj = 0; jj < sShells.GetSize(); ++jj)
        {
            SmShell *pShell = sShells[jj];
            SmExtent3d sBox;
            if (pShell->IsFaceuseShell() && pShell->IsClosed() && pShell->CalculateBoundingBox(sBox) == SM_SUCCESS)
            {
                SmVector3d sSize = sBox.GetSize();
                vCandidates.emplace_back(sSize.x * sSize.y * sSize.z, pShell);
            }
        }
        // stable: equal boxes keep the arrays' order
        std::stable_sort(vCandidates.begin(), vCandidates.end(),
                         [](const std::pair<double, SmShell *> &a, const std::pair<double, SmShell *> &b) { return a.first > b.first; });
        for (const std::pair<double, SmShell *> &rCandidate : vCandidates)
        {
            bool bEncloses = true;
            for (ULONG jj = 0; jj < sShells.GetSize() && bEncloses; ++jj)
            {
                bEncloses = sShells[jj] == rCandidate.second || sShells[jj]->IsInnerShellOf(rCandidate.second) == TRUE;
            }
            if (bEncloses)
            {
                if (rCandidate.second != sShells[0])
                {
                    pRegion->SetList(rCandidate.second); // the circular list now starts at it
                }
                break;
            }
        }
    }
}

SmStatus BrepFromFlat(const FlatBrep &F, const SmContext &crSmContext, SmBrep *&prSmBrep_ii)
{
    prSmBrep_ii = NULL;

    Shape S;
    if (!Validate(F, S))
    {
        return SM_ERR_INVALID_INPUT;
    }

    SmBrepData *pBrepData_ii = NULL;
    SmStatus buildStatus = FlatToBrepData(F, S, crSmContext, pBrepData_ii);
    if (buildStatus != SM_SUCCESS || pBrepData_ii == NULL)
    {
        delete pBrepData_ii;
        return buildStatus != SM_SUCCESS ? buildStatus : SM_ERR;
    }
    SmObjDelete sCleanBrepData(pBrepData_ii);

    // FlatToBrepData restored the packed intersection tolerance as this topology zone tolerance.
    SmZoneTol3d sZoneTol3d = pBrepData_ii->m_sZoneTol3d;

    SmBrep *pBrep_ii = new (crSmContext) SmBrep();
    SmObjDelete sCleanBrep(pBrep_ii);

    // Kludge. No heal for B-reps written by SMLib (USD: customData source == "SMLib").
    ULONG lHealerVersion = (!F.source.empty() && F.source[0] == 1) ? SM_HEALER_VERSION : SM_HEALER_VERSION - 2;
    pBrepData_ii->SetHealerVersion(lHealerVersion);

    // Hard coding SmZoneTol3d here to avoid inflation of tolerance for slender breps, which
    // lead to issues in tessellation.
    pBrepData_ii->m_sZoneTol3d = sZoneTol3d;

    // For performance on tessellation we set the following options FALSE
    // UVTrimCurves will be generated ad hoc
    // Conversion to analytics can be completed as a separate heal step
    SmTArray<SmAttribute *> sAttributes; // no attributes are carried
    SmBoolean bCreateUVTrimCurves = FALSE;
    SmBoolean bConvertToAnalytics = FALSE;
    SmStatus topologyStatus = pBrep_ii->MakeTopologyFromData(
        pBrepData_ii,        // in : target BrepData object
        sAttributes,         // in : list of all attributes to be assigned to Brep entities
        bCreateUVTrimCurves, // in : default TRUE = call SmFace::CreateUVTrimCurves for each face
        bConvertToAnalytics, // in : default TRUE = replace every Surface and Curve with a derived analytic type when appropriate.
        TRUE,                // in : when reading breps set TRUE (use original tol)
        FALSE,               // in : bCheckUVTrimCurves
        FALSE,               // in : bAlwaysHeal
        TRUE                 // in : bHealerIsEnabled
    );
    if (topologyStatus != SM_SUCCESS)
    {
        return topologyStatus;
    }
    PutEnclosingShellsFirst(*pBrep_ii);

    prSmBrep_ii = pBrep_ii;
    sCleanBrep.Clear();
    return SM_SUCCESS;
}

/* ---- Text form ------------------------------------------------------------------------------- */

bool WriteText(const FlatBrep &F, const std::string &path)
{
    FILE *f = std::fopen(path.c_str(), "w");
    if (f == nullptr)
    {
        return false;
    }
    for (const Field &fd : Fields())
    {
        if (fd.i != nullptr && !(F.*(fd.i)).empty())
        {
            const std::vector<int64_t> &v = F.*(fd.i);
            std::fprintf(f, "%s i64 %zu", fd.name, v.size());
            for (int64_t x : v)
            {
                std::fprintf(f, " %lld", static_cast<long long>(x));
            }
            std::fprintf(f, "\n");
        }
        else if (fd.f != nullptr && !(F.*(fd.f)).empty())
        {
            const std::vector<double> &v = F.*(fd.f);
            std::fprintf(f, "%s f64 %zu", fd.name, v.size());
            for (double x : v)
            {
                std::fprintf(f, " %.17g", x);
            }
            std::fprintf(f, "\n");
        }
    }
    return std::fclose(f) == 0;
}

bool ReadText(const std::string &path, FlatBrep &F)
{
    F = FlatBrep();
    std::ifstream in(path);
    if (!in)
    {
        return false;
    }
    std::string name, type, token;
    while (in >> name)
    {
        long long len = -1;
        if (!(in >> type >> len) || len < 0 || len > SMLIB_FLAT_MAX_LEN)
        {
            return false;
        }
        const Field *fd = nullptr;
        for (const Field &x : Fields())
        {
            if (name == x.name)
            {
                fd = &x;
            }
        }
        if (fd == nullptr || (type == "i64") != (fd->i != nullptr) || (type != "i64" && type != "f64"))
        {
            return false;
        }
        for (long long k = 0; k < len; ++k)
        {
            if (!(in >> token))
            {
                return false;
            }
            char *end = nullptr;
            if (fd->i != nullptr)
            {
                long long x = std::strtoll(token.c_str(), &end, 10);
                if (*end != '\0')
                {
                    return false;
                }
                (F.*(fd->i)).push_back(x);
            }
            else
            {
                double x = std::strtod(token.c_str(), &end);
                if (*end != '\0')
                {
                    return false;
                }
                (F.*(fd->f)).push_back(x);
            }
        }
    }
    return in.eof();
}

} // namespace smlib_flat_impl
