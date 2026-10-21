// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
// SPDX-License-Identifier: Apache-2.0

/*******************************************************************//**
* FILE NAME --- SmAdvSurfaceIntersector.cpp
* PURPOSE: This file contains surface intersector methods.
**********************************************************************/

#include "StdAfx.h"

#include <SmAdvSurfaceIntersector.h>
#include <SmGeomUtility.h>
#include <SmConic.h>
#include <algorithm>
#include <cmath>

#ifdef SM_DEBUG_CODE
#include <SmFace.h>
#include <SmBrep.h>
#endif // SM_DEBUG_CODE



// SolidCraft change: interior tangential contact points.
//
// FindInteriorCurves() only traces crossing curves: node pairs whose surface
// normals are (nearly) parallel are skipped, so a point where the two surfaces
// merely touch inside both patches (an isolated tangential intersection point,
// e.g. a paraboloid apex resting on a plane) was never reported.  The helpers
// below solve for such points directly (Yang, Jia & Yan 2023, "Topology
// Guaranteed B-Spline Surface/Surface Intersection", treat isolated tangential
// points as the zeros of S1(u,v) - S2(s,t) at which the normals are parallel):
//
//   F1 = (Sa - Sb).Sb_s     F2 = (Sa - Sb).Sb_t     (Sb(s,t) is the foot of Sa(u,v) on Sb)
//   F3 = Sa_u.(Sb_s x Sb_t) F4 = Sa_v.(Sb_s x Sb_t) (normals parallel)
//
// solved by Newton's method with the exact Jacobian; a root is a contact point
// only if additionally |Sa - Sb| is within tolerance.  The system is regular at
// a nondegenerate (isolated) contact, so convergence is quadratic there.
namespace {

struct SmTpEval
{
  double P[3], Pu[3], Pv[3], Puu[3], Puv[3], Pvv[3];
};

inline double tp_dot(const double a[3], const double b[3]) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
inline void   tp_cross(const double a[3], const double b[3], double r[3])
{
  r[0] = a[1]*b[2]-a[2]*b[1]; r[1] = a[2]*b[0]-a[0]*b[2]; r[2] = a[0]*b[1]-a[1]*b[0];
}
inline void   tp_set(double r[3], const SmVector3d & v) { r[0] = v.x; r[1] = v.y; r[2] = v.z; }

SmStatus tp_Evaluate(const SmSurface & crSrf, const SmPoint2d & crUV, SmTpEval & rE)
{
  SmVector3d aD[9];
  SER(crSrf.Evaluate(crUV, 2, 2, TRUE, TRUE, TRUE, aD, FALSE));
  tp_set(rE.P,   aD[SM_S2(0,0)]);
  tp_set(rE.Pu,  aD[SM_S2(1,0)]);
  tp_set(rE.Pv,  aD[SM_S2(0,1)]);
  tp_set(rE.Puu, aD[SM_S2(2,0)]);
  tp_set(rE.Puv, aD[SM_S2(1,1)]);
  tp_set(rE.Pvv, aD[SM_S2(0,2)]);
  return SM_SUCCESS;
}

// Solve the 4x4 system A x = b (partial pivoting).  FALSE when singular.
SmBoolean tp_Solve4(double A[4][4], double b[4], double x[4])
{
  int i, j, k;
  for (k=0; k<4; k++)
    {
      int iPiv = k;
      for (i=k+1; i<4; i++) if (fabs(A[i][k]) > fabs(A[iPiv][k])) iPiv = i;
      if (fabs(A[iPiv][k]) < 1.0e-300) return FALSE;
      if (iPiv != k)
        {
          for (j=0; j<4; j++) { double t = A[k][j]; A[k][j] = A[iPiv][j]; A[iPiv][j] = t; }
          double t = b[k]; b[k] = b[iPiv]; b[iPiv] = t;
        }
      for (i=k+1; i<4; i++)
        {
          double f = A[i][k]/A[k][k];
          for (j=k; j<4; j++) A[i][j] -= f*A[k][j];
          b[i] -= f*b[k];
        }
    }
  for (i=3; i>=0; i--)
    {
      double s = b[i];
      for (j=i+1; j<4; j++) s -= A[i][j]*x[j];
      x[i] = s/A[i][i];
    }
  // reject a numerically singular system (relative pivot test)
  double dMax = 0.0, dMin = SM_BIG_DOUBLE;
  for (i=0; i<4; i++) { dMax = std::max(dMax, fabs(A[i][i])); dMin = std::min(dMin, fabs(A[i][i])); }
  return (dMin > 1.0e-13*dMax) ? TRUE : FALSE;
}

// Newton iteration for a tangential contact of Sa (u,v) and Sb (s,t), seeded at
// rUVa, rUVb.  On success the parameters are updated, rdGap is |Sa - Sb| and
// rbConverged tells whether the iteration converged (independently of the gap).
SmStatus tp_Newton(const SmSurface & crSa, const SmExtent2d & crDomA,
                   const SmSurface & crSb, const SmExtent2d & crDomB,
                   const SmExtent2d & crSeedA, const SmExtent2d & crSeedB,
                   SmPoint2d & rUVa, SmPoint2d & rUVb,
                   SmBoolean & rbConverged, double & rdGap)
{
  // a seed only has to find a contact near its own patch pair: give up when
  // the iteration leaves the seed patches grown by half their size each side
  SmExtent2d sLocA = crSeedA, sLocB = crSeedB;
  SmPoint2d sGrowA = crSeedA.GetSize(), sGrowB = crSeedB.GetSize();
  sLocA.SetMinMax(crSeedA.GetMin().x - 0.5*sGrowA.x, crSeedA.GetMin().y - 0.5*sGrowA.y,
                  crSeedA.GetMax().x + 0.5*sGrowA.x, crSeedA.GetMax().y + 0.5*sGrowA.y);
  sLocB.SetMinMax(crSeedB.GetMin().x - 0.5*sGrowB.x, crSeedB.GetMin().y - 0.5*sGrowB.y,
                  crSeedB.GetMax().x + 0.5*sGrowB.x, crSeedB.GetMax().y + 0.5*sGrowB.y);

  rbConverged = FALSE;
  rdGap = SM_BIG_DOUBLE;
  SmPoint2d sSizeA = crDomA.GetSize(), sSizeB = crDomB.GetSize();
  double adMaxStep[4] = { 0.25*sSizeA.x, 0.25*sSizeA.y, 0.25*sSizeB.x, 0.25*sSizeB.y };
  double adPar[4] = { rUVa.x, rUVa.y, rUVb.x, rUVb.y };

  // a few Gauss-Newton projection steps first so (s,t) starts at the foot of Sa(u,v)
  for (int it=0; it<8; it++)
    {
      SmTpEval A, B;
      SER(tp_Evaluate(crSa, SmPoint2d(adPar[0],adPar[1]), A));
      SER(tp_Evaluate(crSb, SmPoint2d(adPar[2],adPar[3]), B));
      double D[3] = { A.P[0]-B.P[0], A.P[1]-B.P[1], A.P[2]-B.P[2] };
      double a11 = tp_dot(B.Pu,B.Pu), a12 = tp_dot(B.Pu,B.Pv), a22 = tp_dot(B.Pv,B.Pv);
      double r1 = tp_dot(D,B.Pu), r2 = tp_dot(D,B.Pv);
      double det = a11*a22 - a12*a12;
      if (fabs(det) < 1.0e-300) break;
      double ds = ( a22*r1 - a12*r2)/det;
      double dt = (-a12*r1 + a11*r2)/det;
      ds = std::max(-adMaxStep[2], std::min(adMaxStep[2], ds));
      dt = std::max(-adMaxStep[3], std::min(adMaxStep[3], dt));
      SmPoint2d sNew = crDomB.ClampPoint2d(SmPoint2d(adPar[2]+ds, adPar[3]+dt));
      double dMove = fabs(sNew.x-adPar[2])/std::max(sSizeB.x,SM_EFF_ZERO) + fabs(sNew.y-adPar[3])/std::max(sSizeB.y,SM_EFF_ZERO);
      adPar[2] = sNew.x; adPar[3] = sNew.y;
      if (dMove < 1.0e-12) break;
    }

  for (int it=0; it<25; it++)
    {
      SmTpEval A, B;
      SER(tp_Evaluate(crSa, SmPoint2d(adPar[0],adPar[1]), A));
      SER(tp_Evaluate(crSb, SmPoint2d(adPar[2],adPar[3]), B));
      double D[3] = { A.P[0]-B.P[0], A.P[1]-B.P[1], A.P[2]-B.P[2] };
      double M[3], Ms[3], Mt[3], t1[3], t2[3];
      tp_cross(B.Pu, B.Pv, M);
      tp_cross(B.Puu, B.Pv, t1); tp_cross(B.Pu, B.Puv, t2);
      Ms[0] = t1[0]+t2[0]; Ms[1] = t1[1]+t2[1]; Ms[2] = t1[2]+t2[2];
      tp_cross(B.Puv, B.Pv, t1); tp_cross(B.Pu, B.Pvv, t2);
      Mt[0] = t1[0]+t2[0]; Mt[1] = t1[1]+t2[1]; Mt[2] = t1[2]+t2[2];

      double F[4] = { tp_dot(D,B.Pu), tp_dot(D,B.Pv), tp_dot(A.Pu,M), tp_dot(A.Pv,M) };
      double J[4][4] = {
        { tp_dot(A.Pu,B.Pu), tp_dot(A.Pv,B.Pu), -tp_dot(B.Pu,B.Pu)+tp_dot(D,B.Puu), -tp_dot(B.Pv,B.Pu)+tp_dot(D,B.Puv) },
        { tp_dot(A.Pu,B.Pv), tp_dot(A.Pv,B.Pv), -tp_dot(B.Pu,B.Pv)+tp_dot(D,B.Puv), -tp_dot(B.Pv,B.Pv)+tp_dot(D,B.Pvv) },
        { tp_dot(A.Puu,M),   tp_dot(A.Puv,M),    tp_dot(A.Pu,Ms),                     tp_dot(A.Pu,Mt) },
        { tp_dot(A.Puv,M),   tp_dot(A.Pvv,M),    tp_dot(A.Pv,Ms),                     tp_dot(A.Pv,Mt) } };
      double b[4] = { -F[0], -F[1], -F[2], -F[3] };
      double dx[4];
      if (!tp_Solve4(J, b, dx)) return SM_SUCCESS;

      double dRel = 0.0;
      for (int k=0; k<4; k++)
        {
          dx[k] = std::max(-adMaxStep[k], std::min(adMaxStep[k], dx[k]));
          double dSize = (k<2) ? (k==0 ? sSizeA.x : sSizeA.y) : (k==2 ? sSizeB.x : sSizeB.y);
          dRel = std::max(dRel, fabs(dx[k])/std::max(dSize,SM_EFF_ZERO));
        }
      SmPoint2d sA = crDomA.ClampPoint2d(SmPoint2d(adPar[0]+dx[0], adPar[1]+dx[1]));
      SmPoint2d sB = crDomB.ClampPoint2d(SmPoint2d(adPar[2]+dx[2], adPar[3]+dx[3]));
      adPar[0] = sA.x; adPar[1] = sA.y; adPar[2] = sB.x; adPar[3] = sB.y;
      if (dRel < 1.0e-13) { rbConverged = TRUE; break; }
      if (!sLocA.ContainsPoint2d(sA) || !sLocB.ContainsPoint2d(sB)) return SM_SUCCESS;
    }
  if (!rbConverged) return SM_SUCCESS;

  rUVa = SmPoint2d(adPar[0],adPar[1]);
  rUVb = SmPoint2d(adPar[2],adPar[3]);
  SmPoint3d sPa, sPb;
  SER(crSa.EvaluatePoint(rUVa, sPa));
  SER(crSb.EvaluatePoint(rUVb, sPb));
  rdGap = sPa.DistanceBetween(sPb);
  return SM_SUCCESS;
}

} // end anonymous namespace

/*******************************************************************//**
PURPOSE: SolidCraft change: find isolated tangential contact points in the
   interior of both surfaces, seeded by the node pairs FindInteriorCurves()
   could not trace (near parallel normals / no local intersection), and add
   each one as a degenerate SM_TC_TOUCHING curve.

NOTES: A point is reported only when Newton converges to a parameter pair
   strictly inside both domains, the surfaces meet there within the
   approximation tolerance, their normals are parallel, and the point is not
   on an intersection curve already found.
***********************************************************************/
static SmStatus sm_FindTangentContactPoints
  (const SmContext          & crContext,
   const SmSurface          & crSrf1,
   const SmExtent2d         & crDom1,
   const SmSurface          & crSrf2,
   const SmExtent2d         & crDom2,
   double                     dTol3d,
   const SmTArray<SmExtent2d> & crSeedDom1,
   const SmTArray<SmExtent2d> & crSeedDom2,
   const SmSurfaceIntersector & crSSI,
   SmTArray<SmCurve*>         & r3DCurves,
   SmTArray<SmCurve*>         & rSurface1UVCurves,
   SmTArray<SmCurve*>         & rSurface2UVCurves,
   SmTArray<SmTsectCurveType> & rCurveTypes,
   SmTArray<double>           & rDeviations)
{
  SmTArray<SmPoint3d> sTried;   // contact candidates already examined (accepted or not)
  ULONG lNumCurvesIn = r3DCurves.GetSize();
  // A true contact converges to a numerically zero gap.  A small nonzero gap at
  // a parallel-normal point means a near miss or a small intersection loop of
  // radius ~sqrt(gap/curvature), which must not be reported as a point; the
  // gap limit keeps such loops below the approximation tolerance.
  const double dGapTol = 1.0e-7*dTol3d;
  auto GapOk = [&](double dGap, const SmPoint3d & sP) -> SmBoolean
    {
      double dScale = 1.0 + fabs(sP.x) + fabs(sP.y) + fabs(sP.z);
      return dGap <= std::max(dGapTol, 1.0e-14*dScale) ? TRUE : FALSE;
    };

  // Accept the contact candidate (sUV1,sUV2) when the normals are parallel and it
  // is new; add it as a degenerate SM_TC_TOUCHING curve.
  auto Accept = [&](const SmPoint2d & sUV1, const SmPoint2d & sUV2, double dGap) -> SmStatus
    {
      SmPoint3d sP1, sP2;
      SER(crSrf1.EvaluatePoint(sUV1, sP1));
      SER(crSrf2.EvaluatePoint(sUV2, sP2));
      SmPoint3d sPnt = (sP1 + sP2) / 2.0;
      if (!GapOk(dGap, sPnt)) return SM_SUCCESS;

      for (ULONG k=0; k<sTried.GetSize(); k++)
        { if (sTried[k].DistanceBetween(sPnt) < 100.0*dTol3d) return SM_SUCCESS; }
      sTried.Add(sPnt);

      // normals must be parallel (EvaluateNormal also handles poles)
      SmVector3d sN1, sN2;
      SER(crSrf1.EvaluateNormal(sUV1, TRUE, TRUE, sN1));
      SER(crSrf2.EvaluateNormal(sUV2, TRUE, TRUE, sN2));
      if (sN1.LengthSquared() < SM_EFF_ZERO_SQ || sN2.LengthSquared() < SM_EFF_ZERO_SQ) return SM_SUCCESS;
      double dAng;
      SER(sN1.AngleBetween(sN2, dAng));
      if (dAng > SM_PI/2.0) dAng = SM_PI - dAng;
      if (dAng > 1.0e-4) return SM_SUCCESS;

      // skip points on curves already found
      if (crSSI.IsPointOnCurve(sPnt, r3DCurves)) return SM_SUCCESS;
      for (ULONG k=0; k<lNumCurvesIn; k++)
        {
          SmCurve * pCurve = r3DCurves[k];
          SmExtent1d sIvl = pCurve->GetNaturalInterval();
          SmBoolean bOk = FALSE; double dPar, dDist;
          SER(pCurve->DropPoint(sIvl, sPnt, NULL, 100.0*dTol3d, NULL, bOk, dPar, dDist, SM_SO_INTERSECT));
          if (bOk && dDist < 100.0*dTol3d) return SM_SUCCESS;
        }

      SmBSplineCurve * pNew = NULL;
      SER(SmBSplineCurve::CreateDegenerateCurve(crContext, 3, sPnt, pNew));
      r3DCurves.Add(pNew);
      pNew = NULL;
      SER(SmBSplineCurve::CreateDegenerateCurve(crContext, 2, SmPoint3d(sUV1.x, sUV1.y, 0.0), pNew));
      rSurface1UVCurves.Add(pNew);
      pNew = NULL;
      SER(SmBSplineCurve::CreateDegenerateCurve(crContext, 2, SmPoint3d(sUV2.x, sUV2.y, 0.0), pNew));
      rSurface2UVCurves.Add(pNew);
      rDeviations.Add(dGap);
      rCurveTypes.Add(SM_TC_TOUCHING);
      return SM_SUCCESS;
    };

  // Poles: the Newton system above degenerates where a surface's partials
  // vanish, so test each collapsed domain edge (pole) of either surface
  // directly: a pole lying on the other surface with parallel normals is a
  // contact point (e.g. two surfaces of revolution sharing an axis point).
  for (int iSrf=0; iSrf<2; iSrf++)
    {
      const SmSurface  & rSrf  = iSrf == 0 ? crSrf1 : crSrf2;
      const SmExtent2d & rDom  = iSrf == 0 ? crDom1 : crDom2;
      const SmSurface  & rOth  = iSrf == 0 ? crSrf2 : crSrf1;
      const SmExtent2d & rODom = iSrf == 0 ? crDom2 : crDom1;
      for (int iEdge=0; iEdge<4; iEdge++)
        {
          SmPoint2d aUV[3];
          for (int k=0; k<3; k++)
            {
              double f = 0.25*(k+1);
              aUV[k] = iEdge == 0 ? rDom.Evaluate(f, 0.0) : iEdge == 1 ? rDom.Evaluate(f, 1.0)
                     : iEdge == 2 ? rDom.Evaluate(0.0, f) : rDom.Evaluate(1.0, f);
            }
          SmPoint3d aP[3];
          for (int k=0; k<3; k++) SER(rSrf.EvaluatePoint(aUV[k], aP[k]));
          if (   !GapOk(aP[0].DistanceBetween(aP[1]), aP[1])
              || !GapOk(aP[1].DistanceBetween(aP[2]), aP[1])) continue;   // not a pole
          SmBoolean bOk = FALSE, bMulti = FALSE;
          SmPoint2d sOthUV;
          double dGap = SM_BIG_DOUBLE;
          SER(rOth.DropPoint(aP[1], rODom, NULL, bOk, sOthUV, dGap, bMulti));
          if (!bOk) continue;
          if (iSrf == 0) { SER(Accept(aUV[1], sOthUV, dGap)); }
          else           { SER(Accept(sOthUV, aUV[1], dGap)); }
        }
    }

  for (ULONG i=0; i<crSeedDom1.GetSize(); i++)
    {
      const SmExtent2d & rD1 = crSeedDom1[i];
      const SmExtent2d & rD2 = crSeedDom2[i];
      static const double sadFr[5][2] = { {0.5,0.5}, {0.2,0.2}, {0.8,0.2}, {0.2,0.8}, {0.8,0.8} };
      for (int iSeed=0; iSeed<5; iSeed++)
      for (int iRole=0; iRole<2; iRole++)
        {
          // seed on one patch, start the other at its patch center
          SmPoint2d sUV1, sUV2;
          SmBoolean bConverged = FALSE;
          double dGap = SM_BIG_DOUBLE;
          if (iRole == 0)
            {
              sUV1 = rD1.Evaluate(sadFr[iSeed][0], sadFr[iSeed][1]);
              sUV2 = rD2.Evaluate(0.5, 0.5);
              SER(tp_Newton(crSrf1, crDom1, crSrf2, crDom2, rD1, rD2, sUV1, sUV2, bConverged, dGap));
            }
          else
            {
              sUV2 = rD2.Evaluate(sadFr[iSeed][0], sadFr[iSeed][1]);
              sUV1 = rD1.Evaluate(0.5, 0.5);
              SER(tp_Newton(crSrf2, crDom2, crSrf1, crDom1, rD2, rD1, sUV2, sUV1, bConverged, dGap));
            }
          if (!bConverged) continue;
          SER(Accept(sUV1, sUV2, dGap));
        }
    }
  return SM_SUCCESS;
}

/*******************************************************************//**
PURPOSE: Construct a surface intersector object.

NOTES: 
***********************************************************************/
SmAdvSurfaceIntersector::SmAdvSurfaceIntersector
 (const SmSurface  & crSurface1, 
  const SmExtent2d & crUVDomain1,
  const SmSurface  & crSurface2,
  const SmExtent2d & crUVDomain2,
  SmBoolean  bFromFilletIntersector)      // in : TRUE = Constructing a FilletIntersector object
                                          //      FALSE= not
                                          //      default:[FALSE]
: SmSurfaceIntersector(crSurface1,
                       crUVDomain1,
                       crSurface2,
                       crUVDomain2,
                       bFromFilletIntersector),
  m_dClosestValue(SM_BIG_DOUBLE),
  m_bSeperateGaussMaps(TRUE)
{
  m_apNodes[0] = NULL ;
  m_apNodes[1] = NULL ;

} // end SmAdvSurfaceIntersector::SmAdvSurfaceIntersector constructor

/*******************************************************************//**
PURPOSE: add vector crSingVec to array aSingDirections when it
   is not already within the array and increment rlSingCount.

NOTES: 
***********************************************************************/
static void sm_AddSingDir
  (const SmVector3d & crSingVec,       // in : target vector to add to aSingDirections
   ULONG            & rlSingCount,     // i/o: in/out size of array aSingDirections
   SmVector3d       * aSingDirections) // i/o: array of singular direction vectors
                                       //      to be augmented with 
{
  ULONG ii ;
  SmBoolean bFound = FALSE;

  // see if crSingVec is already in aSingDirections array
  for(ii=0;ii<rlSingCount;ii++) 
    { if (aSingDirections[ii].Dot(crSingVec) > 0.9) 
        { bFound = TRUE ;
          break ;
        }
    }

  // when crSingVec is not in aSingDirections - add it
  if (!bFound) 
    {
      // add crSingVec to aSingDirections and increment counter
      aSingDirections[rlSingCount] = crSingVec ;
      rlSingCount++ ;
    }

} // end static void sm_AddSingDir

/*******************************************************************//**
PURPOSE: Compute TsectPnt->3DCurve->Tangent direction when TsectPnt is at
  a singular point (the two intersecting surfaces are parallel at this point)
  store the result in 
        rTsectPnt.CrvDeriv     = tangent direction to pursue.
    and rTsectPnt.m_ePointType = oneof SM_IP_TANGENT_POINT, - surfaces touch at point                  
                                       SM_IP_TANGENT_CURVE, - point on tangent curve                   
                                       SM_IP_SINGULARITY,   - point at nexus of 4 or 6 intersect curves
                                       SM_IP_COINCIDENCE    - matching 2nd fundamental forms           

NOTES: 
  Assumes: rTsectPnt Surface values are set 
              including surfaceNormal  in rTsectPnt.SrfNorm
                        1st derivs     in rTsectPnt.SrfDu, SrfDv
                        2nd derivs     in rTsectPnt.SrfDuu, SrfDuv, SrfDvv

  When two surfaces intersect at a point where the surfaces are parallel
  to one another (their surface normals are parallel), the intersection
  between the surfaces in the neighborhood of the point 
  fall into one of a few categories and can be distinguished by a study of the
  surface curvatures at the point.  The cases include:

  SM_IP_TANGENT_POINT: The surfaces kiss at this point and intersect at no other nearby points
                       (example: a hill touching a plane)
  SM_IP_TANGENT_CURVE: Point is on a tangent line and has 2 directions of continued intersections
                       in opposite directions
                       (example: a cylinder touching a plane)
  SM_IP_COINCIDENCE  : The two surfaces have the same shape in the neighbor hood of
                       this point and there are intersection points near this point
                       in all directions.
                       (examples: two coincident planes, two coincident hills, two coincident saddles)
  SM_IP_SINGULARITY  : the surfaces are tangent at only this point but has 4 or 6 distinct
                       directions in which more intersecting points can be found
                       (examples: a sadle intersecting a plane and a torus tangent to a cylinder along a circle)
***********************************************************************/
SmStatus SmAdvSurfaceIntersector::ComputeSingularityPoint
  (SmPoint2d [2],                // NOT USED: in : surf1 and surf2 xSectPoint UV values
   SmTsectPnt & rTsectPnt,       // out: xSectPoint Type, 
                                 //                 3DCurve position & tangent, 
                                 //                 3DCurve position & tangent projected to each surface,
                                 //             and Surface position & derivatives  
   SmTsectPnt *pOptPreviousPnt)  // in : Last intersection point on curve being stepped out, NULL to ignore
                                 //      When supplied, used to handle singularity cases.
                                 //      NOTE: only one vector in this is used: CrvDeriv(), the
                                 //      first derivative of the intersection curve; if everything
                                 //      else is unset, that's ok.  

{

// macro to flip output CUrveDir when it happens to be off by 180 degrees
#define ALIGN_DIR_WITH_PREVIOUS_POINT \
  if(pPrev)                           \
   { rTsectPnt.CrvDeriv().AngleBetween(pPrev->CrvDeriv(), dAngleRadians) ;                  \
     if(smos_Fabs(dAngleRadians - SM_PI) < smos_Min( SM_PI / 30.0, 5.0 * m_dThisAngTolRad)) \
       { rTsectPnt.CrvDeriv() = -rTsectPnt.CrvDeriv() ;                                     \
       }                                                                                    \
   }

  // init curve direction when possible - NULL when given a start point
  if(pOptPreviousPnt)
    {
      rTsectPnt.CrvDeriv() = pOptPreviousPnt->CrvDeriv() ;
    }

  // classify surface shape at given points:
  //   oneof: SM_LS_PLANAR,       - U and V isoParam Curves are planar at a point                             
  //          SM_LS_ELLIPTICAL,   - U and V isoParam Curves bend in same direction at a point (a bowl or hill)
  //          SM_LS_HYPERBOLIC,   - U and V isoparam Curves bend in different directions at a point (a saddle)
  //          SM_LS_CYLINDRICAL   - only one of U or V isopara Curve is planar at a point                     
  SmLocalSurfaceType eSurf1Type = smsurf_ComputeSurfacePointType(rTsectPnt.SrfNorm(0),
                                                                 rTsectPnt.SrfDuu(0),
                                                                 rTsectPnt.SrfDvv(0));
  SmLocalSurfaceType eSurf2Type = smsurf_ComputeSurfacePointType(rTsectPnt.SrfNorm(1),
                                                                 rTsectPnt.SrfDuu(1),
                                                                 rTsectPnt.SrfDvv(1));

  // get planar surface index (if any)
  ULONG lSrf =   (eSurf1Type == SM_LS_PLANAR) ? 0
               : (eSurf2Type == SM_LS_PLANAR) ? 1
               : 99 ;

  // locals for multiple intersection directions for singular points
  ULONG      ii, lSingDir = 0 ;
  SmVector3d aSingDirections[8] ; // we have experienced examples with 8 vectors, that is of course a bug - max should be 4

  // some cases project the previous point's CurveDir() to the current point's tangency plane.
  // Work in a scratch variable to prevent corrupting the previous point's CurveDir() value.
  double      dAngleRadians = 0.0;
  SmTsectPnt  sPreviousPnt ; 
  SmTsectPnt *pPrev =  (pOptPreviousPnt || m_vCurvePoints.GetLastNode())
                      ? &sPreviousPnt
                      :  NULL ;  

  // Project pPrev CrvDeriv() vector into the rTsectPnt tangency plane
  if (pPrev)  /* TRUE when (pOptPreviousPnt || m_vCurvePoints.GetLastNode()) */
    {
      // set sPreviousPnt values
      sPreviousPnt.CrvDeriv() =   pOptPreviousPnt 
                                ? pOptPreviousPnt->CrvDeriv()
                                : m_vCurvePoints.GetLastNode()->CrvDeriv() ;

      // project pPrev CrvDeriv() vector into rTsectPnt tangency plane
      SmPoint3d sEndPt = rTsectPnt.CrvPos() + pPrev->CrvDeriv();
      SmPoint3d sPtOnPlane;
      smgu_PointProjectToPlane(sEndPt,
                               rTsectPnt.CrvPos(),
                               rTsectPnt.SrfNorm(0),
                               sPtOnPlane);
      SmVector3d sDir = sPtOnPlane - rTsectPnt.CrvPos();

      // when projection worked - update the pPrev->CrvDeriv direction to that
      if ( sDir.Unitize() == SM_SUCCESS )
        { pPrev->CrvDeriv() = sDir; }

    } // end project pPrev CrvDeriv into this rTsectPnt Tangency plane

  // note: only when pPrev != NULL will we align output CrvDirev to previous Pnt directions

  // Process planar case -
  if(eSurf1Type == SM_LS_PLANAR || eSurf2Type == SM_LS_PLANAR)
    { 
      // get other surface index
      ULONG lOther = 1 - lSrf;

      // get otherSurf curvature properties at point
      double     dGauss, dNormal, dPrinK1, dPrinK2;
      SmVector3d sEFG, sLMN, sPrinKV1, sPrinKV2;
      //      SmVector3d aEvalMat[3][3];
      //      rTsectPnt.CopySrfEvalMatrix( lOther, aEvalMat ) ;

      SER(smsurf_EvaluateGeometric(// aEvalMat,              // in : parametric evaluation                            
                                   rTsectPnt.m_vSurfacePV[lOther],      // in : parametric evaluation
                                   dGauss, dNormal,                     // out: gaussian and mean curvatures                     
                                   dPrinK1, dPrinK2,                    // out: principal curvatures                             
                                   sEFG, sLMN,                          // out: 1st and 2nd fundamental forms                    
                                   sPrinKV1,sPrinKV2)) ;                // out: tangent vectors in principal curvature directions

      // If the normals do not agree then reverse the curvatures
      SmBoolean bOppositeNormals = FALSE;
      if(rTsectPnt.SrfNorm(lOther).Dot(rTsectPnt.SrfNorm(lSrf)) < SM_EFF_ZERO) 
        {
          bOppositeNormals = TRUE;
          sPrinKV1 = - sPrinKV1;
          sPrinKV2 = - sPrinKV2;
        }

      // case: planar surfacePoint is kissing a cylindrical surfacePoint in 1st principal curvature direction
      //   (2 possible directions for 3DCurve->Tangent)
      if (   smos_Fabs(dPrinK1) < SM_EFF_ZERO_SQRT 
          && smos_Fabs(dPrinK2) > SM_EFF_ZERO_SQRT) 
        {
          rTsectPnt.m_ePointType = SM_IP_TANGENT_CURVE;
          rTsectPnt.CrvDeriv()   = sPrinKV1 * m_dCurveTraceDirection;
          ALIGN_DIR_WITH_PREVIOUS_POINT ;
          return SM_SUCCESS;
        }

      // case: planar surfacePoint is kissing a cylindrical surfacePoint in 2nd principal curvature direction
      //   (2 possible directions for 3DCurve->Tangent)
      else if (   smos_Fabs(dPrinK2) < SM_EFF_ZERO_SQRT 
               && smos_Fabs(dPrinK1) > SM_EFF_ZERO_SQRT) 
        {
          rTsectPnt.m_ePointType = SM_IP_TANGENT_CURVE;
          rTsectPnt.CrvDeriv()   = sPrinKV2 * m_dCurveTraceDirection;
          ALIGN_DIR_WITH_PREVIOUS_POINT ;
          return SM_SUCCESS;
        }

      // case: planar surfacePoint is intersecting a saddle surfacePoint 
      //   (4 possible directions for 3DCurve->Tangent)
      else if (dPrinK1 * dPrinK2 < 0.0) 
        {
          // Here is a hyperbolic point and we need to add 
          // The point to the start point list if it is not
          // already there.

          // set PointType = SM_IP_SINGULARITY
          rTsectPnt.m_ePointType = SM_IP_SINGULARITY ;

          // Let us compute the vectors by evaluating the hyperbola
          // K1*x^2 + K2*y^2 = 1 at a large Y value and computing the
          // X.  This should give us a vector in terms of the principle
          // curvature vectors which defines one of the asymptotic 
          // directions of the hyperbola.
          double dY   = 10000.0 ;
          double dTmp = (1-dPrinK1*dY*dY) / dPrinK2 ;
          if (dTmp < 0.0) dTmp = 0.0 ;

          double dX   = smos_Sqrt( dTmp );

          SmVector3d sAsympDir1   = dY * sPrinKV1 + dX * sPrinKV2; 
          SmVector3d sAsympDirIn1 = dY * sPrinKV1 - dX * sPrinKV2;

          // If the normals oppose then switch directions
          if (bOppositeNormals) 
            {
              SmVector3d sTmp = sAsympDir1;
              sAsympDir1      = sAsympDirIn1;
              sAsympDirIn1    = sTmp;
            }

          // unitize the assymptotic directions
          sAsympDir1.Unitize() ;
          sAsympDirIn1.Unitize() ;

          // the other start directions are in the negative assymptotic directions
          SmVector3d sAsympDir2   = - sAsympDir1 ;
          SmVector3d sAsympDirIn2 = - sAsympDirIn1 ;

          // Reverse vectors if tracing backwards
          lSingDir = 4 ;
          aSingDirections[0] = sAsympDir1   * m_dCurveTraceDirection ;
          aSingDirections[1] = sAsympDir2   * m_dCurveTraceDirection ;
          aSingDirections[2] = sAsympDirIn1 * m_dCurveTraceDirection ;
          aSingDirections[3] = sAsympDirIn2 * m_dCurveTraceDirection ;

#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe = FALSE;
          if (bDebugMe) {
              smgfx_SetColor(1,0,0);
              SmPoint3d sSurfPnt = rTsectPnt.CrvPos();
              sSurfPnt.Draw();
              aSingDirections[0].Draw(&sSurfPnt);
              aSingDirections[1].Draw(&sSurfPnt);
              aSingDirections[2].Draw(&sSurfPnt);
              aSingDirections[3].Draw(&sSurfPnt);
              sm_GraphicsLoop();
              m_cpSurface[lOther]->DrawUV(20,10);
              sm_GraphicsLoop();
          }
#endif
        }
      else // case: plane is kissing a hill or bowl
        {
          // GWC: logic error that probably does not matter.
          //      The case of Plane/Plane is not being handled explicitly and ends up
          //      erroneously in this branch. If the case were explicit,
          //      I'd expect that branch return a SM_IP_COINCIDENCE which would
          //      cause the methods calling this one to cull this point from the 
          //      startPoints list, i.e. it would not cause an intersection 
          //      curve to be generated because those methods
          //      seem to ignore coincident area cases leaving those intersections
          //      to be found at a higher level by the coincidence checker.  
          //      The error made here of setting PointType for these cases to SM_IP_TANGENT_POINT
          //      might generate a degenerate point intersection case which
          //      would probably be found and removed as being part of the
          //      coincident area solution.  If that combining does not take
          //      place, then this error here might generate an extra point solution
          //      within a coincident region.  Basically a harmless extra solution
          //      included with the all the other solutions.
          // ToFix: Consider handling the Plane/Plane case explicitly.

          // Surfaces touch at a point but there are no local intersections
          // between the surfaces. 
          rTsectPnt.m_ePointType = SM_IP_TANGENT_POINT;
          return SM_SUCCESS;
       }
    } // end at least one planar point branch

  else // No planar points branch
    {
      SmTArray<SmConic*> sAboveCurves1, sBelowCurves1;
      SmTArray<SmConic*> sAboveCurves2, sBelowCurves2;
      SmObjsDelete<SmConic*> sClean1(&sAboveCurves1);
      SmObjsDelete<SmConic*> sClean2(&sAboveCurves2);
      SmObjsDelete<SmConic*> sClean3(&sBelowCurves1);
      SmObjsDelete<SmConic*> sClean4(&sBelowCurves2);
      
      // Dupin indicatrex curves for 1st surface
      double     dPrinK1 = 0.0,  dPrinK2 = 0.0;
      SmVector3d sPrinKV1, sPrinKV2, sEFG, sLMN ;

      //      SmVector3d aEvalMat[3][3];
      //      rTsectPnt.CopySrfEvalMatrix( 0, aEvalMat );
      SER(MakeDupin(// aEvalMat,                     // in : surfacePoint parametric evaluation matrix
                    rTsectPnt.m_vSurfacePV[0],    // in : surfacePoint parametric evaluation matrix
                    rTsectPnt.SrfNorm(0),         // in : reference SurfaceNormal
                    sAboveCurves1, sBelowCurves1, // out: Dupin Indicatrix curves (ellipses, lines, or hyperbolas)
                    sPrinKV1,      sPrinKV2,      // out: surface tangent vectors in principal curvature directions
                    dPrinK1,       dPrinK2,       // out: principal curvatures
                    sEFG,          sLMN)) ;       // out: first and 2nd fundamental forms
  
      // Compute values of second surface based upon the normal of the first
      double     dPrinK1_2 = 0.0,  dPrinK2_2 = 0.0;
      SmVector3d sPrinKV1_2, sPrinKV2_2, sEFG_2, sLMN_2 ;
      //      rTsectPnt.CopySrfEvalMatrix( 1, aEvalMat );
      SER(MakeDupin(// aEvalMat,                     // in : surfacePoint parametric evaluation matrix
                    rTsectPnt.m_vSurfacePV[1],    // in : surfacePoint parametric evaluation matrix
                    rTsectPnt.SrfNorm(0),         // in : reference SurfaceNormal                                  
                    sAboveCurves2, sBelowCurves2, // out: Dupin Indicatrix curves (ellipses, lines, or hyperbolas) 
                    sPrinKV1_2,    sPrinKV2_2,    // out: surface tangent vectors in principal curvature directions
                    dPrinK1_2,     dPrinK2_2,     // out: principal curvatures
                    sEFG_2,        sLMN_2)) ;     // out: first and 2nd fundamental forms                                    
      
      // if Surf1 has an Umbilical (Spherical) surfacePoint where all curvatures are the same
      //   if Surf2 has same curvature singular point either COINCIDENT or TANGENT_CURVE
      double dScale = SM_EFF_ZERO * (1.0 + smos_Fabs(dPrinK1));
      if (smos_Fabs(dPrinK1-dPrinK2) < dScale) 
        {
          // if Surf2->k1 = Surf1->k1
          if (smos_Fabs(dPrinK1-dPrinK1_2) < SM_EFF_ZERO_SQRT) 
            {
              // and Surf2  k2 is the same - coincident surfacePoints
              if (smos_Fabs(dPrinK2-dPrinK2_2) < SM_EFF_ZERO_SQRT) 
                {
                  rTsectPnt.m_ePointType = SM_IP_COINCIDENCE;
                  return SM_SUCCESS;
                }
              else // only tangent in k1 direction 
                {
                  rTsectPnt.m_ePointType = SM_IP_TANGENT_CURVE;
                  rTsectPnt.CrvDeriv()   = sPrinKV1_2 * m_dCurveTraceDirection;
                  ALIGN_DIR_WITH_PREVIOUS_POINT ;
                  return SM_SUCCESS;
                }
            }

          // else if Surf2->k2 = Surf1->k2 
          if (smos_Fabs(dPrinK2-dPrinK2_2) < SM_EFF_ZERO_SQRT) 
            {
              // only tangent in k2 direction 
              rTsectPnt.m_ePointType = SM_IP_TANGENT_CURVE;
              rTsectPnt.CrvDeriv()   = sPrinKV2_2 * m_dCurveTraceDirection;
              ALIGN_DIR_WITH_PREVIOUS_POINT ;
              return SM_SUCCESS;
            }            
        } // end Surf1 Umbilical (spherical) surfacePoint check

      // arrive here if Surf1 is not Umbilical or Surf2 has different curvature values

      // if Surf2 has an Umbilical (Spherical) surfacePoint where all curvatures are the same
      //   if Surf1 has same curvature singular point is either COINCIDENT or TANGENT_CURVE
      dScale = SM_EFF_ZERO * (1.0 + smos_Fabs(dPrinK1_2));
      if (smos_Fabs(dPrinK1_2-dPrinK2_2) < dScale) 
        {
          // if Surf1->k1 == Surf2->k1 singular Point is either COINCIDENT or TANGENT_CURVE
          if (smos_Fabs(dPrinK1-dPrinK1_2) < SM_EFF_ZERO_SQRT) 
            {
              if (smos_Fabs(dPrinK2-dPrinK2_2) < SM_EFF_ZERO_SQRT) 
                {
                  rTsectPnt.m_ePointType = SM_IP_COINCIDENCE;
                  return SM_SUCCESS;
                }
              else 
                {
                  rTsectPnt.m_ePointType = SM_IP_TANGENT_CURVE;
                  rTsectPnt.CrvDeriv()   = sPrinKV1 * m_dCurveTraceDirection;
                  ALIGN_DIR_WITH_PREVIOUS_POINT ;
                  return SM_SUCCESS;
                }
            }
          
          // if Surf1->k2 == Surf2->k1 singular Point is TANGENT_CURVE
          if (smos_Fabs(dPrinK2-dPrinK2_2) < SM_EFF_ZERO_SQRT) 
            {
              rTsectPnt.m_ePointType = SM_IP_TANGENT_CURVE;
              rTsectPnt.CrvDeriv()   = sPrinKV2 * m_dCurveTraceDirection;
              ALIGN_DIR_WITH_PREVIOUS_POINT ;
              return SM_SUCCESS;
            }            
      } // end Surf2 Umbilical (spherical) surfacePoint check
      
      // arrive here when Surf1 and Surf2 SurfacePoints are not Spherical or are
      // Spherical without a matching curvature value at the other surfacePoint.

      // Second Do a Curvature Check to see if any principle curvatures
      // match in both direction and orientation.

      // when surf1->k1_Direction = surf2->k1_Direction
      SmVector3d sCross = sPrinKV1 * sPrinKV1_2;
      if (sCross.Length() < SM_EFF_ZERO_SQRT*100.0) 
        {
          // when surf1->k1 == surf2->k1
          if (smos_Fabs(dPrinK1-dPrinK1_2) < SM_EFF_ZERO_SQRT) 
            {
              // and surf1->k2 == surf2->k2
              if (smos_Fabs(dPrinK2-dPrinK2_2) < SM_EFF_ZERO_SQRT) 
                {
                  lSingDir = 4;
                  aSingDirections[0] =  sPrinKV1;
                  aSingDirections[1] = -sPrinKV1;
                  aSingDirections[2] =  sPrinKV2;
                  aSingDirections[3] = -sPrinKV2;
                  rTsectPnt.CrvDeriv() = aSingDirections[0];
                  double dAngleMin=0.0;

                  //      if (pPrev)  // The direction must be in the tangent plane of this rTSectPnt
                  //        {
                  //          SmPoint3d sEndPt = rTsectPnt.CrvPos() + pPrev->CrvDeriv();
                  //          SmPoint3d sPtOnPlane;
                  //          smgu_PointProjectToPlane(sEndPt,
                  //                                   rTsectPnt.CrvPos(),
                  //                                   rTsectPnt.SrfNorm(0),
                  //                                   sPtOnPlane);
                  //          SmVector3d sDir = sPtOnPlane - rTsectPnt.CrvPos();
                  //          SER(sDir.Unitize());
                  //          pPrev->CrvDeriv() = sDir;
                  //        }
                  //      if (pPrev == NULL) pPrev = m_vCurvePoints.GetLastNode();

                  // when there is a prevPoint or a LastNode of m_vCurvePoints
                  if (pPrev) 
                    {
                      // select the direction most like the previous point direction
                      // this includes checks of both +/- sPrinKV1 and +/- sPrinKV2
                      SER(pPrev->CrvDeriv().AngleBetween(aSingDirections[0],dAngleMin));
                      for (ULONG kk=1; kk<lSingDir; kk++) 
                        {
                          double dAngleTest = 0.0;
                          SER(pPrev->CrvDeriv().AngleBetween(aSingDirections[kk],dAngleTest));
                          if (dAngleTest < dAngleMin) 
                            {
                              dAngleMin = dAngleTest;
                              rTsectPnt.CrvDeriv() = aSingDirections[kk];
                            }
                        }
                      rTsectPnt.m_ePointType = SM_IP_SINGULARITY;
                    }
                  else // no PreviousPoint to help select a direction
                    {
                      // select the direction with the smallest radius of curvature
                      rTsectPnt.CrvDeriv() =   (smos_Fabs(dPrinK2) < smos_Fabs(dPrinK1))
                                             ? sPrinKV2 * m_dCurveTraceDirection
                                             : sPrinKV1 * m_dCurveTraceDirection ;
                      rTsectPnt.m_ePointType = SM_IP_COINCIDENCE;
                    }
                  
                  return SM_SUCCESS;
                } // end surf1->k2 == surf2->k2 branch
              else // surf1->k2 != surf2->k2 branch
                {
                  rTsectPnt.m_ePointType = SM_IP_TANGENT_CURVE;
                  rTsectPnt.CrvDeriv()   = sPrinKV1 * m_dCurveTraceDirection;
                  ALIGN_DIR_WITH_PREVIOUS_POINT ;
                  return SM_SUCCESS;
                }
            } // end surf1->k1 == surf2->k1 check

          // when surf1->k2 == surf2->k2 while surf1->k1 != surf2->k1
          if (smos_Fabs(dPrinK2-dPrinK2_2) < SM_EFF_ZERO_SQRT) 
            {
              rTsectPnt.m_ePointType = SM_IP_TANGENT_CURVE;
              rTsectPnt.CrvDeriv()   = sPrinKV2 * m_dCurveTraceDirection;
              ALIGN_DIR_WITH_PREVIOUS_POINT ;
              return SM_SUCCESS;
            } // end surf1->k2 == surf2->k2 check
        } // end surf1->k1_Direction = surf2->k1_Direction check
      
      // when surf1->k1_Direction = +/- surf2->k2_Direction
      SmVector3d sCross2 = sPrinKV1 * sPrinKV2_2;
      if (sCross2.Length() < SM_EFF_ZERO_SQRT*100.0) 
        {
          // when surf1->k1 == surf2->k2
          if (smos_Fabs(dPrinK1-dPrinK2_2) < SM_EFF_ZERO_SQRT) 
            {
              // and when surf1->k2 == durf2->k1
              if (smos_Fabs(dPrinK2-dPrinK1_2) < SM_EFF_ZERO_SQRT) 
                {
                  lSingDir = 4;
                  aSingDirections[0] =  sPrinKV1;
                  aSingDirections[1] = -sPrinKV1;
                  aSingDirections[2] =  sPrinKV2;
                  aSingDirections[3] = -sPrinKV2;
                  rTsectPnt.CrvDeriv() = aSingDirections[0];
                  double dAngleMin=0.0;

                  // when there is a prevPoint or a LastNode of m_vCurvePoints
                  if (pPrev) 
                    {
                      // select the direction most like the previous point direction
                      // this includes checks of both +/- sPrinKV1 and +/- sPrinKV2
                      SER(pPrev->CrvDeriv().AngleBetween(aSingDirections[0],dAngleMin));
                      for (ULONG kk=1; kk<lSingDir; kk++) 
                        {
                          double dAngleTest = 0.0;
                          SER(pPrev->CrvDeriv().AngleBetween(aSingDirections[kk],dAngleTest));
                          if (dAngleTest < dAngleMin) 
                            {
                              dAngleMin = dAngleTest;
                              rTsectPnt.CrvDeriv() = aSingDirections[kk];
                            }
                        }
                      rTsectPnt.m_ePointType = SM_IP_SINGULARITY;
                    }
                  else // no PreviousPoint to help select a direction
                    {
                      // select the direction with the smallest radius of curvature
                      rTsectPnt.CrvDeriv() =   (smos_Fabs(dPrinK2) < smos_Fabs(dPrinK1))
                                             ? sPrinKV2 * m_dCurveTraceDirection
                                             : sPrinKV1 * m_dCurveTraceDirection ;
                      rTsectPnt.m_ePointType = SM_IP_COINCIDENCE;
                    }

                  return SM_SUCCESS;
                } // end surf1->k2 == surf2->k1 branch
              else // surf1->k2 != durf2->k1 while surf1->k1 == surf2->k2
                {
                  rTsectPnt.m_ePointType = SM_IP_TANGENT_CURVE;
                  rTsectPnt.CrvDeriv()   = sPrinKV1 * m_dCurveTraceDirection;
                  ALIGN_DIR_WITH_PREVIOUS_POINT ;
                  return SM_SUCCESS;
                }
            } // end surf1->k1 == surf2->k2 check

          // when surf1->k2 == surf2->k1 while surf1->k1 != surf2->k2
          if (smos_Fabs(dPrinK2-dPrinK1_2) < SM_EFF_ZERO_SQRT) 
            {
              rTsectPnt.m_ePointType = SM_IP_TANGENT_CURVE;
              rTsectPnt.CrvDeriv()   = sPrinKV2 * m_dCurveTraceDirection;
              ALIGN_DIR_WITH_PREVIOUS_POINT ;
              return SM_SUCCESS;
            }            
        } // end surf1->k1_Direction = surf2->k2_Direction check
      
      // arrive here when surfaces do not have matching principle curvatures.
      // Next: intersect properly oriented surface Dupin indicatrix curves
      //       to estimate the behavior of the intersection curves radiating
      //       from this nexus point. These intersections should accurately
      //       predict the number of intersections and estimate their directions
      //       as they leave the singularity point.
      // Build: lSingDir        = number of xSect curves radiating from this nexus point
      //        aSingDirections = the directions for each xSect curve radiating from this nexus point (4 for planar/hyperbolic case)

      // angle between surface's natural coordinate sysstems relative to one another.      
      double dAngle = 0.0;
      SER( rTsectPnt.SrfNorm(0).CCWAngleBetween( sPrinKV1, sPrinKV1_2, dAngle));
      
      // Rotate curves from surface 2 into position on surface 1
      //   MakeDupin creates the Dupin indicatrix curves in local coordinates.
      //   Rotate the Surface2 Dupin curves to the relative orientation between surfs 1 and 2.
      for(ii=0;ii<sAboveCurves2.GetSize();ii++) { SER(sAboveCurves2[ii]->Transform2D(dAngle,SmPoint2d(0,0))) ; }
      for(ii=0;ii<sBelowCurves2.GetSize();ii++) { SER(sBelowCurves2[ii]->Transform2D(dAngle,SmPoint2d(0,0))) ; }
      
      // for every surf1->AboveCurves - intersect with surf2->AboveCurve
      SmSolutionArray sSols;
      { for(ii=0;ii<sAboveCurves1.GetSize();ii++)
          {
            // for every surf2->AboveCurves
            for(ULONG jj=0;jj<sAboveCurves2.GetSize();jj++)
              {
                // intersect the dupin indicatrix curves
                SER(sAboveCurves1[ii]->IntersectConic(sAboveCurves1[ii]->GetNaturalInterval(),
                                                     *sAboveCurves2[jj],
                                                      sAboveCurves2[jj]->GetNaturalInterval(),
                                                      SM_EFF_ZERO_SQRT*100.0,
                                                      sSols));
              
                // for every dupin indicatrix curve/curve xSect solution
                for(ULONG kk=0;kk<sSols.GetSize();kk++)
                  {
                    // build the surfaceTangent direction from the curve/curve xSect solution
                    double dT = sSols[kk].m_vStart[0] ;

                    // Get xSect Point in Surface1 natural coordinate system
                    SmPoint3d sP1 ;
                    SER(sAboveCurves1[ii]->EvaluatePoint(dT,sP1)) ;
                    SER(sP1.Unitize()) ;

                    // Get 3D direction for vector = [sP1 - NaturalOrigin] = sP1 
                    SmVector3d sSingDir = sP1.x * sPrinKV1 + sP1.y * sPrinKV2 ;
                    SER(sSingDir.Unitize()) ;

                    // add unique sSingDir vector to aSingDirections 
                    sm_AddSingDir(sSingDir, lSingDir, aSingDirections) ;

                  } // end iter every Dupin Indicatrix curve/curve xSect solution
              } // end iter every surf2->AboveCurves
          } // end iter every surf1->AboveCurves 
      } // end find AboveCurve Dupin Indicatrix curve/curve intersections block

      // for every surf1->BelowCurves - intersect with Surf2 Below Curve 
      { for(ii=0;ii<sBelowCurves1.GetSize();ii++)
          {
            // for every surf2->AboveCurves
            for(ULONG jj=0;jj<sBelowCurves2.GetSize();jj++)
              {
                // intersect the dupin indicatrix curves
                SER(sBelowCurves1[ii]->IntersectConic(sBelowCurves1[ii]->GetNaturalInterval(),
                                                     *sBelowCurves2[jj],
                                                      sBelowCurves2[jj]->GetNaturalInterval(),
                                                      SM_EFF_ZERO_SQRT*100.0,
                                                      sSols)) ;
              
                // for every dupin indicatrix curve/curve xSect solution
                for(ULONG kk=0;kk<sSols.GetSize();kk++)
                  {
                    // build the surfaceTangent direction from the curve/curve xSect solution
                    double dT = sSols[kk].m_vStart[0] ;
                    
                    // Get xSect Point in Surface1 natural coordinate system
                    SmPoint3d sP1 ;
                    SER(sBelowCurves1[ii]->EvaluatePoint(dT,sP1)) ;
                    SER(sP1.Unitize()) ;
                    
                    // Get 3D direction for vector = [sP1 - NaturalOrigin] = sP1 
                    SmVector3d sSingDir = sP1.x * sPrinKV1 + sP1.y * sPrinKV2 ;
                    SER(sSingDir.Unitize()) ;
                    
                    // add unique sSingDir vector to aSingDirections 
                    sm_AddSingDir(sSingDir,lSingDir,aSingDirections) ;

                  } // end iter every Dupin Indicatrix curve/curve xSect solution
              } // end iter every surf2->BelowCurves
          } // end iter every surf1->BelowCurves 
      } // end find BelowCurve Dupin Indicatrix curve/curve intersections block
    } // end No planar points branch

  // arrive here for cases not already handled: 
  //     1.) all             planar/hyperbolic(saddle) surf pairs
  //     2.) nonSpecial case nonPlanar/nonPlanar       surf pairs
  // with:
  //   lSingDir        = number of xSect curves radiating from this nexus point
  //   aSingDirections = the directions for each xSect curve radiating from this nexus point (4 for planar/hyperbolic case)
  //   pPrev           = contains PreviousPnt.CrvDeriv() when pOptPreviousPnt != NULL or m_vCurvePoints has a LastNode,
  //                     else NULL when no previous point is available.
  //                     When (pPrev == NULL) we are processing a xSectCurve Start Point.

  // next: 1. when lSindDir == 0 SM_IP_TANGENT_POINT,
  //       2. when lSingDir == 2 SM_IP_TANGENT_CURVE, (two sides of a tangentCurve xSect) pick best direction, and return
  //       3. when lSingDir  > 2 SM_IP_SINGULARITY, cull all starting directions that point into both surfaces.
  //          3a. for through points (pPrev != NULL) pick best continuing direction from remaining curves
  //          3b. for start points   (pPrev == NULL) pick best start direction
  //              GWC: the 3b. part of the algorithm is weak - it should save all possible start directions in some manner. 

  // Special cases which return before arriving here
  // set rTsectPnt.m_ePointType and rTSectPnt.CrvDeriv() and return SM_SUCCESS
  //       rTsectPnt.m_ePointType set for cases planar/cylindrical                    SM_IP_TANGENT_CURVE, 
  //                                            planar/elliptical(bowl)               SM_IP_TANGENT_POINT,
  //                                            planar/planar                         SM_IP_TANGENT_POINT (gwc:should it be SM_IP_COINCIDENT?)
  //                                            sphere/matching Sphere                SM_IP_COINCIDENT
  //                                            Sphere/a matching principle curvature SM_IS_TANGENT_CURVE
  //              matching principle curvatures k1_1_dir == k1_2_dir 
  //                                        and k1_1 == k1_2 and k2_1 == k2_2         SM_IP_SINGULARITY (with pPrev)
  //                                                                                  SM_IP_COINCIDENCE (without pPrev)
  //                                        or  k1_1 == k1_2 and k2_1 != k2_2         SM_IP_TANGENT_CURVE
  //                                        or  k1_1 != k1_2 and k2_1 == k2_2         SM_IP_TANGENT_CURVE
  //              matching principle curvatures k1_1_dir == k2_2_dir 
  //                                        and k1_1 == k2_2 and k2_1 == k1_2         SM_IP_SINGULARITY (with pPrev)
  //                                                                                  SM_IP_COINCIDENCE (without pPrev)
  //                                        or  k1_1 == k2_2 and k2_1 != k1_2         SM_IP_TANGENT_CURVE
  //                                        or  k1_1 != k2_2 and k2_1 == k1_2         SM_IP_TANGENT_CURVE
  //
  // GWC: could get some nonPlanar/nonPlanar cases  
  //                            elliptical/elliptical (nested bowls)                   SM_IP_TANGENT_POINT
  //                            elliptical/elliptical bowl/hill configuration          SM_IP_TANGENT_POINT
  //                            elliptical/cylindrical (nested bowl in large cylinder) SM_IP_TANGENT_POINT
  //                            elliptical/cylinderical cylinder on hill               SM_IP_TANGENT_POINT
  //                            cylinder/cylinder axis aligned                         SM_IP_TANGENT_CURVE or SM_IP_COINCIDENT

  // no singular directions branch
  if (lSingDir == 0) 
    { 
      rTsectPnt.m_ePointType = SM_IP_TANGENT_POINT;
      return SM_SUCCESS;
    }  
  
  // when number of singular directions is 2
  else if (lSingDir == 2) 
    {
      rTsectPnt.m_ePointType = SM_IP_TANGENT_CURVE ;

      // set sing direct = 1st singDirection found
      rTsectPnt.CrvDeriv() = aSingDirections[0];

      // when there is a PreviousPoint - select best matching direction
      if (pPrev) 
        {
          // select singDirect = direction most parallel to previousPoint->3DCurve direction
          double dAngleMin = 0.0, dAngleMin2 = 0.0;
          SER(pPrev->CrvDeriv().AngleBetween(aSingDirections[0],dAngleMin));
          SER(pPrev->CrvDeriv().AngleBetween(aSingDirections[1],dAngleMin2));
          if (dAngleMin2 < dAngleMin) 
            {
              rTsectPnt.CrvDeriv() = aSingDirections[1];
            }
        }
      
      ALIGN_DIR_WITH_PREVIOUS_POINT ;
      return SM_SUCCESS;
    } // end 2 singular direction branch

  // 3 or more singular directions branch
  else if ( lSingDir > 2 )
    { 
      rTsectPnt.m_ePointType = SM_IP_SINGULARITY ;

      // locals
      SmBoolean bAtStartPoint = (m_vCurvePoints.GetLastNode() == NULL) ;
      
      // If we are starting a curve - use sAsympDir1 and put sAsympDir2 onto start point list.
      if ( bAtStartPoint )
        {
          // seek singularity direction which is inside of the domains.
          ULONG      lNumInside=0 ;
          ULONG      lInsideIndex[8] ;
          SmVector2d aSingUVSurf1[8], aSingUVSurf2[8]; // must be the same size as aSingDirections

          // drop SingDirectons into Surf1 parameter space
          SER(smsurf_DropVectors(rTsectPnt.SrfDu(0),   // in : Surface U_dir tangent                                 
                                 rTsectPnt.SrfDv(0),   // in : Surface V_dir tangent                                 
                                 lSingDir,             // in : Number of vectors to drop                             
                                 aSingDirections,      // in : vectors to drop, sized:[lNumVectorsToDrop]            
                                 aSingUVSurf1)) ;      // out: resulting UVVectors where                             
                                                       //        DropVec3D = u*rDU + v*rDV                           
                                                       //        DropVec3D = vector generated by droping input vec to
                                                       //                  rDU/rDV plane                             
                                                       //        when rDU or rDV are very close to zero - they are   
                                                       //         replaced by unitized vectors.                      
                                                       //      sized:[lNumVectorsToDrop]                             

          // drop SingDirections into Surf2 parameter space
          SER(smsurf_DropVectors(rTsectPnt.SrfDu(1),   // in : Surface U_dir tangent                     
                                 rTsectPnt.SrfDv(1),   // in : Surface V_dir tangent                     
                                 lSingDir,             // in : Number of vectors to drop                 
                                 aSingDirections,      // in : vectors to drop, sized:[lNumVectorsToDrop]
                                 aSingUVSurf2)) ;      // out: resulting UVVectors where                 

          // init output 3DCurveDirection = 1st SingDirection Vector
          rTsectPnt.CrvDeriv() = aSingDirections[0] ;

          // for every SingDirection vector - find those directions which point inside surf1 and surf2
          SmPoint2d sDelta1, sDelta2 ;
          for(ii=0;ii<lSingDir;ii++)
            {
              // when small step from xSectPoint in singDirection is inside Surf1
              sDelta1 = rTsectPnt.UVPos(0) + SM_EFF_ZERO_SQRT * aSingUVSurf1[ii];
              if (m_vUVDomain[0].ContainsPoint2d(sDelta1,SM_EFF_ZERO_SQRT/100000.0)) 
                {
                  // when small step from xSectPoint in singDirection is inside Surf2
                  sDelta2 = rTsectPnt.UVPos(1) + SM_EFF_ZERO_SQRT * aSingUVSurf2[ii];
                  if (m_vUVDomain[1].ContainsPoint2d(sDelta2,SM_EFF_ZERO_SQRT/100000.0)) 
                    {
                      // found a singular direction inside both domains
                      lInsideIndex[lNumInside] = ii;
                      lNumInside++;

                    } // end inside surf2 check
                } // end inside surf1 check
            } // end iter every singDirection vector

          // arrive here when:
          // lSingDir         = total  number of xSectCurve start directions from this nexus point
          // lNumInside       = culled number of XSectCurve start directions that point inside both Surf1 and Surf2 domains
          // lInsideIndex[ii] = index of iith inside pointing direction in arrays aSingDirections, aSingUVSurf1, and aSingUVSurf2 
          //                    sized:[lNumInside]
          // aSingDirections  = array of expected xSect Curve directions radiating from this nexus point, sized:[lSingDir]
          // aSingUVSurf1     = associated array of Surf1 UV directions, sized:[lSingDir] 
          // aSingUVSurf2     = associated array of Surf2 UV directions, sized:[lSingDir] 

          // Decide which of the aSingDirections to use.
          if ( pPrev != NULL ) 
            {
              if ( lNumInside > 0) 
                {
                  // Find the singDirection closest to Prev->3DCurve_Direction.
                  double dAngleMin = 0.0, dAngleTest = 0.0;
                  SER(pPrev->CrvDeriv().AngleBetween(aSingDirections[lInsideIndex[0]], dAngleMin)) ;
                  for(ii=1;ii<lNumInside;ii++)
                    {
                      SER(pPrev->CrvDeriv().AngleBetween(aSingDirections[lInsideIndex[ii]], dAngleTest)) ;
                      if (dAngleTest < dAngleMin)  
                        { 
                          rTsectPnt.CrvDeriv() = aSingDirections[lInsideIndex[ii]] ;
                          dAngleMin = dAngleTest ;

                        } // end found a keeper check
                    } // end iter every other singDirection
                } // end need to find the best singDirection check
              else
                {
                  // None points inside both domains.  Just use the prev direction.
                  rTsectPnt.CrvDeriv() = pPrev->CrvDeriv();
                }
            } // end pPrev not Null
          else
            {
              // No prev.
              if ( lNumInside == 1 )
                {
                  // Easy, just use it.  (pPrev must be Null.)
                  rTsectPnt.CrvDeriv() = aSingDirections[lInsideIndex[0]];
                }
              else if ( lNumInside > 1 )
                {
                  // Pick a best direction?  For now, just use the first one.
                  rTsectPnt.CrvDeriv() = aSingDirections[lInsideIndex[0]];

                  // GWC_NOTE: GWC_NEED_TO_REMEMBER_THAT_MULTIPLE_CURVES_CAN_START_FROM_THIS_POINT_HERE 
                }
              else
                {
                  // No direction was inside both.
                  // Just leave it as it is and hope for the best.
                  // We're probably at a surface boundary starting an xSect curve that leaves the domain.
                  // Expect to generate a degenerate point intersection if we come through this branch. 
                }
            } // end pPrev Null or not.

        } // end this is a StartPoint branch

      else // is NOT a startPoint branch 
        {
          // If we are continuing the curve, take the sXect Curve direction vector
          //  which is closest match for the last point direction vector.
          rTsectPnt.CrvDeriv() = aSingDirections[0];
          double dAngleMin = 0.0, dAngleTest = 0.0;
          SER(pPrev->CrvDeriv().AngleBetween(aSingDirections[0],dAngleMin));
          for (ULONG kk=1; kk<lSingDir; kk++) 
            {
              // direction test
              SER(pPrev->CrvDeriv().AngleBetween(aSingDirections[kk],dAngleTest));
              if (dAngleTest < dAngleMin) 
                {
                  dAngleMin = dAngleTest;
                  rTsectPnt.CrvDeriv() = aSingDirections[kk];
                }

              // neg direction test - not needed if aSingDirections already contains neg directions
              SER(pPrev->CrvDeriv().AngleBetween(-aSingDirections[kk],dAngleTest));
              if (dAngleTest < dAngleMin) 
                {
                  dAngleMin = dAngleTest;
                  rTsectPnt.CrvDeriv() = -aSingDirections[kk];
                }
            }

          // check for a likely problem direction - if we come into this branch
          //  we might be picking the wrong xSect direction going through this 
          //  singularity point.  GWC: To debug draw the surfaces and all the direction vectors
          //  the pPrev SmTSectPoint and highlight the angle to see what's going wrong.
          if (dAngleMin > 20.0 * SM_PI/180.0) 
            {
              // Usually this gets rejected by other tests at a higher level
              // No need to make an error here
//                SE(SM_ERR);  // Possible error if angle is too large (> 30 degrees)

#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe = FALSE ;
              // draw 
              if(bDebugMe)
                {
                  SmPoint3d sCrvPos    = rTsectPnt.CrvPos() ; 
                  const SmSurface *pSurface1 = m_cpSurface[0] ;
                  const SmSurface *pSurface2 = m_cpSurface[1] ;
                  SmFace    *pFace1    = (SmFace *)pSurface1->GetFace() ;
                  SmFace    *pFace2    = (SmFace *)pSurface2->GetFace() ;
                  SmBrep    *pBrep1    = pFace1 ? pFace1->GetBrep() : NULL ;
                  SmBrep    *pBrep2    = pFace2 ? pFace2->GetBrep() : NULL ;
                  SmTArray<SmTsectPnt*> sTsectPnts ;
                  m_vCurvePoints.GetAllNodes(sTsectPnts) ;

                  smgfx_Erase() ;
                  smgfx_SetLook(1,2, 0,0,1) ; if(pBrep1) pBrep1->Draw(TRUE) ; sm_GraphicsLoop() ;
                  smgfx_SetLook(1,2, 0,1,0) ; if(pBrep2) pBrep2->Draw(TRUE) ; sm_GraphicsLoop() ;
                  smgfx_SetLook(1,2, 0,0,1) ; if(pSurface1) pSurface1->DrawUV(8,8,FALSE,NULL,TRUE) ; sm_GraphicsLoop() ;
                  smgfx_SetLook(1,2, 0,1,0) ; if(pSurface2) pSurface2->DrawUV(8,8,FALSE,NULL,TRUE) ; sm_GraphicsLoop() ;
                  smgfx_SetLook(1,2, 0,0,0) ; if(pFace1) pFace1->Draw(SM_DM_CROSSHATCH) ; sm_GraphicsLoop() ;
                  smgfx_SetLook(1,2, 0,0,0) ; if(pFace2) pFace2->Draw(SM_DM_CROSSHATCH) ; sm_GraphicsLoop() ;

                  smgfx_SetLook(3,4, 1,0,0) ; if(pPrev) pPrev->Draw() ; sm_GraphicsLoop() ;
                  for(ii=0;ii<lSingDir;ii++)
                    { smgfx_SetLook(3,4, .2,0,1) ; sCrvPos.Draw(&aSingDirections[ii]) ; sm_GraphicsLoop() ; }

                  for(ii=0;ii<sTsectPnts.GetSize();ii++)
                    { smgfx_SetLook(3,4, 1,0,1) ; sTsectPnts[ii]->Draw() ; sm_GraphicsLoop() ; }
                  sm_GraphicsLoop() ;
                }
#endif // SM_DEBUG_CODE

            }
        } // end NOT a startPoint branch

      return SM_SUCCESS;

    }  // end singularity point with 3 or more singDirections branch
  else 
    {
      // an unexpected branch 
      // lSingDir will equal 1 which is not a valid result.
      // lSingDir counts the number of xSect curves that radiate from this intersection singularity point.
      // Those are expected to come in matching pairs for xSect curves that run
      //  through this point.  Some of those directions will point out of the boundaries
      //  of the target surfaces and will be culled, but that culling has yet to take place.
      // To Debug - figure out why only 1 xSect direction was found radiating from this nexus point.
      return SM_ERR;
    }

#undef ALIGN_DIR_WITH_PREVIOUS_POINT

} // end SmAdvSurfaceIntersector::ComputeSingularityPoint

/*******************************************************************//**
PURPOSE: Create the Dupin inticatrix curves corresponding to the 
         point on a surface of the given intersection point.

NOTES: 
    Planar     SurfacePoints                   get no Dupin Curves
    Elliptical SurfacePoints (bowls and hills) get 1 elliptical curve (either 1 above or below curve)
    Parabloic  SurfacePoints (cylinders)       get 2 line curves (either both above or both below curves)
                                                  (lines start at point and move in opposite directions)
    Hyperbolic SurfacePoints (saddles)         get 4 hyberbola curves bounded by two asymptopes
                                                  (2 above curves and 2 below curves)

    The curves are placed into crAboveDupinCurves 
      when curvature is moving in direction of given Surface Normal
    else the curves are placed into crBelowDupinCurves.
***********************************************************************/
SmStatus SmAdvSurfaceIntersector::MakeDupin
  (SmVector3d           aSurfaceEval[3][3],    // in : surface point evaluation matrix
                                               //      with positions, tangents, and 2nd parametric derivatives
   const SmVector3d   & crReferenceNormal,     // in : Reference Surface Normal (either parallel or opposing Surface Normal)
   SmTArray<SmConic*> & crAboveDupinCurves,    // out: Dupin indicatrix curves associated with positive principal curvature
                                               //      The curve generated by intersecing surface at this point with a plane
   SmTArray<SmConic*> & crBelowDupinCurves,    // out: Dupin indicatrix curves associated with negative principal curvature
                                               //      The curve generated by intersecing surface at this point with a plane
   SmVector3d         & rPrinK1Vec,            // out: Surface Tangent vector in 1st principal curvature direction
   SmVector3d         & rPrinK2Vec,            // out: Surface Tangent vector in 2nd principal curvature direction
   double             & rdPrinK1,              // out: 1st principal curvature
   double             & rdPrinK2,              // out: 2nd principal curvature
   SmVector3d         & rEFG,                  // out: 1st fundamental form
   SmVector3d         & rLMN)                  // out: 2nd fundamental form
  const
{
  // init outputs
  crAboveDupinCurves.ReSet();
  crBelowDupinCurves.ReSet();

  // locals
  double dGauss, dNormal;
  SER(smsurf_EvaluateGeometric(aSurfaceEval,            // in : surface point evaluation matrix 
                               dGauss,dNormal,          // out: gaussian and mean curvature
                               rdPrinK1,rdPrinK2,       // out: principal curvatures
                               rEFG,rLMN,               // out: 1st and 2nd fundamental forms
                               rPrinK1Vec,rPrinK2Vec)); // out: Surface Tangent vectors in 1st and 2nd principal curvature directions
  
  // negate curvatures when surface normal is opposite reference normal
  if (crReferenceNormal.Dot(aSurfaceEval[2][2]) < 0.0) 
    {
      rdPrinK1 = - rdPrinK1;
      rdPrinK2 = - rdPrinK2;
    }

  // no work - surfacePoint is planar
  if (   smos_Fabs(rdPrinK1) < SM_EFF_ZERO 
      && smos_Fabs(rdPrinK2) < SM_EFF_ZERO) 
    {
      // We have a planar point just return here
      return SM_SUCCESS;
    }

  // If one principle curvature is less than 1000 times less than the other 
  // treat the case as a cylinder and zero the small curvature value out. 
  if (smos_Fabs(rdPrinK1) < smos_Fabs(rdPrinK2) / 1000.00) { rdPrinK1 = 0.0; }
  if (smos_Fabs(rdPrinK2) < smos_Fabs(rdPrinK1) / 1000.00) { rdPrinK2 = 0.0; }

  // when the surfacePoint is cylindrical
  if (   smos_Fabs(rdPrinK1) < SM_EFF_ZERO
      || smos_Fabs(rdPrinK2) < SM_EFF_ZERO) 
    {
      // We have a parabolic point - create two parallel lines

      double dPrinK = rdPrinK1;
      // 90 degree rotation because curves
      // are parallel to the local Y axis.
      double dAngle = SM_PI/2.0;
      if (smos_Fabs(dPrinK) < SM_EFF_ZERO) 
        {
          dPrinK = rdPrinK2;
          dAngle = 0.0;  
        }

      // We have a parabolic point. Create two lines leaving this point 
      // in opposite directions on appropriate side of the surface.
      SmConic *pLine1 = new(*m_cpContext) SmConic(SM_CT_LINE,1.0/smos_Sqrt(smos_Fabs(dPrinK)),0.0);
      NER(pLine1);
      SER(pLine1->Transform2D(dAngle,SmPoint2d(0,0)));

      SmConic *pLine2 = new(*m_cpContext) SmConic(SM_CT_LINE,1.0/smos_Sqrt(smos_Fabs(dPrinK)),0.0);
      NER(pLine2);
      SER(pLine2->Transform2D(SM_PI+dAngle,SmPoint2d(0,0)));

      if (dPrinK > 0.0) 
        {
          crAboveDupinCurves.Add(pLine1);
          crAboveDupinCurves.Add(pLine2);
        }
      else 
        {
          crBelowDupinCurves.Add(pLine1);
          crBelowDupinCurves.Add(pLine2);
        }        
      return SM_SUCCESS;

    } // end cylindrical surfacePoint branch

  // surface is bending up or down from point - add ellipse to DupinCurves array
  else if (rdPrinK1 * rdPrinK2 > 0.0) 
    {
      // We have an ellipse
      SmConic *pEllipse1 = new(*m_cpContext) SmConic(SM_CT_ELLIPSE,
                                                     1.0/smos_Sqrt(smos_Fabs(rdPrinK1)),
                                                     1.0/smos_Sqrt(smos_Fabs(rdPrinK2)));
      NER(pEllipse1);
      if (rdPrinK1 > 0.0) { crAboveDupinCurves.Add(pEllipse1) ; }
      else                { crBelowDupinCurves.Add(pEllipse1) ; }
                                  
    }

  // surfacePoint is a saddle - add 4 hyberbolas to DupinCurves array
  else if (rdPrinK1 * rdPrinK2 < 0.0) 
    {
      // We have a hyperbolic point
      // First create hyperbola on right side
      SmConic *pHyp = new(*m_cpContext) SmConic(SM_CT_HYPERBOLA,
                                                1.0/smos_Sqrt(smos_Fabs(rdPrinK1)),
                                                1.0/smos_Sqrt(smos_Fabs(rdPrinK2))) ;
      NER(pHyp);
      if (rdPrinK1 > 0.0) { crAboveDupinCurves.Add(pHyp) ; }
      else                { crBelowDupinCurves.Add(pHyp) ; } 

      // Now do left side hyperbola
      pHyp = new(*m_cpContext) SmConic(SM_CT_HYPERBOLA,
                                       1.0/smos_Sqrt(smos_Fabs(rdPrinK1)),
                                       1.0/smos_Sqrt(smos_Fabs(rdPrinK2))) ;
      NER(pHyp);
      SER(pHyp->Transform2D(SM_PI,SmPoint2d(0,0)));
      if (rdPrinK1 > 0.0) { crAboveDupinCurves.Add(pHyp) ; }
      else                { crBelowDupinCurves.Add(pHyp) ; }        

      // Now do top hyperbola
      pHyp = new(*m_cpContext) SmConic(SM_CT_HYPERBOLA,
                                       1.0/smos_Sqrt(smos_Fabs(rdPrinK2)),
                                       1.0/smos_Sqrt(smos_Fabs(rdPrinK1))) ;
      NER(pHyp);
      SER(pHyp->Transform2D(SM_PI/2.0,SmPoint2d(0,0)));
      if (rdPrinK2 > 0.0) { crAboveDupinCurves.Add(pHyp) ; }
      else                { crBelowDupinCurves.Add(pHyp) ; }        

      // Now do bottom hyperbola
      pHyp = new(*m_cpContext) SmConic(SM_CT_HYPERBOLA,
                                       1.0/smos_Sqrt(smos_Fabs(rdPrinK2)),
                                       1.0/smos_Sqrt(smos_Fabs(rdPrinK2))) ;
      NER(pHyp);
      SER(pHyp->Transform2D(SM_PI+SM_PI/2.0,SmPoint2d(0,0))) ;
      if (rdPrinK2 > 0.0) { crAboveDupinCurves.Add(pHyp) ; }
      else                { crBelowDupinCurves.Add(pHyp) ; }        
    
    } // end Saddle SurfacePoint branch

  return SM_SUCCESS;

} // end SmAdvSurfaceIntersector::MakeDupin

/*******************************************************************//**
PURPOSE: Find interior intersection curves.

NOTES: 
***********************************************************************/
SmStatus SmAdvSurfaceIntersector::FindInteriorCurves
  (SmTArray<SmCurve*>         & r3DCurves,           // out: 
   SmTArray<SmCurve*>         & rSurface1UVCurves,   // out: 
   SmTArray<SmCurve*>         & rSurface2UVCurves,   // out: 
   SmTArray<SmTsectCurveType> & rCurveTypes,         // out: 
   SmTArray<double>           & rDeviations)         // out: 
{
  SmSolution aData[130];
  SmSolutionArray sSolutions(130,aData);

  // Find node pairs in the tree where intersections might exist
  SER(SolveIt(SM_SO_INTERSECT, SM_SR_ALL, m_dThisApproxTol3d,
              SM_BIG_DOUBLE, NULL, sSolutions));

  // SolidCraft change: node pairs that cannot start a crossing curve (near
  // parallel normals, no local intersection) seed the tangential contact
  // point search below.
  SmTArray<SmExtent2d> sTangentSeeds1, sTangentSeeds2;

  // SolidCraft change: 3D curves traced from low-angle start points and then
  // rejected; later low-angle start points on them are not traced again.
  SmTArray<SmCurve*> sRejectedCurves;
  SmObjsDelete<SmCurve*> sCleanRejected(&sRejectedCurves);

  // SolidCraft change: smallest normal angle worth tracing from. On a
  // transversal loop of radius r between surfaces of relative curvature
  // radius R the normals differ by about r/R. Loops of radius above the
  // tolerance on surfaces curved at model scale stay above tol/modelSize,
  // while tangent curves and tangent points (angle at round-off level, or
  // tol-sized offsets from a contact point) fall below it.
  double dMinLoopAngleRad = 0.0;
  {
    SmExtent3d sBox1, sBox2;
    SER(m_cpSurface[0]->CalculateBoundingBox(m_vUVDomain[0],&sBox1));
    SER(m_cpSurface[1]->CalculateBoundingBox(m_vUVDomain[1],&sBox2));
    double dModelSize = smos_Max(sBox1.GetSize().Length(),sBox2.GetSize().Length());
    dMinLoopAngleRad = (dModelSize > m_dThisApproxTol3d) ? 2.0*m_dThisApproxTol3d/dModelSize : SM_PI;
  }

  // Now for each node pair do a test to see if produces a unique start point
  // Note that our first pass at this is naive and slow.
  for (ULONG i=0; i<sSolutions.GetSize(); i++) 
    {
      SmSolution & rSol   = sSolutions[i] ;

      SmPoint2d       sUVGuess1(rSol.m_vStart[0],rSol.m_vStart[1]) ;
      SmPoint2d       sUVGuess2(rSol.m_vStart[2],rSol.m_vStart[3]) ;
      SmTreeNode    * pNode1     = rSol.m_apNodes[0] ; 
      SmTreeNode    * pNode2     = rSol.m_apNodes[1] ; 
      SM_ASSERT(pNode1->m_eAuxDataType == SM_AD_BEZIER_SURFACE) ;
      SM_ASSERT(pNode2->m_eAuxDataType == SM_AD_BEZIER_SURFACE) ;
      SmBezierPatch   * pBezPatch1 = (SmBezierPatch*)pNode1->m_pData; NER(pBezPatch1) ;
      SmBezierPatch   * pBezPatch2 = (SmBezierPatch*)pNode2->m_pData; NER(pBezPatch2) ;
      const SmSurface & sBSS1 = *(SmSurface *)pNode1->GetOwnerObject() ; 
      const SmSurface & sBSS2 = *(SmSurface *)pNode2->GetOwnerObject() ; 

      SM_ASSERT(pNode1->GetOwnerObject() == m_cpSurface[0]) ;
      SM_ASSERT(pNode2->GetOwnerObject() == m_cpSurface[1]) ;
      SM_ASSERT(pNode1->GetOwnerObject() == pBezPatch1->mBA_pSurface) ;
      SM_ASSERT(pNode2->GetOwnerObject() == pBezPatch2->mBA_pSurface) ;

      // Skip cases where surface normals are very close
      // SolidCraft change: a low-angle start point is no longer skipped
      // outright. A small closed loop between nearly tangent surfaces (a plane
      // dipping into a paraboloid or a sphere) has normals within a fraction of
      // a degree everywhere along it, so the angle tests below dropped every
      // start point on it and the loop was lost. Such a start point is now
      // traced on probation: the result is kept only when it is a closed,
      // transversal loop larger than the tolerance (see below); anything else
      // (tangent curves, coincident regions, tangent points) is discarded as
      // before.
      SmBoolean bLowAngleStart = FALSE;
      {
        SmVector3d sNorm1, sNorm2;
        SER(sBSS1.EvaluateNormal(sUVGuess1,TRUE,TRUE,sNorm1));
        SER(sBSS2.EvaluateNormal(sUVGuess2,TRUE,TRUE,sNorm2));

        double dAngleRad;
        SER(sNorm1.AngleBetween(sNorm2,dAngleRad));
        if (dAngleRad > SM_PI/2.0) dAngleRad = SM_PI - dAngleRad;
        if (SM_RAD2DEG(dAngleRad) < 5.0) 
          { // SolidCraft change: a seed for the tangent point search and FindTangentCurves(),
            // and a probation start for a small loop
            sTangentSeeds1.Add(pBezPatch1->GetUVDomain());
            sTangentSeeds2.Add(pBezPatch2->GetUVDomain());
            m_vTangentSeeds.Add(sUVGuess1); m_vTangentSeeds.Add(sUVGuess2);
            bLowAngleStart = TRUE;
          }
      }

#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe = FALSE;
      if (bDebugMe) 
        {
          SmPoint3d sP1, sP2;
          SER(sBSS1.EvaluatePoint(sUVGuess1,sP1));
          SER(sBSS2.EvaluatePoint(sUVGuess2,sP2));
          smgfx_SetColor(1,0,0);
          sP1.Draw();
          smgfx_SetColor(1,0,1);
          sP2.Draw();
          smgfx_SetColor(0,1,0);
          SmVector3d sVec = sP2 - sP1;
          sVec.Draw(&sP1);
          sm_GraphicsLoop();
        }
#endif
      // First let's do a fast test to see if the intersection of
      // these two nodes is close to a 3D curve.  If we pass this
      // test then we can go on to a stricter test.
      //
      // No, this test can exclude legitimate points, and thereby cause
      // intersection curves to be missed.  [090702]

//      SmExtent3d sIntBBox = pNode1->m_sBBox;
//      sIntBBox.ExpandAbsolute(m_dThisApproxTol3d);
//      SER(sIntBBox.Intersect(pNode2->m_sBBox,sIntBBox));
//  
//      SmBoolean bIsNearCurve = FALSE;
//  
//      for (ULONG j=0; j<r3DCurves.GetSize(); j++) 
//        {
//          SmCurve *pCurve = r3DCurves[j];
//          SmCurveCache *pCC = (SmCurveCache*)
//              SmCacheMgr::GetOrCreateObjectCache(SM_OC_CURVE,pCurve); 
//          NER(pCC);
//          SmTree *pTree = pCC->GetTree();
//          if (pTree->IntersectsBox(sIntBBox)) 
//            {
//              bIsNearCurve = TRUE;
//              break;
//            }
//        }
//  
//      if (bIsNearCurve) continue;

      // Don't do the following 'Drop points back and forth'.
      // Algorithmically, it's part of a very slow method, and in
      // this case, it can result in values that are way off,
      // because we're working with very small pieces of surfaces
      // that may not intersect.  [bd, 05Jan06; regressions:
      // Two Tori, and my_shell_demo iters 0, 3, 6]
#if 0  // (skip the flip)
      { // Drop points back and forth to get better guess points
        SmPoint3d sP1, sP2;
        SER(sBSS2.EvaluatePoint(sUVGuess2,sP2));
        SmPoint3d sDropPoints[1];
        sDropPoints[0] = sP2;
        SmVector2d sUVPoints[1];
        SER(sBSS1.DropPointsToTangentPlane(sUVGuess1,TRUE,TRUE,1,sDropPoints,sUVPoints));
        sUVGuess1 = m_vUVDomain[0].ClampPoint2d(sUVPoints[0]);

        SER(sBSS1.EvaluatePoint(sUVGuess1,sP1));
        sDropPoints[0] = sP1;
        SER(sBSS2.DropPointsToTangentPlane(sUVGuess2,TRUE,TRUE,1,sDropPoints,sUVPoints));
        sUVGuess2 = m_vUVDomain[1].ClampPoint2d(sUVPoints[0]);
      }
#endif  // (skip the flip)

// GWC:DEBUG. V6.5, 6/18/04. BACKED OUT A 6.0 to 6.5 CHANGE. 
// change = conditional skip [CAUSES FAILURE IN GINA0148 BUG]
//        if (!sUVDomain1.ContainsPoint2d(sUVGuess1) ||
//            !sUVDomain2.ContainsPoint2d(sUVGuess2) ) {
//            continue;
//        }
      // Get an intersection curve start point and convert it to 3D
      SmSolution sSol;
      SmBoolean bFoundAnswer;
// GWC:DEBUG. V6.5, 6/18/04. BACKED OUT A 6.0 to 6.5 CHANGE. 
// GWC:DEBUG REMOVE [CAUSES FAILURE IN GINA0148 BUG]
//        SER(m_cpSurface[0]->LocalSurfaceIntersect(sUVDomain1,
//            *m_cpSurface[1],sUVDomain2,m_dThisApproxTol3d/100.0,
//            sUVGuess1,sUVGuess2,
//            NULL,NULL,FALSE,
//            bFoundAnswer, sSol));


// GWC:DEBUG RESTORED 6.0 VERSION
      SER(m_cpSurface[0]->LocalSurfaceIntersect(m_vUVDomain[0],
                                               *m_cpSurface[1],
                                                m_vUVDomain[1],
                                                m_dThisApproxTol3d/100.0,
                                                sUVGuess1,
                                                sUVGuess2,
                                                NULL,
                                                NULL,
                                                FALSE,
                                                bFoundAnswer, 
                                                sSol));
// end DEBUG
      // for nodes that do not intersect - skip to next node pair
      if (!bFoundAnswer) 
        {
          sTangentSeeds1.Add(pBezPatch1->GetUVDomain());
          sTangentSeeds2.Add(pBezPatch2->GetUVDomain());
          continue ;
        }

      // Get solution 3d and UV points point
      SmPoint3d sSurfPnt ;
      SmPoint2d sUV (sSol.m_vStart[0], sSol.m_vStart[1]) ;
      SmPoint2d sUV2(sSol.m_vStart[2], sSol.m_vStart[3]) ;
      SER(m_cpSurface[0]->EvaluatePoint(sUV, sSurfPnt) );

      // Skip cases where surface normals are very close
        {
          // solution normals
          SmVector3d sNorm1, sNorm2 ;
          SER(sBSS1.EvaluateNormal(sUV,TRUE,TRUE,sNorm1)) ;
          SER(sBSS2.EvaluateNormal(sUV2,TRUE,TRUE,sNorm2)) ;

          // angle between normals
          double dAngleRad ;
          SER(sNorm1.AngleBetween(sNorm2,dAngleRad)) ;
          if (dAngleRad > SM_PI/2.0) dAngleRad = SM_PI - dAngleRad ;

          // probation for this node pair when angle between normals is small;
          // skip it when the angle is too small for a loop above tolerance
          if (SM_RAD2DEG(dAngleRad) < 2.0) 
            { // SolidCraft change: a seed for the tangent point search and FindTangentCurves()
              if (!bLowAngleStart)
                {
                  sTangentSeeds1.Add(pBezPatch1->GetUVDomain());
                  sTangentSeeds2.Add(pBezPatch2->GetUVDomain());
                }
              m_vTangentSeeds.Add(sUV); m_vTangentSeeds.Add(sUV2);
              bLowAngleStart = TRUE;
            }
          if (bLowAngleStart && dAngleRad < dMinLoopAngleRad)
            {
              continue;
            }
        }

      // skip start points that are already on xsect curves
      SmBoolean bPointIsOnCurve = IsPointOnCurve(sSurfPnt,r3DCurves);
      if ( bPointIsOnCurve ) 
        { continue; }
      if ( bLowAngleStart && IsPointOnCurve(sSurfPnt,sRejectedCurves) )
        { continue; }
          
      // If we made it this far then we have a point which should serve as
      // a good start point.  Load it up and trace it.

      // make a SmTSectPnt for the intersection solution
      SmPoint2d sUVPts[2] ;
      sUVPts[0] = sUV ;
      sUVPts[1] = SmPoint2d(sSol.m_vStart[2],sSol.m_vStart[3]) ;
      SmTsectPnt * pStartTSP = (SmTsectPnt*)m_vTSPntMgr.GetNewElement() ;
      m_dCurveTraceDirection = 1.0 ;
      if(SM_SUCCESS != ComputePointValues(sUVPts, *pStartTSP, NULL)) 
        { continue ;}

      // init trace params
      double dInputTraceDirection = m_dCurveTraceDirection ;
      m_pStartPoint = pStartTSP ;

      // do the trace
      SER(TraceCurve(*pStartTSP,r3DCurves)) ;

      // when appropriate - trace in the other direction as well
      if(   !m_bCurveIsClosed
         && !pStartTSP->IsBounded()
         &&  dInputTraceDirection == m_dCurveTraceDirection)
        {
          SER(ReverseCurveDirection());
          m_dCurveTraceDirection *= -1.0;
          SER(TraceCurve(*pStartTSP,r3DCurves));

          // a cusp will be made if TraceCurve reversed directions on the last call
          SM_ASSERT_MSG(dInputTraceDirection != m_dCurveTraceDirection,
                        _T("SmAdvSurfaceIntersector::FindInteriorCurves: Assumed TraceCurve() would not reverse trace direction - probably building a cusp")) ;
          
          // restore the curve direction
          m_dCurveTraceDirection *= 1.0;
          SER(ReverseCurveDirection());
        }

      // SolidCraft change: a low-angle start point only yields output when it
      // traced a closed loop that is not a single point within tolerance.
      SmBoolean bKeepCurve = TRUE;
      if (bLowAngleStart)
        {
          bKeepCurve = FALSE;
          if (m_bCurveIsClosed)
            {
              SmTsectPnt *aPntData[130];
              SmTArray<SmTsectPnt*> sLoopPnts(130,aPntData);
              m_vCurvePoints.GetAllNodes(sLoopPnts);
              if (sLoopPnts.GetSize() > 2)
                {
                  SmExtent3d sLoopBox(sLoopPnts[0]->CrvPos());
                  for (ULONG lPnt=1; lPnt<sLoopPnts.GetSize(); lPnt++)
                    { sLoopBox.Union(SmExtent3d(sLoopPnts[lPnt]->CrvPos()),sLoopBox); }
                  bKeepCurve = (sLoopBox.GetSize().Length() > 2.0*m_dThisApproxTol3d) ? TRUE : FALSE;
                }
            }
        }

      // turn m_vCurvePoints into an intersection curve
      ULONG lNumCurves = r3DCurves.GetSize();
      SER(FlushCurve(r3DCurves,rSurface1UVCurves,rSurface2UVCurves,rCurveTypes,rDeviations));
      m_pStartPoint = NULL;

      // drop a probationary curve that is not a transversal loop, and
      // remember it and its start point so nearby start points are not
      // traced again
      SmBoolean bReject =    bLowAngleStart
                          && (   r3DCurves.GetSize() == lNumCurves
                              || !bKeepCurve
                              || rCurveTypes.GetLast() == SM_TC_TANGENT
                              || rCurveTypes.GetLast() == SM_TC_COINCIDENT
                              || rCurveTypes.GetLast() == SM_TC_TOUCHING);
      if (bReject)
        {
          while (r3DCurves.GetSize() > lNumCurves)
            {
              sRejectedCurves.Add(r3DCurves.GetLast());
              r3DCurves.RemoveLast();
              if (rSurface1UVCurves.GetSize() > lNumCurves)
                { SmObjDelete sDelUV1(rSurface1UVCurves.GetLast()); rSurface1UVCurves.RemoveLast(); }
              if (rSurface2UVCurves.GetSize() > lNumCurves)
                { SmObjDelete sDelUV2(rSurface2UVCurves.GetLast()); rSurface2UVCurves.RemoveLast(); }
              if (rCurveTypes.GetSize() > lNumCurves) rCurveTypes.RemoveLast();
              if (rDeviations.GetSize() > lNumCurves) rDeviations.RemoveLast();
            }
          SmBSplineCurve *pStartPnt = NULL;
          SER(SmBSplineCurve::CreateDegenerateCurve(*m_cpContext,3,sSurfPnt,pStartPnt));
          sRejectedCurves.Add(pStartPnt);
        }
    
    } // end iter every node/node solution result

  // SolidCraft change: isolated tangential contact points
  if (!m_bDoingSelfIntersection)
    {
      SER(sm_FindTangentContactPoints(*m_cpContext,
                                      *m_cpSurface[0], m_vUVDomain[0],
                                      *m_cpSurface[1], m_vUVDomain[1],
                                      m_dThisApproxTol3d,
                                      sTangentSeeds1, sTangentSeeds2, *this,
                                      r3DCurves, rSurface1UVCurves, rSurface2UVCurves,
                                      rCurveTypes, rDeviations));
    }

  // sll done
  return SM_SUCCESS;

} // end SmAdvSurfaceIntersector::FindInteriorCurves

/*******************************************************************//**
PURPOSE: This is a special solve branch which does some really 
   fancy stuff with the gauss map of two surfaces to determine
   when it might be necessary to try to find interior intersection loops.
   Basically the way it works is to find only one start point after 
   seperability is achieved.  

NOTES: This is a recursive method which traverses down the
   tree (subdivides) until it can do a local solve.
***********************************************************************/
SmStatus SmAdvSurfaceIntersector::SolveBranch
  (SmTreeNode * apBranch[SM_GS_MAX_TREES])
{
  m_bSeperateGaussMaps = TRUE;

  if (!BranchMayContainAnswers(apBranch)) 
    { return SM_SUCCESS; }

  // no work - SolverOperation != INTERSECT - pass the call along to the regular SmSurfaceIntersector::SolveBranch()
  if (m_eSolverOperation != SM_SO_INTERSECT) 
    {
      SER(SmSurfaceIntersector::SolveBranch(apBranch));
      return SM_SUCCESS;
    }

  // when intersecting
  if (m_eSolverOperation == SM_SO_INTERSECT) 
    {
      SmBoolean bSaveAnswer = FALSE;
      
      SmBoolean bAchievedSeperability = FALSE;
      if (!m_bSeperateGaussMaps) 
        {
          SER(BranchIsGaussSeparate(apBranch,bAchievedSeperability));
        }
      
      if (bAchievedSeperability) 
        {
          // If we achieve seperability at the top nodes
          // of both trees than there is no need to do interior
          // curve processing.
          SmTreeNode * pSurfaceNode1 = apBranch[0];
          SmTreeNode * pSurfaceNode2 = apBranch[1];
          if(   pSurfaceNode1->m_pParent == NULL 
             && pSurfaceNode2->m_pParent == NULL) 
            {
              return SM_SUCCESS;
            }
          
          SmTemporaryChangeValue<SmBoolean> sChange(m_bSeperateGaussMaps,TRUE);
          m_dClosestValue = SM_BIG_DOUBLE;
          bSaveAnswer     = TRUE;
          
          if (AreReadyForLocalSolve(apBranch)) 
            {
              SmBoolean bNeedsMoreSubdivision;
              SER(LocalSolve(apBranch,bNeedsMoreSubdivision));
            }
          else 
            {
              SmTreeNode* sChild1[SM_GS_MAX_TREES];
              SmTreeNode* sChild2[SM_GS_MAX_TREES];
              SER(Subdivide(apBranch,sChild1,sChild2));
              SER(SolveBranch(sChild1));
              SER(SolveBranch(sChild2));
            }
        } // end BranchIsGaussSeperate branch
      else // Branch is not GaussSeparate
        {
          if (AreReadyForLocalSolve(apBranch)) 
            {
              SmTemporaryChangeValue<SmBoolean> sChange(m_bSeperateGaussMaps,TRUE);
              m_dClosestValue = SM_BIG_DOUBLE;
              // If we are unable to achieve seperability at the leaf level
              // than we save the answer any how.
              bSaveAnswer = TRUE;
              SmBoolean bNeedsMoreSubdivision;
              SER(LocalSolve(apBranch,bNeedsMoreSubdivision));
            }
          else 
            {
              SmTreeNode* sChild1[SM_GS_MAX_TREES];
              SmTreeNode* sChild2[SM_GS_MAX_TREES];
              SER(Subdivide(apBranch,sChild1,sChild2));
              SER(SolveBranch(sChild1));
              SER(SolveBranch(sChild2));
            }
        }
      
      // If m_dClosestValue < SM_BIG_DOUBLE then some of the leaf nodes
      // intersected -- therefore save the best one as an answer.
      if (bSaveAnswer && m_dClosestValue < SM_BIG_DOUBLE) 
        {
          // Save the midpoint of the leaf nodes which intersect
          
          SmSolution sSolution;
          sSolution.m_eSolutionType = SM_ST_SINGLE_VALUE;
          sSolution.m_lNumVariables = 4;
          sSolution.m_vStart[0]     = m_sUVS[0].x;
          sSolution.m_vStart[1]     = m_sUVS[0].y;
          sSolution.m_vStart[2]     = m_sUVS[1].x;
          sSolution.m_vStart[3]     = m_sUVS[1].y;
          sSolution.m_apNodes[0]    = m_apNodes[0];
          sSolution.m_apNodes[1]    = m_apNodes[1];
          m_pSolutions->Add(sSolution);
        }
    } // end When intersecting check

  // all done
  return SM_SUCCESS;

} // end SmAdvSurfaceIntersector::SolveBranch

/*******************************************************************//**
PURPOSE: Local solver for the SmAdvSurfaceIntersector.  Basically
    all it does is to compare the midpoints of the two patches to
    see if they are closer than a previous answer after it does a 
    tight seperability test.

NOTES: m_eSolverOperation == SM_SO_INTERSECT - take a NR step and update values:
              m_dClosestValue = dDist;
              m_sUVS[0]       = sUV1;
              m_sUVS[1]       = sUV2;
              m_apNodes[0]    = apBranch[0];
              m_apNodes[1]    = apBranch[1];

       m_eSolverOperation == SM_SO_INTERSECTION_TEST Get union of all intersection BBoxes of branch solutions that intersect
       m_pSolutions[i].m_lNumVariables = 7;
                       m_eSolutionType = SM_ST_SINGLE_VALUE;
                       m_vStart[0] = sResult.GetMin().x; // \ 
                       m_vStart[1] = sResult.GetMin().y; // |
                       m_vStart[2] = sResult.GetMin().z; //  \ xyz XSecting BBox union containing possible solution
                       m_vStart[3] = sResult.GetMax().x; //  /
                       m_vStart[4] = sResult.GetMax().y; // |
                       m_vStart[5] = sResult.GetMax().z; // /
                       m_vStart[6] = 1;                  // number of XSecting BBoxs unioned together
                       m_vEnd[0] = sUVDomain1.GetMin().x; // \
                       m_vEnd[1] = sUVDomain1.GetMin().y; //  \ uv XSecting BBox union containing possible solution for Surface1
                       m_vEnd[2] = sUVDomain1.GetMax().x; //  /
                       m_vEnd[3] = sUVDomain1.GetMax().y; // /
                       m_vEnd[4] = sUVDomain2.GetMin().x; // \
                       m_vEnd[5] = sUVDomain2.GetMin().y; //  \ uv XSecting BBox union containing possible solution for Surface2
                       m_vEnd[6] = sUVDomain2.GetMax().x; //  /
                       m_vEnd[7] = sUVDomain2.GetMax().y; // /
***********************************************************************/
SmStatus SmAdvSurfaceIntersector::LocalSolve
  (SmTreeNode * apBranch[SM_GS_MAX_TREES],  // in : candidate nodes for solution
   SmBoolean & rbNeedsMoreSubdivision)      // out: TRUE = 
{
  // init output
  rbNeedsMoreSubdivision = FALSE;

  // SurfaceNode 1 and 2 locals:
  SmTreeNode      * pSurfaceNode1 = apBranch[0];
  SmTreeNode      * pSurfaceNode2 = apBranch[1];
  const SmSurface & sBSP1         = *(SmSurface *)pSurfaceNode1->GetOwnerObject() ; 
  const SmSurface & sBSP2         = *(SmSurface *)pSurfaceNode2->GetOwnerObject() ;
  SmBezierPatch   * pBezPatch1    = (SmBezierPatch*)pSurfaceNode1->m_pData; NER(pBezPatch1);
  SmBezierPatch   * pBezPatch2    = (SmBezierPatch*)pSurfaceNode2->m_pData; NER(pBezPatch2);
  SmExtent2d        sUVDomain1    = pBezPatch1->GetUVDomain();
  SmExtent2d        sUVDomain2    = pBezPatch2->GetUVDomain();

  SM_ASSERT(pSurfaceNode1->m_eAuxDataType == SM_AD_BEZIER_SURFACE);
  SM_ASSERT(pSurfaceNode2->m_eAuxDataType == SM_AD_BEZIER_SURFACE);

  //      SM_ASSERT(pSurfaceNode1->GetOwnerObject()->IsKindOf(SmBSplineSurface_TYPE)) ;
  //      SM_ASSERT(pSurfaceNode2->GetOwnerObject()->IsKindOf(SmBSplineSurface_TYPE)) ;
  SM_ASSERT(pSurfaceNode1->GetOwnerObject() == m_cpSurface[0]) ;
  SM_ASSERT(pSurfaceNode2->GetOwnerObject() == m_cpSurface[1]) ;
  SM_ASSERT(pSurfaceNode1->GetOwnerObject() == pBezPatch1->mBA_pSurface) ;
  SM_ASSERT(pSurfaceNode2->GetOwnerObject() == pBezPatch2->mBA_pSurface) ;

  // Do interval tests and reject if disjoint or adjust
  // if sUVDomain is not completely inside of m_crUVDomain
  if (sUVDomain1.AreDisjoint(m_vUVDomain[0])) { return SM_SUCCESS; }
  if (sUVDomain2.AreDisjoint(m_vUVDomain[1])) { return SM_SUCCESS; }

  // If we are doing self intersection then reject case where domains 
  // are not disjoint.
  if (m_bDoingSelfIntersection) 
    {
      if (!sUVDomain1.AreDisjoint(sUVDomain2)) 
        { return SM_SUCCESS; }
    }

  // when looking for intersections
  if (m_eSolverOperation == SM_SO_INTERSECT) 
    {
      SmPseudoBox sTmpPBox0 = pBezPatch1->GetPseudoBox();
      sTmpPBox0.ExpandAbsolute( m_d3dTolerance );
      if (!sTmpPBox0.AreDisjoint(pBezPatch2->GetPseudoBox())) 
        {
          SmVector3d sSize1 = pSurfaceNode1->m_sBBox.GetSize();
          SmVector3d sSize2 = pSurfaceNode2->m_sBBox.GetSize();
          // If the relative sizes of the boxes of these patches
          // are not close than drop point from one to the other
          // to get starting points.
          SmPoint3d sP1, sP2;
          SmPoint2d sUV1 = sUVDomain1.Evaluate(0.5,0.5);
          SmPoint2d sUV2 = sUVDomain2.Evaluate(0.5,0.5);
          SER(sBSP1.EvaluatePoint(sUV1,sP1));
          SER(sBSP2.EvaluatePoint(sUV2,sP2));
          double dDist = sP1.DistanceBetween(sP2);

          // If the sizes are very different - try taking one Newton step
          // to get a closer guess point on the big surface.
          if (sSize1.Length() > 4.0*sSize2.Length()) 
            {
              SmPoint3d sDropPoints[1];
              sDropPoints[0] = sP2;
              SmVector2d sUVPoints[1];
              SER(sBSP1.DropPointsToTangentPlane(sUV1,TRUE,TRUE,1,sDropPoints,sUVPoints));
              SmPoint2d sUV1New = sUVDomain1.ClampPoint2d(sUVPoints[0]);
              SmPoint3d sP1New;
              SER(sBSP1.EvaluatePoint(sUV1New,sP1New));
              double dNewDist = sP1New.DistanceBetween(sP2);
              if (dNewDist < dDist) 
                {
                  sUV1 = sUV1New;
                  dDist = dNewDist;
                }
            }
          else if (sSize2.Length() > 4.0*sSize1.Length()) 
            {
              SmPoint3d sDropPoints[1];
              sDropPoints[0] = sP1;
              SmVector2d sUVPoints[1];
              SER(sBSP2.DropPointsToTangentPlane(sUV2,TRUE,TRUE,1,sDropPoints,sUVPoints));
              SmPoint2d sUV2New = sUVDomain2.ClampPoint2d(sUVPoints[0]);
              SmPoint3d sP2New;
              SER(sBSP2.EvaluatePoint(sUV2New,sP2New));
              double dNewDist = sP2New.DistanceBetween(sP1);
              if (dNewDist < dDist) 
                {
                  sUV2 = sUV2New;
                  dDist = dNewDist;
                }
            }

          if (dDist < m_dClosestValue) 
            {
              m_dClosestValue = dDist;
              m_sUVS[0] = sUV1;
              m_sUVS[1] = sUV2;
              m_apNodes[0] = apBranch[0];
              m_apNodes[1] = apBranch[1];
            }
        }
    } // end m_eSolverOperation == SM_SO_INTERSECTION branch

  // SM_SO_INTERSECTION_TEST branch
  else if (m_eSolverOperation == SM_SO_INTERSECTION_TEST) 
    {
      // get Patch1 PseudoBox
      SmPseudoBox sTmpPBox0 = pBezPatch1->GetPseudoBox();
      SmPseudoBox sTmpPBox1 = pBezPatch2->GetPseudoBox();
      sTmpPBox0.ExpandAbsolute( m_d3dTolerance );

      // when Patch1 and Patch2 PseudoBoxes intersect
      if (!sTmpPBox0.AreDisjoint(sTmpPBox1)) 
        {
          // Update the bounding box information in the solution.
          SmExtent3d sBBox = pSurfaceNode1->m_sBBox;
          sBBox.ExpandAbsolute( m_d3dTolerance );
          SmExtent3d sResult;
          sBBox.Intersect( pSurfaceNode2->m_sBBox, sResult );

          // Add intersection to solution: 
          //   Union 3D and UV Solution BBoxes, inc solution count
          
          // When this is the first solution
          if (m_pSolutions->GetSize() == 0) 
            {
              SmSolution sSol;
              sSol.m_apNodes[0] = apBranch[0];
              sSol.m_apNodes[1] = apBranch[1];
              sSol.m_lNumVariables = 7;
              sSol.m_eSolutionType = SM_ST_SINGLE_VALUE;
              sSol.m_vStart[0] = sResult.GetMin().x;
              sSol.m_vStart[1] = sResult.GetMin().y;
              sSol.m_vStart[2] = sResult.GetMin().z;
              sSol.m_vStart[3] = sResult.GetMax().x;
              sSol.m_vStart[4] = sResult.GetMax().y;
              sSol.m_vStart[5] = sResult.GetMax().z;
              sSol.m_vStart[6] = 1;
              sSol.m_vEnd[0] = sUVDomain1.GetMin().x;
              sSol.m_vEnd[1] = sUVDomain1.GetMin().y;
              sSol.m_vEnd[2] = sUVDomain1.GetMax().x;
              sSol.m_vEnd[3] = sUVDomain1.GetMax().y;
              sSol.m_vEnd[4] = sUVDomain2.GetMin().x;
              sSol.m_vEnd[5] = sUVDomain2.GetMin().y;
              sSol.m_vEnd[6] = sUVDomain2.GetMax().x;
              sSol.m_vEnd[7] = sUVDomain2.GetMax().y;
              m_pSolutions->Add(sSol);
            } // end 1st solution branch
          else // add this solution to the previous solutions
            {
              SmSolution & rSol = (*m_pSolutions)[0];
              rSol.m_vStart[6] = rSol.m_vStart[6] + 1.0;  // Number of contributors to intersection

              SmExtent3d sOldBox(SmPoint3d(rSol.m_vStart[0],rSol.m_vStart[1],rSol.m_vStart[2]),
                                 SmPoint3d(rSol.m_vStart[3],rSol.m_vStart[4],rSol.m_vStart[5]));

              sResult.Union(sOldBox,sResult);

              rSol.m_vStart[0] = sResult.GetMin().x;
              rSol.m_vStart[1] = sResult.GetMin().y;
              rSol.m_vStart[2] = sResult.GetMin().z;
              rSol.m_vStart[3] = sResult.GetMax().x;
              rSol.m_vStart[4] = sResult.GetMax().y;
              rSol.m_vStart[5] = sResult.GetMax().z;

              SmExtent2d sOldUV1(SmPoint2d(rSol.m_vEnd[0],rSol.m_vEnd[1]),
                                 SmPoint2d(rSol.m_vEnd[2],rSol.m_vEnd[3]));

              SmExtent2d sOldUV2(SmPoint2d(rSol.m_vEnd[4],rSol.m_vEnd[5]),
                                 SmPoint2d(rSol.m_vEnd[6],rSol.m_vEnd[7]));

              sOldUV1.Union(sUVDomain1,sUVDomain1);
              sOldUV2.Union(sUVDomain2,sUVDomain2);
              rSol.m_vEnd[0] = sUVDomain1.GetMin().x;
              rSol.m_vEnd[1] = sUVDomain1.GetMin().y;
              rSol.m_vEnd[2] = sUVDomain1.GetMax().x;
              rSol.m_vEnd[3] = sUVDomain1.GetMax().y;
              rSol.m_vEnd[4] = sUVDomain2.GetMin().x;
              rSol.m_vEnd[5] = sUVDomain2.GetMin().y;
              rSol.m_vEnd[6] = sUVDomain2.GetMax().x;
              rSol.m_vEnd[7] = sUVDomain2.GetMax().y;
            } // end m_pSolutions->GetSize() > 0 branch
        } // end boundnig leafNode bounding boxes intersect check
    } // end m_eSolverOperation == SM_SO_INTERSECTION_TEST check

  // all done
  return SM_SUCCESS;

} // end SmAdvSurfaceIntersector::LocalSolve


/*******************************************************************//**
PURPOSE: Determine if the Gauss Maps of this branch are seperable.

NOTES: 
***********************************************************************/
SmStatus SmAdvSurfaceIntersector::BranchIsGaussSeparate
  (SmTreeNode * apBranch[SM_GS_MAX_TREES], 
   SmBoolean & rbIsSeperable)
{
  SmTreeNode *pSurfaceNode1 = apBranch[0];
  SmTreeNode *pSurfaceNode2 = apBranch[1];
  rbIsSeperable = FALSE;
  if(   (   pSurfaceNode1->m_eAuxDataType == SM_AD_BEZIER_SURFACE
         || pSurfaceNode1->m_eAuxDataType == SM_AD_BEZIER_SURFACE_HEAD)
     && (   pSurfaceNode2->m_eAuxDataType == SM_AD_BEZIER_SURFACE
         || pSurfaceNode2->m_eAuxDataType == SM_AD_BEZIER_SURFACE_HEAD) ) 
    {
      SmBezierPatch *pBezPatch1 = (SmBezierPatch*)pSurfaceNode1->m_pData; NER(pBezPatch1);
      SmBezierPatch *pBezPatch2 = (SmBezierPatch*)pSurfaceNode2->m_pData; NER(pBezPatch2);

      if(pBezPatch1->GetPolarBox().AreDisjoint(pBezPatch2->GetPolarBox())) 
        {
          rbIsSeperable = TRUE;
        }
    }                                                   
  return SM_SUCCESS;

} // end SmAdvSurfaceIntersector::BranchIsGaussSeparate

/*******************************************************************//**
PURPOSE: 

NOTES: 
***********************************************************************/
static ULONG sm_FindLargestNonLeaf
  (ULONG ,                                // lNumTrees
   SmTreeNode * apBranch[SM_GS_MAX_TREES])
{
    SmTreeNode *pSurfaceNode1 = apBranch[0];
    SmTreeNode *pSurfaceNode2 = apBranch[1];

    SmBoolean bLeaf1 = FALSE;
    SmBoolean bLeaf2 = FALSE;
    if (pSurfaceNode1->m_pChild1 == NULL) {
        bLeaf1 = TRUE;
    }
    else if (pSurfaceNode1->m_eAuxDataType == SM_AD_BEZIER_SURFACE &&
        pSurfaceNode1->m_pChild1->m_eAuxDataType != SM_AD_BEZIER_SURFACE) {
        bLeaf1 = TRUE;
    }

    if (pSurfaceNode2->m_pChild1 == NULL) {
        bLeaf2 = TRUE;
    }
    else if (pSurfaceNode2->m_eAuxDataType == SM_AD_BEZIER_SURFACE &&
        pSurfaceNode2->m_pChild1->m_eAuxDataType != SM_AD_BEZIER_SURFACE) {
        bLeaf2 = TRUE;
    }

    // If either are leaf nodes than return the other node to be split
    if (bLeaf1) {
        return 1;
    }
    if (bLeaf2) {
        return 0;
    }

//    ULONG lFound = 0;
//    if (pSurfaceNode2->GetPolarBox().HasLargerSpan(pSurfaceNode1->GetPolarBox())) {
//        lFound = 1;
//    }
//    return lFound;
    return 0;

} // end sm_FindLargestNonLeaf


/*******************************************************************//**
PURPOSE: Find which node should be subdivided and see if we need to
    swap the traversal order to achieve more optimal minimization and
    maximization.

NOTES: 
***********************************************************************/
SmStatus SmAdvSurfaceIntersector::FindBestSubdivisionNode
  (SmTreeNode * apBranch[SM_GS_MAX_TREES],
   ULONG & rlBestSubdivisionIndex,
   SmBoolean & rbSwapTraversalOrder) 
 const
{
    rbSwapTraversalOrder = FALSE;
    rlBestSubdivisionIndex = sm_FindLargestNonLeaf(m_lNumTrees,apBranch);

    return SM_SUCCESS;

} // end SmAdvSurfaceIntersector::FindBestSubdivisionNode

/*******************************************************************//**
PURPOSE:
NOTES:
***********************************************************************/
SmBoolean SmAdvSurfaceIntersector::IsKindOf( SM_TYPE t ) const
{
  return ((SmAdvSurfaceIntersector_TYPE == t) ? TRUE : SmSurfaceIntersector::IsKindOf( (t) ));
}


/*******************************************************************//**
PURPOSE: Pretty print   

NOTES:
***********************************************************************/
void SmAdvSurfaceIntersector::Dump(void) const
{
  // TCHAR sBuff[SM_TBLOCK_SIZE];

  // start
  smos_WriteBuffer(_T("\nBegin SmAdvSurfaceIntersector::Dump()")) ;

  // dump base
  SmSurfaceIntersector::Dump() ;

  // pretty print values

  // end
  smos_WriteBuffer(_T("\nEnd SmAdvSurfaceIntersector::Dump()\n")) ;

} // end SmAdvSurfaceIntersector::Dump

