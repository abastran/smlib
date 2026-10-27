// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

/**********************************************************************//**
* FILE NAME --- SmFilletCliffOverflow.cpp
* PURPOSE: SolidCraft change. The cliff overflow of a convex blend: a blend
*   whose rails leave the faces they roll on across cliffs (an opening's rim,
*   the body's outline) is built as if those faces went on past the cliffs,
*   and the material it removes is subtracted from the body.
*
* The executive computes the fillet as usual on the support faces' surfaces
* (CreateFilletBrep, m_bCliffOverflow): m_pFilletBrep is the sheet of fillet
* and corner faces, its boundary on the surfaces of the faces the blends roll
* on and end on. SubtractRemovalRegion
*   1. cuts the fillet off at the ends its Nx1 corners left uncut (the rails
*      missing their side edges there) by the end face's surface continued
*      (CutFilletAtOpenEnd);
*   2. closes the sheet into the region the blends remove (CloseRemovalRegion):
*      on each target face the sheet's boundary runs on, a face on that face's
*      surface, continued exactly where it needs to be, bounded by the sheet's
*      boundary there, the edges blended, and the side edges (continued along
*      their lines past their ends when they must be) that join the two;
*   3. subtracts that region from the target with SMLib's manifold boolean.
* The result is the body blended as if the cliffs were not there, then cut by
* what the cliffs cut away: subtraction only removes material, so the blend
* never bridges an opening.
*
* Exactness: a plane continues exactly (the same plane, a larger domain); any
* other surface is used within its own domain, and a needed continuation of
* one is refused (SM_FILERR_NOT_EXACT). A side edge continues past its end only
* when it is a straight line (SM_FILERR_NOT_EXACT otherwise). An end the
* fillet runs past is cut only by a planar end face. Configurations the
* construction does not cover (a non-planar end face run past, a closed rail
* on a face, a sheet boundary point where it does not pass from one face to
* another along a side edge) stay what they were before this change: too
* large for the faces, SM_FILERR_NO_ROOM, the reason in the error's comment.
**********************************************************************/

#include "StdAfx.h"

#include <SmFilletExecutive.h>
#include <SmFilletGeom.h>
#include <SmFilletSolver.h>
#include <SmFilletStandardSolver.h>
#include <SmBrep.h>
#include <SmFace.h>
#include <SmEdge.h>
#include <SmVertex.h>
#include <SmLine.h>
#include <SmPlane.h>
#include <SmSurface.h>
#include <SmCurve.h>
#include <SmRegion.h>
#include <SmPointClass.h>
#include <SmTol.h>

#include <algorithm>
#include <map>
#include <set>
#include <vector>

namespace
{

// A curve of a closing face's loop: owned until the face takes it, run from crFrom to crTo.
struct sm_LoopCurve
{
  SmCurve    * m_pCurve = NULL;
  SmOrientType m_eOrient = SM_OT_SAME;
};

// The curves of one loop of a closing face, deleted unless handed on.
struct sm_Loop
{
  std::vector<sm_LoopCurve> m_vCurves;
  ~sm_Loop() { for (sm_LoopCurve & rC : m_vCurves) { delete rC.m_pCurve; } }
};

// The two ends of a curve over its natural interval.
SmStatus sm_CurveEnds(const SmCurve & crCurve, SmPoint3d & rStart, SmPoint3d & rEnd)
{
  SmExtent1d sIvl = crCurve.GetNaturalInterval();
  SER( crCurve.EvaluatePoint( sIvl.GetMin(), rStart ));
  SER( crCurve.EvaluatePoint( sIvl.GetMax(), rEnd ));
  return SM_SUCCESS;
}

// Adds pCurve (taken over) to rLoop, oriented to run from crFrom to crTo.
SmStatus sm_AddCurve(sm_Loop & rLoop, SmCurve * pCurve, const SmPoint3d & crFrom, const SmPoint3d & crTo)
{
  NER(pCurve);
  sm_LoopCurve sC;
  sC.m_pCurve = pCurve;
  rLoop.m_vCurves.push_back(sC);
  SmPoint3d sStart, sEnd;
  SER( sm_CurveEnds( *pCurve, sStart, sEnd ));
  double dSame     = sStart.DistanceBetween( crFrom ) + sEnd.DistanceBetween( crTo );
  double dOpposite = sStart.DistanceBetween( crTo )   + sEnd.DistanceBetween( crFrom );
  rLoop.m_vCurves.back().m_eOrient = (dSame <= dOpposite) ? SM_OT_SAME : SM_OT_OPPOSITE;
  return SM_SUCCESS;
}

// A copy of the edge's curve over the edge's interval.
SmStatus sm_EdgeCurve(const SmContext & crContext, const SmEdge & crEdge, SmCurve *& rpCurve)
{
  rpCurve = NULL;
  const SmCurve * pCurve = crEdge.GetCurve();   NER(pCurve);
  SER( pCurve->Copy( crContext, rpCurve ));     NER(rpCurve);
  SmExtent1d sIvl = crEdge.GetInterval();
  if ( !( sIvl == rpCurve->GetNaturalInterval() ))
    { SER( rpCurve->Trim( sIvl )); }
  return SM_SUCCESS;
}

// The distance from crPnt to the surface of crFace continued, as far as it continues exactly:
// a plane's to the whole plane, any other surface's to its own domain.
SmStatus sm_DistanceToSurface(const SmFace & crFace, const SmPoint3d & crPnt, double & rdDist)
{
  rdDist = SM_BIG_DOUBLE;
  const SmSurface * pSurface = crFace.GetSurface();   NER(pSurface);
  if ( pSurface->IsKindOf( SmPlane_TYPE ))
    {
      const SmAxis2Placement & rFrame = ((const SmPlane*)pSurface)->GetPosition();
      rdDist = smos_Fabs( ( crPnt - rFrame.GetOrigin() ).Dot( rFrame.GetZAxis() ));
      return SM_SUCCESS;
    }
  SmBoolean bDropped = FALSE, bMulti = FALSE;
  SmPoint2d sUV;
  double    dGap = SM_BIG_DOUBLE;
  SER( pSurface->DropPoint( crPnt, pSurface->GetNaturalUVDomain(), NULL, bDropped, sUV, dGap, bMulti ));
  if ( bDropped )
    { rdDist = dGap; }
  return SM_SUCCESS;
}

// Points along a curve, its ends among them.
SmStatus sm_Samples(const SmCurve & crCurve, ULONG lCount, std::vector<SmPoint3d> & rPoints)
{
  SmExtent1d sIvl = crCurve.GetNaturalInterval();
  for (ULONG i=0; i<=lCount; i++)
    {
      SmPoint3d sPnt;
      SER( crCurve.EvaluatePoint( sIvl.Evaluate( (double)i / lCount ), sPnt ));
      rPoints.push_back( sPnt );
    }
  return SM_SUCCESS;
}

// A surface that is crFace's surface, whose domain holds every point of crPoints (all on it):
// a plane the same plane over a domain that covers them, any other surface a copy of itself
// when its own domain does (rbExact FALSE when it does not).
SmStatus sm_CoveringSurface(const SmContext & crContext, const SmFace & crFace, const std::vector<SmPoint3d> & crPoints,
                            double dTol, SmSurface *& rpSurface, SmBoolean & rbExact)
{
  rpSurface = NULL;
  rbExact   = TRUE;
  const SmSurface * pSurface = crFace.GetSurface();   NER(pSurface);
  if ( pSurface->IsKindOf( SmPlane_TYPE ))
    {
      // the plane through the same origin along the same axes, at unit scale: the same points
      const SmAxis2Placement & rFrame = ((const SmPlane*)pSurface)->GetPosition();
      SmPoint3d  sOrigin = rFrame.GetOrigin();
      SmVector3d sX = rFrame.GetXAxis(), sY = rFrame.GetYAxis();
      SmExtent2d sBox;
      SmBoolean  bFirst = TRUE;
      for (const SmPoint3d & rP : crPoints)
        {
          SmVector3d sD = rP - sOrigin;
          SmPoint2d  sUV( sD.Dot( sX ), sD.Dot( sY ));
          if ( bFirst ) { sBox = SmExtent2d( sUV, sUV ); bFirst = FALSE; }
          else          { sBox.AddPoint2d( sUV ); }
        }
      if ( bFirst )
        { SER( SM_ERR ); }
      double dMargin = 0.01 * sBox.GetSize().Length() + 10.0 * dTol;
      SmExtent2d sDomain( sBox.GetMin().x - dMargin, sBox.GetMin().y - dMargin,
                          sBox.GetMax().x + dMargin, sBox.GetMax().y + dMargin );
      rpSurface = new (crContext) SmPlane( sOrigin, sX, sY, SmVector2d( 1.0, 1.0 ), sDomain, &crContext );
      NER(rpSurface);
      return SM_SUCCESS;
    }
  // Within its own domain, or not exactly.
  for (const SmPoint3d & rP : crPoints)
    {
      double dDist = SM_BIG_DOUBLE;
      SER( sm_DistanceToSurface( crFace, rP, dDist ));
      if ( dDist > dTol )
        {
          rbExact = FALSE;
          return SM_SUCCESS;
        }
    }
  SER( pSurface->Copy( crContext, rpSurface ));      NER(rpSurface);
  return SM_SUCCESS;
}

// Whether two points coincide to dTol.
SmBoolean sm_Same(const SmPoint3d & crA, const SmPoint3d & crB, double dTol)
{
  return crA.DistanceBetween( crB ) <= dTol ? TRUE : FALSE;
}

} // namespace

/*******************************************************************//**
PURPOSE: SolidCraft change. Whether a variable round on planar faces holds
    its radius law, between the filleted edge's ends.

NOTES: For SubtractRemovalRegion. A rolling ball of radius r between two
    planes meeting at the interior angle a stands off the edge by
    r (1 / sin(a/2) - 1) where its section (the circle about the ball's
    centre through its two contacts) crosses the wedge's bisector. The fillet
    surface is traced on the faces' continuations; a trace that went wrong
    past them (clamped, or run into the law's continuation where it shrinks
    to nothing) gives a surface that does not, though it is a valid face
    (SmFilletExecutive.cpp, sm_ExtendTraceToItsRails). That point must lie
    on the surface, at seven points along the edge, to twenty times the
    solver's tolerance and a ten-thousandth of the radius; of the two angles the planes make, and the bisector
    either way, the one that fits. Rounds of other solvers, on faces that
    are not both planes, or along an edge that is not straight, pass.
***********************************************************************/
static SmStatus sm_HoldsItsLaw
  (SmFilletSolver * pSolver,   // in : the solver
   SmBoolean      & rbHolds)   // out: FALSE = the round's surface misses its law
{
  rbHolds = TRUE;
  SmVariableRadiusFS * pVRFS = SM_CAST_PTR( SmVariableRadiusFS, pSolver );
  if ( pVRFS == NULL )
    { return SM_SUCCESS; }
  SmFilletSurfaceGenerator * pFSG = pSolver->GetFilletSurfaceGenerator();
  if ( pFSG == NULL || pFSG->GetFilletSurfaceGeneratorType() != SM_FSG_CIRCULAR )
    { return SM_SUCCESS; }
  SmVector3d asNormal[2];
  for (ULONG k=0; k<2; k++)
    {
      SmOffsetSurface * pOffset = pSolver->GetSurface(k);                 NER(pOffset);
      const SmSurface * pBase   = pOffset->GetBaseSurface();              NER(pBase);
      SmPoint3d sAt;
      if ( !pBase->IsPlanar( SM_EFF_ZERO, &sAt, &asNormal[k] ))
        { return SM_SUCCESS; }
    }
  double dPhi = acos( smos_Max( -1.0, smos_Min( 1.0, asNormal[0].Dot( asNormal[1] ))));
  SmTArray<SmFilletGeom*> sFGs;
  pSolver->GetFilletGeoms( sFGs );
  if ( sFGs.GetSize() != 1 || sFGs[0] == NULL || sFGs[0]->GetFilletSurface() == NULL )
    { return SM_SUCCESS; }
  const SmSurface * pFillet = sFGs[0]->GetFilletSurface();
  SmEdge * pEdge = pSolver->GetEdgeuse(0)->GetEdge();                     NER(pEdge);
  NER(pEdge->GetCurve());
  if ( !pEdge->GetCurve()->IsLinear( SM_EFF_ZERO ))
    { return SM_SUCCESS; }
  double dTol = 20.0 * pSolver->GetThisApproxTol3d();
  for (ULONG i=1; i<8 && rbHolds; i++)
    {
      double    dT = pEdge->GetInterval().Evaluate( i / 8.0 );
      SmPoint3d sPnt;
      SER( pEdge->GetCurve()->EvaluatePoint( dT, sPnt ));
      double dRadius = smos_Fabs( pVRFS->GetFilletRadius( dT ));
      // The section's middle: along the wedge's bisector, r (1 / sin(a/2) - 1) off the edge.
      // The planes' normals are both outward (or inward), the bisector along their sum and the
      // angle pi less theirs, or one of each, the bisector along their difference.
      SmBoolean bFits = FALSE;
      for (ULONG c=0; c<4 && !bFits; c++)
        {
          SmVector3d sBisector = ( c < 2 ) ? asNormal[0] + asNormal[1] : asNormal[0] - asNormal[1];
          double     dAngle    = ( c < 2 ) ? SM_PI - dPhi : dPhi;
          double     dLen      = sBisector.Length();
          if ( dLen <= SM_EFF_ZERO_SQRT || dAngle <= SM_EFF_ZERO_SQRT )
            { continue; }
          double    dStandOff = dRadius * ( 1.0 / sin( dAngle / 2.0 ) - 1.0 );
          SmPoint3d sMiddle   = sPnt + sBisector * ( (( c % 2 ) ? -1.0 : 1.0 ) * dStandOff / dLen );
          SmBoolean bDropped = FALSE, bMulti = FALSE;
          SmPoint2d sUV;
          double    dGap = SM_BIG_DOUBLE;
          SER( pFillet->DropPoint( sMiddle, pFillet->GetNaturalUVDomain(), NULL, bDropped, sUV, dGap, bMulti ));
          // the solver's tolerance, and its circular sections' relative accuracy (a non-rational
          // cubic across the round, within about 2e-5 of the radius)
          bFits = ( bDropped && dGap <= dTol + 1.0e-4 * dRadius ) ? TRUE : FALSE;
        }
      rbHolds = bFits;
    }
  return SM_SUCCESS;
} // end sm_HoldsItsLaw

/*******************************************************************//**
PURPOSE: SolidCraft change. Note an end of a fillet that its Nx1 corner left
    uncut in a cliff overflow (SmFilletNx1Corner::CalcCornerEdgeGeom).

NOTES: pEndFace's surface, continued, cuts the fillet off there
    (CutFilletAtOpenEnd).
***********************************************************************/
SmStatus SmFilletExecutive::NoteOpenEnd
  (SmFilletSolver * pSolver,   // in : the solver of the filleted edge
   SmVertex       * pVertex,   // in : the edge's vertex at that end
   SmFace         * pEndFace)  // in : the end face there
{
  NER(pSolver);
  NER(pVertex);
  NER(pEndFace);
  SmFilletOpenEnd sEnd;
  sEnd.m_pSolver  = pSolver;
  sEnd.m_pVertex  = pVertex;
  sEnd.m_pEndFace = pEndFace;
  m_vOpenEnds.push_back( sEnd );
  return SM_SUCCESS;
} // end SmFilletExecutive::NoteOpenEnd

/*******************************************************************//**
PURPOSE: SolidCraft change. Cut the fillet off at an end its Nx1 corner left
    uncut, by the end face's surface continued.

NOTES: The fillet face of the end's solver runs on past the end, its rails
    past their side edges. The end face, planar, continues as its plane: the
    plane sections the fillet face (SmSurface::CreatePlanarSectionCurves,
    exact for the fillet's analytic surfaces), its section nearest the
    edge's end (the end rib) splits it (SmBrep::MergeCurvesOnSurface), and of
    the two pieces the one beyond the
    end goes: the piece kept holds the fillet's point nearest the middle of
    the filleted edge. A rib that does not split the face in two is
    SM_FILERR_NO_ROOM, as is an end face that is not planar.
***********************************************************************/
SmStatus SmFilletExecutive::CutFilletAtOpenEnd
  (const SmFilletOpenEnd & crEnd)   // in : the end left uncut
{
  NER(crEnd.m_pSolver);
  NER(crEnd.m_pEndFace);
  SmFilletSolver * pSolver = crEnd.m_pSolver;
  double           dTol    = pSolver->GetThisApproxTol3d();

  const SmSurface * pEndSurface = crEnd.m_pEndFace->GetSurface();   NER(pEndSurface);
  SmPoint3d  sAt;
  SmVector3d sNormal;
  if ( !pEndSurface->IsPlanar( SM_EFF_ZERO, &sAt, &sNormal ))
    {
      GetFilletErrorInfo()->NoteError( SM_FILERR_NO_ROOM, _T("Fillet Error: a cliff overflow runs past a fillet end whose end face is not planar") );
      return SM_ERR;
    }
  SmPlane sCut( sAt, sNormal, &m_crContext );

  // the fillet's point nearest the middle of the filleted edge
  SmEdge * pFilEdge = pSolver->GetEdgeuse(0)->GetEdge();             NER(pFilEdge);
  SmPoint3d sEdgeMid;
  SER( pFilEdge->GetCurve()->EvaluatePoint( pFilEdge->GetInterval().Evaluate( 0.5 ), sEdgeMid ));

  SmTArray<SmFilletGeom*> sFGs;
  pSolver->GetFilletGeoms( sFGs );
  for (ULONG k=0; k<sFGs.GetSize(); k++)
    {
      SmFilletGeom * pFG = sFGs[k];
      if ( pFG == NULL || pFG->GetFilletSurface() == NULL )
        { continue; }
      SmFace * pFace = SM_CAST_PTR( SmFace, pFG->GetFilletSurface()->GetOwner() );
      if ( pFace == NULL || pFace->GetBrep() != m_pFilletBrep )
        { continue; }
      SmSurface * pSurface = pFace->GetSurface();                      NER(pSurface);

      SmBoolean bDropped = FALSE, bMulti = FALSE;
      SmPoint2d sMidUV;
      double    dGap = SM_BIG_DOUBLE;
      SER( pSurface->DropPoint( sEdgeMid, pFace->GetUVDomain(), NULL, bDropped, sMidUV, dGap, bMulti ));
      if ( !bDropped )
        { SER( SM_ERR ); }
      SmPoint3d sKeep;
      SER( pSurface->EvaluatePoint( sMidUV, sKeep ));

      // section the face by the end's plane, and split it there
      SmApproxTol3d       sApprox( dTol );
      SmTArray<SmCurve*>  s3DCurves, sUVCurves;
      SmObjsDelete<SmCurve*> sClean3D( &s3DCurves ), sCleanUV( &sUVCurves );
      SER( pSurface->CreatePlanarSectionCurves( m_crContext, pFace->GetUVDomain(), sAt, sNormal, &sApprox, NULL, &s3DCurves, &sUVCurves ));
      if ( s3DCurves.GetSize() == 0 )
        { continue; }  // this face does not reach the end
      // The end rib: the section nearest the edge's end (a fillet traced far past its end
      // can cross the plane again out there, where it is cut off anyway).
      ULONG  lRib = 0;
      double dNearest = SM_BIG_DOUBLE;
      for (ULONG c=0; c<s3DCurves.GetSize(); c++)
        {
          SmBoolean bOnCurve = FALSE;
          double    dParam = 0.0, dDist = SM_BIG_DOUBLE;
          SER( s3DCurves[c]->DropPoint( s3DCurves[c]->GetNaturalInterval(), crEnd.m_pVertex->GetPoint(), NULL, SM_BIG_DOUBLE, NULL, bOnCurve, dParam, dDist ));
          if ( bOnCurve && dDist < dNearest )
            {
              dNearest = dDist;
              lRib     = c;
            }
        }
      SmTArray<SmCurve*> sRib3D, sRibUV;
      sRib3D.Add( s3DCurves[lRib] );
      if ( sUVCurves.GetSize() == s3DCurves.GetSize() )
        { sRibUV.Add( sUVCurves[lRib] ); }
      SmTArray<SmFace*> sNewFaces;
      SmTArray<SmEdge*> sOnEdges;
      SER( m_pFilletBrep->MergeCurvesOnSurface( *pSurface, dTol, sRib3D,
                                                ( sRibUV.GetSize() == 1 ) ? &sRibUV : NULL,
                                                sNewFaces, sOnEdges, FALSE, &sCut ));
      if ( sNewFaces.GetSize() == 0 )
        { continue; }  // the section only touches it
      if ( sNewFaces.GetSize() != 1 )
        {
          GetFilletErrorInfo()->NoteError( SM_FILERR_NO_ROOM, _T("Fillet Error: the end face's plane does not cut a fillet that overflows its end in two") );
          return SM_ERR;
        }
      // keep the piece holding the fillet's middle, delete the one beyond the end
      SmFace * apPieces[2] = { pFace, sNewFaces[0] };
      SmBoolean abHolds[2] = { FALSE, FALSE };
      for (ULONG p=0; p<2; p++)
        {
          SmPointClassification sClass;
          SER( apPieces[p]->Point3DClassify( sKeep, SmTol::GetZoneTol3d( apPieces[p] ), TRUE, sClass ));
          abHolds[p] = ( sClass.GetPointClass() != SM_PC_UNKNOWN ) ? TRUE : FALSE;
        }
      if ( abHolds[0] == abHolds[1] )
        {
          GetFilletErrorInfo()->NoteError( SM_FILERR_NO_ROOM, _T("Fillet Error: the end face's plane does not part a fillet that overflows its end from its middle") );
          return SM_ERR;
        }
      SmFace * pBeyond = abHolds[0] ? apPieces[1] : apPieces[0];
      if ( abHolds[1] )
        {
          // the kept piece is the fillet surface's face from here on
          pFG->SetFilletSurface( apPieces[1]->GetSurface() );
        }
      SER( m_pFilletBrep->DeleteFace( pBeyond, TRUE ));
    }
  return SM_SUCCESS;
} // end SmFilletExecutive::CutFilletAtOpenEnd

/*******************************************************************//**
PURPOSE: SolidCraft change. Close m_pFilletBrep, the sheet of fillet and
    corner faces, into the region the blends remove: a face on the surface of
    each target face its boundary runs on.

NOTES: Each lamina edge of the sheet lies on the surface (continued) of one
    target face round the filleted edges' vertices: its host. Consecutive
    lamina edges on one host form a chain. Where a chain ends the sheet's
    boundary passes onto the next host X, along the side edge s of the
    chain's host F that F shares with X and that starts at a filleted edge's
    vertex v: the chain's end lies on s, or on s's line past s's other end.
    The closing face on F's surface runs along the chain, back along s to v,
    along the filleted edges of F from v to the other end's vertex, and out
    along that end's side edge to the chain's start. A face of a box end (a
    rib between two side edges at one vertex) has no filleted edge between.
***********************************************************************/
SmStatus SmFilletExecutive::CloseRemovalRegion()
{
  NER(m_pFilletBrep);
  NER(m_pTargetBrep);

  // the filleted edges, their vertices, and the faces round them
  double dTol = 0.0;
  std::set<SmEdge*>   sBlended;
  std::set<SmVertex*> sBlendedVerts;
  std::set<SmFace*>   sHostSet;
  std::vector<SmFace*> vHosts;
  for (ULONG i=0; i<m_vFilletSolvers.GetSize(); i++)
    {
      SmFilletSolver * pSolver = m_vFilletSolvers[i];                NER(pSolver);
      dTol = smos_Max( dTol, pSolver->GetThisApproxTol3d() );
      SmEdge * pEdge = pSolver->GetEdgeuse(0)->GetEdge();            NER(pEdge);
      sBlended.insert( pEdge );
      for (SmVertex * pV : { pEdge->GetStartVertex(), pEdge->GetEndVertex() })
        {
          NER(pV);
          sBlendedVerts.insert( pV );
          SmTArray<SmFace*> sFaces;
          pV->GetFaces( sFaces );
          for (ULONG f=0; f<sFaces.GetSize(); f++)
            {
              if ( sHostSet.insert( sFaces[f] ).second )
                { vHosts.push_back( sFaces[f] ); }
            }
        }
    }
  double dOnTol = 10.0 * dTol;

  // each lamina edge of the sheet and its host
  SmTArray<SmEdge*> sSheetEdges;
  m_pFilletBrep->GetEdges( sSheetEdges );
  std::vector<SmEdge*>                          vLamina;
  std::map<SmEdge*, SmFace*>                    sHostOf;
  std::map<SmVertex*, std::vector<SmEdge*> >    sLaminaAt;
  for (ULONG i=0; i<sSheetEdges.GetSize(); i++)
    {
      SmEdge * pEdge = sSheetEdges[i];
      if ( !pEdge->IsLamina() )
        { continue; }
      SmCurve * pCurve = NULL;
      SER( sm_EdgeCurve( m_crContext, *pEdge, pCurve ));
      std::vector<SmPoint3d> vPoints;
      SmStatus eStat = sm_Samples( *pCurve, 4, vPoints );
      delete pCurve;
      SER( eStat );
      SmFace * pHost = NULL;
      for (SmFace * pFace : vHosts)
        {
          SmBoolean bOn = TRUE;
          for (size_t p=0; p<vPoints.size() && bOn; p++)
            {
              double dDist = SM_BIG_DOUBLE;
              SER( sm_DistanceToSurface( *pFace, vPoints[p], dDist ));
              bOn = ( dDist <= dOnTol ) ? TRUE : FALSE;
            }
          if ( !bOn )
            { continue; }
          if ( pHost != NULL )
            {
              GetFilletErrorInfo()->NoteError( SM_FILERR_NO_ROOM, _T("Fillet Error: a cliff overflow's fillet boundary lies on the surfaces of two faces") );
              return SM_ERR;
            }
          pHost = pFace;
        }
      if ( pHost == NULL )
        {
          GetFilletErrorInfo()->NoteError( SM_FILERR_NOT_EXACT, _T("Fillet Error: a cliff overflow's fillet boundary lies on no face's exact continuation") );
          return SM_ERR;
        }
      vLamina.push_back( pEdge );
      sHostOf[pEdge] = pHost;
      sLaminaAt[pEdge->GetStartVertex()].push_back( pEdge );
      sLaminaAt[pEdge->GetEndVertex()].push_back( pEdge );
    }

  // the chains of lamina edges on one host, each from its start to its end vertex
  struct sm_Chain
    {
      SmFace                 * m_pHost;
      std::vector<SmEdge*>     m_vEdges;
      std::vector<SmVertex*>   m_vVerts;  // m_vEdges.size() + 1
    };
  std::vector<sm_Chain> vChains;
  std::set<SmEdge*>     sDone;
  for (SmEdge * pFirst : vLamina)
    {
      if ( sDone.count( pFirst ))
        { continue; }
      SmFace * pHost = sHostOf[pFirst];
      // the next lamina edge on the same host at pV, other than pFrom (NULL at a chain's end)
      auto fNext = [&]( SmVertex * pV, SmEdge * pFrom, SmEdge *& rpNext ) -> SmStatus
        {
          rpNext = NULL;
          for (SmEdge * pE : sLaminaAt[pV])
            {
              if ( pE == pFrom || sHostOf[pE] != pHost )
                { continue; }
              if ( rpNext != NULL )
                { SER( SM_ERR ); }  // a branch of the boundary
              rpNext = pE;
            }
          return SM_SUCCESS;
        };
      sm_Chain sChain;
      sChain.m_pHost = pHost;
      std::vector<SmEdge*>   vBack;
      std::vector<SmVertex*> vBackVerts;
      // forward from pFirst's end vertex
      sChain.m_vEdges.push_back( pFirst );
      sChain.m_vVerts.push_back( pFirst->GetStartVertex() );
      sChain.m_vVerts.push_back( pFirst->GetEndVertex() );
      sDone.insert( pFirst );
      SmBoolean bClosed = FALSE;
      for (;;)
        {
          SmEdge * pNext = NULL;
          if ( fNext( sChain.m_vVerts.back(), sChain.m_vEdges.back(), pNext ) != SM_SUCCESS )
            {
              GetFilletErrorInfo()->NoteError( SM_FILERR_NO_ROOM, _T("Fillet Error: a cliff overflow's fillet boundary branches on a face") );
              return SM_ERR;
            }
          if ( pNext == NULL )
            { break; }
          if ( pNext == pFirst )
            { bClosed = TRUE; break; }
          sChain.m_vEdges.push_back( pNext );
          sChain.m_vVerts.push_back( pNext->GetOtherVertex( sChain.m_vVerts.back() ));
          sDone.insert( pNext );
        }
      if ( bClosed )
        {
          GetFilletErrorInfo()->NoteError( SM_FILERR_NO_ROOM, _T("Fillet Error: a cliff overflow's fillet boundary closes on one face") );
          return SM_ERR;
        }
      // backward from pFirst's start vertex
      SmVertex * pAt   = sChain.m_vVerts.front();
      SmEdge   * pFrom = pFirst;
      for (;;)
        {
          SmEdge * pNext = NULL;
          if ( fNext( pAt, pFrom, pNext ) != SM_SUCCESS )
            {
              GetFilletErrorInfo()->NoteError( SM_FILERR_NO_ROOM, _T("Fillet Error: a cliff overflow's fillet boundary branches on a face") );
              return SM_ERR;
            }
          if ( pNext == NULL )
            { break; }
          SmVertex * pOther = pNext->GetOtherVertex( pAt );
          sChain.m_vEdges.insert( sChain.m_vEdges.begin(), pNext );
          sChain.m_vVerts.insert( sChain.m_vVerts.begin(), pOther );
          sDone.insert( pNext );
          pAt   = pOther;
          pFrom = pNext;
        }
      vChains.push_back( sChain );
    }

  // the closing face of each chain
  for (const sm_Chain & rChain : vChains)
    {
      SmFace * pHost = rChain.m_pHost;
      sm_Loop  sLoop;

      // the chain, from its start to its end
      for (size_t e=0; e<rChain.m_vEdges.size(); e++)
        {
          SmCurve * pCurve = NULL;
          SER( sm_EdgeCurve( m_crContext, *rChain.m_vEdges[e], pCurve ));
          SER( sm_AddCurve( sLoop, pCurve, rChain.m_vVerts[e]->GetPoint(), rChain.m_vVerts[e+1]->GetPoint() ));
        }

      // At each end of the chain, the side edge the sheet's boundary passes along onto the next
      // host, and the filleted edges' vertex it starts at.
      SmVertex * apSideVert[2] = { NULL, NULL };    // [0] at the chain's start, [1] at its end
      SmCurve  * apConnector[2] = { NULL, NULL };   // from the vertex to the chain's end
      SmObjDelete sClean0, sClean1;
      for (ULONG k=0; k<2; k++)
        {
          SmVertex * pEnd  = (k == 0) ? rChain.m_vVerts.front() : rChain.m_vVerts.back();
          SmEdge   * pMine = (k == 0) ? rChain.m_vEdges.front() : rChain.m_vEdges.back();
          SmFace   * pNextHost = NULL;
          for (SmEdge * pE : sLaminaAt[pEnd])
            {
              if ( pE == pMine )
                { continue; }
              if ( pNextHost != NULL || sHostOf[pE] == pHost )
                {
                  GetFilletErrorInfo()->NoteError( SM_FILERR_NO_ROOM, _T("Fillet Error: a cliff overflow's fillet boundary branches on a face") );
                  return SM_ERR;
                }
              pNextHost = sHostOf[pE];
            }
          NER(pNextHost);
          SmPoint3d sP = pEnd->GetPoint();
          SmBoolean bInexact = FALSE;
          SmTArray<SmEdge*> sHostEdges;
          pHost->GetEdges( sHostEdges );
          for (ULONG s=0; s<sHostEdges.GetSize() && apSideVert[k] == NULL; s++)
            {
              SmEdge * pSide = sHostEdges[s];
              if ( sBlended.count( pSide ))
                { continue; }
              SmTArray<SmFace*> sSideFaces;
              pSide->GetFaces( sSideFaces );
              ULONG lAt = 0;
              if ( !sSideFaces.FindElement( pNextHost, lAt ))
                { continue; }
              for (SmVertex * pV : { pSide->GetStartVertex(), pSide->GetEndVertex() })
                {
                  if ( apSideVert[k] != NULL || pV == NULL || !sBlendedVerts.count( pV ))
                    { continue; }
                  SmPoint3d sV = pV->GetPoint();
                  if ( sm_Same( sV, sP, dOnTol ))
                    {
                      apSideVert[k] = pV;  // the chain ends at the vertex itself: no connector
                      continue;
                    }
                  // on the side edge itself
                  const SmCurve * pSideCurve = pSide->GetCurve();       NER(pSideCurve);
                  SmBoolean bDropped = FALSE;
                  double    dParam = 0.0, dDist = SM_BIG_DOUBLE, dParamV = 0.0, dDistV = SM_BIG_DOUBLE;
                  SER( pSideCurve->DropPoint( pSide->GetInterval(), sP, NULL, SM_BIG_DOUBLE, NULL, bDropped, dParam, dDist ));
                  if ( bDropped && dDist <= dOnTol )
                    {
                      SER( pSideCurve->DropPoint( pSide->GetInterval(), sV, NULL, SM_BIG_DOUBLE, NULL, bDropped, dParamV, dDistV ));
                      if ( !bDropped || dDistV > dOnTol )
                        { SER( SM_ERR ); }
                      SmCurve * pPiece = NULL;
                      SER( pSideCurve->Copy( m_crContext, pPiece ));      NER(pPiece);
                      SmExtent1d sPart( smos_Min( dParam, dParamV ), smos_Max( dParam, dParamV ));
                      SmStatus eStat = pPiece->Trim( sPart );
                      if ( eStat != SM_SUCCESS )
                        { delete pPiece; SER( eStat ); }
                      apSideVert[k]  = pV;
                      apConnector[k] = pPiece;
                      continue;
                    }
                  // on its line past its other end
                  SmVertex * pW = pSide->GetOtherVertex( pV );           NER(pW);
                  SmVector3d sAlong = pW->GetPoint() - sV;
                  SmVector3d sOut   = sP - sV;
                  double     dLen   = sAlong.Length();
                  if ( dLen <= dOnTol || sOut.Dot( sAlong ) <= 0.0 )
                    { continue; }
                  SmVector3d sCross = sOut * sAlong;
                  if ( sCross.Length() / dLen > dOnTol || sOut.Length() <= dLen )
                    { continue; }
                  if ( !pSideCurve->IsLinear( SM_EFF_ZERO ))
                    {
                      bInexact = TRUE;  // the chain's end lies along it, but it is no line
                      continue;
                    }
                  apSideVert[k]  = pV;
                  apConnector[k] = new (m_crContext) SmLine( sV, sP, 3, &m_crContext );
                  NER(apConnector[k]);
                }
            }
          if ( apSideVert[k] == NULL )
            {
              GetFilletErrorInfo()->NoteError( bInexact ? SM_FILERR_NOT_EXACT : SM_FILERR_NO_ROOM,
                                               bInexact ? _T("Fillet Error: a cliff overflow needs a side edge continued past its end that is not a line")
                                                        : _T("Fillet Error: a cliff overflow's fillet boundary leaves a face along no side edge of the filleted edges") );
              delete apConnector[0];
              delete apConnector[1];
              return SM_ERR;
            }
        }
      sClean0.SetObj( apConnector[0] );
      sClean1.SetObj( apConnector[1] );

      // from the chain's end back along its side edge to the vertex
      if ( apConnector[1] != NULL )
        {
          sClean1.Clear();
          SER( sm_AddCurve( sLoop, apConnector[1], rChain.m_vVerts.back()->GetPoint(), apSideVert[1]->GetPoint() ));
        }
      // along the filleted edges of the host from there to the other end's vertex
      SmTArray<SmEdge*> sHostEdges;
      pHost->GetEdges( sHostEdges );
      SmVertex * pAt   = apSideVert[1];
      SmEdge   * pFrom = NULL;
      for (ULONG lSteps=0; pAt != apSideVert[0]; lSteps++)
        {
          if ( lSteps > sHostEdges.GetSize() )
            { SER( SM_ERR ); }
          SmEdge * pNext = NULL;
          for (ULONG s=0; s<sHostEdges.GetSize(); s++)
            {
              SmEdge * pE = sHostEdges[s];
              if ( pE == pFrom || !sBlended.count( pE ) || ( pE->GetStartVertex() != pAt && pE->GetEndVertex() != pAt ))
                { continue; }
              if ( pNext != NULL )
                {
                  GetFilletErrorInfo()->NoteError( SM_FILERR_NO_ROOM, _T("Fillet Error: a cliff overflow's closing face meets two filleted edges at a vertex") );
                  return SM_ERR;
                }
              pNext = pE;
            }
          if ( pNext == NULL )
            {
              GetFilletErrorInfo()->NoteError( SM_FILERR_NO_ROOM, _T("Fillet Error: a cliff overflow's closing face is not bounded by filleted edges between its side edges") );
              return SM_ERR;
            }
          SmVertex * pOther = pNext->GetOtherVertex( pAt );            NER(pOther);
          SmCurve  * pCurve = NULL;
          SER( sm_EdgeCurve( m_crContext, *pNext, pCurve ));
          SER( sm_AddCurve( sLoop, pCurve, pAt->GetPoint(), pOther->GetPoint() ));
          pFrom = pNext;
          pAt   = pOther;
        }
      // out along the start's side edge to the chain's start
      if ( apConnector[0] != NULL )
        {
          sClean0.Clear();
          SER( sm_AddCurve( sLoop, apConnector[0], apSideVert[0]->GetPoint(), rChain.m_vVerts.front()->GetPoint() ));
        }

      // the host's surface over the loop, exactly
      std::vector<SmPoint3d> vPoints;
      for (const sm_LoopCurve & rC : sLoop.m_vCurves)
        { SER( sm_Samples( *rC.m_pCurve, 8, vPoints )); }
      SmSurface * pSurface = NULL;
      SmBoolean   bExact = TRUE;
      SER( sm_CoveringSurface( m_crContext, *pHost, vPoints, dOnTol, pSurface, bExact ));
      if ( !bExact )
        {
          GetFilletErrorInfo()->NoteError( SM_FILERR_NOT_EXACT, _T("Fillet Error: a cliff overflow needs a face's surface continued past its domain, which is not a plane") );
          return SM_ERR;
        }
      NER(pSurface);
      SmObjDelete sCleanSurface( pSurface );

      // anticlockwise on the surface: the loop's points in its parameters, by their signed area
      double dArea = 0.0;
      std::vector<SmPoint2d> vUV;
      for (const sm_LoopCurve & rC : sLoop.m_vCurves)
        {
          std::vector<SmPoint3d> vAlong;
          SER( sm_Samples( *rC.m_pCurve, 8, vAlong ));
          if ( rC.m_eOrient == SM_OT_OPPOSITE )
            { std::reverse( vAlong.begin(), vAlong.end() ); }
          for (size_t p=0; p+1<vAlong.size(); p++)
            {
              SmBoolean bDropped = FALSE, bMulti = FALSE;
              SmPoint2d sUV;
              double    dGap = SM_BIG_DOUBLE;
              SER( pSurface->DropPoint( vAlong[p], pSurface->GetNaturalUVDomain(), NULL, bDropped, sUV, dGap, bMulti ));
              if ( !bDropped )
                { SER( SM_ERR ); }
              vUV.push_back( sUV );
            }
        }
      for (size_t p=0; p<vUV.size(); p++)
        {
          const SmPoint2d & rA = vUV[p];
          const SmPoint2d & rB = vUV[(p+1) % vUV.size()];
          dArea += rA.x * rB.y - rB.x * rA.y;
        }
      if ( dArea < 0.0 )
        {
          std::reverse( sLoop.m_vCurves.begin(), sLoop.m_vCurves.end() );
          for (sm_LoopCurve & rC : sLoop.m_vCurves)
            { rC.m_eOrient = ( rC.m_eOrient == SM_OT_SAME ) ? SM_OT_OPPOSITE : SM_OT_SAME; }
        }

      SmTArray<ULONG>        sLoopCounts;
      SmTArray<SmCurve*>     s3DCurves;
      SmTArray<SmOrientType> sOrients;
      SmTArray<SmPoint3d>    sLoopPoints;
      sLoopCounts.Add( (ULONG)sLoop.m_vCurves.size() );
      for (const sm_LoopCurve & rC : sLoop.m_vCurves)
        {
          s3DCurves.Add( rC.m_pCurve );
          sOrients.Add( rC.m_eOrient );
        }
      SmRegion * pNewRegion = NULL;
      SmShell  * pNewShell  = NULL;
      SmFace   * pNewFace   = NULL;
      // the new face owns its curves and its surface
      for (sm_LoopCurve & rC : sLoop.m_vCurves)
        { rC.m_pCurve = NULL; }
      sCleanSurface.Clear();
      SER( m_pFilletBrep->MakeFaceWithCurves( m_pFilletBrep->GetInfiniteRegion(), sLoopCounts, &s3DCurves, NULL, sOrients,
                                              sLoopPoints, pSurface, pSurface->GetNaturalUVDomain(), SM_OT_SAME,
                                              pNewRegion, pNewShell, pNewFace ));
      NER(pNewFace);
    }
  return SM_SUCCESS;
} // end SmFilletExecutive::CloseRemovalRegion

/*******************************************************************//**
PURPOSE: SolidCraft change. In a cliff overflow, subtract the region the
    blends remove from the target: the blends built on their faces'
    continuations, then cut by what lies beyond the cliffs.

NOTES: See the file's head. m_pFilletBrep, closed into that region and
    oriented as a solid, is subtracted from m_pTargetBrep by
    SmMerge::ManifoldBoolean (SM_BO_DIFFERENCE), which modifies the target
    in place; the region is deleted after.
***********************************************************************/
SmStatus SmFilletExecutive::SubtractRemovalRegion()
{
  NER(m_pFilletBrep);
  NER(m_pTargetBrep);

  for (ULONG i=0; i<m_vFilletSolvers.GetSize(); i++)
    {
      SmBoolean bHolds = TRUE;
      SER( sm_HoldsItsLaw( m_vFilletSolvers[i], bHolds ));
      if ( !bHolds )
        {
          GetFilletErrorInfo()->NoteError( SM_FILERR_NO_ROOM, _T("Fillet Error: a variable round traced on its faces' continuations does not hold its radius law") );
          return SM_ERR;
        }
    }
  for (const SmFilletOpenEnd & rEnd : m_vOpenEnds)
    { SER( CutFilletAtOpenEnd( rEnd )); }
  SER( CloseRemovalRegion() );

  // stitch the closing faces to the sheet, and orient the closed shell as a solid
  ULONG  lStitched = 0, lLamina = 0;
  double dVertexGap = 0.0, dEdgeGap = 0.0;
  SER( m_pFilletBrep->StitchFaces( m_pFilletBrep->GetTolerance(), lStitched, lLamina, dVertexGap, dEdgeGap ));
  SmBoolean bMaybeOpen = FALSE;
  if ( lLamina == 0 )
    { SER( m_pFilletBrep->OrientTrimmedSurfaces( TRUE, bMaybeOpen, FALSE )); }
  if ( lLamina != 0 || bMaybeOpen || !m_pFilletBrep->IsManifoldSolid() )
    {
      GetFilletErrorInfo()->NoteError( SM_FILERR_BAD_MERGE, _T("Fillet Error: a cliff overflow's removal region does not close into a solid") );
      return SM_ERR;
    }

  // the fillet edges' attribute stays with the executive's own Breps
  SmTArray<SmEdge*> sEdges;
  m_pFilletBrep->GetEdges( sEdges );
  for (ULONG i=0; i<sEdges.GetSize(); i++)
    {
      if ( sEdges[i]->FindAttribute( SM_AI_ALL_FILLET_EDGES ) == m_pAttributeE )
        { sEdges[i]->RemoveAttribute( m_pAttributeE ); }
    }

  SmBrep * pRegion = m_pFilletBrep;
  m_pFilletBrep = NULL;
  SmObjDelete sDeleteRegion( pRegion );
  SmMerge sMerge( *m_pTargetBrep->GetContext(), m_pTargetBrep, pRegion );
  sMerge.SetKeepOtherBrep( TRUE );
  SmBrep * pResult = NULL;
  SmStatus eStat = m_pTargetBrep->IsManifoldSolid()
                   ? sMerge.ManifoldBoolean( SM_BO_DIFFERENCE, pResult )
                   : sMerge.NonManifoldBoolean( SM_BO_DIFFERENCE, pResult );
  if ( eStat != SM_SUCCESS || pResult != m_pTargetBrep )
    {
      GetFilletErrorInfo()->NoteError( SM_FILERR_BAD_MERGE, _T("Fillet Error: subtracting a cliff overflow's removal region failed") );
      return SM_ERR;
    }
  return SM_SUCCESS;
} // end SmFilletExecutive::SubtractRemovalRegion
