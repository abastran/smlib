// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
// SPDX-License-Identifier: Apache-2.0

/**********************************************************************//**
* FILE NAME --- SmFilletRemoval.cpp
* PURPOSE: Source file for SmFilletRemoval class.
**********************************************************************/

#include "StdAfx.h"

#include <SmFilletRemoval.h>

#ifndef __SMBREP_H__
#include <SmBrep.h>
#endif

#ifndef __SMVERTEX_H__
#include <SmVertex.h>
#endif

#include <SmBSplineSurface.h>
#include <SmBSplineCurve.h>
#include <SmEdgeuse.h>
#include <SmSolutionArray.h>
#include <SmAssertArray.h>

/*******************************************************************//**
PURPOSE: SmFilletRemovel.

NOTES:
  The surfaces of the four faces adjacent to the Fillet Face
  being removed must all be (derived from) SmBSplineSurface,
  and the four edge curves coming into the Fillet Face must
  be (derived from) SmBSplineCurve.

 About the algorithm:
 We first collect all of the topology, and then
 do all the intersections, etc., before we actually
 change anything, in case something fails before
 in those earlier operations.

 About surface extensions:
 The tangent surfaces, and sometimes the end surfaces,
 will have to be extended to fill in where the fillet trimmed them.
 Our method is:
 - In the information-gathering phase, just evaluate them
   outside their domain.
 - In the modification phase, do the actual extensions.
   Use a 'natural' extension (SM_CT_CINFINITY) because that gives
   the same surface that's evaluated in the first phase,
   when evaluating outside the domain.
 Advantages to this are:
 - If anything goes wrong before the modification phase,
   then the surfaces aren't changed;
 - This way we know exactly how far to extend them;
 - Natural extension doesn't introduce any knots or discontinuities,
   and it will behave well enough (and as expected) for the amount
   that they're being extended.

***********************************************************************/


/*******************************************************************//**
PURPOSE: SmFilletRemovel constructor.

NOTES:
***********************************************************************/
SmFilletRemoval::SmFilletRemoval( SmFace *pFaceToDelete, double d3dTolerance )
{
  m_pFaceToDelete = pFaceToDelete;
  m_d3dTol        = d3dTolerance;

  m_aTanEdges[0] = m_aTanEdges[1] = NULL;
  m_aTanFaces[0] = m_aTanFaces[1] = NULL;
  m_aEndEdges[0] = m_aEndEdges[1] = NULL;
  m_aEndFaces[0] = m_aEndFaces[1] = NULL;

  m_pNewEdgeCurve3d = NULL;
  m_aNewEdgeCurveUVs[0] = m_aNewEdgeCurveUVs[1] = NULL;

  m_dNewEdgeDeviation = -1;
}
/*******************************************************************//**
PURPOSE: SmFilletRemovel destructor.

NOTES:
***********************************************************************/
SmFilletRemoval::~SmFilletRemoval()
{
  // Don't delete these: they're in the new topology.
  m_pNewEdgeCurve3d     = NULL;
  m_aNewEdgeCurveUVs[0] = NULL;
  m_aNewEdgeCurveUVs[1] = NULL;
}
/*******************************************************************//**
PURPOSE: SmFilletCornerInfo constructor.

NOTES:
***********************************************************************/
SmFilletCornerInfo::SmFilletCornerInfo()
{
    m_pNewSideEdgeExtension3d  = NULL;
    m_pNewSideEdgeExtEndFaceUV = NULL;
    m_pNewSideEdgeExtTanFaceUV = NULL;
    m_pSideEdge = NULL;
}
/*******************************************************************//**
PURPOSE: SmFilletCornerInfo destructor.

NOTES:
***********************************************************************/
SmFilletCornerInfo::~SmFilletCornerInfo()
{
  // Don't delete these: they're in the new topology.
  m_pNewSideEdgeExtension3d  = NULL;
  m_pNewSideEdgeExtEndFaceUV = NULL;
  m_pNewSideEdgeExtTanFaceUV = NULL;
  m_pSideEdge = NULL;
}


/*******************************************************************//**
PURPOSE: Main method of SmFilletRemovel: delete the Fillet and close up the hole.

NOTES:
***********************************************************************/
SmStatus SmFilletRemoval::RemoveFace()
{
#ifdef SM_DEBUG_CODE
  SmBrep *pBrep = m_pFaceToDelete->GetBrep(); NER( pBrep );
  SmBoolean bDebugMe = FALSE;
if ( bDebugMe ) {
  SM_ASSERT_VALID( pBrep );
  pBrep->ValidatePointers();
  pBrep->Dump();
}
#endif

  // Collect the two tangent edges, two end edges,
  // and the Faces across those edges.

  SER( this->CollectFilletTopology() );

  ULONG lWhichEnd, lWhichSide;

  // For both ends:
  for ( lWhichEnd = 0; lWhichEnd < 2; lWhichEnd++ )
  {
      SER( this->CalcNewCornerVertices( lWhichEnd ) );

      // For both sides at this end:
      for ( lWhichSide = 0; lWhichSide < 2; lWhichSide++ )
      {
          SER( this->FindSideEdgeExtensions( lWhichEnd, lWhichSide ) );
      }
  }

  // Intersect the surfaces of the tangent Faces (which have been extended).
  SER( this->IntersectTangentSurfaces() );


  // Ok, we've collected all of the info without incident.
  // Start digging into the Brep.

  SER( this->DoRemoval() );

#ifdef SM_DEBUG_CODE
if ( bDebugMe ) {
  SM_ASSERT_VALID( pBrep );
  pBrep->ValidatePointers();
  pBrep->Dump();
}
#endif

  return SM_SUCCESS;

} // end Execute


/*******************************************************************//**
PURPOSE: Gather information about the topology of the Fillet Face and its neighbors.

NOTES: Fillet Face must have exactly four edges, with two opposite
   Edges being tangent Edges, and the other two non-tangent.

   Input: just the results of the constructor: m_pFaceToDelete, and dTol.

   Output: fills in:
   - two End Edges and two Tangent (side, rail) Edges:
     m_aTanEdges[2]
     m_aEndEdges[2]
   - The adjacent Faces, across from those four Edges:
     m_aTanFaces[2]
     m_aEndFaces[2]
   - Orientation information:
     m_eRailParam
     m_bEndsSwapped
     m_bSidesSwapped
   On each of the four Corners:
   - the Vertex
     m_pVertex
   - uv values of m_pVertex in three Faces:
     m_sUVFilSurf
     m_sUVEndSurfOld
     m_sUVTanSurfOld
   - the Side Edge and its orientation (the third Edge incident at the Vertex)
     m_pSideEdge
     m_bHighEndOfSideEdge

***********************************************************************/
SmStatus SmFilletRemoval::CollectFilletTopology()
{
  // Identify two tangent Edges and two end Edges.

  // Assumption: GetEdges() returns Edges ordered around the Loop.
  SmTArray< SmEdge* > sEdges;
  m_pFaceToDelete->GetEdges( sEdges );
  if ( sEdges.GetSize() != 4 ) {
      SER( SM_ERR );
  }

  double dTangencyTolDeg = 2.0;  // [for fillet1018.smb, Face[5]; was 0.5]

  if ( sEdges[0]->IsTangentEdge( dTangencyTolDeg ) )
  {
      if (  sEdges[1]->IsTangentEdge( dTangencyTolDeg ) ) SER( SM_ERR );
      if ( !sEdges[2]->IsTangentEdge( dTangencyTolDeg ) ) SER( SM_ERR );
      if (  sEdges[3]->IsTangentEdge( dTangencyTolDeg ) ) SER( SM_ERR );

      this->m_aTanEdges[0] = sEdges[0];
      this->m_aTanEdges[1] = sEdges[2];
      this->m_aEndEdges[0] = sEdges[1];
      this->m_aEndEdges[1] = sEdges[3];
  }
  else if ( sEdges[1]->IsTangentEdge( dTangencyTolDeg ) )
  {
      if (  sEdges[0]->IsTangentEdge( dTangencyTolDeg ) ) SER( SM_ERR );
      if (  sEdges[2]->IsTangentEdge( dTangencyTolDeg ) ) SER( SM_ERR );
      if ( !sEdges[3]->IsTangentEdge( dTangencyTolDeg ) ) SER( SM_ERR );

      this->m_aTanEdges[0] = sEdges[1];
      this->m_aTanEdges[1] = sEdges[3];
      this->m_aEndEdges[0] = sEdges[0];
      this->m_aEndEdges[1] = sEdges[2];
  }
  else
  {
      SER( SM_ERR );
  }

  // Get the four corner vertices.

  ULONG i, lWhichEnd, lWhichSide;

  // First collect both verts from all four Edges: eight total,
  // of which there will be four matching pairs.
  SmVertex * aEndVerts[2][2];
  SmVertex * aTanVerts[2][2];

  SmTArray< SmVertex* > aVerts;
  for ( i = 0; i < 2; i++ )
  {
      m_aTanEdges[i]->GetVertices( aVerts );
      if ( aVerts.GetSize() != 2 )
        { SER( SM_ERR ); }
      aTanVerts[i][0] = aVerts[0];
      aTanVerts[i][1] = aVerts[1];

      m_aEndEdges[i]->GetVertices( aVerts );
      if ( aVerts.GetSize() != 2 )
        { SER( SM_ERR ); }
      aEndVerts[i][0] = aVerts[0];
      aEndVerts[i][1] = aVerts[1];
  }

  // Now sort these out.
  for ( lWhichEnd = 0; lWhichEnd < 2; lWhichEnd++ )
  {
      for ( lWhichSide = 0; lWhichSide < 2; lWhichSide++ )
      {
          // Naming:
          SmFilletCornerInfo &rCorn = m_aCornerInfo[lWhichEnd][lWhichSide];

          if ( aEndVerts[lWhichEnd][0] == aTanVerts[lWhichSide][0] )
          {
              rCorn.m_pVertex = aEndVerts[lWhichEnd][0];
          }
          else if ( aEndVerts[lWhichEnd][0] == aTanVerts[lWhichSide][1] )
          {
              rCorn.m_pVertex = aEndVerts[lWhichEnd][0];
          }
          else if ( aEndVerts[lWhichEnd][1] == aTanVerts[lWhichSide][0] )
          {
              rCorn.m_pVertex = aEndVerts[lWhichEnd][1];
          }
          else if ( aEndVerts[lWhichEnd][1] == aTanVerts[lWhichSide][1] )
          {
              rCorn.m_pVertex = aEndVerts[lWhichEnd][1];
          }
          else
          {
              SER( SM_ERR );
          }
      }
  } // end double loop on two ends / two sides.


  // Find the adjacent Faces.
  // Also, all Surfaces have to be SmBSplineSurfaces.
  // because of CreateExtendedSurface().
  SmTArray< SmFace* >sFaces;
  SmFace *pFilFace = this->m_pFaceToDelete;

  for ( lWhichEnd = 0; lWhichEnd < 2; lWhichEnd++ )
  {
      this->m_aTanEdges[lWhichEnd]->GetFaces( sFaces );
      if ( sFaces.GetSize() != 2 )
        { SER( SM_ERR ); }

      if ( sFaces[0] == pFilFace ) {
          this->m_aTanFaces[lWhichEnd] = sFaces[1];
      }
      else if ( sFaces[1] == pFilFace ) {
          this->m_aTanFaces[lWhichEnd] = sFaces[0];
      }
      else
        { SER( SM_ERR ); }

      SmSurface *pSrf = this->m_aTanFaces[lWhichEnd]->GetSurface();
      if ( pSrf == NULL || ! pSrf->IsKindOf( SmBSplineSurface_TYPE ) )
        { SER( SM_ERR ); }

      this->m_aEndEdges[lWhichEnd]->GetFaces( sFaces );
      if ( sFaces.GetSize() != 2 )
        { SER( SM_ERR ); }

      if ( sFaces[0] == pFilFace ) {
          this->m_aEndFaces[lWhichEnd] = sFaces[1];
      }
      else if ( sFaces[1] == pFilFace ) {
          this->m_aEndFaces[lWhichEnd] = sFaces[0];
      }
      else
        { SER( SM_ERR ); }

      pSrf = this->m_aEndFaces[lWhichEnd]->GetSurface();
      if ( pSrf == NULL || ! pSrf->IsKindOf( SmBSplineSurface_TYPE ) )
        { SER( SM_ERR ); }

  } // end for loop on two ends of fillet

  // Fill in surface uv info for the vertices.
  // There are four Vertices, each is in three Faces.
  // Also find the side edges: the third edge coming into each corner vertex.
  SmVertexuse *pVU;

  for ( lWhichEnd = 0; lWhichEnd < 2; lWhichEnd++ )
  {
      for ( lWhichSide = 0; lWhichSide < 2; lWhichSide++ )
      {
          // For syntactic convenience: rename this corner-info object.
          SmFilletCornerInfo & rCorn = m_aCornerInfo[lWhichEnd][lWhichSide];

          SmVertex *pVtx = rCorn.m_pVertex;

          // Fillet Face:
          pVU = pVtx->GetVertexuseOfFace( m_pFaceToDelete );  NER( pVU );
          SER( pVU->ComputeUVPoint( rCorn.m_sUVFilSurf ));

          // End Face:
          pVU = pVtx->GetVertexuseOfFace( m_aEndFaces[lWhichEnd] );   NER( pVU );
          SER( pVU->ComputeUVPoint( rCorn.m_sUVEndSurfOld ));

          // Tangent Face:
          pVU = pVtx->GetVertexuseOfFace( m_aTanFaces[lWhichSide] );  NER( pVU );
          SER( pVU->ComputeUVPoint( rCorn.m_sUVTanSurfOld ));

          // Also find the side edges.
          pVtx->GetEdges( sEdges );
          if ( sEdges.GetSize() != 3 )
          {
              // More than three edges: we do things differently; see the .h file.
              rCorn.m_pSideEdge = NULL;
          }
          else
          {
              for ( i = 0; i < 3; i++ )
              {
                  SmEdge *pThisEdge = sEdges[i];
                  if ( pThisEdge == m_aEndEdges[ lWhichEnd  ] ) continue;
                  if ( pThisEdge == m_aTanEdges[ lWhichSide ] ) continue;
                  rCorn.m_pSideEdge = pThisEdge;
                  break;
              }
              NER( rCorn.m_pSideEdge );

              // We will also want to know which end of the side edge we're at.
              SmVertex *pSideVtx = rCorn.m_pSideEdge->GetStartVertex();
              if ( pSideVtx == rCorn.m_pVertex )
              {
                  rCorn.m_bHighEndOfSideEdge = FALSE;
              }
              else
              {
                  pSideVtx = rCorn.m_pSideEdge->GetOtherVertex( pSideVtx );
                  if ( pSideVtx == rCorn.m_pVertex )
                  {
                      rCorn.m_bHighEndOfSideEdge = TRUE;
                  }
                  else
                    { SER( SM_ERR ); }  // If this happens it's our own fault.
              }
          } // end if exactly 3 Edges for this vertex

      } // end inner loop
  } // end double loop on four corners, filling in surface uv and side-edge info

  // Now set the orientation of the fillet surface w.r.t. our Ends and Sides:
  // members m_eRailParam, m_bEndsSwapped, and m_bSidesSwapped.
  // Look at the directions of the vectors between corner uv values.
  SmVector2d sVecEndToEnd(     m_aCornerInfo[1][0].m_sUVFilSurf
                             - m_aCornerInfo[0][0].m_sUVFilSurf );
  SmVector2d sVecSideToSide(   m_aCornerInfo[0][1].m_sUVFilSurf
                             - m_aCornerInfo[0][0].m_sUVFilSurf );
  sVecEndToEnd.Unitize();
  sVecSideToSide.Unitize();

  if ( smos_Fabs( sVecEndToEnd.x ) > smos_Fabs( sVecEndToEnd.y ) )
  {
      m_eRailParam = SM_SP_U;
      SM_ASSERT( smos_Fabs( sVecSideToSide.y ) > smos_Fabs( sVecSideToSide.x ));
      m_bEndsSwapped  = sVecEndToEnd.x   < 0;
      m_bSidesSwapped = sVecSideToSide.y < 0;
  }
  else
  {
      m_eRailParam = SM_SP_V;
      SM_ASSERT( smos_Fabs( sVecSideToSide.x ) > smos_Fabs( sVecSideToSide.y ));
      m_bEndsSwapped  = sVecEndToEnd.y   < 0;
      m_bSidesSwapped = sVecSideToSide.x < 0;
  }

  return SM_SUCCESS;

} // end SmFilletRemoval::CollectFilletTopology


/*******************************************************************//**
PURPOSE: Calculate what will be the new corner vertices (one at each end).

NOTES:
   Does a surface/surface/surface intersection on each end,
   between the end Face and the two tangent Faces.

   Input: The results of CollectFilletTopology().

   Output: fills in:
   - two new corner vertex positions, and their new tolerances:
     m_aNewVtxPts[2]
     m_aNewVtxTols[2]
   - The uv's of these points in the two end Faces:
     m_aUVEndSurfNew[2]
   On each of the four Corners:
   - the new uv value of the new vertex position in the Tangent Face:
     m_sUVTanSurfNew

***********************************************************************/
SmStatus SmFilletRemoval::CalcNewCornerVertices( ULONG lWhichEnd )
{
  // Find what will be the new ends of the new Edge:
  // on each end, srf/srf/srf intersect the two tangent Faces and the end Face.

  SmSurface * pTanSurfs[2];
  pTanSurfs[0] = m_aTanFaces[0]->GetSurface();  NER( pTanSurfs[0] );
  pTanSurfs[1] = m_aTanFaces[1]->GetSurface();  NER( pTanSurfs[1] );
  SmExtent2d sTanSurfDomains[2];
  sTanSurfDomains[0] = pTanSurfs[0]->GetNaturalUVDomain();
  sTanSurfDomains[1] = pTanSurfs[1]->GetNaturalUVDomain();
  SmBoolean bFoundAnswer;
  SmSolution sSol;

  SmSurface *pEndSurf  = m_aEndFaces[lWhichEnd]->GetSurface();  NER(pEndSurf);
  SmExtent2d sEndSurfDomain = pEndSurf->GetNaturalUVDomain();

  // Temporarily enable OutOfBounds evaluations for all SmBSplineSurfaces.
  SmBoolean bTmp0 = FALSE, bTmp1 = FALSE, bTmp2 = FALSE ; 
  SmBSplineSurface * pBSplineSurface0 = SM_CAST_PTR(SmBSplineSurface, pTanSurfs[0]) ;
  SmBSplineSurface * pBSplineSurface1 = SM_CAST_PTR(SmBSplineSurface, pTanSurfs[1]) ;
  SmBSplineSurface * pBSplineSurface2 = SM_CAST_PTR(SmBSplineSurface, pEndSurf    ) ;
  SmTemporaryChangeValue<SmBoolean> sClean0(pBSplineSurface0 ? pBSplineSurface0->GetOutOfBoundsEnabled() : bTmp0, TRUE) ;
  SmTemporaryChangeValue<SmBoolean> sClean1(pBSplineSurface1 ? pBSplineSurface1->GetOutOfBoundsEnabled() : bTmp1, TRUE) ;
  SmTemporaryChangeValue<SmBoolean> sClean2(pBSplineSurface2 ? pBSplineSurface2->GetOutOfBoundsEnabled() : bTmp2, TRUE) ;

  // Also expand the domains, for srf/srf/srf intersect.
  sTanSurfDomains[0].ExpandRelative( 2.0 );
  sTanSurfDomains[1].ExpandRelative( 2.0 );
  sEndSurfDomain.    ExpandRelative( 2.0 );

  SmPoint2d sUVEnd = ( m_aCornerInfo[lWhichEnd][0].m_sUVEndSurfOld +
                       m_aCornerInfo[lWhichEnd][1].m_sUVEndSurfOld ) / 2.0;
  SmPoint2d sUVTanSurf0( m_aCornerInfo[lWhichEnd][0].m_sUVTanSurfOld );
  SmPoint2d sUVTanSurf1( m_aCornerInfo[lWhichEnd][1].m_sUVTanSurfOld );

  SmStatus eStat = pEndSurf->LocalSurfaceSurfaceIntersect( sEndSurfDomain,
        *(pTanSurfs[0]), sTanSurfDomains[0], *(pTanSurfs[1]), sTanSurfDomains[1],
        sUVEnd, sUVTanSurf0, sUVTanSurf1,
        bFoundAnswer, sSol
  );

  if ( eStat != SM_SUCCESS || ! bFoundAnswer )
    { SER( SM_ERR ); }

  // Store the new vtx position and the three surface uv's.
  // Surface uv's first.
  sUVEnd.Set     ( sSol.m_vStart[0], sSol.m_vStart[1] );
  sUVTanSurf0.Set( sSol.m_vStart[2], sSol.m_vStart[3] );
  sUVTanSurf1.Set( sSol.m_vStart[4], sSol.m_vStart[5] );

  m_aUVEndSurfNew[lWhichEnd] = sUVEnd;

  m_aCornerInfo[lWhichEnd][0].m_sUVTanSurfNew = sUVTanSurf0;
  m_aCornerInfo[lWhichEnd][1].m_sUVTanSurfNew = sUVTanSurf1;

  // Evaluate all three surfaces, check 3d distances,
  // and use the average position.

  SmPoint3d sVtxPtEnd, sVtxPtTan0, sVtxPtTan1;
  pEndSurf->EvaluatePoint    ( sUVEnd,      sVtxPtEnd  );
  pTanSurfs[0]->EvaluatePoint( sUVTanSurf0, sVtxPtTan0 );
  pTanSurfs[1]->EvaluatePoint( sUVTanSurf1, sVtxPtTan1 );

  double dDist00 = sVtxPtEnd.DistanceBetween ( sVtxPtTan0 );
  double dDist01 = sVtxPtEnd.DistanceBetween ( sVtxPtTan1 );
  double dDist11 = sVtxPtTan0.DistanceBetween( sVtxPtTan1 );
  double dTol = m_d3dTol;
  if (  dDist00 > dTol || dDist01 > dTol || dDist11 > dTol )
    { SER( SM_ERR ); }

  SmPoint3d sNewVtxPt = ( sVtxPtEnd + sVtxPtTan0 + sVtxPtTan1 ) / 3.0;

  m_aNewVtxPts [ lWhichEnd ] = sNewVtxPt;
  m_aNewVtxTols[ lWhichEnd ] = smos_3Max( dDist00, dDist01, dDist11 );

  return SM_SUCCESS;

} // end SmFilletRemoval::CalcNewCornerVertices


/*******************************************************************//**
PURPOSE: Find the extension curves for the four side Edges.

NOTES:
   This loops over the four corners.

   First checks whether the curve of the side Edge is already long enough.
   If it is, then just calculates the parameter of the new Vertex on it,
   for updating the side Edge's domain.
   If the side Edge does need extension, calculates the extension curve,
   which is the intersection of the end Face and the tangent Face (at
   this corner) running from the original Vertex to the new Vertex.

   Input: The results of CollectFilletTopology();
                     and CalcNewCornerVertices().

   Output: fills in:
   On each of the four Corners:
   - Either: The new extension curve pieces, one 3d and two uv curves:
       m_pNewSideEdgeExtension3d
       m_pNewSideEdgeExtEndFaceUV
       m_pNewSideEdgeExtTanFaceUV
   - Or: the parameter of the new corner Vertex on the side Edge:
     - m_dNewSideEdgeParam

   Which of the two is set is signalled by the curve-extension pointers:
   if m_pNewSideEdgeExtension3d is Null, then m_dNewSideEdgeParam is set.

   Note: if there are more than three Edges incident on a corner vertex,
   then no Edges get extended.  In that case, the end Edge will not
   be deleted, but its geometry will be replaced with the intersection
   curve of the end Face and the tangent Face.  So in that case, just
   store the intersection, don't try extending.

***********************************************************************/
SmStatus SmFilletRemoval::FindSideEdgeExtensions( ULONG lWhichEnd, ULONG lWhichSide )
{
  // The curves of the side edges coming in to the new corner vertex
  // might not be long enough.  Check whether they are, and if not,
  // intersect the last bit and join the curves.
  // We will have to update the domains of the Edges accordingly.
  // Drop the new vertex point to each curve.  This will tell us
  // whether the curve is long enough, and if so, give us the
  // parameter value on it, for updating the Edge domain.

  // We'll use the parameter value at the appropriate end of the side edge
  // as a guess parameter.
  double    dTol      = m_d3dTol;
  double    dCrvDist  = 2*dTol;
  double    dCrvParam = 0;
  SmBoolean bSuccess  = FALSE;
  SmStatus  eStat     = SM_ERR;

  SmBSplineCurve  *pIntCrv3d = NULL, *pIntCrvUVEnd = NULL, *pIntCrvUVTan = NULL;
  SmTsectCurveType eCurveType;
  double           dDeviation = 0.0;

  // For syntactic convenience: rename this corner-info object:
  SmFilletCornerInfo & rCorn = m_aCornerInfo[lWhichEnd][lWhichSide];

  SmEdge *pSideEdge = rCorn.m_pSideEdge;

  if ( pSideEdge != NULL )
  {
      double dCrvGuessParam =
          rCorn.m_bHighEndOfSideEdge ? pSideEdge->GetInterval().GetMax()
                                     : pSideEdge->GetInterval().GetMin();

      SmCurve *pCurve = pSideEdge->GetCurve();
      eStat = pCurve->DropPoint(pCurve->GetNaturalInterval(), // in : target curve allowed domain
                                m_aNewVtxPts[lWhichEnd],      // in : Point to drop to curve
                                NULL,                         // in : For Snapping, NULL=NoSnap, Vec pointing to CrvInside from desired EndPt for snapping.
                                                              //      when DropPt is within dDistTol in cpVecToCrvInside dir of EndPt, Snap to Ivl EndPt.
                                                              //      Vec handles ambiguities and makes sure that curves are not just touching at the ends.
                                dTol*10,                      // in : Used for SnapDist, and SM_SO_INTERSECT sXSectTol3d distance.
                                                              //      sXSectTol3d = dDistTol ;      SM_SO_INTERSECT 
                                                              //      sXSectTol3d = SM_BIG_DOUBLE ; SM_SO_MINIMIZE, SM_SO_NORMALIZE
                                                              //      sSnapTol3d  = dDistTol ;      SM_SO_MINIMIZE, SM_SO_NORMALIZE, SM_SO_INTERSECT    
                                &dCrvGuessParam,              // in : If given, call cheap LocalPointSolve() else call expensive GlobalPointSolve()
                                bSuccess,                     // out: TRUE = found a drop point
                                dCrvParam,                    // out: found drop curve param
                                dCrvDist) ;                   // out: found drop distance
                                                              // in : SM_SO_MINIMIZE = allow nonNormal drops near endPoints
                                                              //      SM_SO_INTERSECT= like Minimize but point must be within dDistTol
                                                              //      SM_SO_NORMALIZE= exclude nonNormal drops near endPoints
                                                              //      default:[SM_SO_MINIMIZE] to preserve original behavior
  }


  if ( eStat == SM_SUCCESS && bSuccess && dCrvDist <= dTol )
  {
      // Good, the curve is long enough, just update the Edge domain.
      rCorn.m_dNewSideEdgeParam = dCrvParam;

      // Just to make sure:
      rCorn.m_pNewSideEdgeExtension3d  = NULL;

      SmBSplineCurve *pUVCurve         = NULL;
      SmBSplineCurve *pUVCurveExtended = NULL;
      SmEdgeuse *pEdgeuse;
      SmExtent1d sSideEdgeExtent = pSideEdge->GetInterval();
      sSideEdgeExtent.AddValue( dCrvParam );

      const SmContext *pContext = pSideEdge->GetContext();

      // Do End-face uv curve first, then Tan-face uv curve.
      pEdgeuse = pSideEdge->GetEdgeuseOfFace( m_aEndFaces[ lWhichEnd ] );
      if(pEdgeuse) pUVCurve = pEdgeuse->GetUVTrimCurve();

      if ( pUVCurve != NULL )
      {
          SmExtent1d sUVExtent = pUVCurve->GetNaturalInterval();
          if ( ! sSideEdgeExtent.IsContainedBy( sUVExtent, SM_EFF_ZERO ) )
          {
              pUVCurve->CreateExtendedCurve( *pContext,
                      rCorn.m_dNewSideEdgeParam, SM_CT_CINFINITY,
                      pUVCurveExtended, TRUE );

              // Check whether it's good enough.
              if ( pUVCurveExtended != NULL )
              {
                  SmPoint3d sUVEval;
                  pUVCurveExtended->EvaluatePoint( dCrvParam, sUVEval );
                  if ( sUVEval.CloserThan( pSideEdge->GetTolerance(),
                          m_aUVEndSurfNew[ lWhichEnd ] ) )
                  {
                      // Good enough.  Swap out the uv curve for the extended one.
                      // This call deletes the old one, if 3rd arg is TRUE.
                      pEdgeuse->SetUVTrimCurve( pUVCurveExtended, 0, TRUE) ; // TRUE = delete existing UVTrimCurve
                  }
                  else
                  {
                      // At this point, just clear out the uv curve,
                      // it will be re-created sooner or later.
                      pEdgeuse->SetUVTrimCurve( NULL, 0, TRUE) ; // TRUE = delete existing UVTrimCurve
                      delete pUVCurveExtended; pUVCurveExtended = NULL;
                  }
              } // end if we created an extended curve for the uv trim curve.
          } // end if the side edge's uv trim curve was too short.
      } // end if the side edge had a uv trim curve in the End face.

      // Now do the uv-curve check again, for the Tan-face uv curve.
      pEdgeuse = pSideEdge->GetEdgeuseOfFace( m_aTanFaces[ lWhichSide ] );
      if(pEdgeuse) pUVCurve = pEdgeuse->GetUVTrimCurve();

      if ( pUVCurve != NULL )
      {
          SmExtent1d sUVExtent = pUVCurve->GetNaturalInterval();
          if ( ! sSideEdgeExtent.IsContainedBy( sUVExtent, SM_EFF_ZERO ) )
          {
              pUVCurve->CreateExtendedCurve( *pContext,
                      rCorn.m_dNewSideEdgeParam, SM_CT_CINFINITY,
                      pUVCurveExtended, TRUE );

              // Check whether it's good enough.
              if ( pUVCurveExtended != NULL )
              {
                  SmPoint3d sUVEval;
                  pUVCurveExtended->EvaluatePoint( dCrvParam, sUVEval );
                  if ( sUVEval.CloserThan( pSideEdge->GetTolerance(),
                          rCorn.m_sUVTanSurfNew ) )
                  {
                      // Good enough.  Swap out the uv curve for the extended one.
                      // This call deletes the old one, if 3rd arg is TRUE.
                      pEdgeuse->SetUVTrimCurve( pUVCurveExtended, 0, TRUE) ; // TRUE = delete existing UVTrimCurve
                  }
                  else
                  {
                      // At this point, just clear out the uv curve,
                      // it will be re-created sooner or later.
                      pEdgeuse->SetUVTrimCurve( NULL, 0, TRUE) ; // TRUE = delete existing UVTrimCurve
                      delete pUVCurveExtended; pUVCurveExtended = NULL;
                  }
              } // end if we created an extended curve for the uv trim curve.
          } // end if the side edge's uv trim curve was too short.
      } // end if the side edge had a uv trim curve in the End face.
  } // end if the existing side edge was long enough to reach the new int point.

  else
  {
      // We have to intersect the end part and append to the edge curve.

      SmTArray< SmPoint2d > sSurface0UVPoints;
      SmTArray< SmPoint2d > sSurface1UVPoints;

      sSurface0UVPoints.Add( rCorn.m_sUVEndSurfOld );
      sSurface0UVPoints.Add( m_aUVEndSurfNew[ lWhichEnd ] );

      sSurface1UVPoints.Add( rCorn.m_sUVTanSurfOld );
      sSurface1UVPoints.Add( rCorn.m_sUVTanSurfNew );

      // Definitely want to give it a start direction.
      // Just start-to-end should be fine, it won't curve much
      // over the radius of a fillet.
      SmVector3d sStartDir(
          m_aNewVtxPts[ lWhichEnd ] - rCorn.m_pVertex->GetPoint() );

      SmSurface *pEndSurf = m_aEndFaces[ lWhichEnd  ]->GetSurface();  NER(pEndSurf);
      SmSurface *pTanSurf = m_aTanFaces[ lWhichSide ]->GetSurface();  NER(pTanSurf);
      SmExtent2d sUVDomainEnd( pEndSurf->GetNaturalUVDomain() );
      SmExtent2d sUVDomainTan( pTanSurf->GetNaturalUVDomain() );
      // Temporarily enable OutOfBounds evaluations for all SmBSplineSurfaces.
      SmBoolean bTmp0 = FALSE, bTmp1 = FALSE ; 
      SmBSplineSurface * pBSplineSurface0 = SM_CAST_PTR(SmBSplineSurface, pEndSurf) ;
      SmBSplineSurface * pBSplineSurface1 = SM_CAST_PTR(SmBSplineSurface, pTanSurf) ;
      SmTemporaryChangeValue<SmBoolean> sClean0(pBSplineSurface0 ? pBSplineSurface0->GetOutOfBoundsEnabled() : bTmp0, TRUE) ;
      SmTemporaryChangeValue<SmBoolean> sClean1(pBSplineSurface1 ? pBSplineSurface1->GetOutOfBoundsEnabled() : bTmp1, TRUE) ;

      // Also expand the domains, for srf/srf/srf intersect.
      sUVDomainEnd.ExpandRelative( 2.0 );
      sUVDomainTan.ExpandRelative( 2.0 );

      const SmContext *pContext = pEndSurf->GetContext();

      eStat = pEndSurf->PointBasedSurfaceIntersect( *pContext,
                                                     sUVDomainEnd, 
                                                     *pTanSurf, 
                                                     sUVDomainTan,
                                                     sSurface0UVPoints, 
                                                     sSurface1UVPoints,
                                                     FALSE, 
                                                     FALSE, 
                                                     &sStartDir, 
                                                     NULL,
                                                     SM_CAST_APPROXTOL3D_PTR(&dTol), 
                                                     NULL,
                                                     pIntCrv3d, 
                                                     pIntCrvUVEnd, 
                                                     pIntCrvUVTan, 
                                                     eCurveType, 
                                                     dDeviation
                                                 );

      if ( eStat != SM_SUCCESS || pIntCrv3d == NULL )
        { SER( SM_ERR ); }

      // Save the curves in the kit.
      rCorn.m_pNewSideEdgeExtension3d  = pIntCrv3d;
      rCorn.m_pNewSideEdgeExtEndFaceUV = pIntCrvUVEnd;
      rCorn.m_pNewSideEdgeExtTanFaceUV = pIntCrvUVTan;

  } // end else: had to intersect end portion of side edge curve.


  return SM_SUCCESS;

} // end SmFilletRemoval::FindSideEdgeExtensions


/*******************************************************************//**
PURPOSE: Find what will be the new Edge curve.

NOTES:
   Intersect the surfaces of the two tangent Faces.

   Input: The results of CollectFilletTopology();
                         CalcNewCornerVertices();
                     and FindSideEdgeExtensions().

   Output: fills in:
   - The new intersection curve (3d and two uv curves) and its deviation:
     m_pNewEdgeCurve3d
     m_aNewEdgeCurveUVs[2]
     m_dNewEdgeDeviation

***********************************************************************/
SmStatus SmFilletRemoval::IntersectTangentSurfaces()
{
  const SmContext *pContext = m_aTanEdges[0]->GetContext();

  // Collect the uv points for point-based surface intersect.
  SmTArray<SmPoint2d> sSurface0UVPoints(2);
  SmTArray<SmPoint2d> sSurface1UVPoints(2);
  sSurface0UVPoints.Add( m_aCornerInfo[0][0].m_sUVTanSurfNew );
  sSurface0UVPoints.Add( m_aCornerInfo[1][0].m_sUVTanSurfNew );
  sSurface1UVPoints.Add( m_aCornerInfo[0][1].m_sUVTanSurfNew );
  sSurface1UVPoints.Add( m_aCornerInfo[1][1].m_sUVTanSurfNew );

  SmBSplineCurve  *pIntCrv3d = NULL, *pIntCrvUV0 = NULL, *pIntCrvUV1 = NULL;
  SmTsectCurveType eCurveType;
  double           dDeviation = 0.0;

  // Get a start direction.
  SmVector3d sStartDir;
  SmPoint3d  sPt;
  SmVector3d sSu, sSv;
  SmPoint2d  sUV = m_aCornerInfo[0][0].m_sUVFilSurf;

  SER( m_pFaceToDelete->GetSurface()->Evaluate1stDerivatives( sUV, TRUE, TRUE,
      sPt, sSu, sSv ));

  if ( m_eRailParam == SM_SP_U )
  {
      sStartDir = m_bEndsSwapped ? -sSu : sSu;
  }
  else
  {
      sStartDir = m_bEndsSwapped ? -sSv : sSv;
  }

  SmSurface *pTanSurf0 = m_aTanFaces[0]->GetSurface();
  SmSurface *pTanSurf1 = m_aTanFaces[1]->GetSurface();
  SmExtent2d sUVDomain1( pTanSurf0->GetNaturalUVDomain() );
  SmExtent2d sUVDomain2( pTanSurf1->GetNaturalUVDomain() );

  // Temporarily enable OutOfBounds evaluations for all SmBSplineSurfaces.
  SmBoolean bTmp0 = FALSE, bTmp1 = FALSE ; 
  SmBSplineSurface * pBSplineSurface0 = SM_CAST_PTR(SmBSplineSurface, pTanSurf0) ;
  SmBSplineSurface * pBSplineSurface1 = SM_CAST_PTR(SmBSplineSurface, pTanSurf1) ;
  SmTemporaryChangeValue<SmBoolean> sClean0(pBSplineSurface0 ? pBSplineSurface0->GetOutOfBoundsEnabled() : bTmp0, TRUE) ;
  SmTemporaryChangeValue<SmBoolean> sClean1(pBSplineSurface1 ? pBSplineSurface1->GetOutOfBoundsEnabled() : bTmp1, TRUE) ;

  // Also expand the domains, for srf/srf/srf intersect.
  sUVDomain1.ExpandRelative( 2.0 );
  sUVDomain2.ExpandRelative( 2.0 );

  SmStatus eStat = pTanSurf0->PointBasedSurfaceIntersect( *pContext,
                                                           sUVDomain1, 
                                                           *pTanSurf1, 
                                                           sUVDomain2,
                                                           sSurface0UVPoints, 
                                                           sSurface1UVPoints,
                                                           FALSE, 
                                                           FALSE, 
                                                           &sStartDir, 
                                                           NULL,
                                                           SM_CAST_APPROXTOL3D_PTR(&m_d3dTol), 
                                                           NULL,
                                                           pIntCrv3d, 
                                                           pIntCrvUV0, 
                                                           pIntCrvUV1, 
                                                           eCurveType, 
                                                           dDeviation);

  if ( eStat != SM_SUCCESS || pIntCrv3d == NULL )
    { SER( SM_ERR ); }

  m_pNewEdgeCurve3d     = pIntCrv3d;
  m_aNewEdgeCurveUVs[0] = pIntCrvUV0;
  m_aNewEdgeCurveUVs[1] = pIntCrvUV1;
  m_dNewEdgeDeviation   = dDeviation;

  return SM_SUCCESS;

} // end SmFilletRemoval::IntersectTangentSurfaces

/*******************************************************************//**
PURPOSE: Move a Vertex-Edge combination from one Face to an adjacent Face.

NOTES:
 This subroutine just does the topology.

 Another way to look at it: The other Edge coming into pVtxToMove that is
 also in pFromFace -- move that Edge off of pVtxToMove, to the other
 end of pEdgeToMove.

 This is a rather specialized topology-changing routine,
 which is why it is implemented here instead of SmBrep or something.
***********************************************************************/
SmStatus SmFaceRemoval::MoveVtxAndEdgeToAdjacentFace( SmVertex *pVtxToMove,
            SmEdge   *pEdgeToMove,
            SmFace   *pFromFace
          )
{
  // Grab this before we change anything:
  SmVertex *pOppVtx = pEdgeToMove->GetOtherVertex( pVtxToMove );
  NER( pOppVtx );

  // We work just with the circular linked lists of Edgeuses, plus the Vertexuses.

  // We need to move two Edgeuses from FromFace to ToFace.
  // (The two on pEdgeToMove that are in pFromFace.)
  // Start with the one that is oriented in the direction ( pVtxToMove -> pEdgeToMove ),
  // i.e., whose Vertexuse corresponds to pVtxToMove.
  SmEdgeuse *pEuToMove = pEdgeToMove->GetEdgeuseOfFace( pFromFace );
  if(pEuToMove == NULL) return SM_ERR;  

  if (pEuToMove && pEuToMove->GetVertexuse()->GetVertex() != pVtxToMove )
    { pEuToMove = pEuToMove->GetMate(); }

  // Make sure that this EU is not the one in the list that the EU's Loopuse points to.
  if (pEuToMove && pEuToMove->GetLoopuse()->GetEUorVU() == pEuToMove )
  {
      // Just slide it around the loop by one.
      pEuToMove->GetLoopuse()->SetEdgeuse( pEuToMove->GetCCWEdgeuse() );
  }

  // Get the EUs in the adjacent face where pEuToMove should go.
  // 'Before' means CW.
  SmEdgeuse *pEuInsertBefore = pEuToMove->m_pCW->GetRadial();
  SmEdgeuse *pEuInsertAfter  = pEuInsertBefore->m_pCW;

  // Grab the vertexuse to move, before we change anything.
  SmVertexuse *pVuToMove = pEuInsertBefore->GetVertexuse();

  // Unhitch pEuToMove from its doubly-linked EU list.
  pEuToMove->m_pCCW->m_pCW = pEuToMove->m_pCW;
  pEuToMove->m_pCW->m_pCCW = pEuToMove->m_pCCW;

  // Insert between Before and After.
  pEuInsertAfter->m_pCCW = pEuToMove;
  pEuInsertBefore->m_pCW = pEuToMove;
  pEuToMove->m_pCW  = pEuInsertAfter;
  pEuToMove->m_pCCW = pEuInsertBefore;

  // Change pEuToMove's Loopuse.
  pEuToMove->m_pSorLU = pEuInsertAfter->m_pSorLU;

  // And move the Vertexuse.

  // And move the Vertexuse.
  pVtxToMove->Remove ( pVuToMove );
  pOppVtx->PostInsert( pVuToMove );  // (pre/post? does it matter?)



  // Now repeat for the Mates.  CW and CCW get switched.
  pEuToMove = pEuToMove->GetMate();

  // Make sure that this EU is not the one in the list that the EU's Loopuse points to.
  if ( pEuToMove->GetLoopuse()->GetEUorVU() == pEuToMove )
  {
      // Just slide it around the loop by one.
      pEuToMove->GetLoopuse()->SetEdgeuse( pEuToMove->GetCCWEdgeuse() );
  }

  // Get the EUs in the adjacent face where pEuToMove should go.
  // 'Before' means CW.
  pEuInsertAfter  = pEuToMove->m_pCCW->GetRadial();
  pEuInsertBefore = pEuInsertAfter->m_pCCW;

  // Grab the vertexuse to move, before we change anything.
  pVuToMove = pEuToMove->m_pCCW->GetVertexuse();

  // Unhitch pEuToMove from its doubly-linked EU list.
  pEuToMove->m_pCCW->m_pCW = pEuToMove->m_pCW;
  pEuToMove->m_pCW->m_pCCW = pEuToMove->m_pCCW;

  // Insert between Before and After.
  pEuInsertAfter->m_pCCW = pEuToMove;
  pEuInsertBefore->m_pCW = pEuToMove;
  pEuToMove->m_pCW  = pEuInsertAfter;
  pEuToMove->m_pCCW = pEuInsertBefore;

  // Change pEuToMove's Loopuse.
  pEuToMove->m_pSorLU = pEuInsertAfter->m_pSorLU;

  // And move the Vertexuse.
  pVtxToMove->Remove ( pVuToMove );
  pOppVtx->PostInsert( pVuToMove );  // (pre/post? does it matter?)

  return SM_SUCCESS;

} // end MoveVtxAndEdgeToAdjacentFace


/*******************************************************************//**
PURPOSE: Do the topology editing: remove the Fillet Face and patch the hole.

NOTES:
   Input: This object must be completely filled:
          The results of CollectFilletTopology();
                         CalcNewCornerVertices();
                         FindSideEdgeExtensions();
                     and IntersectTangentSurfaces().

   Output: the Brep is modified.

***********************************************************************/
SmStatus SmFilletRemoval::DoRemoval()
{
  SmBrep * pBrep = m_pFaceToDelete->GetBrep();    NER( pBrep );
  pBrep->Notify(SM_NO_PRE_EDIT, pBrep, NULL, NULL);

  // Before we can delete the Fillet Face, we have to check a complication.
  // If any corner had more than three Edges incident on the vertex,
  // then that End Edge will not be deleted,
  // it must be moved from the Fillet Face to the Tangent Face on the
  // side of the corner in question.
  // This must be done before deleting the Fillet Face, because that would
  // lose the two Edgeuses of the End Edge in the Fillet Face;
  // instead we'll use those Edgeuses to move the End Edge out of the
  // Fillet Face and into the Tangent Face.

  if ( m_aCornerInfo[0][0].NoSideEdge() || m_aCornerInfo[0][1].NoSideEdge() )
  {
      // Here we keep the edge, but it needs to be changed.
      // It was the edge between the Fillet Face and the End Face,
      // and it will now be between the Tangent Face and the End Face.
      // The geometry for that intersection is stored in the Corner's
      // three EdgeExtention curves (one 3d, two uv curves).

      ULONG lWhichSide = ( m_aCornerInfo[0][0].NoSideEdge() ) ? 0 : 1;

      SmFilletCornerInfo & rCorn = m_aCornerInfo[0][lWhichSide];

      MoveVtxAndEdgeToAdjacentFace(
              rCorn.m_pVertex,  // Vertex to move
              m_aEndEdges[0],   // Edge to move
              m_pFaceToDelete   // From Face
      );

      // Now change its geometry.
      SmCurve *pOldCrv = m_aEndEdges[0]->GetCurve();       NER( pOldCrv );
      SmCurve *pNewCrv = rCorn.m_pNewSideEdgeExtension3d;  NER( pNewCrv );

      // Check whether the new geometry goes in the right direction.
      SmVector3d sPtDer0[2], sPtDer1[2];
      pOldCrv->Evaluate( m_aEndEdges[0]->GetInterval().Evaluate(0.5), 1, FALSE, sPtDer0 );
      pNewCrv->Evaluate( pNewCrv->GetNaturalInterval().Evaluate(0.5), 1, FALSE, sPtDer1 );
      if ( sPtDer0[1].Dot( sPtDer1[1] ) < 0.0 )
      {
          SmExtent1d sDummy(0,1);
          pNewCrv->ReverseParameterization( sDummy, sDummy );
          rCorn.m_pNewSideEdgeExtEndFaceUV->ReverseParameterization( sDummy, sDummy );
          rCorn.m_pNewSideEdgeExtTanFaceUV->ReverseParameterization( sDummy, sDummy );
      }

      // Set m_aEndEdges[0]->Curve
      m_aEndEdges[0]->Replace3DCurve( rCorn.m_pNewSideEdgeExtension3d ); // note: works for Edge and CEdge, deletes Edge->Curve and UVTrimCurves

      // Set Edgeuse->UVTrimCurve = pUVTrimCurve for Edgeuse connecting thisEdge to TargetFace
      m_aEndEdges[0]->SetTrimCurve( m_aEndFaces[0],  
                                    rCorn.m_pNewSideEdgeExtEndFaceUV,
                                    m_dNewEdgeDeviation, 
                                    FALSE );    // FALSE = do not validate geometry here, they are in an inconsistent state.
      m_aEndEdges[0]->SetTrimCurve( m_aTanFaces[lWhichSide], 
                                    rCorn.m_pNewSideEdgeExtTanFaceUV,
                                    m_dNewEdgeDeviation, 
                                    FALSE );    // FALSE = do not validate geometry here, they are in an inconsistent state.

  } // end no-side-edge check, end [0]

  if ( m_aCornerInfo[1][0].NoSideEdge() || m_aCornerInfo[1][1].NoSideEdge() )
  {
      // Same for other end.

      ULONG lWhichSide = ( m_aCornerInfo[1][0].NoSideEdge() ) ? 0 : 1;

      SmFilletCornerInfo & rCorn = m_aCornerInfo[1][lWhichSide];

      MoveVtxAndEdgeToAdjacentFace(
              rCorn.m_pVertex,  // Vertex to move
              m_aEndEdges[1],   // Edge to move
              m_pFaceToDelete   // From Face
      );

      // Now change its geometry.
      SmCurve *pOldCrv = m_aEndEdges[1]->GetCurve();       NER( pOldCrv );
      SmCurve *pNewCrv = rCorn.m_pNewSideEdgeExtension3d;  NER( pNewCrv );

      // Check whether the new geometry goes in the right direction.
      SmVector3d sPtDer0[2], sPtDer1[2];
      pOldCrv->Evaluate( m_aEndEdges[1]->GetInterval().Evaluate(0.5), 1, FALSE, sPtDer0 );
      pNewCrv->Evaluate( pNewCrv->GetNaturalInterval().Evaluate(0.5), 1, FALSE, sPtDer1 );
      if ( sPtDer0[1].Dot( sPtDer1[1] ) < 0.0 )
      {
          SmExtent1d sDummy(0,1);
          pNewCrv->ReverseParameterization( sDummy, sDummy );
          rCorn.m_pNewSideEdgeExtEndFaceUV->ReverseParameterization( sDummy, sDummy );
          rCorn.m_pNewSideEdgeExtTanFaceUV->ReverseParameterization( sDummy, sDummy );
      }

      // Set m_aEndEdges[1]->Curve
      m_aEndEdges[1]->Replace3DCurve( rCorn.m_pNewSideEdgeExtension3d ); // note: works for Edge and CEdge, deletes Edge->Curve and UVTrimCurves

      // Set Edgeuse->UVTrimCurve = pUVTrimCurve for Edgeuse connecting thisEdge to TargetFace
      m_aEndEdges[1]->SetTrimCurve( m_aEndFaces[1],  
                                    rCorn.m_pNewSideEdgeExtEndFaceUV,
                                    m_dNewEdgeDeviation, 
                                    FALSE );    // FALSE = do not validate geometry here, they are in an inconsistent state.
      m_aEndEdges[1]->SetTrimCurve( m_aTanFaces[lWhichSide], 
                                    rCorn.m_pNewSideEdgeExtTanFaceUV,
                                    m_dNewEdgeDeviation, 
                                    FALSE );    // FALSE = do not validate geometry here, they are in an inconsistent state.

  } // end no-side-edge check, end [1]

  // Now we can delete the fillet Face.
  // Two FALSE flags:
  //  don't delete Edges, etc.;
  //  and don't combine shells/regions.
  pBrep->DeleteFace( m_pFaceToDelete, FALSE, FALSE );


  // Loop over both End faces and Side (tangent) faces,
  // extending them and their surfaces to include any new uv points.
  ULONG i;
  SmFace           * pFace = NULL;
  SmBSplineSurface * pSrf  = NULL;
  SmExtent2d         sFaceDomain;
  SmStatus           eStat;

  for ( i = 0; i < 2; i++ )
  {
      // Update the end faces.
      pFace = m_aEndFaces[ i ];
      sFaceDomain = pFace->GetUVDomain();
      sFaceDomain.AddPoint2d( m_aUVEndSurfNew[ i ] );
      pFace->SetUVDomain( sFaceDomain );

      pSrf = SM_CAST_PTR( SmBSplineSurface, pFace->GetSurface() );
      if ( pSrf != NULL )
      {
          eStat = pSrf->ExpandBoundarySpan( sFaceDomain.GetMin().x, SM_SP_U );
          eStat = pSrf->ExpandBoundarySpan( sFaceDomain.GetMax().x, SM_SP_U );
          eStat = pSrf->ExpandBoundarySpan( sFaceDomain.GetMin().y, SM_SP_V );
          eStat = pSrf->ExpandBoundarySpan( sFaceDomain.GetMax().y, SM_SP_V );
          eStat = pSrf->UpdateAnalyticalDomain( sFaceDomain );
      }

      // Now do the side (tangent) faces.  These have two new uv points.
      pFace = m_aTanFaces[ i ];
      sFaceDomain = pFace->GetUVDomain();
      sFaceDomain.AddPoint2d( m_aCornerInfo[0][i].m_sUVTanSurfNew );
      sFaceDomain.AddPoint2d( m_aCornerInfo[1][i].m_sUVTanSurfNew );
      pFace->SetUVDomain( sFaceDomain );

      // Now update its surface.
      pSrf = SM_CAST_PTR( SmBSplineSurface, pFace->GetSurface() );
      if ( pSrf != NULL )
      {
          eStat = pSrf->ExpandBoundarySpan( sFaceDomain.GetMin().x, SM_SP_U );
          eStat = pSrf->ExpandBoundarySpan( sFaceDomain.GetMax().x, SM_SP_U );
          eStat = pSrf->ExpandBoundarySpan( sFaceDomain.GetMin().y, SM_SP_V );
          eStat = pSrf->ExpandBoundarySpan( sFaceDomain.GetMax().y, SM_SP_V );
          eStat = pSrf->UpdateAnalyticalDomain( sFaceDomain );
      }
  } // end loop on two End faces and two Tangent faces.


  // Now, for both ends:
  //   Set the new vertex positions,
  //   and extend the two side edges.

  ULONG lWhichEnd, lWhichSide;

  for ( lWhichEnd = 0; lWhichEnd < 2; lWhichEnd++ )
  {
      // Set the vertex positions,
      // and adjust their tolerances (only if bigger: TRUE arg).
      SmPoint3d sNewVtxPt = m_aNewVtxPts [ lWhichEnd ];
      double    dNewTol   = m_aNewVtxTols[ lWhichEnd ];

      // But leave the corner alone if we didn't record a side edge
      // (more than three vertices at that vertex).
      if ( !m_aCornerInfo[lWhichEnd][0].NoSideEdge() )
      {
          m_aCornerInfo[lWhichEnd][0].m_pVertex->SetPoint( sNewVtxPt );
#ifdef SM_USE_NEWTOL
          SM_NEWTOL_LINE m_aCornerInfo[lWhichEnd][0].m_pVertex->ClearLocalZoneTol3d( );
#else // SM_USE_OLDTOL
          SM_OLDTOL_LINE m_aCornerInfo[lWhichEnd][0].m_pVertex->SetTolerance( dNewTol, TRUE );
#endif // SM_USE_OLDTOL
      }

      if ( !m_aCornerInfo[lWhichEnd][1].NoSideEdge() )
      {
          m_aCornerInfo[lWhichEnd][1].m_pVertex->SetPoint( sNewVtxPt );
#ifdef SM_USE_NEWTOL
          SM_NEWTOL_LINE m_aCornerInfo[lWhichEnd][1].m_pVertex->ClearLocalZoneTol3d( );
#else // SM_USE_OLDTOL
          SM_OLDTOL_LINE m_aCornerInfo[lWhichEnd][1].m_pVertex->SetTolerance( dNewTol, TRUE );
#endif // SM_USE_OLDTOL
      }


      // For the side edges: append the last bit of intersection if necessary.
      for ( lWhichSide = 0; lWhichSide < 2; lWhichSide++ )
      {
          SmFilletCornerInfo &rCorn = m_aCornerInfo[lWhichEnd][lWhichSide];

          if ( rCorn.NoSideEdge() )
            { continue; } // more than three edges, don't extend anything.

          SmEdge *pSideEdge = rCorn.m_pSideEdge;  NER( pSideEdge );

          SmExtent1d sNewSideEdgeDomain = pSideEdge->GetInterval();

          if ( rCorn.m_pNewSideEdgeExtension3d == NULL )
          {
              // Curve was long enough: just update the Edge's domain.
              sNewSideEdgeDomain.AddValue( rCorn.m_dNewSideEdgeParam );

              pSideEdge->SetInterval( sNewSideEdgeDomain );

              // And the two uv curves: in End face, and Tangent face.
              SmSurface *pSurf = m_aEndFaces[ lWhichEnd ]->GetSurface();
              SmBSplineCurve *pUVCurve = pSideEdge->GetUVTrimCurveOfSurface( pSurf );
              if ( pUVCurve != NULL )
              {
                  //SmBoolean bBool = pUVCurve->IsKindOf( SmBSplineCurve_TYPE );
                  //SM_TYPE tType = pUVCurve->GetType();
              }
          }
          else
          {
              // Append this curve piece to the side edge's curve.
              // First, we had better make sure that the side edge curve
              // doesn't extend past the original vertex position.

              SmCurve *pCurve = pSideEdge->GetCurve();  NER( pCurve );

              // pCurve will have to be a SmBSplineCurve, for JoinWith().
              SmBSplineCurve *pBSplCurve = SM_CAST_PTR( SmBSplineCurve, pCurve );
              NER( pBSplCurve );

              SmExtent1d sCurveDomain = pBSplCurve->GetNaturalInterval();
              if ( rCorn.m_bHighEndOfSideEdge )
              {
                  sCurveDomain.SetMinMax( sCurveDomain.GetMin(),
                                          sNewSideEdgeDomain.GetMax() );
              }
              else
              {
                  sCurveDomain.SetMinMax( sNewSideEdgeDomain.GetMin(),
                                          sCurveDomain.GetMax() );
              }
              pBSplCurve->Trim( sCurveDomain ); // may snap sIvl by tol to existing knots

              // Now append the new intersection curve.
              SER( pBSplCurve->JoinWith( rCorn.m_bHighEndOfSideEdge ? 1 : 0,
                                    rCorn.m_pNewSideEdgeExtension3d, 0 ));

              // Done with this:
              if ( rCorn.m_pNewSideEdgeExtension3d != NULL )
              {
                  delete rCorn.m_pNewSideEdgeExtension3d;
                  rCorn.m_pNewSideEdgeExtension3d = NULL;
              }

              // Find the SideEdge's new domain, and update the Edge.
              // Refresh this, after trim and join.
              // Notes: JoinWith() can completely change the parameterization,
              // and, the far end of the Edge (the vertex) is not necessarily
              // at the end of the curve.
              // So: Update the Edge's domain:
              //   This end: same as this end of curve.
              //   Far end: have to drop the far vertex to the curve.

              // Some locals:
              double t0, t1;
              SmBoolean bSuccess;
              double dGuessParam, dDist;
              SmPoint3d sVtxPt;
              sCurveDomain = pBSplCurve->GetNaturalInterval();
              if ( rCorn.m_bHighEndOfSideEdge )
              {
                  t1 = sCurveDomain.GetMax();

                  // Set t0 at the far end.
                  sVtxPt = pSideEdge->GetStartVertex()->GetPoint();
                  dGuessParam = sNewSideEdgeDomain.GetMin();
                  dGuessParam = sCurveDomain.ClampValue( dGuessParam );
                  eStat = pBSplCurve->DropPoint(sCurveDomain,  // in : target curve allowed domain
                                                sVtxPt,        // in : Point to drop to curve
                                                NULL,          // in : For Snapping, NULL=NoSnap, Vec pointing to CrvInside from desired EndPt for snapping.
                                                               //      when DropPt is within dDistTol in cpVecToCrvInside dir of EndPt, Snap to Ivl EndPt.
                                                               //      Vec handles ambiguities and makes sure that curves are not just touching at the ends.
                                                m_d3dTol,      // in : Used for SnapDist, and SM_SO_INTERSECT sXSectTol3d distance.
                                                               //      sXSectTol3d = dDistTol ;      SM_SO_INTERSECT 
                                                               //      sXSectTol3d = SM_BIG_DOUBLE ; SM_SO_MINIMIZE, SM_SO_NORMALIZE
                                                               //      sSnapTol3d  = dDistTol ;      SM_SO_MINIMIZE, SM_SO_NORMALIZE, SM_SO_INTERSECT    
                                                &dGuessParam,  // in : If given, call cheap LocalPointSolve() else call expensive GlobalPointSolve()
                                                bSuccess,      // out: TRUE = found a drop point
                                                t0,            // out: found drop curve param
                                                dDist) ;       // out: found drop distance
                                                               // in : SM_SO_MINIMIZE = allow nonNormal drops near endPoints
                                                               //      SM_SO_INTERSECT= like Minimize but point must be within dDistTol
                                                               //      SM_SO_NORMALIZE= exclude nonNormal drops near endPoints
                                                               //      default:[SM_SO_MINIMIZE] to preserve original behavior
                  if ( eStat != SM_SUCCESS || !bSuccess || dDist > m_d3dTol )
                    { SER( SM_ERR ); }
              }
              else // we're at the low end of Edge
              {
                  t0 = sCurveDomain.GetMin();

                  // Set t1 at the far end.
                  SmVertex *pNearVtx = pSideEdge->GetStartVertex();
                  sVtxPt = pSideEdge->GetOtherVertex( pNearVtx )->GetPoint();
                  dGuessParam = sNewSideEdgeDomain.GetMax();
                  dGuessParam = sCurveDomain.ClampValue( dGuessParam );
                  eStat = pBSplCurve->DropPoint(sCurveDomain, // in : target curve allowed domain
                                                sVtxPt,       // in : Point to drop to curve
                                                NULL,         // in : For Snapping, NULL=NoSnap, Vec pointing to CrvInside from desired EndPt for snapping.
                                                              //      when DropPt is within dDistTol in cpVecToCrvInside dir of EndPt, Snap to Ivl EndPt.
                                                              //      Vec handles ambiguities and makes sure that curves are not just touching at the ends.
                                                m_d3dTol,     // in : Used for SnapDist, and SM_SO_INTERSECT sXSectTol3d distance.
                                                              //      sXSectTol3d = dDistTol ;      SM_SO_INTERSECT 
                                                              //      sXSectTol3d = SM_BIG_DOUBLE ; SM_SO_MINIMIZE, SM_SO_NORMALIZE
                                                              //      sSnapTol3d  = dDistTol ;      SM_SO_MINIMIZE, SM_SO_NORMALIZE, SM_SO_INTERSECT    
                                                &dGuessParam, // in : If given, call cheap LocalPointSolve() else call expensive GlobalPointSolve()
                                                bSuccess,     // out: TRUE = found a drop point
                                                t1,           // out: found drop curve param
                                                dDist) ;      // out: found drop distance
                                                              // in : SM_SO_MINIMIZE = allow nonNormal drops near endPoints
                                                              //      SM_SO_INTERSECT= like Minimize but point must be within dDistTol
                                                              //      SM_SO_NORMALIZE= exclude nonNormal drops near endPoints
                                                              //      default:[SM_SO_MINIMIZE] to preserve original behavior
                  if ( eStat != SM_SUCCESS || !bSuccess || dDist > m_d3dTol )
                    { SER( SM_ERR ); }

              } // end if (high or low end of side edge)

              sNewSideEdgeDomain.SetMinMax( t0, t1 );
              pSideEdge->SetInterval( sNewSideEdgeDomain );

              // ... and then we have to do the uv curves also.
              // These domains should be the same as the curve.
              if ( rCorn.m_pNewSideEdgeExtEndFaceUV != NULL )
              {
                  // This gets the same-oriented EU in the Face:
                  SmEdgeuse *pSideEU = pSideEdge->GetEdgeuseOfFace( m_aEndFaces[ lWhichEnd ] );
                  if(pSideEU)
                  {
                    SmBSplineCurve *pSideUVCurve = pSideEU->GetUVTrimCurve();

                    if(pSideUVCurve != NULL)
                    {
                      SER( pSideUVCurve->JoinWith( rCorn.m_bHighEndOfSideEdge ? 1 : 0,
                           rCorn.m_pNewSideEdgeExtEndFaceUV, 0 ) );
                      pSideUVCurve->EditParameterization( sCurveDomain );

                      // Done with this:
                      if(rCorn.m_pNewSideEdgeExtEndFaceUV != NULL)
                      {
                        delete rCorn.m_pNewSideEdgeExtEndFaceUV;
                        rCorn.m_pNewSideEdgeExtEndFaceUV = NULL;
                      }
                    }
                  } // end if pSideEU
              } // end if end-face extension uv curve not null, so append it.

              if ( rCorn.m_pNewSideEdgeExtTanFaceUV != NULL )
              {
                  // This gets the same-oriented EU in the Face:
                  SmEdgeuse *pSideEU = pSideEdge->GetEdgeuseOfFace( m_aTanFaces[lWhichSide] );
                  if(pSideEU)
                  {
                    SmBSplineCurve *pSideUVCurve = pSideEU->GetUVTrimCurve();
                    if(pSideUVCurve != NULL)
                    {
                      SER( pSideUVCurve->JoinWith( rCorn.m_bHighEndOfSideEdge ? 1 : 0,
                           rCorn.m_pNewSideEdgeExtTanFaceUV, 0 ) );
                      pSideUVCurve->EditParameterization( sCurveDomain );

                      // Done with this:
                      if(rCorn.m_pNewSideEdgeExtTanFaceUV != NULL)
                      {
                        delete rCorn.m_pNewSideEdgeExtTanFaceUV;
                        rCorn.m_pNewSideEdgeExtTanFaceUV = NULL;
                      }
                    }
                  } // end if pSideEU
              } // end if tangent-face extension uv curve not null, so append it.
          } // end if side-edge extension 3d curve not null, so append it.
      } // end loop on two Sides at this end, updating side edges.
  } // end loop on both Ends.


  // Install the intersection curve of the two tangent faces
  // into the two tangent edges -- they will be glued.
  // First check the relative orientation of the two tangent edges: needed for Glue.
  SmCurve  * pTanCurve0 = m_aTanEdges[0]->GetCurve();
  SmCurve  * pTanCurve1 = m_aTanEdges[1]->GetCurve();
  SmExtent1d sDomain0   = m_aTanEdges[0]->GetInterval();
  SmExtent1d sDomain1   = m_aTanEdges[1]->GetInterval();
  SmVector3d aPV0[2];
  SmVector3d aPV1[2];

  pTanCurve0->Evaluate( sDomain0.Evaluate( 0.5 ), 1, TRUE, aPV0 );
  pTanCurve1->Evaluate( sDomain1.Evaluate( 0.5 ), 1, TRUE, aPV1 );

  SmOrientType eOrient =  (aPV0[1].Dot(aPV1[1]) > 0.0)
                        ? SM_OT_SAME
                        : SM_OT_OPPOSITE;

  // Note, don't install the new curve in both edges.
  // One edge will be deleted, which deletes its curve.
  // Set the second Edge's curve to NULL, and make sure
  // that we pass it second into GlueEdges().
  // And just to make sure of that, rename the Edges:
  SmEdge *pSurviveEdge = m_aTanEdges[0];  
  SmEdge *pDeleteEdge  = m_aTanEdges[1];  

  // We have to get the direction of the new Edge curve right:
  // it has to run in the same direction as the Edge.
  // We just evaluated the original curve (pTanCurve0), in aPV0.
  // Evaluate the new curve and check directions.

  sDomain0 = m_pNewEdgeCurve3d->GetNaturalInterval();
  m_pNewEdgeCurve3d->Evaluate( sDomain0.Evaluate(0.5), 1, TRUE, aPV1 );
  if ( aPV0[1].Dot( aPV1[1] ) < 0 )
  {
      m_pNewEdgeCurve3d->ReverseParameterization( sDomain0, sDomain1 );
      if ( m_aNewEdgeCurveUVs[0] != NULL )
      {
          sDomain0 = m_aNewEdgeCurveUVs[0]->GetNaturalInterval();
          m_aNewEdgeCurveUVs[0]->ReverseParameterization( sDomain0, sDomain1 );
      }
      if ( m_aNewEdgeCurveUVs[1] != NULL )
      {
          sDomain0 = m_aNewEdgeCurveUVs[1]->GetNaturalInterval();
          m_aNewEdgeCurveUVs[1]->ReverseParameterization( sDomain0, sDomain1 );
      }
  }

  pSurviveEdge->SetCurve( m_pNewEdgeCurve3d, TRUE ); // TRUE = delete preExisting Edge->Curve
                                                     // side effect: delete current pSurviveEdge->UVTrimCurve
  m_pNewEdgeCurve3d->SetOwner( pSurviveEdge );
  pDeleteEdge ->SetCurve(NULL, TRUE ) ; // TRUE = delete preExisting Edge->Curve
                                        // side effect: delete current pSurviveEdge->UVTrimCurve

  // Delete the edge curves that were there before.
  // delete pTanCurve0; pTanCurve0 = NULL; // already done
  // delete pTanCurve1; pTanCurve1 = NULL; // already done
  pTanCurve0 = NULL; // done with pTanCurve0, it's been deleted
  pTanCurve1 = NULL; // done with pTanCurve1, it's been deleted

  // and update the Edges' intervals.
  // Would be just the domain of the new int curve: it went vertex-to-vertex.
  sDomain0 = m_pNewEdgeCurve3d->GetNaturalInterval();
  pSurviveEdge->SetInterval( sDomain0 );
  pDeleteEdge ->SetInterval( sDomain0 );


  // Now we can delete the two end Edges, and glue the vertices,
  // unless either vertex at either end has more than three vertices
  // (in which case we didn't record a side edge).
  if ( ! ( m_aCornerInfo[0][0].NoSideEdge() || m_aCornerInfo[0][1].NoSideEdge() ) )
  {
      SER( pBrep->DeleteEdge( m_aEndEdges[0] ));

      SER( pBrep->GlueVertices( m_aCornerInfo[0][0].m_pVertex,
                                m_aCornerInfo[0][1].m_pVertex ));
  }
  if ( ! ( m_aCornerInfo[1][0].NoSideEdge() || m_aCornerInfo[1][1].NoSideEdge() ) )
  {
      SER( pBrep->DeleteEdge( m_aEndEdges[1] ));

      SER( pBrep->GlueVertices( m_aCornerInfo[1][0].m_pVertex,
                                m_aCornerInfo[1][1].m_pVertex ));
  }

  // Now glue the two tangent edges, which will become the new (non-tangent) Edge.
  SER( pBrep->GlueEdges( pSurviveEdge, eOrient, 0.0, pDeleteEdge ));

  // Set Edgeuse->UVTrimCurve = pUVTrimCurve for Edgeuse connecting thisEdge to TargetFace
  pSurviveEdge->SetTrimCurve( m_aTanFaces[0], 
                              m_aNewEdgeCurveUVs[0],
                              m_dNewEdgeDeviation, 
                              FALSE );     // FALSE = do not validate geometry here, they are in an inconsistent state.

  // Set Edgeuse->UVTrimCurve = pUVTrimCurve for Edgeuse connecting thisEdge to TargetFace
  pSurviveEdge->SetTrimCurve( m_aTanFaces[1], 
                              m_aNewEdgeCurveUVs[1],
                              m_dNewEdgeDeviation, 
                              TRUE );     // TRUE = validate geometry here. On the last one its possible

  // And that should do it.
  return SM_SUCCESS;

} // end SmFilletRemoval::DoRemoval

/*******************************************************************//**
PURPOSE: SolidCraft change: SmRingFilletRemoval constructor.

NOTES:
***********************************************************************/
SmRingFilletRemoval::SmRingFilletRemoval( SmFace *pFaceToDelete, double d3dTolerance )
{
  m_pFaceToDelete = pFaceToDelete;
  m_d3dTol        = d3dTolerance;
  m_pSeam         = NULL;
  m_pNewEdgeCurve = NULL;
  m_bNotExact     = FALSE;
  for ( ULONG ii = 0; ii < 2; ii++ )
  {
      m_aRails[ii]      = NULL;
      m_aTanFaces[ii]   = NULL;
      m_aVertices[ii]   = NULL;
      m_aSideEdges[ii]  = NULL;
      m_aSideAtHigh[ii] = FALSE;
      m_aSideParams[ii] = 0.0;
      m_aSideCurves[ii] = NULL;
      for ( ULONG jj = 0; jj < 2; jj++ )
      {
          m_aSideUVs[ii][jj]      = NULL;
          m_aSideUVOwners[ii][jj] = NULL;
      }
  }
}

/*******************************************************************//**
PURPOSE: SolidCraft change: SmRingFilletRemoval destructor.

NOTES: The new curves are owned here until DoRemoval() hands them to the
   Brep (and clears these pointers).
***********************************************************************/
SmRingFilletRemoval::~SmRingFilletRemoval()
{
  delete m_pNewEdgeCurve;
  m_pNewEdgeCurve = NULL;
  for ( ULONG ii = 0; ii < 2; ii++ )
  {
      delete m_aSideCurves[ii];
      m_aSideCurves[ii] = NULL;
      for ( ULONG jj = 0; jj < 2; jj++ )
      {
          delete m_aSideUVs[ii][jj];
          m_aSideUVs[ii][jj] = NULL;
      }
  }
}

/*******************************************************************//**
PURPOSE: SolidCraft change: TRUE when pFace is bounded as a ring fillet.

NOTES: Two closed tangent manifold edges, and at most one edge whose only
   face is pFace (its seam). Nothing about the surface is checked: that is
   SmFace::IdentifyFillet()'s job.
***********************************************************************/
SmBoolean SmRingFilletRemoval::IsRing( const SmFace *pFace, double dTangencyTolDeg )
{
  if ( pFace == NULL )
    { return FALSE; }
  SmTArray<SmEdge*> sEdges;
  pFace->GetEdges( sEdges );
  ULONG lRails = 0, lSeams = 0;
  for ( ULONG ii = 0; ii < sEdges.GetSize(); ii++ )
  {
      SmEdge *pE = sEdges[ii];
      SmTArray<SmFace*> sFaces;
      pE->GetFaces( sFaces );
      if ( sFaces.GetSize() == 1 && sFaces[0] == pFace )
        { lSeams++; continue; }
      if ( pE->IsClosed() && pE->IsManifold() && pE->IsTangentEdge( dTangencyTolDeg ) )
        { lRails++; continue; }
      return FALSE;
  }
  return lRails == 2 && lSeams <= 1;
}

/*******************************************************************//**
PURPOSE: SolidCraft change: remove the ring fillet and close the hole with
   the sharp edge.

NOTES: See the class notes.
***********************************************************************/
SmStatus SmRingFilletRemoval::RemoveFace()
{
  SER( this->CollectRingTopology() );
  SER( this->IntersectTangentSurfaces() );
  SER( this->FindSideEdgeExtension( 0 ) );
  SER( this->FindSideEdgeExtension( 1 ) );

  // Everything is known: change the Brep.
  SER( this->DoRemoval() );
  return SM_SUCCESS;

} // end SmRingFilletRemoval::RemoveFace

/*******************************************************************//**
PURPOSE: SolidCraft change: collect the rails, the tangent faces, the
   rails' vertices, the fillet's seam and the tangent faces' seams at the
   vertices.

NOTES:
***********************************************************************/
SmStatus SmRingFilletRemoval::CollectRingTopology()
{
  const double dTangencyTolDeg = 2.0; // as SmFilletRemoval::CollectFilletTopology
  SmFace *pFilFace = m_pFaceToDelete;
  NER( pFilFace );
  if ( !IsRing( pFilFace, dTangencyTolDeg ) )
    { SER( SM_ERR ); }

  SmTArray<SmEdge*> sEdges;
  pFilFace->GetEdges( sEdges );
  ULONG lRail = 0;
  for ( ULONG ii = 0; ii < sEdges.GetSize(); ii++ )
  {
      SmTArray<SmFace*> sFaces;
      sEdges[ii]->GetFaces( sFaces );
      if ( sFaces.GetSize() == 1 && sFaces[0] == pFilFace )
        { m_pSeam = sEdges[ii]; continue; }
      if ( sFaces.GetSize() != 2 || lRail > 1 )
        { SER( SM_ERR ); }
      m_aRails[lRail]    = sEdges[ii];
      m_aTanFaces[lRail] = ( sFaces[0] == pFilFace ) ? sFaces[1] : sFaces[0];
      lRail++;
  }
  if ( lRail != 2 || m_aTanFaces[0] == m_aTanFaces[1] || m_aTanFaces[0] == pFilFace || m_aTanFaces[1] == pFilFace )
    { SER( SM_ERR ); }
  for ( ULONG ii = 0; ii < 2; ii++ )
  {
      SmSurface *pSrf = m_aTanFaces[ii]->GetSurface();
      if ( pSrf == NULL || !pSrf->IsKindOf( SmBSplineSurface_TYPE ) )
        { SER( SM_ERR ); }
      m_aVertices[ii] = m_aRails[ii]->GetStartVertex();
      NER( m_aVertices[ii] );
  }

  // The fillet's seam joins the two rails' vertices.
  if ( m_pSeam != NULL )
  {
      SmVertex *pV0 = m_pSeam->GetStartVertex();
      SmVertex *pV1 = m_pSeam->GetOtherVertex( pV0 );
      if ( !( ( pV0 == m_aVertices[0] && pV1 == m_aVertices[1] ) || ( pV0 == m_aVertices[1] && pV1 == m_aVertices[0] ) ) )
        { SER( SM_ERR ); }
  }

  // At each rail's vertex, besides the rail and the fillet's seam: at most one seam of the
  // tangent face (an edge whose only face it is).
  for ( ULONG ii = 0; ii < 2; ii++ )
  {
      SmTArray<SmEdge*> sAtVertex;
      m_aVertices[ii]->GetEdges( sAtVertex );
      for ( ULONG jj = 0; jj < sAtVertex.GetSize(); jj++ )
      {
          SmEdge *pE = sAtVertex[jj];
          if ( pE == m_aRails[ii] || pE == m_pSeam )
            { continue; }
          SmTArray<SmFace*> sFaces;
          pE->GetFaces( sFaces );
          if ( m_aSideEdges[ii] != NULL || sFaces.GetSize() != 1 || sFaces[0] != m_aTanFaces[ii] || pE->IsClosed() )
            { SER( SM_ERR ); }
          m_aSideEdges[ii]  = pE;
          m_aSideAtHigh[ii] = ( pE->GetStartVertex() != m_aVertices[ii] );
      }
  }
  return SM_SUCCESS;

} // end SmRingFilletRemoval::CollectRingTopology

/*******************************************************************//**
PURPOSE: SolidCraft change: the sharp edge: where the two tangent faces'
   surfaces meet, the closed curve nearest to the rails, running as rail 0
   and starting at the new vertex.

NOTES: The new vertex lies where the sharp edge crosses a tangent face's
   seam (its seam edge is to run on to it), or where the curve starts when
   neither tangent face has a seam at the rail.
***********************************************************************/
SmStatus SmRingFilletRemoval::IntersectTangentSurfaces()
{
  const SmContext *pContext = m_pFaceToDelete->GetContext();
  NER( pContext );
  SmSurface *pSurf0 = m_aTanFaces[0]->GetSurface();  NER( pSurf0 );
  SmSurface *pSurf1 = m_aTanFaces[1]->GetSurface();  NER( pSurf1 );

  SmBoolean bUseSurfaceEdges[2] = { TRUE, TRUE };
  SmApproxTol3d sSSITol = m_d3dTol / 4.0;
  double dAngTolRad = SM_DEG2RAD( 20.0 );
  SmTArray<SmCurve*> s3DCurves;          SmObjsDelete<SmCurve*> sClean3D( &s3DCurves );
  SmTArray<SmCurve*> sUVCurves0;         SmObjsDelete<SmCurve*> sCleanUV0( &sUVCurves0 );
  SmTArray<SmCurve*> sUVCurves1;         SmObjsDelete<SmCurve*> sCleanUV1( &sUVCurves1 );
  SER( pSurf0->GlobalSurfaceIntersect( *pContext, pSurf0->GetNaturalUVDomain(), *pSurf1, pSurf1->GetNaturalUVDomain(),
                                       bUseSurfaceEdges, &sSSITol, &dAngTolRad, &s3DCurves, &sUVCurves0, &sUVCurves1,
                                       NULL, NULL ));

  // The rails' samples, and how far apart the rails are: a sharp edge lies within a few
  // chords of its fillet's rails (r tan(phi/2) from each, the chord being 2 r sin(phi/2),
  // for faces turning by phi < 170 degrees: within 6 chords).
  const ULONG lNumSamples = 8;
  SmPoint3d sRailPts[2][lNumSamples];
  double dChord = 0.0;
  for ( ULONG ii = 0; ii < 2; ii++ )
  {
      for ( ULONG kk = 0; kk < lNumSamples; kk++ )
      {
          SER( m_aRails[ii]->GetPrimaryEdgeuse()->NormalizedEvaluate( ( kk + 0.5 ) / lNumSamples, FALSE, sRailPts[ii][kk] ));
      }
  }
  for ( ULONG kk = 0; kk < lNumSamples; kk++ )
  {
      SmBoolean bFound = FALSE;
      double dParam = 0.0, dDist = 0.0;
      SmCurve *pRail1 = m_aRails[1]->GetCurve();  NER( pRail1 );
      SER( pRail1->DropPoint( m_aRails[1]->GetInterval(), sRailPts[0][kk], NULL, m_d3dTol, NULL, bFound, dParam, dDist ));
      if ( !bFound )
        { SER( SM_ERR ); }
      dChord = smos_Max( dChord, dDist );
  }
  double dReach = 6.0 * dChord;

  // The closed curve every rail sample lies within reach of, the nearest.
  ULONG  lFound = s3DCurves.GetSize();
  double dBest  = SM_BIG_DOUBLE;
  for ( ULONG cc = 0; cc < s3DCurves.GetSize(); cc++ )
  {
      SmCurve *pC = s3DCurves[cc];
      if ( pC == NULL || pC->IsDegenerate( m_d3dTol ) || !pC->IsClosed( pC->GetNaturalInterval(), m_d3dTol ) )
        { continue; }
      double dWorst = 0.0;
      for ( ULONG ii = 0; ii < 2 && dWorst <= dReach; ii++ )
      {
          for ( ULONG kk = 0; kk < lNumSamples && dWorst <= dReach; kk++ )
          {
              SmBoolean bFound = FALSE;
              double dParam = 0.0, dDist = 0.0;
              SER( pC->DropPoint( pC->GetNaturalInterval(), sRailPts[ii][kk], NULL, m_d3dTol, NULL, bFound, dParam, dDist ));
              dWorst = bFound ? smos_Max( dWorst, dDist ) : SM_BIG_DOUBLE;
          }
      }
      if ( dWorst <= dReach && dWorst < dBest )
      {
          dBest  = dWorst;
          lFound = cc;
      }
  }
  if ( lFound >= s3DCurves.GetSize() )
  {
      // The surfaces don't meet near the fillet within their domains: the sharp edge would lie
      // on a continuation of one of them.
      m_bNotExact = TRUE;
      SER( SM_ERR );
  }
  SmBSplineCurve *pSharp = SM_CAST_PTR( SmBSplineCurve, s3DCurves[lFound] );
  NER( pSharp );
  s3DCurves[lFound] = NULL; // ours now
  SmObjDelete sCleanSharp( pSharp );

  // Run as rail 0.
  {
      SmPoint3d sPt;
      SmVector3d sPV[2], sQV[2];
      SmCurve *pRail0 = m_aRails[0]->GetCurve();  NER( pRail0 );
      double dMid = m_aRails[0]->GetInterval().Evaluate( 0.5 );
      SER( pRail0->Evaluate( dMid, 1, TRUE, sPV ));
      SmBoolean bFound = FALSE;
      double dParam = 0.0, dDist = 0.0;
      SER( pSharp->DropPoint( pSharp->GetNaturalInterval(), sPV[0], NULL, m_d3dTol, NULL, bFound, dParam, dDist ));
      if ( !bFound )
        { SER( SM_ERR ); }
      SER( pSharp->Evaluate( dParam, 1, TRUE, sQV ));
      // the rail's own direction: its edge's orientation on the curve
      if ( sPV[1].Dot( sQV[1] ) < 0.0 )
      {
          SmExtent1d sOld = pSharp->GetNaturalInterval(), sNew;
          SER( pSharp->ReverseParameterization( sOld, sNew ));
      }
  }

  // The new vertex: where the sharp edge crosses a tangent face's seam at its rail's vertex.
  SmExtent1d sIvl = pSharp->GetNaturalInterval();
  double dStartAt = sIvl.GetMin();
  SmBoolean bHaveVertex = FALSE;
  for ( ULONG ii = 0; ii < 2; ii++ )
  {
      SmEdge *pSide = m_aSideEdges[ii];
      if ( pSide == NULL )
        { continue; }
      SmBSplineCurve *pSideCurve = SM_CAST_PTR( SmBSplineCurve, pSide->GetCurve() );  NER( pSideCurve );
      SmSolution      sData[8];
      SmSolutionArray sSolutions( 8, sData );
      SER( pSharp->GlobalCurveIntersect( sIvl, *pSideCurve, pSideCurve->GetNaturalInterval(), m_d3dTol, sSolutions ));
      if ( sSolutions.GetSize() == 0 )
      {
          // The seam's curve stops short of the sharp edge (the fillet trimmed it): a line runs on
          // along itself (its end span continued: the same line), anything else would have to be
          // continued past its curve.
          if ( !pSideCurve->IsLinear( m_d3dTol ) )
          {
              m_bNotExact = TRUE;
              SER( SM_ERR );
          }
          SmExtent1d sSideIvl = pSideCurve->GetNaturalInterval();
          double dEnd = m_aSideAtHigh[ii] ? sSideIvl.GetMax() : sSideIvl.GetMin();
          SmVector3d sPV[2];
          SER( pSideCurve->Evaluate( dEnd, 1, TRUE, sPV ));
          double dSpeed = sPV[1].Length();
          if ( !( dSpeed > SM_EFF_ZERO ) )
            { SER( SM_ERR ); }
          double dReachParam = ( dReach + 10.0 * m_d3dTol ) / dSpeed;
          double dFar = m_aSideAtHigh[ii] ? dEnd + dReachParam : dEnd - dReachParam;
          SmBSplineCurve *pLonger = NULL;
          SER( pSideCurve->CreateExtendedCurve( *pContext, dFar, SM_CT_CINFINITY, pLonger, TRUE ));
          NER( pLonger );
          m_aSideCurves[ii] = pLonger;
          pSideCurve = pLonger;
          SER( pSharp->GlobalCurveIntersect( sIvl, *pSideCurve, pSideCurve->GetNaturalInterval(), m_d3dTol, sSolutions ));
          if ( sSolutions.GetSize() == 0 )
            { SER( SM_ERR ); }
      }
      // One crossing (a closed curve starting there meets it at both its ends).
      SmPoint3d sAt;
      SER( pSharp->EvaluatePoint( sSolutions[0].m_vStart[0], sAt ));
      for ( ULONG kk = 1; kk < sSolutions.GetSize(); kk++ )
      {
          SmPoint3d sOther;
          SER( pSharp->EvaluatePoint( sSolutions[kk].m_vStart[0], sOther ));
          if ( sOther.DistanceBetween( sAt ) > m_d3dTol )
            { SER( SM_ERR ); }
      }
      if ( bHaveVertex )
      {
          // Both tangent faces have a seam there: they must cross the sharp edge together.
          if ( sAt.DistanceBetween( m_sNewVertex ) > m_d3dTol )
            { SER( SM_ERR ); }
      }
      else
      {
          m_sNewVertex = sAt;
          dStartAt     = sSolutions[0].m_vStart[0];
          bHaveVertex  = TRUE;
      }
      m_aSideParams[ii] = sSolutions[0].m_vStart[1];
  }
  if ( !bHaveVertex )
  {
      SER( pSharp->EvaluatePoint( sIvl.GetMin(), m_sNewVertex ));
  }

  // Start the closed curve at the new vertex: the piece after it, then the piece before it.
  double dSnap = SM_EFF_ZERO * ( 1.0 + sIvl.GetMaxDimension() );
  if ( dStartAt > sIvl.GetMin() + dSnap && dStartAt < sIvl.GetMax() - dSnap )
  {
      SmCurve *pCopy = NULL;
      SER( pSharp->Copy( *pContext, pCopy ));
      SmBSplineCurve *pBefore = SM_CAST_PTR( SmBSplineCurve, pCopy );
      SmObjDelete sCleanBefore( pCopy );
      NER( pBefore );
      SmExtent1d sAfterIvl( dStartAt, sIvl.GetMax() );
      SmExtent1d sBeforeIvl( sIvl.GetMin(), dStartAt );
      SER( pSharp->Trim( sAfterIvl ));
      SER( pBefore->Trim( sBeforeIvl ));
      double dGap = m_d3dTol;
      SER( pSharp->JoinWith( 1, pBefore, 0, &dGap ));
  }
  sCleanSharp.Clear();
  m_pNewEdgeCurve = pSharp;
  return SM_SUCCESS;

} // end SmRingFilletRemoval::IntersectTangentSurfaces

/*******************************************************************//**
PURPOSE: SolidCraft change: a tangent face's seam at a rail's vertex runs on
   to the new vertex: on its own curve (found by IntersectTangentSurfaces())
   and on its UV trims, each continued along its own line.

NOTES: A UV trim that does not reach a point of the face's surface at the
   new vertex fails the removal: it would have to be rebuilt.
***********************************************************************/
SmStatus SmRingFilletRemoval::FindSideEdgeExtension( ULONG lWhichSide )
{
  SmEdge *pSide = m_aSideEdges[lWhichSide];
  if ( pSide == NULL )
    { return SM_SUCCESS; }
  const SmContext *pContext = pSide->GetContext();  NER( pContext );
  SmSurface *pSurf = m_aTanFaces[lWhichSide]->GetSurface();  NER( pSurf );
  double dNew = m_aSideParams[lWhichSide];

  SmTArray<SmEdgeuse*> sUses;
  SER( pSide->GetEdgeuses( sUses ));
  ULONG lPair = 0;
  for ( ULONG ii = 0; ii < sUses.GetSize(); ii++ )
  {
      SmEdgeuse *pOwner = NULL;
      SmBSplineCurve *pUV = sUses[ii]->GetUVTrimCurve( &pOwner );
      if ( pUV == NULL || pOwner != sUses[ii] )
        { continue; }
      if ( lPair > 1 )
        { SER( SM_ERR ); }
      SmBSplineCurve *pExtended = NULL;
      if ( pUV->GetNaturalInterval().ContainsValue( dNew ) )
      {
          SmCurve *pCopy = NULL;
          SER( pUV->Copy( *pContext, pCopy ));
          pExtended = SM_CAST_PTR( SmBSplineCurve, pCopy );
      }
      else
      {
          SER( pUV->CreateExtendedCurve( *pContext, dNew, SM_CT_CINFINITY, pExtended, TRUE ));
      }
      NER( pExtended );
      m_aSideUVs[lWhichSide][lPair]      = pExtended;
      m_aSideUVOwners[lWhichSide][lPair] = sUses[ii];
      lPair++;

      // It must reach the face's surface at the new vertex.
      SmPoint3d sUVAt, sAt;
      SER( pExtended->EvaluatePoint( dNew, sUVAt ));
      SER( pSurf->EvaluatePoint( SmPoint2d( sUVAt[0], sUVAt[1] ), sAt ));
      if ( sAt.DistanceBetween( m_sNewVertex ) > m_d3dTol )
        { SER( SM_ERR ); }
  }
  if ( lPair == 0 )
    { SER( SM_ERR ); }
  return SM_SUCCESS;

} // end SmRingFilletRemoval::FindSideEdgeExtension

/*******************************************************************//**
PURPOSE: SolidCraft change: change the Brep: the fillet goes, its rails
   become the sharp edge.

NOTES:
***********************************************************************/
SmStatus SmRingFilletRemoval::DoRemoval()
{
  SmBrep *pBrep = m_pFaceToDelete->GetBrep();  NER( pBrep );
  NER( m_pNewEdgeCurve );
  pBrep->Notify( SM_NO_PRE_EDIT, pBrep, NULL, NULL );

  // The fillet goes; its seam (a wire edge then) with it.
  SER( pBrep->DeleteFace( m_pFaceToDelete, TRUE, FALSE ));
  m_pFaceToDelete = NULL;
  m_pSeam = NULL;

  // The tangent faces' seams run on to the new vertex.
  for ( ULONG ii = 0; ii < 2; ii++ )
  {
      SmEdge *pSide = m_aSideEdges[ii];
      if ( pSide == NULL )
        { continue; }
      SmExtent1d sSideIvl = pSide->GetInterval();
      if ( m_aSideAtHigh[ii] )
        { sSideIvl.SetMinMax( sSideIvl.GetMin(), m_aSideParams[ii] ); }
      else
        { sSideIvl.SetMinMax( m_aSideParams[ii], sSideIvl.GetMax() ); }
      if ( m_aSideCurves[ii] != NULL )
      {
          // Its line continued (the old part keeps its parameters); its UV trims go with the
          // old curve and are replaced below.
          pSide->SetCurve( m_aSideCurves[ii], TRUE );
          m_aSideCurves[ii]->SetOwner( pSide );
          m_aSideCurves[ii] = NULL; // the edge's now
      }
      pSide->SetInterval( sSideIvl );
      for ( ULONG jj = 0; jj < 2; jj++ )
      {
          if ( m_aSideUVOwners[ii][jj] == NULL )
            { continue; }
          m_aSideUVOwners[ii][jj]->SetUVTrimCurve( m_aSideUVs[ii][jj], 0.0, TRUE );
          m_aSideUVs[ii][jj] = NULL; // the edgeuse's now
      }
  }

  // Both vertices to the new one, both rails onto the sharp edge's curve.
  m_aVertices[0]->SetPoint( m_sNewVertex );
  m_aVertices[1]->SetPoint( m_sNewVertex );

  SmEdge *pSurviveEdge = m_aRails[0];
  SmEdge *pDeleteEdge  = m_aRails[1];
  SmVector3d sPV0[2], sPV1[2];
  SER( pSurviveEdge->GetCurve()->Evaluate( pSurviveEdge->GetInterval().Evaluate( 0.5 ), 1, TRUE, sPV0 ));
  SER( pDeleteEdge ->GetCurve()->Evaluate( pDeleteEdge ->GetInterval().Evaluate( 0.5 ), 1, TRUE, sPV1 ));
  SmOrientType eOrient = ( sPV0[1].Dot( sPV1[1] ) > 0.0 ) ? SM_OT_SAME : SM_OT_OPPOSITE;

  SmExtent1d sNewIvl = m_pNewEdgeCurve->GetNaturalInterval();
  pSurviveEdge->SetCurve( m_pNewEdgeCurve, TRUE ); // deletes the rail's curve and UV trims
  m_pNewEdgeCurve->SetOwner( pSurviveEdge );
  m_pNewEdgeCurve = NULL; // the edge's now
  pSurviveEdge->SetInterval( sNewIvl );
  pDeleteEdge->SetCurve( NULL, TRUE );
  pDeleteEdge->SetInterval( sNewIvl );

  // The vertex with a seam on it survives (a seam's end is that vertex).
  SmVertex *pSurviveVtx = ( m_aSideEdges[0] != NULL || m_aSideEdges[1] == NULL ) ? m_aVertices[0] : m_aVertices[1];
  SmVertex *pDeleteVtx  = ( pSurviveVtx == m_aVertices[0] ) ? m_aVertices[1] : m_aVertices[0];
  SER( pBrep->GlueVertices( pSurviveVtx, pDeleteVtx ));
  m_aVertices[0] = m_aVertices[1] = pSurviveVtx;

  SER( pBrep->GlueEdges( pSurviveEdge, eOrient, 0.0, pDeleteEdge ));
  m_aRails[1] = pSurviveEdge;

  // Deleting the fillet split each shell it alone held together (the boss's top from the
  // rest) into one for each part; joined again along the sharp edge they are one shell.
  SER( pBrep->JoinShellsAtVertex( pSurviveVtx ));

  // UV trims of the sharp edge on both tangent faces, and their domains to take them in.
  SmTArray<SmEdgeuse*> sUses;
  SER( pSurviveEdge->GetEdgeuses( sUses ));
  for ( ULONG ii = 0; ii < sUses.GetSize(); ii++ )
  {
      SmBSplineCurve *pUV = NULL;
      SER( sUses[ii]->GetOrCreateUVTrimCurve( pUV ));
      NER( pUV );
  }
  for ( ULONG ii = 0; ii < 2; ii++ )
  {
      SmFace *pFace = m_aTanFaces[ii];
      SmEdgeuse *pEU = pSurviveEdge->GetEdgeuseOfFace( pFace );  NER( pEU );
      SmBSplineCurve *pUV = pEU->GetUVTrimCurve();  NER( pUV );
      SmExtent2d sFaceDomain = pFace->GetUVDomain();
      SmExtent1d sUVIvl = pUV->GetNaturalInterval();
      for ( ULONG kk = 0; kk <= 16; kk++ )
      {
          SmPoint3d sUV;
          SER( pUV->EvaluatePoint( sUVIvl.Evaluate( kk / 16.0 ), sUV ));
          sFaceDomain.AddPoint2d( SmPoint2d( sUV[0], sUV[1] ));
      }
      // The face lies within its surface's domain (nothing was continued).
      SmExtent2d sSurfDomain = pFace->GetSurface()->GetNaturalUVDomain(), sInSurface;
      SER( sFaceDomain.Intersect( sSurfDomain, sInSurface ));
      pFace->SetUVDomain( sInSurface );
  }
  return SM_SUCCESS;

} // end SmRingFilletRemoval::DoRemoval

/*******************************************************************//**
PURPOSE: Delete a fillet face and fill in the hole.

NOTES: Deletes the fillet face and extends and intersects
     the two tangent adjacent faces to create a sharp edge.

  The fillet face must have exactly four edges, with one pair
  of opposite edges being tangent edges, and the other two
  non-tangent.

  The surfaces of the four faces adjacent to the Fillet Face
  being removed must all be (derived from) SmBSplineSurface,
  and the four edge curves coming into the Fillet Face must
  be (derived from) SmBSplineCurve.

  It sometimes happens that a fillet does not meet the criteria
  for removal, but after removing an adjacent fillet, it does.
  For example, two mitered fillets of different radii will result
  in the larger fillet having five edges; removing the smaller
  fillet will leave the larger one with four edges.  It is therefore
  advisable to run this more than once, perhaps in a loop on all
  fillet faces until no more are removed.

***********************************************************************/
SmStatus SmBrep::DeleteFilletFace( SmFace *pFilletFace, double * pd3dTol )
{
  NER( pFilletFace );

  // check state - m_bEditingEnabled must be on
  if ( !m_bEditingEnabled ) 
    SER_MSG(SM_ERR, _T("Tried to edit Brep with m_bEditingEnabled == FALSE"));

  double d3dTol = ( pd3dTol != NULL ) ? *pd3dTol : static_cast<double>(pFilletFace->GetTolerance());

  // SolidCraft change: a ring fillet (closed rails, no ends) has its own removal.
  if ( SmRingFilletRemoval::IsRing( pFilletFace, 2.0 ) )
  {
      SmRingFilletRemoval sRingRemoval( pFilletFace, d3dTol );
      return sRingRemoval.RemoveFace();
  }

  SmFilletRemoval sRemoval( pFilletFace, d3dTol );

  return sRemoval.RemoveFace();

} // end DeleteFilletFace




