// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
// SPDX-License-Identifier: Apache-2.0

/**********************************************************************//**
* FILE NAME --- SmOffsetExecutive.cpp 
* PURPOSE: Source code file for SmOffsetExecutive object.
**********************************************************************/

#include "StdAfx.h"
#include <algorithm>
#include <unordered_set>
#include <utility>
#include <vector>

#include <SmOffsetExecutive.h>
#include <SmSurfaceCache.h>
#include <SmMerge.h>
#include <SmStitch.h>
#include <SmGeomUtility.h>
#include <SmOffsetSurface.h>
#include <SmTopologySolver.h>
#include <SmPrimitiveCreation.h>
#include <SmPlane.h>
#include <SmTrimmingTools.h>
#include <SmAssertArray.h>
#include <SmLoop.h>

#ifdef SM_DEBUG_CODE
#include <SmTopologyTraverser.h>
#endif // SM_DEBUG_CODE


//#define SM_VALIDATE_TOPOLOGY 1

/*******************************************************************//**
PURPOSE: Constructor for the SmOffsetExecutive.

NOTES: 
***********************************************************************/
SmOffsetExecutive::SmOffsetExecutive
  (const SmContext & crContext,
   SmOffsetOperationType eOffsetOperation,
   SmOffsetGeometryCreation & rGeometryCreation,
   SmBrep *pOriginalBrep)
 : m_crContext(crContext), 
   m_eOffsetOperation(eOffsetOperation), 
   m_rGeometryCreation(rGeometryCreation), 
   m_pOriginalBrep(pOriginalBrep),
   m_pOffsetBrep(NULL), 
   m_bInset(FALSE),
   m_bExtendConvexEdges(FALSE),
   m_bMergeResults(FALSE),
   m_bCreateOffsetSolid(FALSE), 
   m_pShellFaces(NULL),
   m_vFUToF(crContext), 
   m_vEUToEU(crContext),
   m_vEToF(crContext), 
   m_vVToF(crContext), 
   m_vBlendToFaces(crContext), 
   m_vSurfExtSurf(*(new (crContext) SmMapPtrToPtr<SmSurface, SmSurface>(&crContext)))
  { }

/*******************************************************************//**
PURPOSE: Destructor for the SmOffsetExecutive

NOTES: 
***********************************************************************/
SmOffsetExecutive::~SmOffsetExecutive()
{
  if (m_pShellFaces) { delete m_pShellFaces; m_pShellFaces = NULL ; }
  SM_ASSERT(&m_vSurfExtSurf != NULL) ; delete &m_vSurfExtSurf ; 

  // These lists can contain pointers into m_pOffsetBrep, which will be deleted.
  // Get the pointers out of there before they become stale.  [B453]
  m_vFUToF.RemoveAll();
  m_vEUToEU.RemoveAll();
  m_vEToF.RemoveAll();
  m_vVToF.RemoveAll();
  m_vBlendToFaces.RemoveAll();

  if (m_pOffsetBrep) { delete m_pOffsetBrep; m_pOffsetBrep = NULL ; }
}


//------------------------------------------
// User interface methods
//------------------------------------------

/*******************************************************************//**
PURPOSE: Static High Level Function to Offset a Brep.

NOTES: 
   Creates a new Brep that is offset from the original.

   In the description of the arguments below, 'convex' and 'concave' refer
   to the direction of the offset.  A convex offset is one where, at an edge,
   the edge is convex in the direction of the offset.  For example, with a
   cube, when offsetting outward, all edges have convex offsets, and when
   offsetting inwards, all offsets are concave.  With a convex offset, the
   two faces meeting at an edge will be offset away from each other, leaving
   a gap that must be filled.  The argument 'bDoExtendedOffset' tells how to
   fill the gap.  With a concave offset, the two offset faces will intersect
   each other; 'bDoSelfInt' tells what to do about that.
***********************************************************************/
SmStatus SmOffsetExecutive::OffsetBrep
  (const SmContext & crContext,
   SmBrep *pBrep,               // in : Will be consumed in creation of offset brep.
   double dOffsetDistance,      // in : Offset distance (negative for insets).
   SmBoolean bDoExtendedOffset, // in : For convex offsets, to close the gap:
                                //      TRUE = Extend and intersect the two offset faces.
                                //      FALSE= Create a fillet-like face to join the faces smoothly.
   SmBoolean bDoSelfInt,        // in : For concave offsets, to handle the self-intersection:
                                //      TRUE = intersect and trim the intersecting offsets.
                                //      FALSE= Do nothing: leave the overlapping faces.
   SmBrep *& rpOffset)          // out: result of the offset operation
{
    rpOffset = NULL;
    
    double dOffsetRadius = smos_Fabs(dOffsetDistance);     
    double dTolerance    = pBrep->GetTolerance();

    SmOffsetGeometryCreation   sOGC(dOffsetRadius, dTolerance);
    SmOffsetExecutive          sExec(crContext, SM_OO_SOLID_OFFSET, sOGC, pBrep);
    if (dOffsetDistance < 0.0)
      { sExec.SetInset(TRUE); }
    sExec.SetExtendConvexEdges(bDoExtendedOffset);
    sExec.SetMergeResults(bDoSelfInt);

    if (sExec.DoSolidOffset(rpOffset) != SM_SUCCESS)  // note: increments unlocked mark value
      { SER(SM_ERR); }

    return SM_SUCCESS;

} // end SmOffsetExecutive::OffsetBrep

/*******************************************************************//**
PURPOSE: Static High Level Function to Shell a Brep.

NOTES:
   Shelling creates an offset of the given Brep and, if the input is a
   manifold body, performs a Boolean Difference operation to create a
   hollow body whose inner and outer shells are the original Brep and
   the offset.

   Faces should be oriented correctly.  Call  pBrep->StitchAndOrient();
   and pBrep->ShrinkGeometry() if needed, before calling this routine,

   See the description of 'convex' and 'concave' in the header notes
   for SmOffsetExecutive::OffsetBrep(), above.
***********************************************************************/
SmStatus SmOffsetExecutive::ShellBrep
  (const SmContext & crContext,               // in : context for new object construction
   SmBrep          * pBrepToShell,            // in : Will be either deleted or passed
                                              //      back in rpShellBrep - just assume it is deleted and
                                              //      assume that rpShellBrep is a brand new brep.
   double            dOffsetDistance,         // in : Offset distance (negative for insets).
   SmBoolean         bDoExtendedOffset,       // in : For convex offsets, to close the gap:
                                              //      TRUE = Extend and intersect the two offset faces.
                                              //      FALSE= Create a fillet-like face to join the faces smoothly.
   SmBoolean         bDoSelfInt,              // in : For concave offsets, to handle the self-intersection:
                                              //      TRUE = intersect and trim the intersecting offsets.
                                              //      FALSE= Do nothing: leave the overlapping faces.
   SmBoolean         bCreateOffsetSolid,      // in : If the input Brep has lamina edges, i.e.,
                                              //      it is some kind of sheet body:
                                              //      TRUE = connect the lamina edges and their
                                              //      offsets with ruled surfaces, to close the
                                              //      result into a solid body.
                                              //      It assumes that we have an open object with 
                                              //      extended offsetting.
                                              //      FALSE= Don't: leave the offsets separate.
   const SmTArray<SmFace*> & crFacesToShell,  // in : A subset of pBrepToShell faces that will
                                              //      be copied to output without an offset
                                              //      If any faces are specified, the resulting Brep
                                              //      will not be a shell, it will just be the
                                              //      offset body.
   SmBrep                 *& rpShellBrep)     // out: result of the shelling operation
{
  SmBoolean bInputOwnershipTransferred = FALSE;
  return ShellBrep(crContext, pBrepToShell, dOffsetDistance,
                   bDoExtendedOffset, bDoSelfInt,
                   bCreateOffsetSolid, crFacesToShell,
                   rpShellBrep, bInputOwnershipTransferred);
}

/*******************************************************************//**
PURPOSE: ShellBrep overload that reports whether ownership of the input
         was transferred to the operation.

NOTES:
   Geometry and legacy input lifetime behavior are identical to ShellBrep().
   rbInputOwnershipTransferred is TRUE when pBrepToShell was deleted or was
   returned as rpShellBrep.  When FALSE, pBrepToShell remains a separate
   allocation owned by the caller, including on early failures.
***********************************************************************/
SmStatus SmOffsetExecutive::ShellBrep
  (const SmContext & crContext,
   SmBrep          * pBrepToShell,
   double            dOffsetDistance,
   SmBoolean         bDoExtendedOffset,
   SmBoolean         bDoSelfInt,
   SmBoolean         bCreateOffsetSolid,
   const SmTArray<SmFace*> & crFacesToShell,
   SmBrep                 *& rpShellBrep,
   SmBoolean              & rbInputOwnershipTransferred)
{
  // init output
  rpShellBrep = NULL;
  rbInputOwnershipTransferred = FALSE;
  
  // locals
  double  dOffsetRadius = smos_Fabs(dOffsetDistance); 
  double  dTolerance    = pBrepToShell->GetTolerance();
  SmBrep *pOffsetBrep   = NULL;
  SmBrep *pBrepShell    = NULL;
  SmTArray<SmFace*> sShellFaces, sShellFaces1;
  SmBoolean bInputIsManifoldSolid = pBrepToShell->IsManifoldSolid();
  
  // Error if offset radius is within intersection tolerance
  if (dOffsetRadius < 2.*dTolerance)
    { SER(SM_ERR_INVALID_INPUT); }  

  {
    // set up offset executive
    SmOffsetGeometryCreation sOGC( dOffsetRadius, dTolerance);
    SmOffsetExecutive sExec(crContext,SM_OO_SOLID_OFFSET,sOGC,pBrepToShell);

    // negative OffsetDistances are used for insets
    if (dOffsetDistance < 0.0) { sExec.SetInset(TRUE); }

    // set OffsetExecutive parameters
    sExec.SetShellFaces(crFacesToShell);
    sExec.SetExtendConvexEdges(bDoExtendedOffset);
    sExec.SetMergeResults(bDoSelfInt);
    sExec.SetCreateOffsetSolid(bCreateOffsetSolid);
    
    // create the offset
    SmStatus eStat = sExec.DoSolidOffset( pBrepShell );  // note: increments unlocked mark value
    if ( eStat != SM_SUCCESS )
      { SER(eStat); }
#ifdef SM_DEBUG_CODE
    SmBoolean bDebugMe = FALSE;
    if (bDebugMe)
    {
        pBrepToShell->Dump();
        pBrepShell->Dump();
        smgfx_Erase() ;
        smgfx_SetLook(2, 2, 0, 0, 1) ; if (pBrepToShell) pBrepToShell->Draw(TRUE) ; sm_GraphicsLoop() ;
        smgfx_SetLook(3, 2, 0, 1, 0) ; if (pBrepShell) pBrepShell->Draw(TRUE) ; sm_GraphicsLoop() ;
        sm_GraphicsLoop() ;
    }
#endif // SM_DEBUG_CODE

    // watch for memory leaks
    pOffsetBrep = sExec.GetOffsetBrep();
    if (pOffsetBrep) { delete pOffsetBrep; pOffsetBrep = NULL ; }
  }

  // An inset of a solid whose offset leaves nothing (a wall as thick as the body's inradius or
  // thicker) has no cavity to take away: the difference would hand back the input unhollowed.
  // An offset that is not a closed solid bounds no cavity either.
  if (dOffsetDistance < 0.0 && bInputIsManifoldSolid)
    {
      SmStatus eNoCavity = SM_SUCCESS;
      if (pBrepShell == NULL || pBrepShell->GetNumFaces() == 0)
        { eNoCavity = SM_ERR_INVALID_INPUT; }
      else if (!pBrepShell->IsManifoldSolid())
        { eNoCavity = SM_ERR; }
      if (eNoCavity != SM_SUCCESS)
        {
          delete pBrepShell;
          SER(eNoCavity);
        }
    }

  // when offset result is not manifold
  if (!pBrepShell->IsManifoldSolid()) 
    {
      // and not CreateOffsetSolid
      if (!bCreateOffsetSolid) 
        {
          // return the nonManifold result
          rpShellBrep = pBrepShell;
          return SM_SUCCESS;   // was SM_ERR: ??
        }
    }
  
  double dAngleTol = 20.0*SM_PI/180.0;

  // Replace all SmOffsetSurfaces with approximate SmBSplineSurfaces
  SER( ReplaceImplicitOffsets( pBrepShell, dTolerance ));
  // Changed tolerance arg from 10.0*dTolerance to improve downstream tolerance consistency [B577]
#ifdef SM_DEBUG_CODE
  SmBoolean bDebugMe = FALSE;
  if (bDebugMe)
  {
      pBrepToShell->Dump();
      pBrepShell->Dump();
      smgfx_Erase() ;
      smgfx_SetLook(2, 2, 0, 0, 1) ; if (pBrepToShell) pBrepToShell->Draw(TRUE) ; sm_GraphicsLoop() ;
      smgfx_SetLook(3, 2, 0, 1, 0) ; if (pBrepShell) pBrepShell->Draw(TRUE) ; sm_GraphicsLoop() ;
      sm_GraphicsLoop() ;
  }
#endif // SM_DEBUG_CODE
  // 
  // Only sheet bodies need the original/offset merge.
  if (bCreateOffsetSolid && !bInputIsManifoldSolid)
    {
#ifdef SM_DEBUG_CODE
      if(bDebugMe)
        {
          pBrepToShell->Dump();
          pBrepShell->Dump();
          smgfx_Erase() ;
          smgfx_SetLook(2,2, 0,0,1) ; if(pBrepToShell) pBrepToShell->Draw(TRUE) ; sm_GraphicsLoop() ;
          smgfx_SetLook(3,2, 0,1,0) ; if(pBrepShell  ) pBrepShell  ->Draw(TRUE) ; sm_GraphicsLoop() ;
          sm_GraphicsLoop() ;
        }
#endif // SM_DEBUG_CODE

      // combine topology without any intersection or coincidence checking
      SER(pBrepShell->MergeBrep(*pBrepToShell));
      SM_ASSERT(pBrepToShell != NULL) ; delete pBrepToShell ; pBrepToShell = NULL ;
      rbInputOwnershipTransferred = TRUE;

#ifdef SM_VALIDATE_TOPOLOGY
      pBrepShell->ValidatePointers();
#endif // SM_VALIDATE_TOPOLOGY

      // Stitch lamina edges.
      // for 4 progressively larger tolerances (or up to half the offset distance)
      ULONG lStitched, ii, lLamina=1;
      double dMaxVGap = 0.0, dMaxEGap = 0.0;
      double dStitchTol3d = dTolerance / 10.0;
      for( ii=0; ii<4 && lLamina >0; ii++ )
        {
          SER_DELETE(pBrepShell->StitchFaces( dStitchTol3d, lStitched, lLamina, dMaxVGap, dMaxEGap ),
                     pBrepShell);

          dStitchTol3d *= 10.0;
          if ( dStitchTol3d > smos_Fabs( dOffsetDistance ) / 2.0 ) { break; }
        }

#ifdef SM_VALIDATE_TOPOLOGY
      pBrepShell->ValidatePointers();
#endif // SM_VALIDATE_TOPOLOGY

      rpShellBrep = pBrepShell;
    }

  else if (dOffsetDistance < 0.0) 
    {
      SmBrep *pResult = NULL;
      SmMerge sMerge(crContext,pBrepToShell,pBrepShell,dTolerance,dAngleTol);
      SmStatus eStat = sMerge.ManifoldBoolean(SM_BO_DIFFERENCE, pResult);
      if (eStat != SM_SUCCESS)
        { SER(eStat); }
      rbInputOwnershipTransferred = TRUE;
      rpShellBrep = pResult;
    }
  else 
    {
      SmBrep *pResult = NULL;
      SmMerge sMerge(crContext,pBrepShell,pBrepToShell,dTolerance,dAngleTol);
      SmStatus eStat = sMerge.ManifoldBoolean(SM_BO_DIFFERENCE, pResult);
      if (eStat != SM_SUCCESS)
        { SER(eStat); }
      rbInputOwnershipTransferred = TRUE;
      rpShellBrep = pResult;
    }
#ifdef SM_DEBUG_CODE
  if (bDebugMe)
  {
      rpShellBrep->Dump();
      smgfx_Erase() ;
      smgfx_SetLook(2, 2, 0, 0, 1) ; if (rpShellBrep) rpShellBrep->Draw(TRUE) ; sm_GraphicsLoop() ;
      sm_GraphicsLoop() ;
  }
#endif // SM_DEBUG_CODE
  return SM_SUCCESS;

} // end SmOffsetExecutive::ShellBrep

/*******************************************************************//**
PURPOSE: This is the main top level method to offset a solid brep.

NOTES:
   Increments unlocked mark value.

METHOD ---
  1. Every face is copied either in its current position or in an
     offset position and added to rpResult
  2. A cap face is created and added to rpResult for every edge
     where a gap was created in rpResult where a pair of radial
     faces were offset into a convex sector.
  3. A cap face is created and added to rpResult for every vertex
     where a gap was created in rpResult by offsetting all the 
     faces attached to the vertex.  A cap face is only needed
     at a vertex when every edge attached to the vertex needed
     a cap face.
***********************************************************************/
SmStatus SmOffsetExecutive::DoSolidOffset
  (SmBrep *& rpResult,             // out: Brep containing offset surfaces for all faces
                                   //      and edges and vertices as needed to fill gaps.
   SmBoolean bSkipStitchMerge)        // in : TRUE = pResult contains all offset Surfaces not stitched into a manifold Brep.
                                   //               surface intersections will not be marked.
                                   //      FALSE= pResult contains All offset surfaces booleaned together.
                                   //               surfaces will be trimmed at all surf/surf intersections.
                                   //      default:[FALSE]
{
#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe   = FALSE;
SmBoolean bDrawFaces = FALSE;
  if(bDebugMe)
    {
      if(bDebugMe && bDrawFaces) 
        { bDebugMe = TRUE ;
          bDrawFaces = TRUE ; 
        }
      else
        { bDebugMe = FALSE ;
          bDrawFaces = FALSE ; 
        }
      Dump() ;
      SM_DUMP_AND_ASSERT_VALID(m_pOriginalBrep) ;

      smgfx_Erase() ;
      smgfx_SetLook(1,2, 0,0,1) ; if(m_pOriginalBrep) m_pOriginalBrep->DrawUV(TRUE) ; sm_GraphicsLoop() ;
      smgfx_SetLook(1,2, 0,1,0) ; if(m_pOffsetBrep) m_pOffsetBrep->DrawUV(TRUE) ; sm_GraphicsLoop() ;
      sm_GraphicsLoop() ;
    }
#endif // SM_DEBUG_CODE

  // pass the call along for extending convex edges 
  if (m_bExtendConvexEdges) 
    {
      SER(DoExtendedOffset(rpResult));   // note: increments unlocked mark value
      return SM_SUCCESS;
    }

  NER(m_pOriginalBrep);

  SmTemporaryChangeValue< SmBoolean > sDoingBool(SM_CONST_CAST(SmContext*, &m_crContext)->GetDoingBooleanRef(), TRUE);

  // allocate an empty new m_pOffsetBrep to contain intermediate results
  m_pOffsetBrep = new (m_crContext) SmBrep();
  // No:   [B453]
  //SmObjDelete sCleanOffBrep( m_pOffsetBrep );

  SM_OLDTOL_LINE m_pOffsetBrep->SetTolerance(m_rGeometryCreation.GetThisApproxTol3d());

// Remove Composites
// //cbi_CEdge13:
// #ifdef SM_NO_COMPOSITES
//   m_pOffsetBrep->m_bMakeComposites = FALSE;
// #else
//   m_pOffsetBrep->m_bMakeComposites = TRUE;
// #endif // SM_NO_COMPOSITES

  m_pOffsetBrep->m_bEditingEnabled = TRUE;

  // First step is to offset faces

  // face offset locals
  SmTArray<SmFace*>           sFaces;
  m_pOriginalBrep->GetFaces(sFaces);
  SmSurface*                  sOffsetSurface = NULL;
  SmSSIData                   sIntData;
  SmTArray<SmCurve*>          s3DCurves;
  SmTArray<SmBSplineCurve*>   sUVCurves1;
  SmTArray<SmOrientType>      sOrients;
  SmTArray<SmPoint3d>         sLoopPoints;
  SmTArray<ULONG>             sCurveLoops;

  // for every OriginBrep face - build a copy or an offset copy in m_pOffsetBrep
  for (ULONG ii=0; ii<sFaces.GetSize(); ii++) 
    {
      SmFace *pOrigFace = sFaces[ii];
      ULONG lIndex;

      // when given a shell face list - just copy the shell faces into 
      if (   m_pShellFaces 
          && m_pShellFaces->FindElement(pOrigFace,lIndex)) 
        {
          // Just copy this face instead of offsetting it

          // get ordered copies of all Face->Edge->Curves and Face->Edge->Vertices
          SER(pOrigFace->CreateCurvesFromFace(m_crContext,   // in : context for new object construction
                                              SM_OT_SAME,    // in : use Faceuse whose orientation is the same as eOrientation
                                              sCurveLoops,   // out: number of edges in each edge loop
                                              &s3DCurves,    // out: optional copies of all 3D loop curves - note: one or both of pOpt3DCurves and
                                              &sUVCurves1,   // out: optional copies of all 2D loop curves -         pOptUVCurves must be nonNULL
                                              sOrients,      // out: Curve orientation in relation to its loop for each curve
                                              sLoopPoints)); // out: point position of each of face's vertex loops
          SmRegion  *pNewRegion = NULL;
          SmShell   *pNewShell = NULL;
          SmFace    *pNewFace = NULL;
          SmSurface *pOrigSurface = pOrigFace->GetSurface(); NER(pOrigSurface);
          SmSurface *pSurface = NULL;

          // Copy pOrigSurface for new face, when possible as an analytic surface
          SER(pOrigSurface->CopyAndAddAnalytics(m_crContext,pSurface));

          // Make a copy of this face in m_pOffsetBrep->InfiniteRegion
          SER(m_pOffsetBrep->MakeFaceWithCurves
                   (m_pOffsetBrep->GetInfiniteRegion(), // in : region to contain new topology objects
                    sCurveLoops,                        // in : 1 entry per loop, value = loop edge count, 1st entry=outer loop
                    &s3DCurves,                         // in : opt ordered 3d trimming curves assigned to loops per sLoopEUCounts
                    &sUVCurves1,                        // in : opt ordered 2d trimming curves assigned to loops per sLoopEUCounts
                    sOrients,                           // in : associated orients for each trimming curve, SM_OT_SAME or SM_OT_OPPOSITE
                    sLoopPoints,                        // in : Point positions to build SmVertex VertexLoops
                    pSurface,                           // in : new face->Surface
                    pOrigFace->GetUVDomain(),           // in : domain of Surface used by face
                    SM_OT_SAME,                         // in : Surface orient, oneof SM_OT_SAME or SM_OT_OPPOSITE
                    pNewRegion,                         // out: New region if any. NULL when building trimmed surfaces, may be NotNULL for solids.
                    pNewShell,                          // out: New shell if any.  Trimmed surfaces always create a new shell.
                    pNewFace));                         // out: the new face

          // relate InfiniteRegion bounding OriginalBrep->faceuse to OffsetBrep->new Face
          // relate OriginalBrep->faceuse->edgeuses with m_pOffsetBrep->new Face edgeuses
          SmFaceuse *pFU = pOrigFace->GetUpwardFaceuse();
          if (!pFU->GetShell()->GetRegion()->IsVoid()) 
            {
              pFU = pFU->GetMate();
            }

          // map InfiniteRegion bounding pFU to NewFace and map their edgeuses to one another
          SER( m_vFUToF.RelatePair(pFU,pNewFace) );
          SER(MapEdgeuses(pFU,pNewFace));

#ifdef SM_DEBUG_CODE
          // draw 
          if(bDebugMe)
            {
              Dump() ;
              SmExtent2d sSurfaceDomain = pFU->GetFace()->GetSurface()->GetNaturalUVDomain() ;
              pFU->GetFace()->GetSurface()->Dump() ;
              if(sOffsetSurface)
              {
                SmExtent2d sOffsetDomain = sOffsetSurface->GetNaturalUVDomain();
                sOffsetSurface->Dump();
              }
              smgfx_Erase() ;
              smgfx_SetLook(1,2, 0,0,1) ; m_pOriginalBrep->Draw(TRUE) ; sm_GraphicsLoop() ;
              smgfx_SetLook(2,3, 0,1,0) ; m_pOffsetBrep->Draw(TRUE) ; sm_GraphicsLoop() ;
              smgfx_SetLook(1,2, 0,1,1) ; pOrigFace->GetSurface()->DrawUV(3,3) ; sm_GraphicsLoop() ;
              smgfx_SetLook(1,2, 0,0,0) ; pOrigFace->Draw(SM_DM_CROSSHATCH) ; sm_GraphicsLoop() ;
              smgfx_SetLook(1,2, 1,0,1) ; pNewFace->GetSurface()->DrawUV(7,7) ; sm_GraphicsLoop() ;
              smgfx_SetLook(1,2, 1,.5,0); pNewFace->Draw(SM_DM_CROSSHATCH, 5,5) ; sm_GraphicsLoop() ;
              smgfx_SetLook(3,4, 0,1,1) ; Draw() ; sm_GraphicsLoop() ;
              sm_GraphicsLoop() ;
            }
#endif // SM_DEBUG_CODE

        } // end pOrigFace is just copied branch
      else // pOrigFace is to be offset and copied branch
        {
          // m_bInset == FALSE: get faceuse bounding infinite region
          // m_bInset == TRUE : get faceuse bounding interior region
          SmFaceuse *pFU = pOrigFace->GetUpwardFaceuse();
          if (pFU->GetShell()->GetRegion() != m_pOriginalBrep->GetInfiniteRegion()) 
            { pFU = pFU->GetMate(); }
          if (m_bInset) { pFU = pFU->GetMate(); }

          // build offset surfaces from faceuse - assume out surface(s) domains map may be trimmed
          // but not scaled from the input surface domains so that OutSurface->UVPts = InSurface->UVPts.
          //  - don't check for self-intersections.
          //    Not all virtual functions are the same (they should be)
          //    class SmSurfOfExtrusion - identifies and removes self-intersections
          //                              no matter the bSkipSelfIntersection boolean flag.
          //    class SmBSplineSurface  - generates self intersecting surfaces when
          //                              bSkipSelfIntersection == TRUE.
          //  - when bSkipSelfIntersection == FALSE single surface offsets are checked for self intersections
          sOffsetSurface = NULL;
          SER_DELETE(m_rGeometryCreation.FaceuseOffset(m_crContext,       // in : context for new obj construction
                                                pFU,               // in : Faceuse->Surface to copy and offset
                                                TRUE,              // in : TRUE =
                                                                   //      FALSE=
                                                sOffsetSurface,    // out: New Surface (more than 1 for self-XSect cases)
                                                sIntData),       // out: List of SelfXSect curves (when bSkipSelfIntersection == FALSE)
                                                sOffsetSurface) ;  // JLMCC hunting memory leaks
          if ( sOffsetSurface == NULL )
            { continue; }

          // build a face in m_pOffsetBrep made with sOffsetSurface and pFU->Edge->UVTrimCurves. 
          //      map pFU to NewFace in m_vFUToF
          //  and map all pFU->edgeuses to NewFace->Edgeuses in m_vEUToEU
          //  When OffsetSurface->Domain has been trimmed (due to singularities) the
          //  pFU->Edge->UVTrimCurves are trimmed to their intersections with the smaller domain.
          SER(BuildOffsetOfFaceuse(pFU, sOffsetSurface, sIntData));

#ifdef SM_DEBUG_CODE
          // draw 
           if(bDebugMe)
            {
              ULONG di = 0;

              Dump() ;
              SmExtent2d sSurfaceDomain = pFU->GetFace()->GetSurface()->GetNaturalUVDomain() ;
              SmExtent2d sOffsetDomain  = (sOffsetSurface) ? sOffsetSurface->GetNaturalUVDomain() : SmExtent2d(0,1,0,1) ;
              pFU->GetFace()->GetSurface()->Dump() ;
              
              SmBoolean bOk = SM_DUMP_AND_ASSERT_VALID(sOffsetSurface) ;
              if(!bOk)
              { // a place for a breakpoint
                di++; di--;
              }
             

              smgfx_Erase() ;
              smgfx_SetLook(1,2, 0,0,1) ; m_pOriginalBrep->Draw(TRUE) ; sm_GraphicsLoop() ;
              smgfx_SetLook(1,2, 0,1,0) ; m_pOffsetBrep->Draw(TRUE) ; sm_GraphicsLoop() ;
              smgfx_SetLook(1,2, 0,1,1) ; pOrigFace->GetSurface()->DrawUV() ; sm_GraphicsLoop() ;
              smgfx_SetLook(1,2, 0,0,0) ; pOrigFace->Draw(SM_DM_CROSSHATCH) ; sm_GraphicsLoop() ;
              
              SmFace    *pNewFace = (SmFace *)sOffsetSurface->GetFace() ;
              smgfx_SetLook(1,2, 1,1,0) ; sOffsetSurface->DrawUV() ; sm_GraphicsLoop() ;
              smgfx_SetLook(1,2, 0,0,0) ; pNewFace->Draw(SM_DM_CROSSHATCH) ; sm_GraphicsLoop() ;
       
              smgfx_SetLook(3,4, 0,1,1) ; Draw() ; sm_GraphicsLoop() ;
              sm_GraphicsLoop() ;
            }
#endif // SM_DEBUG_CODE

        } // end pOrigFace is to be offset and copied branch
    } // end iter every m_pOriginalBrep face building a copy into the m_pOffsetBrep

  // remove curve and surface portions not being used by edges and faces.
  // This may help with self-intersection problems.

  SER(m_pOffsetBrep->ShrinkGeometry());

#ifdef SM_DEBUG_CODE
  // draw and dump input and offset breps after offsetting all faces
  if (bDebugMe) 
    {
      SM_ASSERT_VALID(m_pOriginalBrep) ;
      SM_ASSERT_VALID(m_pOffsetBrep) ;

      ULONG kk ;
      SmTArray<SmFace *> sOriginalFaces ;
      SmTArray<SmFace *> sOffestFaces ;
      m_pOriginalBrep->GetFaces(sOriginalFaces) ;
      m_pOffsetBrep->GetFaces(sOffestFaces);

      smgfx_Erase();
      smgfx_SetLook(1,2, 0,0,1) ; m_pOriginalBrep->Draw(TRUE); sm_GraphicsLoop();
      smgfx_SetLook(2,3, 0,1,0) ; m_pOffsetBrep->Draw(TRUE); sm_GraphicsLoop();
      if(bDrawFaces)
        {
          smgfx_SetLook(1,2, 0,1,1) ; for(kk=0;kk<sOriginalFaces.GetSize();kk++)
                                        { SmFace    *pFace    = sOriginalFaces[kk] ;
                                          SmSurface *pSurface = pFace->GetSurface() ;
                                          pSurface->DrawUV() ; sm_GraphicsLoop() ;
                                        }
          smgfx_SetLook(1,2, 1,1,0) ; for(kk=0;kk<sOffestFaces.GetSize();kk++)
                                        { SmFace    *pFace    = sOffestFaces[kk] ;
                                          SmSurface *pSurface = pFace->GetSurface() ;
                                          pSurface->DrawUV() ; sm_GraphicsLoop() ;
                                        }
        }
      smgfx_SetLook(3,4, 0,1,1) ; Draw() ; sm_GraphicsLoop() ;
      sm_GraphicsLoop();
    }
#endif // SM_DEBUG_CODE

  // for every edge make and add faces to fill gaps created by offsetting radial
  // face pairs into convex regions - the needed faces will be
  //  like circular fillet surfaces

  // Get a list of mapped edgeuses in a standard order
  SmEdgeuse * sEEUsData[65];  
  SmEdgeuse * sSecondsData[65];  
  SmEdge    * sEdgesData[256]; 
   
  SmTArray<SmEdgeuse*> sOffsetEdgeuses;
  SmTArray<SmEdgeuse*> sEEUs   (65,  (SmEdgeuse**)sEEUsData) ;
  SmTArray<SmEdgeuse*> sSeconds(65,  (SmEdgeuse**)sSecondsData) ;
  SmTArray<SmEdge*>    sEdges  (256, (SmEdge**)   sEdgesData) ;

  // for every edge - build list of Brep->edgeuses that have been mapped to offset edgeuses
  m_pOriginalBrep->GetEdges(sEdges);
  for (ULONG ii=0; ii<sEdges.GetSize(); ii++) 
    {
      SmEdge *pE = sEdges[ii];
      SER( pE->GetEdgeuses(sEEUs) );
      for (ULONG ieu=0; ieu<sEEUs.GetSize(); ieu++) 
        {
          SmEdgeuse *pEU = sEEUs[ieu];

          // when edgeuse has been mapped to offset edgeuses
          m_vEUToEU.GetAllSeconds(sSeconds);
          if (sSeconds.GetSize() > 0) 
            {
              // add it to the list
              sOffsetEdgeuses.Add(pEU);
            }
        } // end iter every edgeuse
    } // end iter every edge

  SmTArray<SmEdgeuse*> sOffSetEUS, sRadialEUS;
  SmTArray<SmEdge*>    sNewEdges;

  // Generate edge offsets - cap surfaces to fill gaps made by offsetting radial face pairs into a convex sector 

  // for every OffsetBrep->edgeuse - add a CapFace when needed to fill gap between offset raial face pairs

  // increment and lock an unlocked mark
    {
      SmNewMarkAndLock sMarkLock( m_pOffsetBrep->GetContext(), SM_MT_ALLMARKS) ; // increment and lock any unlocked mark
      SmMarkType eMarkType = sMarkLock.GetMarkType() ;

      for (ULONG ii=0; ii<sOffsetEdgeuses.GetSize(); ii++) 
        {
          SmEdgeuse *pEU = sOffsetEdgeuses[ii];
          pEU->Mark(eMarkType);

          // skip marked edgeuses
          SmEdgeuse *pRadial = pEU->GetRadial();
          if (pRadial->IsMarked(eMarkType)) { continue; }

          // use mate when its the offset edge
          ULONG lFound;
          if (    sOffsetEdgeuses.FindElement(pRadial->GetMate(),lFound)
              && !sOffsetEdgeuses.FindElement(pRadial,lFound)) 
            {
              // skip marked edgeuses
              pRadial = pRadial->GetMate();
              if (pRadial->IsMarked(eMarkType)) { continue; }
            }

          // add CapFaces to OffsetBrep to fill gaps between offset faces.
          // Gaps form when radial face pairs are offset into a convex (bigger than 180 sector)
          SE(BuildOffsetOfEdgeuses(pEU,pRadial,eMarkType)); // note: uses without increment eMarkType value

#ifdef SM_DEBUG_CODE
          // draw and dump input and offset breps after offsetting all edges and faces
          if (bDebugMe) 
            {
              SM_ASSERT_VALID(m_pOriginalBrep) ;
              SM_ASSERT_VALID(m_pOffsetBrep) ;

              smgfx_Erase();
              smgfx_SetLook(1,2, 0,0,1) ; m_pOriginalBrep->Draw(TRUE); sm_GraphicsLoop();
              smgfx_SetLook(1,2, 0,1,0) ; m_pOffsetBrep->Draw(TRUE); sm_GraphicsLoop();
              sm_GraphicsLoop();
            }
#endif // SM_DEBUG_CODE

        } // end iter every offset edgeuse
    } // end scope for SmNewMarkAndLock object

#ifdef SM_DEBUG_CODE
  // draw and dump input and offset breps after offsetting all edges and faces
  if (bDebugMe) 
    {
      m_pOriginalBrep->Dump() ;
      m_pOffsetBrep->Dump();

      smgfx_Erase();
      smgfx_SetLook(1,2, 0,0,1) ; m_pOriginalBrep->Draw(TRUE); sm_GraphicsLoop();
      smgfx_SetLook(1,2, 0,1,0) ; m_pOffsetBrep->Draw(TRUE); sm_GraphicsLoop();
      sm_GraphicsLoop();
    }
#endif // SM_DEBUG_CODE


  // For Lamina edges, build a ruled surface between the edge and its offset (if requested).
  if ( m_bCreateOffsetSolid )
    {
      for (ULONG ii=0; ii<sOffsetEdgeuses.GetSize(); ii++)
        {
          SmEdgeuse *pEU = sOffsetEdgeuses[ii];
          if ( pEU->GetEdge()->IsLamina() )
            { SE( BuildSideWallOfEdgeuse( pEU ) ); }
        }
    }

#ifdef SM_VALIDATE_TOPOLOGY
      m_pOffsetBrep->ValidatePointers();
#endif // SM_VALIDATE_TOPOLOGY

//    SER(m_pOffsetBrep->ShrinkGeometry());
#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe2 = FALSE;
  if (bDebugMe2) 
    {
      smgfx_Erase();
      m_pOffsetBrep->Draw();
      sm_GraphicsLoop();
    }
#endif // SM_DEBUG_CODE


  // Generate Vertex Cap Surfaces around corner vertices where face offsets have created holes

  // for every vertex
  SmTArray<SmVertex*> sVertices;
  m_pOriginalBrep->GetVertices(sVertices);
  for (ULONG ii=0; ii<sVertices.GetSize(); ii++)
    {
      SmSurface *pNewSurface;

      // Build a vertex offset CapSurface when needed
      //  - all edges connected to the vertex have an edge offset capSurface
      SmStatus eStat = (BuildOffsetOfVertex(sVertices[ii],s3DCurves,sOrients,pNewSurface));

      if ( eStat != SM_SUCCESS )
        { SER( eStat ); }

      if (!pNewSurface) { continue; }

#ifdef SM_DEBUG_CODE
      // draw and dump input and offset breps after offsetting all edges and faces
      if (bDebugMe) 
        {
          smgfx_Erase();
          smgfx_SetLook(1,2, 0,0,1) ; m_pOriginalBrep->Draw(TRUE); sm_GraphicsLoop();
          smgfx_SetLook(1,2, 0,1,0) ; m_pOffsetBrep->Draw(TRUE); sm_GraphicsLoop();
          smgfx_SetLook(1,2, 1,1,0) ; pNewSurface->DrawUV(6, 6) ; sm_GraphicsLoop() ;
          sm_GraphicsLoop();
        }
#endif // SM_DEBUG_CODE
      SmFace *pNewFace = NULL;
      if (s3DCurves.GetSize() > 0) 
        {
          sCurveLoops.ReSet();
          sCurveLoops.Add(s3DCurves.GetSize());
          SmRegion *pNewRegion = NULL;
          SmShell *pNewShell = NULL;

          // 
          SER(m_pOffsetBrep->MakeFaceWithCurves(m_pOffsetBrep->GetInfiniteRegion(), // in : region to contain new topology objects
                                                sCurveLoops,                        // in : 1 entry per loop, value = loop edge count, 1st entry=outer loop
                                                &s3DCurves,                         // in : opt ordered 3d trimming curves assigned to loops per sLoopEUCounts
                                                NULL,                               // in : opt ordered 2d trimming curves assigned to loops per sLoopEUCounts
                                                sOrients,                           // in : associated orients for each trimming curve, SM_OT_SAME or SM_OT_OPPOSITE
                                                sLoopPoints,                        // in : Point positions to build SmVertex VertexLoops
                                                pNewSurface,                        // in : new face->Surface
                                                pNewSurface->GetNaturalUVDomain(),  // in : domain of Surface used by face
                                                SM_OT_SAME,                         // in : Surface orient, oneof SM_OT_SAME or SM_OT_OPPOSITE
                                                pNewRegion,                         // out: New region if any. NULL when building trimmed surfaces, may be NotNULL for solids.
                                                pNewShell,                          // out: New shell if any.  Trimmed surfaces always create a new shell.
                                                pNewFace));                         // out: the new face
        }
      else //  
        {
          //
          SER(m_pOffsetBrep->CreateFaceFromSurface(pNewSurface,
              pNewSurface->GetNaturalUVDomain(),pNewFace));
        }

      // register vertex/newFace pair
      SER( m_vVToF.RelatePair( sVertices[ii], pNewFace ) );

#ifdef SM_DEBUG_CODE
      // draw originalBrep(blue), offsetBrep(green), newFace(cyan), newSurface(yellow)
      if (bDebugMe) 
        {
          m_pOriginalBrep->Dump() ;
          m_pOffsetBrep->Dump();

          smgfx_Erase();
          smgfx_SetLook(1,2, 0,0,1) ; m_pOriginalBrep->Draw(TRUE); sm_GraphicsLoop();
          smgfx_SetLook(1,2, 0,1,0) ; m_pOffsetBrep->Draw(TRUE); sm_GraphicsLoop();
          smgfx_SetLook(1,2, 0,1,1) ; pNewFace->Draw(SM_DM_CROSSHATCH) ; sm_GraphicsLoop() ;
          smgfx_SetLook(1,2, 1,1,0) ; pNewSurface->DrawUV(6, 6) ; sm_GraphicsLoop() ;
          sm_GraphicsLoop();
        }
#endif // SM_DEBUG_CODE
    } // end iter every originalBrep vertex

  // arrive here after copying all faces (offset or not) and building
  // end faces as needed to fill in the gaps

#ifdef SM_DEBUG_CODE
  // draw and dump input and offset breps after offsetting all vertices, edges and faces
  if (bDebugMe) 
    {
      SmTArray<SmFace*> sFacesDebug ; 
      m_pOffsetBrep->GetFaces( sFacesDebug ) ;

      m_pOriginalBrep->Dump() ;
      m_pOffsetBrep->Dump();

      smgfx_Erase();
      smgfx_SetLook(1,2, 0,0,1) ; m_pOriginalBrep->Draw(TRUE); sm_GraphicsLoop();
      smgfx_SetLook(1,2, 0,1,0) ; m_pOffsetBrep->Draw(TRUE); sm_GraphicsLoop();
      smgfx_SetLook(1,2, 0,1,1) ; for(ULONG kk=0;kk<sFacesDebug.GetSize();kk++)
                                    { SmFace *pFace = sFacesDebug[kk] ;
                                      if(pFace) { pFace->Draw( SM_DM_CROSSHATCH ); sm_GraphicsLoop(); }
                                    }
      sm_GraphicsLoop();
    }
#endif // SM_DEBUG_CODE

  // all done when not stitching and merging
  if (bSkipStitchMerge) 
    {
      rpResult = m_pOffsetBrep;
      //sCleanOffBrep.Clear();
      m_pOffsetBrep = NULL;

      return SM_SUCCESS;
    }

  // piecewise merge - create and populate pResult
  SmBrep *pResult =  NULL;
  // SER( PiecewiseMerge( pResult ) ); // increments unlocked mark values

  {
      SmStatus sErr = (PiecewiseMerge(pResult));
      if (sErr != SM_SUCCESS)
      {
          if (pResult)
          {
              delete pResult;
              pResult = NULL;
          }
          smos_ErrorMessage(sErr, FILE_NAME, LINE_NUMBER, (TCHAR*)0, (TCHAR*)0, FUNC_NAME);
          return (sErr);
      }
  }

  SmObjDelete sClean2( pResult );

  // arrive here when pResult contains all the offset and copied surfaces
  // merged (intersected and inserted) together. It won't be manifold
  // because all the overlapping pieces of the offset surfaces are still in pResult.

#ifdef SM_DEBUG_CODE
  // draw originalBrep(blue), offsetBrep(green), piecewiseMergeResult(yellow)
  if (bDebugMe) 
    {
      ULONG di ;
      m_pOriginalBrep->Dump() ; // input
      if ( m_pOffsetBrep ) { m_pOffsetBrep->Dump(); }    // intermediate structure
      pResult->Dump() ;         // the output - all faces merged - not yet manifold
      SmTArray<SmRegion *> sRegions ;
      pResult->GetRegions(sRegions) ; 

      SM_ASSERT_VALID(pResult) ;

      smgfx_Erase();
      smgfx_SetLook(1,2) ;        if(pResult) pResult->Draw(TRUE); sm_GraphicsLoop();
      smgfx_SetLook(1,2, 0,0,1) ; if(m_pOriginalBrep) m_pOriginalBrep->Draw(TRUE); sm_GraphicsLoop();
      smgfx_SetLook(1,2, 0,1,0) ; if(m_pOffsetBrep) m_pOffsetBrep->Draw(TRUE); sm_GraphicsLoop(); 
      for(di=0;di<sRegions.GetSize();di++)
        { smgfx_Erase() ;             if(sRegions[di]) sRegions[di]->DumpTopology() ; 
          smgfx_SetLook(1,2, 1,0,0) ; if(sRegions[di]) sRegions[di]->Draw() ; sm_GraphicsLoop() ;
          smgfx_SetLook(1,2) ;        if(pResult) pResult->Draw(TRUE); sm_GraphicsLoop() ;
        }
      sm_GraphicsLoop();
    }
#endif // SM_DEBUG_CODE

  // Do some stitching to close those unwanted gaps 
  //    this is a protective step that would not be needed in an ideal world.
  ULONG lStitchedEdges;
  double dVGap, dEGap;
  SmBoolean bProducesSolid;
  SmStatus sRtn = SmStitch::StitchIntoSolid(pResult,bProducesSolid,lStitchedEdges,dVGap,dEGap);
  SM_ASSERT(sRtn == SM_SUCCESS) ;

#ifdef SM_DEBUG_CODE
  // draw originalBrep(blue), offsetBrep(green), piecewiseMergeAndStitchResult(yellow)
  if (bDebugMe) 
    {
      m_pOriginalBrep->Dump() ; // input                                             
      if ( m_pOffsetBrep ) { m_pOffsetBrep->Dump(); } // intermediate structure                            
      pResult->Dump() ;         // the output - all faces merged and stitched so that faces 
                                //              near one another boundaries are now connected                                                     

      SM_ASSERT_VALID(pResult) ;

      smgfx_Erase();
      // draw input, intermediate, and output breps
      smgfx_SetLook(1,2) ;        pResult->Draw(TRUE); sm_GraphicsLoop();
      smgfx_SetLook(1,2, 0,0,1) ; m_pOriginalBrep->Draw(TRUE); sm_GraphicsLoop();
      if ( m_pOffsetBrep ) {
          smgfx_SetLook(1,2, 0,1,0) ; m_pOffsetBrep->Draw(TRUE); sm_GraphicsLoop();
      }

      if(bDrawFaces)
        {
          
          SmTArray<SmFaceuse *> sFaceuses ;
          SmTopologyTraverser   sTraverser ;
          SmBoolean bIsClosed ;
          
          SmTArray<SmFace *>    sFacesDraw; 
          pResult->GetFaces( sFacesDraw ) ;
          SmFaceuse *pFaceuse = sFacesDraw.GetSize() > 8 ? sFacesDraw[8]->GetUpwardFaceuse() : NULL ;
          if(pFaceuse) { SmNewMarkAndLock sMarkLock(pResult->GetContext(), SM_MT_ALLMARKS) ; // increments newly locked mark
                         sTraverser.FaceuseTraversal(pFaceuse,    // in : target faceuse (gets marked)                              
                                                     bIsClosed,   // out: TRUE=faceuse connects to closed faceuse set, FALSE=doesn't
                                                     &sFaceuses,  // out: List of all connected faces (get marked), NULL to ignore.
                                                     sMarkLock) ; // in : specify mark for target objects (not incremented)
                       }
          ULONG kk;
          smgfx_SetLook(1,2, 0,1,0) ; for(kk=0;kk<sFacesDraw.GetSize();kk++)
                                        {
                                          SmFace *pFaceA = sFacesDraw[kk] ;
                                          pFaceA->Draw(SM_DM_CROSSHATCH) ; sm_GraphicsLoop() ;
                                        }

          smgfx_SetLook(1,2, 1,0,0) ; if(pFaceuse) pFaceuse->Draw(SM_DM_CROSSHATCH) ; sm_GraphicsLoop() ;
          smgfx_SetLook(1,2, 0,1,1) ; for(kk=0;kk<sFaceuses.GetSize();kk++)
                                        {
                                          SmFaceuse *pFaceuseA = sFaceuses[kk] ;
                                          smgfx_Erase() ;
                                          smgfx_SetLook(1,2) ; pResult->Draw(TRUE); sm_GraphicsLoop();
                                          smgfx_SetLook(3,5, 0,1,0) ; pFaceuseA->GetFace()->Draw(SM_DM_CROSSHATCH) ; sm_GraphicsLoop() ;
                                          smgfx_SetLook(1,2, 0,1,0) ; pFaceuseA->Draw(SM_DM_FACEUSENEIGHBORS) ; sm_GraphicsLoop() ;
                                          pFaceuseA->GetFace()->Dump() ;
                                          sm_GraphicsLoop() ;
                                        }
        } // end DrawFaces check
      sm_GraphicsLoop();
    } // end if (bDebugMe)
#endif // SM_DEBUG_CODE

  // Here we will need to select regions to keep and
  // send in to make manifold.
  if (pResult->MakeManifold() == SM_SUCCESS) 
    {
      // label every region as inside or outside
      SmTArray<SmRegion*> sRegions;
      pResult->GetRegions(sRegions);
      for (ULONG kkk=0; kkk<sRegions.GetSize(); kkk++) 
        {
          SmRegion *pReg = sRegions[kkk];
#ifdef SM_DEBUG_CODE
          if (bDebugMe)
            {
              smgfx_Erase();
              pReg->Draw(); sm_GraphicsLoop();
              sm_GraphicsLoop();
            }
#endif //SM_DEBUG_CODE
          if (pReg == pResult->GetInfiniteRegion()) { pReg->SetIsVoid(TRUE); }
          else                                      { pReg->SetIsVoid(FALSE); }
        } // end iter every region labeling inside/outside status
    } // end MakeManifold a success check

#ifdef SM_DEBUG_CODE
  // draw originalBrep(blue), offsetBrep(green), piecewiseMergeResult(yellow)
  if (bDebugMe) 
    {
      m_pOriginalBrep->Dump() ; // input
      if ( m_pOffsetBrep )
        { m_pOffsetBrep->Dump(); }    // intermediate structure
      pResult->Dump() ;         // the output - all faces merged - not yet manifold

      SM_ASSERT_VALID(pResult) ;

      smgfx_Erase();
      smgfx_SetLook(1,2) ; pResult->Draw(TRUE); sm_GraphicsLoop();
      smgfx_SetLook(1,2, 0,0,1) ; m_pOriginalBrep->Draw(TRUE); sm_GraphicsLoop();
      if ( m_pOffsetBrep )
        { smgfx_SetLook(1,2, 0,1,0) ; m_pOffsetBrep->Draw(TRUE); sm_GraphicsLoop(); }
      sm_GraphicsLoop();
    }
#endif // SM_DEBUG_CODE

  // keep the result, set output and return
  sClean2.Clear();
  rpResult = pResult;

  // Do this by hand to make sure our pointer is zeroed out.
  //sCleanOffBrep.Clear();
  if ( m_pOffsetBrep != NULL ) { delete m_pOffsetBrep; m_pOffsetBrep = NULL; }

  return SM_SUCCESS;

} // end SmOffsetExecutive::DoSolidOffset

/*******************************************************************//**
PURPOSE: Static helper: order edgeuses by their Brep (faces, faceuses,
   loops, reversed), not by their heap addresses, so that the extended
   offset gives the same result in every process.

NOTES: Edgeuses not found in crBrep go last, in their input order.
***********************************************************************/
static void sm_OrderEdgeusesByBrep
  (const SmBrep         & crBrep,      // in : Brep that owns the edgeuses
   SmTArray<SmEdgeuse*> & rEdgeuses)   // i/o: edgeuses to order
{
  std::unordered_set<SmEdgeuse*> sUnplaced;
  for (ULONG ii=0; ii<rEdgeuses.GetSize(); ii++)
    { sUnplaced.insert(rEdgeuses[ii]); }

  SmTArray<SmEdgeuse*> sBrepOrder, sFaceuseEdgeuses;
  SmTArray<SmFace*> sFaces;
  crBrep.GetFaces(sFaces);
  for (ULONG ii=0; ii<sFaces.GetSize(); ii++)
    {
      SmFaceuse *pFaceuses[2] = { NULL, NULL };
      sFaces[ii]->GetFaceuses(pFaceuses[0], pFaceuses[1]);
      for (SmFaceuse *pFU : pFaceuses)
        {
          if (pFU == NULL) { continue; }
          pFU->GetEdgeuses(sFaceuseEdgeuses);
          sBrepOrder.Append(sFaceuseEdgeuses);
        }
    }

  SmTArray<SmEdgeuse*> sOrdered;
  for (ULONG ii=sBrepOrder.GetSize(); ii-- > 0; )
    {
      if (sUnplaced.erase(sBrepOrder[ii]) > 0) { sOrdered.Add(sBrepOrder[ii]); }
    }
  for (ULONG ii=0; ii<rEdgeuses.GetSize(); ii++)
    {
      if (sUnplaced.count(rEdgeuses[ii]) > 0) { sOrdered.Add(rEdgeuses[ii]); }
    }
  rEdgeuses = sOrdered;
} // end sm_OrderEdgeusesByBrep

/*******************************************************************//**
PURPOSE: This is the main top level method to offset a solid brep
    that extends convex edges instead of making rounded fillets.

NOTES: increments unlocked mark value
       also called by SmBrep::LocalOperation to replace one or more surfaces with others
***********************************************************************/
SmStatus SmOffsetExecutive::DoExtendedOffset
  (SmBrep *& rpResult)                        // out: pointer to new result
{
  NER(m_pOriginalBrep); // Must have one of these to get started

  // allocate a temporary empty new m_pOffsetBrep to contain intermediate results
  m_pOffsetBrep = new (m_crContext) SmBrep();
  // No: we delete this explicitly.  This can lead to a crash.  [B471 B473]
  //SmObjDelete sCleanOffBrep( m_pOffsetBrep );

  SM_OLDTOL_LINE m_pOffsetBrep->SetTolerance(m_rGeometryCreation.GetThisApproxTol3d());
  m_pOffsetBrep->m_bEditingEnabled = TRUE;

// Remove Composites
// //cbi_CEdge14:
// #ifdef SM_NO_COMPOSITES
//   m_pOffsetBrep->m_bMakeComposites = FALSE;
// #else
//   m_pOffsetBrep->m_bMakeComposites = TRUE;
// #endif // SM_NO_COMPOSITES

  // First step is to offset faces
  ULONG ii;
  SmTArray<SmFace*>           sFaces;
  SmSurface*                  sOffsetSurface = NULL;
  SmSSIData                   sIntData;
  SmTArray<SmCurve*>          s3DCurves;
  SmTArray<SmBSplineCurve*>   sUVCurves1;
  SmTArray<SmOrientType>      sOrients;
  SmTArray<SmPoint3d>         sLoopPoints;
  SmTArray<ULONG>             sCurveLoops;
  m_pOriginalBrep->GetFaces(sFaces);

#ifdef SM_DEBUG_CODE
ULONG di ;
// ULONG GWC_SET_NEXT_TWO_LINES_TO_FALSE_BEFORE_RELEASE ;
SmBoolean bDebugMe    = FALSE;
SmBoolean bDebugMeAll = FALSE ; 
  if (bDebugMe || bDebugMeAll) 
    {
      SM_PTR_ARRAY(sShells,   SmShell,   8) ; 
      SM_PTR_ARRAY(sEdges,    SmEdge,   32) ;
      SM_PTR_ARRAY(sVertices, SmVertex, 32) ;
      SM_DUMP_AND_ASSERT_VALID(m_pOriginalBrep) ;
      SM_DUMP_AND_ASSERT_VALID(m_pOffsetBrep) ;
      m_pOriginalBrep->GetEdges(sEdges) ;
      m_pOriginalBrep->GetVertices(sVertices) ;
      m_pOriginalBrep->GetShells(sShells) ;

      smgfx_Erase();
      smgfx_SetLook(1,2, 0,0,1) ; m_pOriginalBrep->Draw(TRUE); sm_GraphicsLoop();
      smgfx_SetLook(3,4, 1,0,1) ; for(di=0;di<sShells.GetSize();di++) { SM_DUMP_AND_ASSERT_VALID(sShells[di]) ;
                                                                        if(sShells[di]) { sShells[di]->Draw(); sm_GraphicsLoop(); }
                                                                        sm_GraphicsLoop();
                                                                      }
      smgfx_SetLook(1,2, 1,0,0) ; for(di=0;di<sFaces.GetSize();di++) { SM_DUMP_AND_ASSERT_VALID(sFaces[di]) ;
                                                                       if(sFaces[di]) { smgfx_SetLook(1,2, 1,0,0) ; sFaces[di]->DrawUV(); sm_GraphicsLoop();
                                                                                        smgfx_SetLook(3,4, 1,0,1) ; sFaces[di]->DrawLoopusesForFaceuse() ; sm_GraphicsLoop() ;
                                                                                        sm_GraphicsLoop();
                                                                                      }
                                                                     }
      smgfx_SetLook(3,4, 0,0,1) ; for(di=0;di<sEdges.GetSize();di++) { SM_DUMP_AND_ASSERT_VALID(sEdges[di]) ;
                                                                       if(sEdges[di]) { sEdges[di]->Draw(); sm_GraphicsLoop(); }
                                                                       sm_GraphicsLoop();
                                                                     }
      smgfx_SetLook(4,5, 1,0,0) ; for(di=0;di<sVertices.GetSize();di++) { SM_DUMP_AND_ASSERT_VALID(sVertices[di]) ;
                                                                          if(sVertices[di]) { sVertices[di]->Draw(); sm_GraphicsLoop(); }
                                                                          sm_GraphicsLoop();
                                                                        }
      sm_GraphicsLoop();
    }
#endif // SM_DEBUG_CODE

  // Offset each Face, unless it's in m_pShellFaces, then just copy it.
  for (ii=0; ii<sFaces.GetSize(); ii++) 
    {
      SmFace *pOrigFace = sFaces[ii];
      ULONG lIndex;

      // when face is one of the given faces to shell
      if (   m_pShellFaces 
          && m_pShellFaces->FindElement(pOrigFace,lIndex)) 
        {
          // Just copy this face instead of offsetting it

          // make a copy of face->Edge->Curves and face->vertexLoop->Points 
          SER(pOrigFace->CreateCurvesFromFace(m_crContext,   // in : context for new object construction
                                              SM_OT_SAME,    // in : use Faceuse whose orientation is the same as eOrientation
                                              sCurveLoops,   // out: number of edges in each edge loop
                                              &s3DCurves,    // out: optional copies of all 3D loop curves - note: one or both of pOpt3DCurves and
                                              &sUVCurves1,   // out: optional copies of all 2D loop curves -         pOptUVCurves must be nonNULL
                                              sOrients,      // out: Curve orientation in relation to its loop for each curve
                                              sLoopPoints)); // out: point position of each of face's vertex loops      
          // locals
          SmRegion  *pNewRegion = NULL;
          SmShell   *pNewShell = NULL;
          SmFace    *pNewFace = NULL;
          SmSurface *pOrigSurface = pOrigFace->GetSurface(); NER(pOrigSurface);
          SmSurface *pSurface = NULL;
          
          // Copy pOrigSurface for new face, when possible as an analytic surface
          SER(pOrigSurface->CopyAndAddAnalytics(m_crContext,pSurface));

          // make a new face which is a copy of the current face in the output Brep 
          SER(m_pOffsetBrep->MakeFaceWithCurves(m_pOffsetBrep->GetInfiniteRegion(),  // in : region to contain new topology objects
                                                sCurveLoops,                         // in : 1 entry per loop, value = loop edge count, 1st entry=outer loop
                                                &s3DCurves,                          // in : opt ordered 3d trimming curves assigned to loops per sLoopEUCounts
                                                &sUVCurves1,                         // in : opt ordered 2d trimming curves assigned to loops per sLoopEUCounts
                                                sOrients,                            // in : associated orients for each trimming curve, SM_OT_SAME or SM_OT_OPPOSITE
                                                sLoopPoints,                         // in : Point positions to build SmVertex VertexLoops
                                                pSurface,                            // in : new face->Surface
                                                pOrigFace->GetUVDomain(),            // in : domain of Surface used by face
                                                SM_OT_SAME,                          // in : Surface orient, oneof SM_OT_SAME or SM_OT_OPPOSITE
                                                pNewRegion,                          // out: New region if any. NULL when building trimmed surfaces, may be NotNULL for solids.
                                                pNewShell,                           // out: New shell if any.  Trimmed surfaces always create a new shell.
                                                pNewFace));                          // out: the new face

#ifdef SM_DEBUG_CODE
          if (bDebugMe || bDebugMeAll) 
            {
              smgfx_Erase();
              smgfx_SetLook(1,2, 0,0,1) ; m_pOriginalBrep->Draw(TRUE); sm_GraphicsLoop();
              smgfx_SetLook(1,2, 1,0,0) ; pOrigFace->DrawUV(); sm_GraphicsLoop();
              smgfx_SetLook(4,5, 1,0,0) ; pOrigSurface->DrawParams(); sm_GraphicsLoop();
              smgfx_SetLook(3,4, 1,0,1) ; pOrigFace->DrawLoopusesForFaceuse() ; sm_GraphicsLoop() ;
              smgfx_SetLook(1,2, 1,0,0) ; pNewFace->DrawUV(); sm_GraphicsLoop();
              smgfx_SetLook(4,5, 1,0,0) ; pSurface->DrawParams(); sm_GraphicsLoop();
              smgfx_SetLook(5,6, 0,1,1) ; pNewFace->DrawLoopusesForFaceuse() ; sm_GraphicsLoop() ;
              sm_GraphicsLoop();
            }

#endif // SM_DEBUG_CODE
          pOrigFace->Notify(SM_NO_COPY, pNewFace, SM_NO_GET_BREP(pNewFace), SM_NO_GET_BREP(pOrigFace)) ; // Basically copied pNewFace from pOrigFace.

          // select original Face->Faceuse to offset 
          //  m_bInset == FALSE - get faceuse bounding infinite region
          //  m_bInset == TRUE  - get faceuse not bounding infinite region
          SmFaceuse *pFU = pOrigFace->GetUpwardFaceuse();
          if (pFU->GetShell()->GetRegion() != m_pOriginalBrep->GetInfiniteRegion()) 
            { pFU = pFU->GetMate(); }
          if (m_bInset)
            { pFU = pFU->GetMate(); }

          // relate this faceuse to the newface
          SER( m_vFUToF.RelatePair(pFU,pNewFace) );
          SER(MapEdgeuses(pFU,pNewFace));

#ifdef SM_DEBUG_CODE
          if (bDebugMe || bDebugMeAll) 
            {
              SM_PTR_ARRAY(sOrigShells,   SmShell,  32) ; SM_PTR_ARRAY(sOffsetShells,   SmShell,  32) ;
              SM_PTR_ARRAY(sOffsetFaces,  SmFace,   32) ;
              SM_PTR_ARRAY(sOrigEdges,    SmEdge,   32) ; SM_PTR_ARRAY(sOffsetEdges,    SmEdge,   32) ;
              SM_PTR_ARRAY(sOrigVertices, SmVertex, 32) ; SM_PTR_ARRAY(sOffsetVertices, SmVertex, 32) ;
              m_pOriginalBrep->GetShells(sOrigShells) ; 
              m_pOriginalBrep->GetEdges(sOrigEdges) ;
              m_pOriginalBrep->GetVertices(sOrigVertices) ;

              m_pOffsetBrep->GetShells(sOffsetShells) ; 
              m_pOffsetBrep->GetFaces(sOffsetFaces) ;
              m_pOffsetBrep->GetEdges(sOffsetEdges) ;
              m_pOffsetBrep->GetVertices(sOffsetVertices) ;

              SM_DUMP_AND_ASSERT_VALID(m_pOriginalBrep) ;
              smgfx_Erase();
              smgfx_SetLook(1,2, 0,0,1) ; m_pOriginalBrep->Draw(TRUE); sm_GraphicsLoop();
              for(di=0;di<sOrigShells.GetSize();di++) { SM_DUMP_AND_ASSERT_VALID(sOrigShells[di]) ;
                                                        if(sOrigShells[di]) { smgfx_SetLook(3,4, 1,0,1) ; sOrigShells[di]->Draw() ; sm_GraphicsLoop() ;
                                                                              sm_GraphicsLoop() ;
                                                                            }
                                                      }
              for(di=0;di<sFaces.GetSize();di++) { SM_DUMP_AND_ASSERT_VALID(sFaces[di]) ;
                                                   if(sFaces[di]) { smgfx_SetLook(1,2, 1,0,0) ; sFaces[di]->DrawUV(); sm_GraphicsLoop();
                                                                    smgfx_SetLook(3,4, 1,0,1) ; sFaces[di]->DrawLoopusesForFaceuse() ; sm_GraphicsLoop() ;
                                                                    sm_GraphicsLoop() ;
                                                                  }
                                                 }
              smgfx_SetLook(3,4, 0,0,1) ; for(di=0;di<sOrigEdges.GetSize();di++) { SM_DUMP_AND_ASSERT_VALID(sOrigEdges[di]) ;
                                                                                   if(sOrigEdges[di]) { sOrigEdges[di]->Draw(); sm_GraphicsLoop(); }
                                                                                   sm_GraphicsLoop() ;
                                                                                 }
              smgfx_SetLook(4,5, 1,0,0) ; for(di=0;di<sOrigVertices.GetSize();di++) { SM_DUMP_AND_ASSERT_VALID(sOrigVertices[di]) ;
                                                                                      if(sOrigVertices[di]) { sOrigVertices[di]->Draw(); sm_GraphicsLoop(); }
                                                                                      sm_GraphicsLoop() ;
                                                                                    }
              sm_GraphicsLoop();

              SM_DUMP_AND_ASSERT_VALID(m_pOffsetBrep) ;
              smgfx_SetLook(1,2, 0,1,0) ; m_pOffsetBrep->Draw(TRUE); sm_GraphicsLoop();
              for(di=0;di<sOffsetShells.GetSize();di++) { SM_DUMP_AND_ASSERT_VALID(sOffsetShells[di]) ;
                                                          if(sOffsetShells[di]) { smgfx_SetLook(3,4, 1,0,1) ; sOffsetShells[di]->Draw() ; sm_GraphicsLoop() ;
                                                                                  sm_GraphicsLoop() ;
                                                                                }
                                                        }
              for(di=0;di<sOffsetFaces.GetSize();di++) { SM_DUMP_AND_ASSERT_VALID(sOffsetFaces[di]) ;
                                                         if(sOffsetFaces[di]) { smgfx_SetLook(1,2, 1,0,0) ; sOffsetFaces[di]->DrawUV(); sm_GraphicsLoop();
                                                                                smgfx_SetLook(3,4, 0,1,1) ; sOffsetFaces[di]->DrawLoopusesForFaceuse() ; sm_GraphicsLoop() ;
                                                                                sm_GraphicsLoop() ;
                                                                              }
                                                       }
              smgfx_SetLook(3,4, 0,0,1) ; for(di=0;di<sOffsetEdges.GetSize();di++) { SM_DUMP_AND_ASSERT_VALID(sOffsetEdges[di]) ;
                                                                                     if(sOffsetEdges[di]) { sOffsetEdges[di]->Draw(); sm_GraphicsLoop(); }
                                                                                     sm_GraphicsLoop() ;
                                                                                   }
              smgfx_SetLook(4,5, 1,0,0) ; for(di=0;di<sOffsetVertices.GetSize();di++) { SM_DUMP_AND_ASSERT_VALID(sOffsetVertices[di]) ;
                                                                                        if(sOffsetVertices[di]) { sOffsetVertices[di]->Draw(); sm_GraphicsLoop(); }
                                                                                        sm_GraphicsLoop() ;
                                                                                      }
              sm_GraphicsLoop();
            }
#endif // SM_DEBUG_CODE

        } // end face is one of the given faces to shell branch

      else // face is not one of the given faces to shell branch: offset it.
        {
          // select original Face->Faceuse to offset 
          //  m_bInset == FALSE - get faceuse pointing into infinite region
          //  m_bInset == TRUE  - get faceuse not pointing into infinite region
          SmFaceuse *pFU = pOrigFace->GetUpwardFaceuse();
          if (!pFU->GetShell()->GetRegion()->IsVoid()) 
            { pFU = pFU->GetMate(); }

          if (m_bInset)
            { pFU = pFU->GetMate(); }

          // create the offset surface for this faceuse - OutSurface(s) domains = InSurface domain, but may be trimmed for singularites.
          SER(m_rGeometryCreation.FaceuseOffset(m_crContext,        // in : context for new obj construction
                                                pFU,                // in : Faceuse->Surface to copy and offset
                                                TRUE,               // in : TRUE =
                                                                    //      FALSE=
                                                sOffsetSurface,     // out: New Surface (more than 1 for self-XSect cases)
                                                sIntData)) ;        // out: List of SelfXSect curves (when bSkipSelfIntersection == FALSE)
                                                                  
          // An inset face whose offset collapses (a convex round thinner than the wall) gets no
          // offset face: the faces beside it are mated across it below, as at a sharp edge.
          if (sOffsetSurface == NULL) 
            {
              if (m_bInset && m_bExtendConvexEdges)
                {
                  m_sCollapsedFaces.Add(pOrigFace) ;
                  continue ;
                }
              SER(SM_ERR);
            }

          // signal self intersection problems
          // But -- we passed True to FaceuseOffset, telling it not to check this...
          ULONG lNumSelfInt = sIntData.m_v3DCurves.GetSize();
          if (   (m_bMergeResults      && lNumSelfInt > 0)  
              || (m_bExtendConvexEdges && lNumSelfInt > 5)) 
            {
              // Unexpected problems found in the offset
              SER(SM_ERR);
            }

          // Build an offset of pFU into m_pOffsetBrep and 
          //  1. relate pFU to newFace in m_vFUToF
          //  2. relate pFU->edges to newFace->edges in m_vEUToEU
          SER(BuildOffsetOfFaceuse(pFU,sOffsetSurface,sIntData));

#ifdef SM_DEBUG_CODE
          if (bDebugMe || bDebugMeAll) 
            {
              SM_PTR_ARRAY(sOrigShells,   SmShell,  32) ; SM_PTR_ARRAY(sOffsetShells,   SmShell,  32) ;
              SM_PTR_ARRAY(sOffsetFaces,  SmFace,   32) ;
              SM_PTR_ARRAY(sOrigEdges,    SmEdge,   32) ; SM_PTR_ARRAY(sOffsetEdges,    SmEdge,   32) ;
              SM_PTR_ARRAY(sOrigVertices, SmVertex, 32) ; SM_PTR_ARRAY(sOffsetVertices, SmVertex, 32) ;

              m_pOriginalBrep->GetShells  (sOrigShells) ;
              m_pOriginalBrep->GetEdges   (sOrigEdges) ;
              m_pOriginalBrep->GetVertices(sOrigVertices) ;

              m_pOffsetBrep->GetShells  (sOffsetShells) ;
              m_pOffsetBrep->GetFaces   (sOffsetFaces) ;
              m_pOffsetBrep->GetEdges   (sOffsetEdges) ;
              m_pOffsetBrep->GetVertices(sOffsetVertices) ;

              SM_DUMP_AND_ASSERT_VALID(m_pOriginalBrep) ;
              smgfx_Erase();
              smgfx_SetLook(1,2, 0,0,1) ; m_pOriginalBrep->Draw(TRUE); sm_GraphicsLoop();
              for(di=0;di<sOrigShells.GetSize();di++) { SM_DUMP_AND_ASSERT_VALID(sOrigShells[di]) ;
                                                        if(sOrigShells[di]) { smgfx_SetLook(3,4, 1,0,1) ; sOrigShells[di]->Draw() ; sm_GraphicsLoop() ;
                                                                              sm_GraphicsLoop() ;
                                                                            }
                                                      }
              for(di=0;di<sFaces.GetSize();di++) { SM_DUMP_AND_ASSERT_VALID(sFaces[di]) ;
                                                   if(sFaces[di]) { smgfx_SetLook(1,2, 1,0,0) ; sFaces[di]->DrawUV(); sm_GraphicsLoop();
                                                                    smgfx_SetLook(3,4, 1,0,1) ; sFaces[di]->DrawLoopusesForFaceuse() ; sm_GraphicsLoop() ;
                                                                    sm_GraphicsLoop() ;
                                                                  }
                                                 }
              smgfx_SetLook(3,4, 0,0,1) ; for(di=0;di<sOrigEdges.GetSize();di++) { SM_DUMP_AND_ASSERT_VALID(sOrigEdges[di]) ;
                                                                                   if(sOrigEdges[di]) { sOrigEdges[di]->Draw(); sm_GraphicsLoop(); }
                                                                                   sm_GraphicsLoop() ;
                                                                                 }
              smgfx_SetLook(4,5, 1,0,0) ; for(di=0;di<sOrigVertices.GetSize();di++) { SM_DUMP_AND_ASSERT_VALID(sOrigVertices[di]) ;
                                                                                      if(sOrigVertices[di]) { sOrigVertices[di]->Draw(); sm_GraphicsLoop(); }
                                                                                      sm_GraphicsLoop() ;
                                                                                    }
              sm_GraphicsLoop();

              SM_DUMP_AND_ASSERT_VALID(m_pOffsetBrep) ;
              smgfx_SetLook(1,2, 0,1,0) ; m_pOffsetBrep->Draw(TRUE); sm_GraphicsLoop();
              for(di=0;di<sOffsetShells.GetSize();di++) { SM_DUMP_AND_ASSERT_VALID(sOffsetShells[di]) ;
                                                          if(sOffsetShells[di]) { smgfx_SetLook(3,4, 1,0,1) ; sOffsetShells[di]->Draw() ; sm_GraphicsLoop() ;
                                                                                  sm_GraphicsLoop() ;
                                                                                }
                                                        }
              for(di=0;di<sOffsetFaces.GetSize();di++) { SM_DUMP_AND_ASSERT_VALID(sOffsetFaces[di]) ;
                                                         if(sOffsetFaces[di]) { smgfx_SetLook(1,2, 1,0,0) ; sOffsetFaces[di]->DrawUV(); sm_GraphicsLoop();
                                                                                smgfx_SetLook(3,4, 0,1,1) ; sOffsetFaces[di]->DrawLoopusesForFaceuse() ; sm_GraphicsLoop() ;
                                                                                sm_GraphicsLoop() ;
                                                                              }
                                                       }
              smgfx_SetLook(3,4, 0,0,1) ; for(di=0;di<sOffsetEdges.GetSize();di++) { SM_DUMP_AND_ASSERT_VALID(sOffsetEdges[di]) ;
                                                                                     if(sOffsetEdges[di]) { sOffsetEdges[di]->Draw(); sm_GraphicsLoop(); }
                                                                                     sm_GraphicsLoop() ;
                                                                                   }
              smgfx_SetLook(4,5, 1,0,0) ; for(di=0;di<sOffsetVertices.GetSize();di++) { SM_DUMP_AND_ASSERT_VALID(sOffsetVertices[di]) ;
                                                                                        if(sOffsetVertices[di]) { sOffsetVertices[di]->Draw(); sm_GraphicsLoop(); }
                                                                                        sm_GraphicsLoop() ;
                                                                                      }
              sm_GraphicsLoop();
            }
#endif // SM_DEBUG_CODE

        } // end face is not one of the given faces to shell branch
    } // end iter every face

  // Now that we have all faces offset.  Let's do some work to replace fillet
  // and blending offsets which are implicit now with a real surface.
//
// Not implemented, so don't bother calling it.
//  for (ULONG ii=0; ii<sFaces.GetSize(); ii++) 
//    {
//      SmFace *pOrigFace = sFaces[ii];
//
//      // skip faces on the Shell Only List
//      ULONG lIndex;
//      if (   m_pShellFaces
//          && m_pShellFaces->FindElement(pOrigFace,lIndex)) 
//        {
//          continue;
//        }
//
//      // when m_bInset == FALSE : pFU = Faceuse bounding infinite region
//      //      m_bInset == TRUE  : pFU = faceuse bounding internal region
//      SmFaceuse *pFU = pOrigFace->GetUpwardFaceuse();
//      if (!pFU->GetShell()->GetRegion()->IsVoid()) {
//          pFU = pFU->GetMate();
//      }
//      if (m_bInset) pFU = pFU->GetMate();
//
//      // not yet implemented 
//      // plan    - Build tolerance free new fillet surfaces
//      // current - no effects
//      SER(BuildFilletOffsetsOfFaceuse(pFU));
//
//    } // end iter every face

  // arrive here after all original faces have been copied (with or without offset)

  SmTArray<SmEdgeuse*> sOffsetEdgeuses;
  SmTArray<SmEdgeuse*> sOffSetEUS, sRadialEUS;
  SmTArray<SmEdge*>    sNewEdges;
  m_vEUToEU.GetAllFirsts(sOffsetEdgeuses);
  sm_OrderEdgeusesByBrep(*m_pOriginalBrep, sOffsetEdgeuses);

#ifdef SM_VALIDATE_TOPOLOGY
      m_pOffsetBrep->ValidatePointers();
#endif // SM_VALIDATE_TOPOLOGY

#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe11 = FALSE;
  if (bDebugMe11 || bDebugMeAll) 
    {
      m_pOriginalBrep->Dump() ; SM_ASSERT_VALID(m_pOriginalBrep) ;
      m_pOffsetBrep->Dump() ;   SM_ASSERT_VALID(m_pOffsetBrep) ;

      smgfx_Erase();
      smgfx_SetLook(1,2, 0,0,1) ; m_pOriginalBrep->Draw(TRUE); sm_GraphicsLoop() ;
      smgfx_SetLook(1,2, 0,1,1) ; m_pOffsetBrep->Draw(TRUE); sm_GraphicsLoop() ;
      sm_GraphicsLoop();
    }
#endif // SM_DEBUG_CODE

  // increment and lock an unlocked mark
  SmNewMarkAndLock sMarkLock( m_pOffsetBrep->GetContext(), SM_MT_ALLMARKS) ; // increment and lock any unlocked mark
  SmMarkType       eMarkType = sMarkLock.GetMarkType() ;

  // Mate the faces beside each collapsed face across it: the offsets of the faces on either side
  // of a collapsed round meet where they intersect, the round's rails becoming one sharp edge.
  // (The rails' ends are trimmed by TrimExtendedEdge, which also looks across collapsed faces.)
  ULONG lUnbridged = 0 ;
  double dCollapseTanTolDeg = GetTangencyTolDegrees();
  for (ii=0; ii<sOffsetEdgeuses.GetSize(); ii++)
    {
      SmEdgeuse *pEU     = sOffsetEdgeuses[ii];
      SmEdgeuse *pRadial = pEU->GetRadial();
      if (   pRadial == NULL
          || pRadial == pEU
          || !IsCollapsedFace(pRadial->GetFace())
          || !pRadial->IsTangentSector(dCollapseTanTolDeg)
          || pEU->IsMarked(eMarkType))
        { continue; }
      SmEdgeuse *pAcross    = GetEdgeuseAcrossCollapsedFace(pRadial, NULL);
      SmEdgeuse *pOffEU     = GetOffsetEdgeuse(pEU);
      SmEdgeuse *pOffAcross = pAcross ? GetOffsetEdgeuse(pAcross) : NULL;
      if (   pOffEU == NULL
          || pOffAcross == NULL
          || m_rGeometryCreation.MateEdges(m_crContext, this, pEU, pAcross, pOffEU, pOffAcross, eMarkType) != SM_SUCCESS)
        {
          lUnbridged++ ;
          continue;
        }
      pEU->Mark(eMarkType);     pEU->GetMate()->Mark(eMarkType);
      pAcross->Mark(eMarkType); pAcross->GetMate()->Mark(eMarkType);
    } // end iter every offset edgeuse beside a collapsed face

  // when extending convex edges - not filleting them 
  // mark offset edges which have tangent sectors
  if (GetExtendConvexEdges()) 
    {
      // for every edgeuse that has been offset
      for (ii=0; ii<sOffsetEdgeuses.GetSize(); ii++) 
        {
          SmEdgeuse *pEU     = sOffsetEdgeuses[ii];
          SmEdgeuse *pRadial = pEU->GetRadial();
          if (IsCollapsedFace(pRadial->GetFace())) { continue; } // mated across it above

          // Mark the EU, and skip it if its radial partner is marked.
          if (pRadial->IsMarked(eMarkType)) { continue; }

          // If Radial is not in the OffsetEdgeuses list but Radial's Mate is,
          // then swap Radial's Mate in for Radial, and skip if that is marked.
          ULONG lFound;
          if (   sOffsetEdgeuses.FindElement(pRadial->GetMate(),lFound)
              && !sOffsetEdgeuses.FindElement(pRadial,lFound)) 
            {
              pRadial = pRadial->GetMate();
              if (pRadial->IsMarked(eMarkType)) { continue; }
            }

#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe6 = FALSE;
static ULONG lCount      = 1 ; lCount++ ;
static ULONG lDebugCount = 0 ;
          if (bDebugMe6 || lCount == lDebugCount || bDebugMeAll) 
            {
              SM_DUMP_AND_ASSERT_VALID(m_pOriginalBrep) ;
              SM_DUMP_AND_ASSERT_VALID(m_pOffsetBrep) ;

              smgfx_Erase();
              smgfx_SetLook(3,5, 1,0,0); pEU->Draw(); sm_GraphicsLoop();
              smgfx_SetLook(1,2, 0,0,0); m_pOffsetBrep->Draw(); sm_GraphicsLoop();
              sm_GraphicsLoop();
            }
#endif // SM_DEBUG_CODE

          // when edgeuse is a tangent sector
          double dTanTolDeg = GetTangencyTolDegrees();
          if (pEU->IsTangentSector(dTanTolDeg)) 
            {
              // 
              SER(BuildOffsetOfEdgeuses(pEU,pRadial,eMarkType)); // note: uses without increment eMarkType value
              pRadial->Mark(eMarkType); 
              pRadial->GetMate()->Mark(eMarkType);
              pEU->Mark(eMarkType); 
              pEU->GetMate()->Mark(eMarkType);
#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe5 = FALSE;
              if (bDebugMe5 || bDebugMeAll) 
                {
                  SM_DUMP_AND_ASSERT_VALID(m_pOriginalBrep) ;
                  SM_DUMP_AND_ASSERT_VALID(m_pOffsetBrep) ;

                  smgfx_Erase();
                  smgfx_SetLook(3,5, 1,0,0); pEU->Draw(); sm_GraphicsLoop();
                  smgfx_SetLook(1,2, 0,0,0); m_pOffsetBrep->Draw(); sm_GraphicsLoop();
                  sm_GraphicsLoop();
                }
#endif // SM_DEBUG_CODE
            }
#ifdef SM_VALIDATE_TOPOLOGY
      m_pOffsetBrep->ValidatePointers();
#endif // SM_VALIDATE_TOPOLOGY
         }
    } // end extending convex edges check

#ifdef SM_VALIDATE_TOPOLOGY
      m_pOffsetBrep->ValidatePointers();
#endif // SM_VALIDATE_TOPOLOGY

  // Generate edge offsets 
  for(ii=0;ii<sOffsetEdgeuses.GetSize();ii++) 
    {
      // Mark the EU, and skip it if its radial partner is marked.
      SmEdgeuse * pEU     = sOffsetEdgeuses[ii];
      SmEdgeuse * pRadial = pEU->GetRadial();
      if (IsCollapsedFace(pRadial->GetFace())) { continue; } // mated across it above
      pEU->Mark(eMarkType);
      if (pRadial->IsMarked(eMarkType)) 
        { continue; }

      // If non-lamina (EU has a mate),
      // and Radial is not in the OffsetEdgeuses list but Radial's Mate is,
      // then swap Radial's Mate in for Radial, and skip if that is marked.
      ULONG lFound;
      if (!pEU->GetEdge()->IsLamina()) 
        {
          if (   sOffsetEdgeuses.FindElement(pRadial->GetMate(),lFound)
              && !sOffsetEdgeuses.FindElement(pRadial,lFound)) 
            {
              pRadial = pRadial->GetMate();
              if (pRadial->IsMarked(eMarkType)) 
                { continue; }
            }
        }

#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe7 = FALSE;
static ULONG lCount      = 1 ; lCount++ ;
static ULONG lDebugCount = 0 ;
      if (bDebugMe7 || lCount == lDebugCount || bDebugMeAll) 
        {
          m_pOffsetBrep->Dump() ;   SM_ASSERT_VALID(m_pOffsetBrep) ;

          smgfx_Erase();
          smgfx_SetLook(3,4, 0,0,1); pEU->Draw(); sm_GraphicsLoop();
          smgfx_SetLook(3,4, 0,1,0); pRadial->Draw(); sm_GraphicsLoop();
          smgfx_SetLook(1,2); m_pOffsetBrep->Draw(TRUE); sm_GraphicsLoop();
          sm_GraphicsLoop();
        }
#endif // SM_DEBUG_CODE

      // Build cap faces to fill in gaps between offset faces
      SER(BuildOffsetOfEdgeuses(pEU,          // in : target edgeuse that has been offset
                                pRadial,      // in : radial mate to target edgeuse that has been offset
                                eMarkType)) ; // in : uses without increment eMarkType value
                                              // note: uses without increment eMarkType value

#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe8 = FALSE;
      if (bDebugMe8 || bDebugMeAll) 
        {
          m_pOffsetBrep->Dump();
          // Note: This Brep can have huge gaps here: we moved an Edge.
          SM_ASSERT_VALID(m_pOffsetBrep) ;

          smgfx_Erase();
          smgfx_SetLook(3,4, 1,0,0); pEU->Draw(); sm_GraphicsLoop();
          smgfx_SetLook(1,2); m_pOffsetBrep->Draw(TRUE); sm_GraphicsLoop();
          sm_GraphicsLoop();
        }
#endif // SM_DEBUG_CODE
    } // end iter every offset edgeuse


  // A local operation has no merge to close what a step wall leaves open at its ends.
  if (m_eOffsetOperation == SM_OO_LOCAL_OPERATION)
    { SER(CloseLoopGaps()); }

  // The offset of an edge along a collapsed face's end (where the round meets the face across its
  // end) now runs from the corner where the mated faces meet back to that corner: squeeze it out.
  for (ii=0; ii<sOffsetEdgeuses.GetSize(); ii++)
    {
      SmEdgeuse *pEU     = sOffsetEdgeuses[ii];
      SmEdgeuse *pRadial = pEU->GetRadial();
      if (   pRadial == NULL
          || pRadial == pEU
          || !IsCollapsedFace(pRadial->GetFace())
          || pRadial->IsTangentSector(dCollapseTanTolDeg))
        { continue; }
      SmEdgeuse *pOffEU = GetOffsetEdgeuse(pEU);
      SmEdge    *pOffEdge = pOffEU ? pOffEU->GetEdge() : NULL;
      SmVertex  *pV1 = pOffEdge ? pOffEdge->GetStartVertex() : NULL;
      SmVertex  *pV2 = pOffEdge ? pOffEdge->GetEndVertex()   : NULL;
      if (   pOffEdge == NULL
          || pV1 == pV2
          || pV1->GetPoint().DistanceBetween(pV2->GetPoint()) > 10.0 * pOffEdge->GetTolerance())
        {
          lUnbridged++ ;
          continue;
        }
      SER(m_vEUToEU.RemovePair(pEU, pOffEU));
      SER(m_pOffsetBrep->SqueezeEdge(pOffEdge, pV1));
    } // end iter every offset edgeuse along a collapsed face's end

  // A collapsed face that could not be bridged (no face across it, or no meeting) leaves the
  // offset open: refuse it rather than build a body with holes.
  if (lUnbridged > 0)
    { SER(SM_ERR_DEGENERATE_SURFACE); }

  SER(RemoveExcessFaces( m_pOffsetBrep ));

#ifdef SM_VALIDATE_TOPOLOGY
      m_pOffsetBrep->ValidatePointers();
#endif // SM_VALIDATE_TOPOLOGY

  // For Lamina edges, build a ruled surface between the edge and its offset.
  if (m_bCreateOffsetSolid)
    {
      for (ULONG kk=0; kk<sOffsetEdgeuses.GetSize(); kk++) 
        {
          SmEdgeuse *pEU = sOffsetEdgeuses[kk];
          if (pEU->GetEdge()->IsLamina()) 
            { SE(BuildSideWallOfEdgeuse(pEU)); }
        }
    }

#ifdef SM_VALIDATE_TOPOLOGY
      m_pOffsetBrep->ValidatePointers();
#endif // SM_VALIDATE_TOPOLOGY

//    SER(m_pOffsetBrep->ShrinkGeometry());
#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe2 = FALSE;
  if (bDebugMe2 || bDebugMeAll) 
    {
      m_pOffsetBrep->Dump();
      smgfx_Erase();
      smgfx_SetLook( 1,2, 0,0,1 ); m_pOffsetBrep->Draw(); sm_GraphicsLoop();
      sm_GraphicsLoop();
    }
#endif // SM_DEBUG_CODE

  // See what we need to do with each face in terms of extensions.
  SmTArray<SmFace*> sOffsetFaces;
  m_pOffsetBrep->GetFaces(sOffsetFaces);
  for (ii=0; ii<sOffsetFaces.GetSize(); ii++) 
    {
       SmFace *pF = sOffsetFaces[ii];
       if (pF->IsMarked(eMarkType)) 
         { continue; }

      if ( m_sClosingFaces.IsIn( pF ) )  // [B549]
        { continue; }

       pF->Mark(eMarkType);
       SmSurface * pOrigSurf   = pF->GetSurface();
       SmSurface * pExtSurf    = NULL ;
       SmExtent2d  sOldDomain  = pF->GetUVDomain();
       SmExtent2d  sCurrDomain = pOrigSurf->GetNaturalUVDomain();
       SmExtent2d  sExtDomain;
       SmExtent2d  sTrimDomain;

       // Fetch or create Original Surf extension
       SER(GetCreateExtendedSurface( pOrigSurf, sExtDomain, pExtSurf ));

       // Temporarily change Face to use the Surface extensions.  Without the extension,  
       // CalculateUVDomainFromUVTrimCurves() has problems with uv points outside the domain.
       // Use SmBrep::ReplaceSurface to properly handle composite surfaces
       SER(m_pOffsetBrep->ReplaceSurface(pF,pExtSurf,UNSURE, NULL, FALSE));

       // pF->SetUVDomain(sExtDomain) ;
       // pF->SetSurface(pExtSurf) ;

       // Pass False: allow CalculateUVDomain to return bigger domain than currently exists. [B359]
       SER(pF->CalculateUVDomainFromUVTrimCurves( sTrimDomain, FALSE ));

       // next: pick the surface and domain size for pF, delete the unused surface

       // When TrimDomain is smaller than CurrDomain - face uses Curr Surface
       if(sTrimDomain.IsContainedBy(sCurrDomain, SM_EFF_ZERO))
         {
           // Restore Face's original Surface and domain.
           SER(m_pOffsetBrep->ReplaceSurface(pF,pOrigSurf,UNSURE, NULL, FALSE));
           // pF->SetUVDomain(sOldDomain);
           // pF->SetSurface(pOrigSurf);

           // Note: this leaves a stale pointer in m_vSurfExtSurf table.
           SM_ASSERT(pExtSurf != NULL) ; delete pExtSurf ; pExtSurf = NULL ;
         }
       // when TrimDomain is bigger than CurrDomain and smaller than ExtDomain - Face uses Extended surface
       else if(   !sTrimDomain.IsContainedBy(sCurrDomain, SM_EFF_ZERO)
               &&  sTrimDomain.IsContainedBy(sExtDomain, SM_EFF_ZERO))
         {
           // Set Face to use extended Surf and Domain
           //   gwc:already done - no need to duplicate these calls here
           // pF->SetUVDomain(sExtDomain) ;
           // pF->SetSurface(pExtSurf) ;

           // Note: this leaves a stale pointer in m_vSurfExtSurf table.
           SM_ASSERT(pOrigSurf != NULL) ; delete pOrigSurf ; pOrigSurf = NULL ;

         } 
       // when TrimDomain is bigger than ExtDomain - signal an error - restore face - delete ExtSurface
       else
         { 
           // Restore Face's original Surface and domain.
           SER(m_pOffsetBrep->ReplaceSurface(pF,pOrigSurf,UNSURE, NULL, FALSE));
           // pF->SetUVDomain(sOldDomain);
           // pF->SetSurface(pOrigSurf);

           // Note: this leaves a stale pointer in m_vSurfExtSurf table.
           SM_ASSERT(pExtSurf != NULL) ; delete pExtSurf ; pExtSurf = NULL ;

           SER(SM_ERR) ; 
         }
       
       // gwc: bug fix - replaced following block with above block   
       // if (!sTrimDomain.IsContainedBy(sCurrDomain, SM_EFF_ZERO)) 
       //   {
       //     if (!sTrimDomain.IsContainedBy(sExtDomain, SM_EFF_ZERO)) 
       //       { SER(SM_ERR); }
       // 
       //     SER(m_pOffsetBrep->ReplaceSurface(pF,pExtSurf,FALSE,NULL,FALSE));
       //     pF->SetUVDomain(sTrimDomain);
       //   }
       // else 
       //   {
       //     // Note: this leaves a stale pointer in m_vSurfExtSurf table.
       //     SM_ASSERT(pExtSurf != NULL) ; delete pExtSurf ; pExtSurf = NULL ;
       //   }
//         pF->FillVertexGaps();  This causes serious problems

    } // end for each face in the offset Brep
  
#ifdef SM_VALIDATE_TOPOLOGY
      m_pOffsetBrep->ValidatePointers();
#endif // SM_VALIDATE_TOPOLOGY

#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe4 = FALSE;
  if (bDebugMe4 || bDebugMeAll) 
    {
      m_pOffsetBrep->ValidatePointers();
      SmTrimmingTools::ClearUVTrimCurves(m_pOffsetBrep);
      SmTrimmingTools::CreateUVTrimCurves(m_pOffsetBrep);
      m_pOffsetBrep->ValidatePointers();

      smgfx_Erase();
      smgfx_SetLook(1,2, 0,0,0); if(m_pOffsetBrep) m_pOffsetBrep->Draw(TRUE); sm_GraphicsLoop();
      sm_GraphicsLoop();

      smgfx_SetUVHatchCount(20,20);
      smgfx_SetDrawingMode(SM_DM_CROSSHATCH);
      smgfx_Erase();
      smgfx_SetLook(1,2, 0,0,1); if(m_pOffsetBrep) m_pOffsetBrep->Draw(TRUE); sm_GraphicsLoop();
      sm_GraphicsLoop();

      smgfx_SetDrawingMode(SM_DM_NORMALS);
      smgfx_Erase();
      smgfx_SetLook(1,2, 0,1,0); if(m_pOffsetBrep) m_pOffsetBrep->Draw(TRUE); sm_GraphicsLoop();
      sm_GraphicsLoop();

      smgfx_SetDrawingMode(SM_DM_WIREFRAME);
    }
#endif // SM_DEBUG_CODE

  // SolidCraft change: an offset face whose loops cross is rebuilt from the stretches of its
  // loops that bound it, before the merge sees it.
  SER(ResolveCrossingLoops());

  SmBrep *pResult;
  SER(PiecewiseMerge(pResult));
  SmObjDelete sClean2(pResult);

  if ( ! m_bCreateOffsetSolid )
    { SER(RemoveExcessFaces(pResult)); }

  // arrive here when pResult contains all the offset and copied surfaces
  // merged (intersected and inserted) together. It won't be manifold
  // because all the overlapping pieces of the offset surfaces are still in pResult.

#ifdef SM_DEBUG_CODE
  // draw originalBrep(blue), offsetBrep(green), piecewiseMergeResult(yellow)
  if (bDebugMe || bDebugMeAll) 
    {
      m_pOriginalBrep->Dump() ; // input
      if ( m_pOffsetBrep ) { m_pOffsetBrep->Dump(); }    // intermediate structure
      pResult->Dump() ;         // the output - all faces merged - not yet manifold

      SM_ASSERT_VALID(pResult) ;

      smgfx_Erase();
      smgfx_SetLook(1,2) ; pResult->Draw(TRUE); sm_GraphicsLoop();
      smgfx_SetLook(1,2, 0,0,1) ; m_pOriginalBrep->Draw(TRUE); sm_GraphicsLoop();
      if ( m_pOffsetBrep )
        { smgfx_SetLook(1,2, 0,1,0) ; m_pOffsetBrep->Draw(TRUE); sm_GraphicsLoop(); }
      sm_GraphicsLoop();
    }
#endif // SM_DEBUG_CODE

  // Do some stitching to close those unwanted gaps
  //    this is a protective step that would not be needed in an ideal world.
  // SolidCraft change: within the merged body's tolerance only. StitchIntoSolid doubles its
  // tolerance until no lamina edge is left (up to a 200th of the body's size), but the lamina
  // edges a merge leaves are the rims of the pieces its intersections cut off (an offset face's
  // part beyond a neighbour's offset), which MakeManifold removes below: escalating, it welded
  // two of them 0.12 mm apart, an arc and the chord between its ends, into one edge of tolerance
  // 0.24 mm, and the offset came apart there.
  {
    SmStitchCallback sStitchCallback;
    SmStitch sStitch(sStitchCallback, pResult->GetTolerance(), TRUE, FALSE);
    sStitch.m_bFastEdgeCompare        = TRUE;
    sStitch.m_bSplitEdgesWithVertices = FALSE;
    sStitch.m_bMakingManifoldSolid    = TRUE;
    ULONG  lStitchedEdges = 0, lLaminaEdges = 0;
    double dVGap = 0.0, dEGap = 0.0;
    SER(sStitch.DoStitching(pResult, NULL, NULL, lStitchedEdges, lLaminaEdges, dVGap, dEGap));
  }

#ifdef SM_DEBUG_CODE
  SmBoolean bDrawFaces = FALSE;
  // draw originalBrep(blue), offsetBrep(green), piecewiseMergeAndStitchResult(yellow)
  if (bDebugMe || bDebugMeAll) 
    {
      m_pOriginalBrep->Dump() ; // input                                             
      if ( m_pOffsetBrep ) { m_pOffsetBrep->Dump(); } // intermediate structure                            
      pResult->Dump() ;         // the output - all faces merged and stitched so that faces 
                                //              near one another boundaries are now connected                                                     

      SM_ASSERT_VALID(pResult) ;

      smgfx_Erase();
      // draw input, intermediate, and output breps
      smgfx_SetLook(1,2) ;        pResult->Draw(TRUE); sm_GraphicsLoop();
      smgfx_SetLook(1,2, 0,0,1) ; m_pOriginalBrep->Draw(TRUE); sm_GraphicsLoop();
      if ( m_pOffsetBrep ) {
          smgfx_SetLook(1,2, 0,1,0) ; m_pOffsetBrep->Draw(TRUE); sm_GraphicsLoop();
      }

      if(bDrawFaces)
        {
          SmTArray<SmFace *>    sFacesDraw ;
          SmTArray<SmFaceuse *> sFaceuses ;
          SmTopologyTraverser   sTraverser ;
          SmBoolean bIsClosed ;
          
          pResult->GetFaces( sFacesDraw ) ;
          SmFaceuse *pFaceuse = sFacesDraw.GetSize() > 8 ? sFacesDraw[8]->GetUpwardFaceuse() : NULL ;
          if(pFaceuse) { SmNewMarkAndLock sMarkLockDebug(pResult->GetContext(), SM_MT_ALLMARKS) ; // increments newly locked mark
                         sTraverser.FaceuseTraversal(pFaceuse,    // in : target faceuse (gets marked)                              
                                                     bIsClosed,   // out: TRUE=faceuse connects to closed faceuse set, FALSE=doesn't
                                                     &sFaceuses,  // out: List of all connected faces (get marked), NULL to ignore.
                                                      sMarkLockDebug ) ; // in : specify mark for target objects (not incremented)
                       }
          ULONG kk;
          smgfx_SetLook(1,2, 0,1,0) ; for(kk=0;kk<sFacesDraw.GetSize();kk++)
                                        {
                                          SmFace *pFaceA = sFacesDraw[kk] ;
                                          pFaceA->Draw(SM_DM_CROSSHATCH) ; sm_GraphicsLoop() ;
                                        }

          smgfx_SetLook(1,2, 1,0,0) ; if(pFaceuse) pFaceuse->Draw(SM_DM_CROSSHATCH) ; sm_GraphicsLoop() ;
          smgfx_SetLook(1,2, 0,1,1) ; for(kk=0;kk<sFaceuses.GetSize();kk++)
                                        {
                                          SmFaceuse *pFaceuseA = sFaceuses[kk] ;
                                          smgfx_Erase() ;
                                          smgfx_SetLook(1,2) ; pResult->Draw(TRUE); sm_GraphicsLoop();
                                          smgfx_SetLook(3,5, 0,1,0) ; pFaceuseA->GetFace()->Draw(SM_DM_CROSSHATCH) ; sm_GraphicsLoop() ;
                                          smgfx_SetLook(1,2, 0,1,0) ; pFaceuseA->Draw(SM_DM_FACEUSENEIGHBORS) ; sm_GraphicsLoop() ;
                                          pFaceuseA->GetFace()->Dump() ;
                                          sm_GraphicsLoop() ;
                                        }
        } // end DrawFaces check
      sm_GraphicsLoop();
    } // end if (bDebugMe)
#endif // SM_DEBUG_CODE

  // Make manifold and set void regions [B550]
  // Here we will need to select regions to keep and
  // send in to make manifold.
  if (pResult->MakeManifold() == SM_SUCCESS) 
    {
      // make sure infinite regions are set correctly
      SmTArray<SmRegion*> sRegions;
      pResult->GetRegions(sRegions);
      if (sRegions.GetSize() > 1) 
        {
          SmTemporaryChangeValue<SmBoolean> sChange(pResult->m_bEditingEnabled,FALSE);
          SER(pResult->FindAndSetInfiniteRegion());
        }
      // label every region as inside or outside
      pResult->GetRegions(sRegions);
      for (ULONG kkk=0; kkk<sRegions.GetSize(); kkk++) 
        {
          SmRegion *pReg = sRegions[kkk];
#ifdef SM_DEBUG_CODE
          if (bDebugMe || bDebugMeAll)
            {
              smgfx_Erase();
              pReg->Draw(); sm_GraphicsLoop();
              sm_GraphicsLoop();
            }
#endif //SM_DEBUG_CODE
          if (pReg == pResult->GetInfiniteRegion()) { pReg->SetIsVoid(TRUE); }
          else                                      { pReg->SetIsVoid(FALSE); }
        } // end iter every region labeling inside/outside status
    } // end MakeManifold a success check

#ifdef SM_DEBUG_CODE
  // draw originalBrep(blue), offsetBrep(green), piecewiseMergeResult(yellow)
  if (bDebugMe || bDebugMeAll) 
    {
      m_pOriginalBrep->Dump() ; // input
      if ( m_pOffsetBrep )
        { m_pOffsetBrep->Dump(); }    // intermediate structure
      pResult->Dump() ;         // the output - all faces merged - not yet manifold

      SM_ASSERT_VALID(pResult) ;

      smgfx_Erase();
      smgfx_SetLook(1,2) ; pResult->Draw(TRUE); sm_GraphicsLoop();
      smgfx_SetLook(1,2, 0,0,1) ; m_pOriginalBrep->Draw(TRUE); sm_GraphicsLoop();
      if ( m_pOffsetBrep )
        { smgfx_SetLook(1,2, 0,1,0) ; m_pOffsetBrep->Draw(TRUE); sm_GraphicsLoop(); }
      sm_GraphicsLoop();
    }
#endif // SM_DEBUG_CODE

  // keep the result, set output and return
  sClean2.Clear();
  rpResult = pResult;

#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe3 = FALSE;
  if (bDebugMe3 || bDebugMeAll) 
    {
      pResult->Dump();
      pResult->ValidatePointers();
      SmTrimmingTools::ClearUVTrimCurves(pResult);
      SmTrimmingTools::CreateUVTrimCurves(pResult);
      pResult->ValidatePointers();

      smgfx_Erase();
      smgfx_SetLook(1,2, 0,0,0); if(pResult) pResult->Draw(); sm_GraphicsLoop();
      sm_GraphicsLoop();

      smgfx_SetUVHatchCount(20,20);
      smgfx_SetDrawingMode(SM_DM_CROSSHATCH);
      smgfx_Erase();
      smgfx_SetLook(1,2, 0,0,1); if(pResult) pResult->Draw(); sm_GraphicsLoop();
      sm_GraphicsLoop();

      smgfx_SetDrawingMode(SM_DM_NORMALS);
      smgfx_Erase();
      smgfx_SetLook(1,2, 0,1,0); if(pResult) pResult->Draw(); sm_GraphicsLoop();
      sm_GraphicsLoop();

      smgfx_SetDrawingMode(SM_DM_WIREFRAME);
    }
#endif // SM_DEBUG_CODE

  // The pointer to m_pOffsetBrep has been passed around
  // to other variables, so do this by hand here.
  //sCleanOffBrep.Clear();
  if ( m_pOffsetBrep != NULL ) { delete m_pOffsetBrep; m_pOffsetBrep = NULL; }

  return SM_SUCCESS;

} // end SmOffsetExecutive::DoExtendedOffset

//------------------------------------------
// Internal methods
//------------------------------------------

/*******************************************************************//**
PURPOSE: Static helper function: Call Notify on an Edge and its Vertices.

NOTES: 
   Attribute propagation (i.e., Notify) is awkward in this operation,
   mainly because copying of entities is not done in a simple manner.
   For example, Faces are copied/offset by creating curves from the face
   ( CreateCurvesFromFace ), then copying the surface, then calling
   MakeFaceWithCurves from that geometry.  Therefore, we have to be
   very careful about calling Notify throughout this operation.  [B180]
***********************************************************************/
static void sm_NotifyEdgeAndVertices
 ( SmEdge *pOldEdge, 
   SmEdge *pNewEdge, 
   SmNotifyOperation eOp )
{
  SM_ASSERT_MSG(eOp == SM_NO_COPY, _T("sm_NotifyEdgeAndVertices only works on eOp == SM_NO_COPY")) ;

  pOldEdge->Notify(eOp, pNewEdge, SM_NO_GET_BREP(pNewEdge), SM_NO_GET_BREP(pOldEdge));
  SmVertex *pOV1 = pOldEdge->GetStartVertex();
  SmVertex *pNV1 = pNewEdge->GetStartVertex();
  pOV1->Notify(eOp, pNV1, SM_NO_GET_BREP(pNV1), SM_NO_GET_BREP(pOV1));

  // Fortunately, Notify checks for existing attributes of the type given
  // and will not attach duplicates.

  SmVertex *pOV2 = pOldEdge->GetOtherVertex( pOV1 );
  if ( pOV2 != pOV1 )
    {
      SmVertex *pNV2 = pNewEdge->GetOtherVertex( pNV1 );
      pOV2->Notify(eOp, pNV2, SM_NO_GET_BREP(pNV2), SM_NO_GET_BREP(pOV2));
    }

} // end sm_NotifyEdgeAndVertices

/*******************************************************************//**
PURPOSE: Replace every SmOffsetSurface type surface in Brep with
            an approximate SmBSplineSurface.

NOTES: 
***********************************************************************/
SmStatus SmOffsetExecutive::ReplaceImplicitOffsets
  (SmBrep * pBrepToFix,             // in : target Brep to modify
   double dThisApproxTol3d)        // in : max allowed deviation in approximation surfaces
{
  // remove curve and surface portions not used by edges and faces
  SER( pBrepToFix->ShrinkGeometry() );

  // locals
  SmTArray<SmFace*>    sFixFaces;
  SmTArray<SmSurface*> sNewSurfaces;
  SmTArray<SmFace*>    sFaces;
  pBrepToFix->GetFaces(sFaces);
  ULONG ii;

  // for every face
  for (ii=0; ii<sFaces.GetSize(); ii++) 
  {
      SmFace          *pF       = sFaces[ii];
      SmSurface       *pSurface = pF->GetSurface();
      SmOffsetSurface *pOff     = SM_CAST_PTR(SmOffsetSurface,pSurface);

      // skip non SmOffsetSurfaces
      if (!pOff) { continue; }

      // prepare to build approximation to exact offset surface
      SmBSplineSurface *pOffNurb = NULL;

      // get SmOffsetSurface base (use extended surface if there is one)
      const SmSurface *pBase =   (pOff->GetExtendedBaseSurface())
                            ? pOff->GetExtendedBaseSurface()
                            : pOff->GetBaseSurface() ;

      // make sure base surface is an SmBSplineSurface
      SmBSplineSurface *pNurb = SM_CAST_PTR(SmBSplineSurface,pBase);
      if (!pNurb) { SER(SM_ERR); }

      // use NLIB function to build offset surface approximation
      SER(pNurb->ApproximateOffsetSurface(*pBrepToFix->GetContext(),pOff->GetOffsetDistance(),
          dThisApproxTol3d,pOffNurb, 0));

      // remember face/surface pair 
      sFixFaces.Add(pF);
      sNewSurfaces.Add(pOffNurb);

  } // end iter every face

  // for every face that had an approximate surface built for an exact offset surface
  for (ii=0; ii<sFixFaces.GetSize(); ii++) 
  {
      SmSurface * pNewSurface = sNewSurfaces[ii];
      SmFace    * pF          = sFixFaces[ii];
      pF->Notify(SM_NO_PRE_EDIT, pF, SM_NO_GET_BREP(pF), NULL);
      SmExtent2d sDomain = pNewSurface->GetNaturalUVDomain();
      pF->SetUVDomain(sDomain);

      // swap the approx surface for the exact offset surface
      SER(pBrepToFix->ReplaceSurface(pF,pNewSurface,TRUE,NULL));
      pF->Notify(SM_NO_POST_EDIT, pF, SM_NO_GET_BREP(pF), NULL);
    
  } // end iter every modified face 

//    SER(pBrepToFix->LocalOperation(sFixFaces,sNewSurfaces));
  return SM_SUCCESS;

} // end SmOffsetExecutive::ReplaceImplicitOffsets

// Forget the pairs of a relation whose second is gone (or, with cpLive, is not
// live). The pointers are only compared, never followed.
template <class FIRST, class SECOND>
static void sm_ForgetSeconds(SmRelation<FIRST,SECOND> & rRelation,
                             const std::unordered_set<SECOND*> & crGone,
                             const std::unordered_set<SECOND*> * cpLive)
{
    SmTArray<SECOND*> sSeconds;
    rRelation.GetAllSeconds(sSeconds);
    for (ULONG i=0; i<sSeconds.GetSize(); i++)
    {
        SECOND *pSecond = sSeconds[i];
        if (crGone.count(pSecond) == 0 && (cpLive == NULL || cpLive->count(pSecond) != 0)) { continue; }
        SmTArray<FIRST*> sFirsts;
        rRelation.GetFirsts(pSecond, sFirsts);
        for (ULONG j=0; j<sFirsts.GetSize(); j++)
        {
            (void)rRelation.RemovePair(sFirsts[j], pSecond); // SolidCraft: each pair was just listed, so it exists
        }
    }
}

// Forget the pairs of a relation whose first is gone.
template <class FIRST, class SECOND>
static void sm_ForgetFirsts(SmRelation<FIRST,SECOND> & rRelation, const std::unordered_set<FIRST*> & crGone)
{
    SmTArray<FIRST*> sFirsts;
    rRelation.GetAllFirsts(sFirsts);
    for (ULONG i=0; i<sFirsts.GetSize(); i++)
    {
        FIRST *pFirst = sFirsts[i];
        if (crGone.count(pFirst) == 0) { continue; }
        SmTArray<SECOND*> sSeconds;
        rRelation.GetSeconds(pFirst, sSeconds);
        for (ULONG j=0; j<sSeconds.GetSize(); j++)
        {
            (void)rRelation.RemovePair(pFirst, sSeconds[j]); // SolidCraft: each pair was just listed, so it exists
        }
    }
}

/*******************************************************************//**
PURPOSE: Remove Excess Faces of the Result by measurement to the original.

NOTES: 
***********************************************************************/
SmStatus SmOffsetExecutive::RemoveExcessFaces(SmBrep * pResult)
{
    // If we are not merging don't do this.
    if (!m_bMergeResults) {
        return SM_SUCCESS; 
    }
    // If shelling this will also not work.
    if (m_pShellFaces && m_pShellFaces->GetSize() > 0) {
        return SM_SUCCESS;
    }


    SmTArray<SmEdge*> sEdges;
    SmTArray<SmFace*> sDeleteFaces;
    SmTArray<SmFace*> sEdgeFaces;
    
    double dMinDistance = m_rGeometryCreation.GetOffsetDistance() * 0.95;
    
    SmSolution sData[16];
    SmSolutionArray sSolutions(16,sData);

    pResult->GetEdges(sEdges);
    for (ULONG i=0; i<sEdges.GetSize(); i++) {
        SmEdge *pE = sEdges[i];
        SmPoint3d sPnt;
        SER( pE->GetPrimaryEdgeuse()->NormalizedEvaluate(0.5,FALSE,sPnt) );  // TRUE = UV Eval, FALSE = 3d Eval
        sSolutions.ReSet();
        m_pOriginalBrep->m_bEditingEnabled = FALSE;
        SER(SmTopologySolver::BrepPointSolve(m_pOriginalBrep,sPnt,SM_SO_MINIMIZE,
            SM_SR_SINGLE,m_pOriginalBrep->GetTolerance(),SM_BIG_DOUBLE,NULL,sSolutions));
        m_pOriginalBrep->m_bEditingEnabled = TRUE;
        if (sSolutions.GetSize() == 1)
        {
            SmSolution &rSol = sSolutions[0];
            if (rSol.m_vStart.m_dSolutionValue < dMinDistance)
            {
                pE->GetFaces(sEdgeFaces);
                SmBoolean bSkip = FALSE;
                for (ULONG kk=0; kk<sEdgeFaces.GetSize(); kk++)
                {
                    SmFace *pF = sEdgeFaces[kk];
                    ULONG lIndex;
                    if (m_pShellFaces && m_pShellFaces->FindElement(pF,lIndex)) {
                        bSkip = TRUE; // This is ok because it is part of a shell face.
                        break;
                    }
                }
                if (bSkip) { continue; }
                for (ULONG j=0; j<sEdgeFaces.GetSize(); j++)
                {
                    SmFace *pF = sEdgeFaces[j];
                    sDeleteFaces.AddUnique(pF);
                }
            }
        }
    }

    // The faces about to be removed and their edgeuses, which RemoveFaces deletes
    std::unordered_set<SmFace*> sGoneFaces;
    std::unordered_set<SmEdgeuse*> sGone;
    for (ULONG i=0; i<sDeleteFaces.GetSize(); i++)
    {
        sGoneFaces.insert(sDeleteFaces[i]);
        SmTArray<SmEdge*> sFaceEdges;
        sDeleteFaces[i]->GetEdges(sFaceEdges);
        for (ULONG j=0; j<sFaceEdges.GetSize(); j++)
        {
            SmTArray<SmEdgeuse*> sUses;
            SER( sFaceEdges[j]->GetEdgeuses(sUses) );
            for (ULONG kk=0; kk<sUses.GetSize(); kk++)
            {
                if (sUses[kk]->GetFace() == sDeleteFaces[i])
                {
                    sGone.insert(sUses[kk]);
                }
            }
        }
    }

    // Clean up those faces we have collected.
    SER( pResult->RemoveFaces(sDeleteFaces) );   // increments unlocked mark value

    // Forget the deleted faces and edgeuses, so later lookups in the relations (lamina side walls,
    // PiecewiseMerge) cannot return them. For m_vEUToEU also forget any edgeuse no result edge has.
    if (!sGoneFaces.empty())
    {
        std::unordered_set<SmEdgeuse*> sLive;
        SmTArray<SmEdge*> sResultEdges;
        pResult->GetEdges(sResultEdges);
        for (ULONG i=0; i<sResultEdges.GetSize(); i++)
        {
            SER( sResultEdges[i]->GetEdgeuses(sLive) );
        }
        sm_ForgetSeconds(m_vEUToEU, sGone, &sLive);
        sm_ForgetSeconds(m_vFUToF, sGoneFaces, (const std::unordered_set<SmFace*> *)NULL);
        sm_ForgetSeconds(m_vEToF, sGoneFaces, (const std::unordered_set<SmFace*> *)NULL);
        sm_ForgetSeconds(m_vVToF, sGoneFaces, (const std::unordered_set<SmFace*> *)NULL);
        sm_ForgetSeconds(m_vBlendToFaces, sGoneFaces, (const std::unordered_set<SmFace*> *)NULL);
        sm_ForgetFirsts(m_vBlendToFaces, sGoneFaces);
    }

    return SM_SUCCESS;

} // end RemoveExcessFaces

/*******************************************************************//**
PURPOSE: Get (with creation if necessary) the extended surface of this face.

NOTES: 
***********************************************************************/
SmStatus SmOffsetExecutive::GetCreateExtendedSurface
 (SmSurface  * pBase,             // in : base surface to be extended 
  SmExtent2d & rUVDomainOfOffset, // out: Extended Surface Domain
  SmSurface *& rpExtSurface)      // out: The extended surface
{
  // check input 
  NER(pBase);

  // Dist to extend
  double dDistanceToExtend = m_rGeometryCreation.GetOffsetDistance() * 5.0;

  // retrieve existing ExtendedSurface by Key if it exists
  SmSurface *pExtendedSurface = (SmSurface*)m_vSurfExtSurf.At(pBase);

  // no work - extended surface already exists
  if (pExtendedSurface)
    {
      // set outputs - and return
      rpExtSurface      = pExtendedSurface;
      rUVDomainOfOffset = pExtendedSurface->GetNaturalUVDomain();
      return SM_SUCCESS;
    }

  // arrive here when an extended surface has to be built

  // init output
  rpExtSurface = NULL;

  // locals
  SmSurface  * pOrigSurf = pBase;
  SmExtent2d   sExtDomain;

  // when OriginalSurf is a BSplineSurface
  if (pOrigSurf != NULL)
    {
      SmSurface * pNewExtSurf = NULL ;
      // create the ExtendedSurface
      SER(pOrigSurf->CreateExtendedSurface(m_crContext,
                                           dDistanceToExtend,
                                           SM_CT_G1R,
                                           pNewExtSurf));
      NER(pNewExtSurf);

      // set output
      rUVDomainOfOffset = pNewExtSurf->GetNaturalUVDomain();
      rpExtSurface      = pNewExtSurf;
    }
  else // OriginalSurf is not a BSplineSurface branch 
    {
      // If we have an offset surface just make it bigger.
      SmOffsetSurface *pOffset = SM_CAST_PTR(SmOffsetSurface,pBase);
      if (pOffset)
        {
          SM_DBG_WARN( _T( "This should be obsolete. Removal pending." ) );
          // The way we extend an SmOffsetSurface is to call EvaluatePoint()
          // at a uv outside the domain: that creates or extends its ExtendedSurface
          // as a side effect.

          // copy the OffsetSurface and set up to allow evaluates to extend the surface
          SmSurface * pSurface;
          SER(pOffset->Copy(m_crContext,pSurface)); NER(pSurface);
          pOffset->Notify(SM_NO_PRE_EDIT, pOffset, SM_NO_GET_OWNER(pOffset), NULL);
          pOffset = SM_CAST_PTR(SmOffsetSurface,pSurface);
          pOffset->SetAllowExtension(TRUE);

          // Calculate uv positions that are outside the domain
          // by a distance corresponding to the input 3d distance.
          SmExtent2d sDomain = pOffset->GetNaturalUVDomain();
          SmVector2d sDomainSize = sDomain.GetSize();
          SmPoint3d sPnt, sDU, sDV;

          // estimate New Min and Max values as linear extensions of the boundary cross tangents
          // Evaluate1stDerivatives() returns non-zero first derivatives.
          SER(pOffset->Evaluate1stDerivatives(sDomain.GetMin(),TRUE,TRUE,sPnt,sDU,sDV));
          SmPoint2d sDiffMin(dDistanceToExtend/sDU.Length(),dDistanceToExtend/sDV.Length());
          if ( smos_Fabs( sDiffMin.x ) > sDomainSize.x * 10.0 )
            { sDiffMin.x = ( sDiffMin.x > 0 ) ? sDomainSize.x * 10.0 : - sDomainSize.x * 10.0; }
          if ( smos_Fabs( sDiffMin.y ) > sDomainSize.y * 10.0 )
            { sDiffMin.y = ( sDiffMin.y > 0 ) ? sDomainSize.y * 10.0 : - sDomainSize.y * 10.0; }

          SmPoint2d sNewMin = sDomain.GetMin() - sDiffMin;

          SER(pOffset->Evaluate1stDerivatives(sDomain.GetMax(),TRUE,TRUE,sPnt,sDU,sDV));
          SmPoint2d sDiffMax(dDistanceToExtend/sDU.Length(),dDistanceToExtend/sDV.Length());
          if ( smos_Fabs( sDiffMax.x ) > sDomainSize.x * 10.0 )
            { sDiffMax.x = ( sDiffMax.x > 0 ) ? sDomainSize.x * 10.0 : - sDomainSize.x * 10.0; }
          if ( smos_Fabs( sDiffMax.y ) > sDomainSize.y * 10.0 )
            { sDiffMax.y = ( sDiffMax.y > 0 ) ? sDomainSize.y * 10.0 : - sDomainSize.y * 10.0; }

          SmPoint2d sNewMax = sDomain.GetMax() + sDiffMax;

          // increase the domain by evaluating at the estimated extensions
          SER(pOffset->EvaluatePoint(sNewMin,sPnt));
          SER(pOffset->EvaluatePoint(sNewMax,sPnt));

          // set outputs
          rUVDomainOfOffset = pOffset->GetNaturalUVDomain();
          rpExtSurface = pOffset;
        }
      else
        {
          SER(SM_ERR); // Don't know what to do with this kind of surface
        }
    }
    
  m_vSurfExtSurf.Insert(pBase,rpExtSurface);

  return SM_SUCCESS;

} // end GetCreateExtendedSurface

/*******************************************************************//**
PURPOSE: Do a smart Piecewise Merge of the offset Brep.

NOTES:

METHOD ---
  if(m_bMergeResults == TRUE)
    call SmMerge::PiecewiseMerge on m_pOffsetBrep after taking the
      time to relate every pair of offset surfaces that should
      already be sharing an edge.
  else
    SmStitch::StitchIntoShell or SmStitch::StitchIntoSolid on m_pOffsetBrep
     
***********************************************************************/
SmStatus SmOffsetExecutive::PiecewiseMerge
  (SmBrep *& rpResult)                     // out: result of PiecewiseMerge or stitching
                                           //      applied m_pOffsetBrep
{
  rpResult = NULL;

  // If we are not merging just do a stitch
  if ( ! m_bMergeResults )
    {
      // Clean up relationships.
      m_vFUToF.RemoveAll();
      m_vEUToEU.RemoveAll();
      m_vEToF.RemoveAll();
      m_vVToF.RemoveAll();
      m_vBlendToFaces.RemoveAll();

      // let rpResult = m_pOffsetBrep
      rpResult      = m_pOffsetBrep;
      m_pOffsetBrep = NULL;

      ULONG     lStitchedEdges;
      double    dVGap, dEGap;
      SmBoolean bProducesSolid;

      if (m_bCreateOffsetSolid) 
        {
          // Tol: StitchIntoShell() adds in the tols of the topology objects,
          // so the Brep's tolerance should be plenty.
          double dStitchTol3d = rpResult->GetTolerance();
          SER( SmStitch::StitchIntoShell( rpResult, TRUE, dStitchTol3d, lStitchedEdges, dVGap, dEGap ) );
        }
      else 
        {
          SER( SmStitch::StitchIntoSolid( rpResult, bProducesSolid, lStitchedEdges, dVGap, dEGap ) );
        }

      SmTArray<SmRegion*> sRegions;
      rpResult->GetRegions(sRegions);
      if (sRegions.GetSize() > 1) 
        {
          SmTemporaryChangeValue<SmBoolean> sChange(rpResult->m_bEditingEnabled,FALSE);
          SER(rpResult->FindAndSetInfiniteRegion());
        }
      return SM_SUCCESS;

    } // end !m_bMergeResults so just do a stitch and return check

  // locals
  SmRelation<SmSurface,SmSurface> sNoSSI(m_crContext);
  SmTArray<SmFaceuse*> sFaceuses;
  SmTArray<SmEdge*>    sEdges;
  SmTArray<SmFace*>    sEOffFaces;
  SmTArray<SmFace*>    sOffFaces;
  SmTArray<SmFace*>    sOffFaces2;
  SmTArray<SmFace*>    sBlended;
  SmTArray<SmFace*>    sTanFaces;
  SmTArray<SmEdgeuse*> sEdgeEUS;
  SmTArray<SmVertex*>  sVertices;
  ULONG ii, jj, kk, ll, mm;

  // Set up non-intersecting pairs into a relationship
  //  this includes FaceOffset/EdgeOffset   faces from face/face->Edge pairs
  //                EdgeOffset/VertexOffset faces from edge/edge->Vertex pairs
  // This is a heuristic to speed things up.  It is possible
  // that this will fail in some cases where surfaces self-intersect

  // for every origBrep faceuse that has been offset or copied
  m_vFUToF.GetAllFirsts(sFaceuses);
  for (ii=0; ii<sFaceuses.GetSize(); ii++) 
    {
      SmFaceuse *pFU = sFaceuses[ii];
      SmFace    *pF  = pFU->GetFace();

      // for every OrigFace->OffsetFace
      m_vFUToF.GetSeconds(pFU,sOffFaces);
      for (jj=0; jj<sOffFaces.GetSize(); jj++) 
        {
          SmFace    *pOffFace = sOffFaces[jj];
          SmSurface *pOffSurf = pOffFace->GetSurface();
          
          // for every OrigFace->Edge
          pF->GetEdges(sEdges);
          for (kk=0; kk<sEdges.GetSize(); kk++) 
            {
              SmEdge *pE = sEdges[kk];

              // for every origEdge offset capFace
              m_vEToF.GetSeconds(pE,sEOffFaces);
              for (ll=0; ll<sEOffFaces.GetSize(); ll++) 
                {
                  SmFace    *pEOffFace = sEOffFaces[ll];
                  if (!pEOffFace) { continue; }
                  SmSurface *pEOffSurf = pEOffFace->GetSurface();

                  // Relate corner blend to blended surfaces
                  m_vBlendToFaces.GetSeconds(pEOffFace, sBlended);
                  for (mm=0; mm<sBlended.GetSize(); mm++) 
                    {
                      SmFace *pBlendedF = sBlended[mm];
                      (void)sNoSSI.RelatePair(pBlendedF->GetSurface(),pEOffSurf); // SolidCraft: a pair this speed-up misses is intersected like any other
                      (void)sNoSSI.RelatePair(pEOffSurf, pBlendedF->GetSurface()); // SolidCraft: a pair this speed-up misses is intersected like any other
                    }

                  // relate faceOffsets to their edgeOffsets
                  (void)sNoSSI.RelatePair(pOffSurf, pEOffSurf); // SolidCraft: a pair this speed-up misses is intersected like any other
                  (void)sNoSSI.RelatePair(pEOffSurf, pOffSurf); // SolidCraft: a pair this speed-up misses is intersected like any other
                  
                  // Now do first vertex
                  SmVertex *pV        = pE->GetVertex();
                  SmFace   *pVOffFace = m_vVToF.GetSecond(pV);
                  if (pVOffFace) 
                    {
                      // relate start VertexOffsets to EdgeOffsets
                      SmSurface *pVOffSurf = pVOffFace->GetSurface();
                      (void)sNoSSI.RelatePair(pEOffSurf,pVOffSurf); // SolidCraft: a pair this speed-up misses is intersected like any other
                      (void)sNoSSI.RelatePair(pVOffSurf,pEOffSurf); // SolidCraft: a pair this speed-up misses is intersected like any other
                    }
                  
                  // Now do second vertex
                  pV        = pE->GetOtherVertex(pV);
                  pVOffFace = m_vVToF.GetSecond(pV);
                  if (pVOffFace) 
                    {
                      // relate end VertexOffsets to EdgeOffsets
                      SmSurface *pVOffSurf = pVOffFace->GetSurface();
                      (void)sNoSSI.RelatePair(pEOffSurf,pVOffSurf); // SolidCraft: a pair this speed-up misses is intersected like any other
                      (void)sNoSSI.RelatePair(pVOffSurf,pEOffSurf); // SolidCraft: a pair this speed-up misses is intersected like any other
                    }
                  
                } // end iter 'll' every origEdge offset capFace

              if (sEOffFaces.GetSize() != 0) 
                {
                  continue;
                }

              // Determine if we have a tangent sector across this edge and
              // eliminate intersections if we do.
              // SolidCraft change: with extended convex edges, every edge without an offset cap
              // face had the offsets of its two faces mated along it (EdgeuseOffset), so those
              // faces already share their intersection edge too.
              double dTanTolDeg = GetTangencyTolDegrees();
              if (pE->IsManifold() && (m_bExtendConvexEdges || pE->GetPrimaryEdgeuse()->IsTangentSector(dTanTolDeg))) 
                {
                  SER( pE->GetEdgeuses(sEdgeEUS) );
                  ULONG ieu, ifac;
                  for (ieu=0; ieu<sEdgeEUS.GetSize(); ieu++) 
                    {
                      SmEdgeuse *pEU = sEdgeEUS[ieu];
                      if (pEU->GetFaceuse() == pFU) 
                        {
                          SmEdgeuse *pEURad = pEU->GetRadial();
                          SmFaceuse *pFURad = pEURad->GetFaceuse();
                          m_vFUToF.GetSeconds(pFURad,sOffFaces2);
                          for (ifac=0; ifac<sOffFaces2.GetSize(); ifac++) 
                            {
                              SmFace *pOffFRad = sOffFaces2[ifac];
                              SmSurface *pOffSurfRad = pOffFRad->GetSurface();
                              (void)sNoSSI.RelatePair(pOffSurf,pOffSurfRad); // SolidCraft: a pair this speed-up misses is intersected like any other
                              (void)sNoSSI.RelatePair(pOffSurfRad,pOffSurf); // SolidCraft: a pair this speed-up misses is intersected like any other
                            }
                        }
                    } // end for each Edgeuse
                } // end if tangent sector
            } // end iter 'kk' every OrigFace->Edge
          
          // get every face/vertex offset face pair
          pF->GetVertices(sVertices);
          for (kk=0; kk<sVertices.GetSize(); kk++) 
            {
              SmVertex *pV = sVertices[kk];
              SmFace *pVOffFace = m_vVToF.GetSecond(pV);
              if (!pVOffFace) { continue; }
              SmSurface *pVOffSurf = pVOffFace->GetSurface();
              (void)sNoSSI.RelatePair(pOffSurf,pVOffSurf); // SolidCraft: a pair this speed-up misses is intersected like any other
              (void)sNoSSI.RelatePair(pVOffSurf,pOffSurf); // SolidCraft: a pair this speed-up misses is intersected like any other

            } // end iter 'kk' every origFace->Vertex
        } // end iter 'jj' every OrigFace->Offset
    } // end iter 'ii' every OrigBrep face that was offset or copied

  // arrive here when sNoSSI contains a relationship for
  // every pair of surfaces that should already share an edge without
  // having to be intersected with one another

  // Clean up relationships.
  m_vFUToF.RemoveAll();
  m_vEUToEU.RemoveAll();
  m_vEToF.RemoveAll();
  m_vVToF.RemoveAll();
  m_vBlendToFaces.RemoveAll();

  // Shrink the geometry - this helps some for self intersections
  SER(m_pOffsetBrep->ShrinkGeometry());
  SM_DUMP_AND_ASSERT2_VALID(m_pOffsetBrep) ;

  // Now Get things ready for piecewise merge
  SmBrep *pMergedOffset = new (m_crContext) SmBrep(); NER(pMergedOffset);
  SM_OLDTOL_LINE pMergedOffset->SetTolerance(m_pOriginalBrep->GetTolerance());

  pMergedOffset->m_bEditingEnabled = TRUE;

// Remove Composites
// //cbi_CEdge11:
// #ifdef SM_NO_COMPOSITES
//   pMergedOffset->m_bMakeComposites = FALSE;
// #else
//   pMergedOffset->m_bMakeComposites = TRUE;
// #endif // SM_NO_COMPOSITES

  SmObjDelete sClean2(pMergedOffset);

#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe9 = FALSE;
  if (bDebugMe9) 
    {
      sClean2.Clear(); 
      rpResult = m_pOffsetBrep;
      m_pOffsetBrep = NULL;
      return SM_SUCCESS;
    }
#endif // SM_DEBUG_CODE

  double dTol = pMergedOffset->GetTolerance();
  double dAngleTol = 20.0*SM_PI/180.0;
  SmMerge sMerge(m_crContext,pMergedOffset,m_pOffsetBrep,dTol,dAngleTol);
  SmBrep *pResult;

#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe = FALSE;
  if (bDebugMe) 
    {
      sNoSSI.Dump() ;
    }
#endif // SM_DEBUG_CODE

  SER(sMerge.PiecewiseMerge(&sNoSSI,TRUE,TRUE,pResult));  // increments unlocked mark values

  // Note, m_pOffsetBrep still exists here.  The caller should clean it up.

  sClean2.Clear();
  rpResult = pResult;

  return SM_SUCCESS;

} // end SmOffsetExecutive::PiecewiseMerge

/*******************************************************************//********
PURPOSE: Sort an array of edgeuse by parameter value.

NOTES:
*****************************************************************************/
static void sm_SortEdgeuses
  (SmTArray<SmEdgeuse*> & rEdgeuses)  // i/o: sort edgeuses by parameter value
{
  // Do a simple bubble sort of edgeuses by minimum interval value
  SmBoolean bDone = FALSE;
  while (!bDone) 
    {
      bDone = TRUE;
      for (ULONG ii=1; ii<rEdgeuses.GetSize(); ii++) 
        {
          SmEdgeuse *pLast = rEdgeuses[ii-1];
          SmEdgeuse *pCurr = rEdgeuses[ii];
          if (pLast->GetEdge()->GetInterval().GetMin() >
              pCurr->GetEdge()->GetInterval().GetMin()) 
            {
              bDone = FALSE;
              rEdgeuses[ii-1] = pCurr;
              rEdgeuses[ii]   = pLast;
            }
        }
    }
} // end sm_SortEdgeuses

/*******************************************************************//********
PURPOSE:

NOTES:
*****************************************************************************/
static SmStatus sm_InsertBreaks
  (const SmTArray<SmEdgeuse*> & crEdgeuses,  // in :
   SmTArray<SmEdgeuse*> & rEUToSplit)        // i/o:
{
  ULONG ii, jj, kk;
  for (kk=0; kk<2; kk++) 
    {
      for (ii=0; ii<crEdgeuses.GetSize(); ii++) 
        {
          SmEdgeuse *pEU = crEdgeuses[ii];
          SmEdge *pEdge = pEU->GetEdge();
          double dParam = pEdge->GetInterval().GetMin();
          if (kk==1)
            { dParam = pEdge->GetInterval().GetMax(); }

          for (jj=0; jj<rEUToSplit.GetSize(); jj++) 
            {
              SmEdgeuse *pEUSplit = rEUToSplit[jj];
              SmEdge *pEdgeSplit = pEUSplit->GetEdge();
              SmBrep *pBrep = pEdgeSplit->GetBrep();
              SmExtent1d sIvl = pEdgeSplit->GetInterval();
              if (sIvl.ContainsValue(dParam) && !sIvl.IsValueOnBoundary(dParam)) 
                {
                  SmEdge *pNewE1 = NULL, *pNewE2 = NULL;
                  SmVertex *pNewVertex = NULL;
                  if (pBrep->MakeVertexSplitEdge(pEdgeSplit,dParam,
                      pNewE1,pNewE2,pNewVertex) != SM_SUCCESS) 
                    {
                      continue; // Can't split if too close to end point
                    }
                  SmEdgeuse *pCCW = pEUSplit->GetCCWEdgeuse();
                  SmEdgeuse *pCW = pEUSplit->GetCWEdgeuse();
                  if (   pCCW->GetEdge() == pNewE1
                      || pCCW->GetEdge() == pNewE2) 
                    {
                      rEUToSplit.Add(pCCW);
                    }
                  if (   pCW->GetEdge() == pNewE1
                      || pCW->GetEdge() == pNewE2) 
                    {
                      rEUToSplit.Add(pCW);
                    }
                  break;
                }
            }
        }
    }

  return SM_SUCCESS;

} // end sm_InsertBreaks

/*******************************************************************//********
PURPOSE:   Given an array of edgeuses (presumed) to be a sequence from a single edge,
              remove all elements less than 10% of the overall length.

NOTES:
*****************************************************************************/
static SmStatus sm_RemoveSmallStuff
  (SmTArray<SmEdgeuse*> & rEdgeuses)   // i/o: set of end-to-end edgeuses on input.
                                       //      All members less than 10% of overall length removed
                                       //      on exit.
{
    // Note that this is just a starting point and is not really the
    // best way to do it.
    SmTArray<SmEdgeuse*> sEdgeuses;
    SmTArray<double> dLengths;
    double dTotalLength = 0.0;
    ULONG ii, kk;

    // for every input edgeuse - get edgeuse->curve lengths
    for (ii=0; ii<rEdgeuses.GetSize(); ii++) 
      {
        SmEdgeuse *pEU    = rEdgeuses[ii];
        SmCurve   *pCurve = pEU->GetEdge()->GetCurve();
        double     dLeng  = pCurve->ApproximateLength(pEU->GetEdge()->GetInterval(),10);
        dLengths.Add(dLeng);
        dTotalLength += dLeng;
      }

    // for every length
    for (kk=0; kk<dLengths.GetSize(); kk++) 
      {
        // save all edgeuses more than 10% of total length
        if (dLengths[kk] > dTotalLength/10.0) 
          {
            sEdgeuses.Add(rEdgeuses[kk]);
          }
      }

    // set output
    rEdgeuses.ReSet();
    rEdgeuses.Append(sEdgeuses);
    return SM_SUCCESS;

} // end sm_RemoveSmallStuff

/*******************************************************************//********
PURPOSE: Given a list of Edgeuses along an Edge, and a list of Radial EUs
   along the same Edge, find any in the first list that do not have a match
   in the second list.

NOTES:
   The Edgeuses in each list all run along the same Edge, covering different intervals.
   For each EU in rEdgeuses, find the corresponding EU in rRadialEdgeuses.
   This is done by matching the mid-parameter if the rEdgeuses EU to the
   interval of the rRadialEdgeuses EU.  If no match is found for any rEdgeuses EU,
   add that EU to rGapEdgeuses.  rEdgeuses is returned with all of its original EUs
   that do have matches, i.e., all those that are not in rGapEdgeuses.
*****************************************************************************/
static void sm_ExtractGaps
  (SmTArray<SmEdgeuse*> & rEdgeuses,        // i/o: Edgeuses along an Edge
   SmTArray<SmEdgeuse*> & rRadialEdgeuses,  // in : Radial Edgeuses along the same Edge
   SmTArray<SmEdgeuse*> & rGapEdgeuses)     // out: Edgeuses from rEdgeuses with no match in rRadialEdgeuses
{
  SmTArray<SmEdgeuse*> sEUFound;

  ULONG ii, jj;
  for (ii=0; ii<rEdgeuses.GetSize(); ii++)
    {
      SmEdgeuse *pEU = rEdgeuses[ii];
      SmBoolean bEUFound = FALSE;
      double dMid = pEU->GetEdge()->GetInterval().Evaluate(0.5);
      for (jj=0; jj<rRadialEdgeuses.GetSize(); jj++) 
        {
          SmExtent1d sIvl = rRadialEdgeuses[jj]->GetEdge()->GetInterval();
          if (sIvl.ContainsValue(dMid)) 
            {
              sEUFound.Add(pEU);
              bEUFound = TRUE;
              break;
            }
        }
      if (!bEUFound) 
        {
          rGapEdgeuses.Add(pEU);
        }
    }

  rEdgeuses.ReSet();
  rEdgeuses.Append(sEUFound);

} // end sm_ExtractGaps

/*******************************************************************//**
PURPOSE: Get the edgeuse corresponding to the offset of the input edgeuse.

NOTES: 
***********************************************************************/
SmEdgeuse* SmOffsetExecutive::GetOffsetEdgeuse
  (SmEdgeuse *cpOrigEU)
{
   SM_PTR_ARRAY(sOffEUS,SmEdgeuse,16);
   m_vEUToEU.GetSeconds(cpOrigEU,sOffEUS);
   if (sOffEUS.GetSize() == 1) {
       return sOffEUS[0];
   }
   return NULL;
} // end GetOffsetEdgeuse

/*******************************************************************//**
PURPOSE: The original edgeuse beyond a collapsed face (an inset face whose
   offset vanished), seen from one of its edgeuses.

NOTES: pEU is an edgeuse of the collapsed face. Across a rail (a tangent
   edge of a four-sided round), the result is the radial of the opposite rail.
   Across an end edge (where the round meets a face across its end), seen
   from pFrom, one of that edge's vertices, it is the radial of the rail at
   the edge's other vertex. NULL when pEU's face did not collapse, or the
   face across collapsed too, or the face has another shape.
***********************************************************************/
SmEdgeuse * SmOffsetExecutive::GetEdgeuseAcrossCollapsedFace
  (const SmEdgeuse *pEU,      // in : edgeuse of a collapsed face
   const SmVertex  *pFrom)    // in : for an end edge, the vertex it is seen from; NULL for a rail
 const
{
  if (pEU == NULL || !IsCollapsedFace(pEU->GetFace())) { return NULL; }
  double dTanTolDeg = m_rGeometryCreation.GetTangencyTolDegrees();
  SmEdgeuse *pNext = NULL;
  if (pEU->IsTangentSector(dTanTolDeg))
    {
      SmTArray<SmEdgeuse*> sLoopEUs;
      pEU->GetLoopuse()->GetEdgeuses(sLoopEUs);
      if (sLoopEUs.GetSize() != 4) { return NULL; }
      pNext = pEU->GetCCWEdgeuse()->GetCCWEdgeuse();
    }
  else
    {
      if (pFrom == NULL) { return NULL; }
      // the end edge's vertex away from pFrom (by position: pFrom may be a vertex of the rail beside it)
      SmVertex *pV0  = pEU->GetEdge()->GetStartVertex();
      SmVertex *pV1  = pEU->GetEdge()->GetEndVertex();
      SmVertex *pFar =   pV0->GetPoint().DistanceBetween(pFrom->GetPoint()) > pV1->GetPoint().DistanceBetween(pFrom->GetPoint())
                       ? pV0
                       : pV1 ;
      SmEdgeuse *aSides[2] = { pEU->GetCWEdgeuse(), pEU->GetCCWEdgeuse() };
      for (ULONG ii=0; ii<2; ii++)
        {
          SmEdge *pE = aSides[ii]->GetEdge();
          if (aSides[ii] != pEU && (pE->GetStartVertex() == pFar || pE->GetEndVertex() == pFar))
            { pNext = aSides[ii]; break; }
        }
    }
  SmEdgeuse *pAcross = pNext ? pNext->GetRadial() : NULL;
  if (pAcross == NULL || pAcross == pNext || IsCollapsedFace(pAcross->GetFace())) { return NULL; }
  return pAcross;

} // end SmOffsetExecutive::GetEdgeuseAcrossCollapsedFace

/*******************************************************************//**
PURPOSE: Syncronize the corresponding offset edgeuses of an edgeuse.
   This is done by splitting edges on one side corresponding to breaks
   on the other side.  After all is said and done there should be two
   lists with edgeuses in each that have the same parametric ranges.

   Get list of mapped edgeuses (those in the offset brep) for both
   the input Edgeuse and its radial mate.

NOTES: Right now it assumes that total range of parameters
   for the corresponding offset edgeuses.
***********************************************************************/
SmStatus SmOffsetExecutive::SyncronizeOffsetEdgeuses
  (SmEdgeuse            * pOrigEdgeuse,          // in : one orig Brep edgeuse to synchronize
   SmEdgeuse            * pOtherEdgeuse,         // in : the orig Brep radial mate to synchronize
   SmTArray<SmEdgeuse*> & rOffsetEdgeuses,       // out: the offset edgeuses mapped to pOrigEdgeuse
   SmTArray<SmEdgeuse*> & rOtherOffsetEdgeuses,  // out: the offset edgeuses mapped to pOtherEdgeuse
   SmTArray<SmEdgeuse*> & rOffGapEdgeuses,       // out: offset Edgeuses with no match in rOtherOffsetEdgeuses
   SmTArray<SmEdgeuse*> & rOtherOffGapEdgeuses)  // out: other offset Edgeuses with no match in rOffsetEdgeuses
{
  // check input
  NER(pOrigEdgeuse);

  // init output
  rOffsetEdgeuses.ReSet();
  rOtherOffsetEdgeuses.ReSet();
  
  // locals
  SmTArray<SmEdgeuse*> sOffEUS, sOtherOffEUS;

  // get offset edgeuses mapped to OrigEdgeuse
  m_vEUToEU.GetSeconds(pOrigEdgeuse,sOffEUS);

  // remove short edgeuses - keep the long ones 
  SER(sm_RemoveSmallStuff(sOffEUS));

  // when OrigEdgeuse has no mapped offset edgeuses
  if (sOffEUS.GetSize() == 0) 
    {
      // get and sort offset edgeuses mapped to pOtherEdgeuse
      m_vEUToEU.GetSeconds(pOtherEdgeuse,rOtherOffsetEdgeuses);
      sm_SortEdgeuses(rOtherOffsetEdgeuses);

      // set both arrays to same size 
      rOffsetEdgeuses.SetSize(rOtherOffsetEdgeuses.GetSize());
      return SM_SUCCESS;
    }

  // get offset edgeuses mapped to pOtherEdgeuse
  m_vEUToEU.GetSeconds(pOtherEdgeuse,sOtherOffEUS);

  // remove short edgeuses - keep the long ones
  SER(sm_RemoveSmallStuff(sOtherOffEUS));
  
  // when pOtherEdgeuse has no mapped offset edgeuses
  if (sOtherOffEUS.GetSize() == 0) 
    {
      // place and sort offset edgeuses mapped to OrigEdgeuse in output
      rOffsetEdgeuses.Append(sOffEUS);
      sm_SortEdgeuses(rOffsetEdgeuses);

      // make both returned arrays the same size
      rOtherOffsetEdgeuses.SetSize(rOffsetEdgeuses.GetSize());
      return SM_SUCCESS;
    }
      
  // Now insert new vertices into edges to split at coresponding breaks
  SER(sm_InsertBreaks(sOffEUS,sOtherOffEUS));
  SER(sm_InsertBreaks(sOtherOffEUS,sOffEUS));

#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe = FALSE;
  if (bDebugMe) 
    {
      SmBrep *pBrep = pOrigEdgeuse->GetBrep() ;

      smgfx_Erase();
      smgfx_SetLook(1,2, 0,0,1) ; if(pBrep) pBrep->Draw(TRUE) ; sm_GraphicsLoop() ;
      smgfx_SetLook(3,4, 1,0,0) ; for (ULONG di=0; di<sOffEUS.GetSize(); di++) 
                                     {
                                       SmEdgeuse *pEU = sOffEUS[di];
                                       pEU->GetEdge()->GetInterval().Dump(); smos_WriteBuffer(_T("\n"));
                                       smgfx_ChangeColor(di!=0?TRUE:FALSE); pEU->GetEdge()->Draw(); sm_GraphicsLoop() ;
                                     }
      smgfx_SetLook(3,4, 1,0,0) ; for (ULONG dj=0; dj<sOtherOffEUS.GetSize(); dj++) 
                                     {
                                       SmEdgeuse *pEU = sOtherOffEUS[dj];
                                       pEU->GetEdge()->GetInterval().Dump(); smos_WriteBuffer(_T("\n"));
                                       smgfx_ChangeColor(dj!=0?TRUE:FALSE); pEU->GetEdge()->Draw(); sm_GraphicsLoop() ;
                                     }
    }
#endif // SM_DEBUG_CODE

  // sort offset edgeuses by parameter value
  sm_SortEdgeuses(sOffEUS);
  sm_SortEdgeuses(sOtherOffEUS);
  
  //
  sm_ExtractGaps(sOffEUS,sOtherOffEUS,rOffGapEdgeuses);
  sm_ExtractGaps(sOtherOffEUS,sOffEUS,rOtherOffGapEdgeuses);

  // set the output
  rOffsetEdgeuses.Append(sOffEUS);
  rOtherOffsetEdgeuses.Append(sOtherOffEUS);

  return SM_SUCCESS;

} // end SmOffsetExecutive::SyncronizeOffsetEdgeuses

/*******************************************************************//**
PURPOSE: Create a ruled surface between a lamina original edge
   and its offset edge.

NOTES: 
   If pEdgeuse or its offset are Null, returns SM_SUCCESS.
***********************************************************************/
SmStatus SmOffsetExecutive::BuildSideWallOfEdgeuse(SmEdgeuse *pEdgeuse)
{
    if ( pEdgeuse == NULL ) { return SM_SUCCESS; }

    SmEdgeuse *pOffEU = GetOffsetEdgeuse( pEdgeuse );
    if ( pOffEU   == NULL ) { return SM_SUCCESS; }

#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe = FALSE;
static ULONG lCount      = 1 ; lCount++ ;
static ULONG lDebugCount = 0 ;
    if (bDebugMe || lCount == lDebugCount) {
        smgfx_Erase();
        smgfx_SetLook(3,5, 1,0,0); pEdgeuse->Draw(); sm_GraphicsLoop();
        smgfx_SetLook(2,4, 1,1,0); pEdgeuse->GetFace()->Draw(); sm_GraphicsLoop();
        smgfx_SetLook(3,5, 0,0,1); pOffEU->Draw(); sm_GraphicsLoop();
        smgfx_SetLook(2,4, 0,1,1); pOffEU->GetFace()->Draw(); sm_GraphicsLoop();
        smgfx_SetLook(1,2, 0,0,0); pEdgeuse->GetEdge()->GetBrep()->Draw(); sm_GraphicsLoop();
        sm_GraphicsLoop();
    }
#endif // SM_DEBUG_CODE

    SmSurface *pNewSurface = NULL;
    SER(SmBrep::CreateEdgeEdgeBlend( m_crContext,    // in : 
                                     0,              // in : 0 = ruled surface
                                     pEdgeuse,       // in : 
                                     pOffEU,         // in : 
                                     1, NULL,        // in : 1 = Curves in same direction, NULL = no blend options
                                     pNewSurface));  // out: 

    SmFace *pNewFace = NULL;
    SER(m_pOffsetBrep->CreateFaceFromSurface(pNewSurface,
        pNewSurface->GetNaturalUVDomain(),pNewFace));

    m_sClosingFaces.Add( pNewFace );  // [B549]

#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe2 = FALSE;
    if (bDebugMe2) {
        // Note, AssertValid() will flag coincident vertices here.
        // Check the new face.
        SmTrimmingTools::CheckFace( pNewFace, pNewFace->GetTolerance() );
        smgfx_Erase();
        smgfx_SetLook(1,2, 0,0,0); m_pOffsetBrep->Draw(); sm_GraphicsLoop();
        smgfx_SetLook(3,5, 0,0,1); pNewFace->DrawUV(5,5); sm_GraphicsLoop();
        smgfx_SetLook(2,4, 0,1,1); pNewSurface->DrawUV(4,4); sm_GraphicsLoop();
        sm_GraphicsLoop();
    }
#endif // SM_DEBUG_CODE

    return SM_SUCCESS;

} // end BuildSideWallOfEdgeuse

/*******************************************************************//**
PURPOSE: Close the gaps a local operation's step walls leave in the loops
   of the faces at the walls' ends.

NOTES: A face beside a step (a side wall of the block whose top piece moved)
   had one vertex where its edge along the moved face met its edge along the
   face that stayed. The moved edge took that vertex along (ReplaceEdgeCurves
   moves an edge's vertices to its new ends), and the other edge still starts
   where it did: the loop's consecutive edgeuses no longer meet at their
   vertex. Each such gap gets the straight edge from the vertex to the edge
   end left behind, which is where the step wall's end lies, and which the
   stitch then joins to it. A gap whose straight edge leaves the face's surface
   is no step end: the operation fails.
***********************************************************************/
SmStatus SmOffsetExecutive::CloseLoopGaps()
{
  NER(m_pOffsetBrep);
  double dApproxTol = m_rGeometryCreation.GetThisApproxTol3d();
  SmTArray<SmFace*> sFaces;
  m_pOffsetBrep->GetFaces(sFaces);
  for (ULONG ii=0; ii<sFaces.GetSize(); ii++)
    {
      SmFace *pFace = sFaces[ii];
      // Each gap closed changes the loop: look again until none is left (each gap is closed once).
      ULONG lGuard = 0;
      SmBoolean bClosed = TRUE;
      while (bClosed)
        {
          bClosed = FALSE;
          if (++lGuard > 64) { SER(SM_ERR); }
          SmFaceuse *pFU = pFace->GetUpwardFaceuse(); NER(pFU);
          SmTArray<SmLoopuse*> sLoopuses;
          pFU->GetLoopuses(sLoopuses);
          for (ULONG jj=0; jj<sLoopuses.GetSize() && !bClosed; jj++)
            {
              SmTArray<SmEdgeuse*> sEUs;
              sLoopuses[jj]->GetEdgeuses(sEUs);
              ULONG lNum = sEUs.GetSize();
              for (ULONG kk=0; kk<lNum && !bClosed; kk++)
                {
                  SmEdgeuse *pEU   = sEUs[kk];
                  SmEdgeuse *pNext = sEUs[(kk+1) % lNum];
                  if (pEU == pNext) { continue; }
                  SmVertex *pV = pNext->GetVertexuse()->GetVertex();
                  if (pEU->GetMate()->GetVertexuse()->GetVertex() != pV) { continue; }
                  double dTol = 10.0 * smos_Max( smos_Max( (double)pEU->GetEdge()->GetTolerance(), (double)pNext->GetEdge()->GetTolerance() ), dApproxTol );
                  SmPoint3d sEnd, sStart;
                  SER(pEU->NormalizedEvaluate(1.0, FALSE, sEnd));
                  SER(pNext->NormalizedEvaluate(0.0, FALSE, sStart));
                  if (sEnd.DistanceBetween(sStart) <= dTol) { continue; }
                  SmPoint3d sV = pV->GetPoint();
                  if (sEnd.DistanceBetween(sV) <= dTol && sStart.DistanceBetween(sV) > dTol)
                    { SER(CloseLoopGap(pFace, pNext, TRUE, sV, dTol)); }
                  else if (sStart.DistanceBetween(sV) <= dTol && sEnd.DistanceBetween(sV) > dTol)
                    { SER(CloseLoopGap(pFace, pEU, FALSE, sV, dTol)); }
                  else
                    { SER(SM_ERR); }
                  bClosed = TRUE;
                } // end iter every pair of consecutive edgeuses
            } // end iter every loop
        } // end until no gap is left
    } // end iter every offset face

  // A step wall's corners where an edge of another face runs past them (the edge between a
  // side face's pieces, one beside the moved face and trimmed to the step, one beside the face
  // that stayed): that edge is split there, so that each part meets its own neighbour. (In a
  // local operation the only cap faces are step walls.)
  SmTArray<SmFace*> sWalls;
  m_vEToF.GetAllSeconds(sWalls);
  SmTArray<SmVertex*> sCorners;
  for (ULONG ii=0; ii<sWalls.GetSize(); ii++)
    {
      SmTArray<SmVertex*> sV;
      sWalls[ii]->GetVertices(sV);
      for (ULONG jj=0; jj<sV.GetSize(); jj++) { sCorners.AddUnique(sV[jj]); }
    }
  for (ULONG ii=0; ii<sCorners.GetSize(); ii++)
    {
      SmPoint3d sP = sCorners[ii]->GetPoint();
      SmTArray<SmEdge*> sEdges;
      m_pOffsetBrep->GetEdges(sEdges);
      for (ULONG jj=0; jj<sEdges.GetSize(); jj++)
        {
          SmEdge *pE = sEdges[jj];
          double dTol = 10.0 * smos_Max( (double)pE->GetTolerance(), dApproxTol );
          if (   pE->GetStartVertex()->GetPoint().DistanceBetween(sP) <= dTol
              || pE->GetEndVertex()->GetPoint().DistanceBetween(sP)   <= dTol)
            { continue; }
          SmBoolean bFound = FALSE;
          double dParam = 0.0, dDist = 0.0;
          SER(pE->ClosestPoint(sP, bFound, dParam, dDist));
          if (!bFound || !(dDist <= dTol)) { continue; }
          SmEdge   *pE1 = NULL;
          SmEdge   *pE2 = NULL;
          SmVertex *pNewV = NULL;
          SER(m_pOffsetBrep->MakeVertexSplitEdge(pE, dParam, pE1, pE2, pNewV));
        } // end iter every edge
    } // end iter every step wall corner
  return SM_SUCCESS;

} // end SmOffsetExecutive::CloseLoopGaps

/*******************************************************************//**
PURPOSE: Close one loop gap: pEdgeuse's start (bAtStart) or end lies away
   from its vertex at crVertexPoint. The edge's curve gets the straight
   segment from the vertex to that end joined on, and is split where they
   join: the segment becomes an edge of its own.
***********************************************************************/
SmStatus SmOffsetExecutive::CloseLoopGap
  (SmFace          * pFace,          // in : face whose loop has the gap
   SmEdgeuse       * pEdgeuse,       // in : edgeuse whose end lies away from its vertex
   SmBoolean         bAtStart,       // in : TRUE = its start, FALSE = its end
   const SmPoint3d & crVertexPoint,  // in : where its vertex is
   double            dTol)           // in : distance within which points meet
{
  NER(pFace);
  NER(pEdgeuse);
  SmEdge    *pEdge = pEdgeuse->GetEdge(); NER(pEdge);
  SmSurface *pSurf = pFace->GetSurface(); NER(pSurf);

  // The edge curve's end at the gap: its start where the edgeuse runs with it and starts there.
  SmBoolean bSame   = (pEdgeuse->GetOrientation() == SM_OT_SAME);
  ULONG     lCrvEnd = (bAtStart == bSame) ? 0 : 1;
  SmBSplineCurve *pCurve = pEdge->CreateTrimmedNURBSCurve(m_crContext);
  SmObjDelete sCleanCurve(pCurve);
  NER(pCurve);
  SmExtent1d sIvl = pCurve->GetNaturalInterval();
  SmPoint3d  sGapEnd;
  SER(pCurve->EvaluatePoint(lCrvEnd == 0 ? sIvl.GetMin() : sIvl.GetMax(), sGapEnd));

  // A vertex that moved along the edge into it (a step down: the moved face's edge pulled the
  // corner in along this one): the edge is trimmed to it.
  {
    SmSolution      sTData[16];
    SmSolutionArray sTSolutions(16, sTData);
    SER(pCurve->GlobalPointSolve(sIvl, SM_SO_MINIMIZE, crVertexPoint, dTol, NULL, NULL, SM_SR_SINGLE, sTSolutions));
    if (sTSolutions.GetSize() > 0 && sTSolutions[0].m_vStart.m_dSolutionValue <= dTol)
      {
        double dAt = sTSolutions[0].m_vStart[0];
        SmExtent1d sKeep = (lCrvEnd == 0) ? SmExtent1d(dAt, sIvl.GetMax()) : SmExtent1d(sIvl.GetMin(), dAt);
        if (!(sKeep.GetLength() > SM_EFF_ZERO * (1.0 + sIvl.GetMaxDimension())))
          { SER(SM_ERR); }
        SER(pCurve->Trim(sKeep));
        sCleanCurve.Clear();
        pEdge->SetCurve(pCurve, TRUE);
        pCurve->SetOwner(pEdge);
        pEdge->SetInterval(pCurve->GetNaturalInterval());
        return SM_SUCCESS;
      }
  }

  // The straight segment must lie on the face's surface (as extended for the offset: the face's
  // own surface ends at its old boundary).
  SmSurface  *pExtSurf = NULL;
  SmExtent2d  sExtDomain;
  SER(GetCreateExtendedSurface(pSurf, sExtDomain, pExtSurf));
  NER(pExtSurf);
  SmSolution      sData[16];
  SmSolutionArray sSolutions(16, sData);
  for (ULONG ii=1; ii<4; ii++)
    {
      double dW = ii / 4.0;
      SmPoint3d sPnt = crVertexPoint + (sGapEnd - crVertexPoint) * dW;
      SER(pExtSurf->GlobalPointSolve(sExtDomain, SM_SO_MINIMIZE, sPnt, dTol, NULL, SM_SR_SINGLE, sSolutions));
      if (sSolutions.GetSize() < 1 || !(sSolutions[0].m_vStart.m_dSolutionValue <= dTol))
        { SER(SM_ERR); }
    }

  // The segment, from the vertex to the gap end, joined at the gap end.
  SmBSplineCurve *pSegment = NULL;
  SER(SmBSplineCurve::CreateLineSegment(m_crContext, 3, crVertexPoint, sGapEnd, pSegment));
  SmObjDelete sCleanSegment(pSegment);
  NER(pSegment);
  SER(pCurve->JoinWith(lCrvEnd, pSegment, 1, &dTol));

  // Where they join: the parameter of the gap end on the joined curve.
  SmExtent1d sJoined = pCurve->GetNaturalInterval();
  SmSolution      sCData[16];
  SmSolutionArray sCSolutions(16, sCData);
  SER(pCurve->GlobalPointSolve(sJoined, SM_SO_MINIMIZE, sGapEnd, dTol, NULL, NULL, SM_SR_SINGLE, sCSolutions));
  if (sCSolutions.GetSize() < 1)
    { SER(SM_ERR); }
  double dJoint = sCSolutions[0].m_vStart[0];

  // The edge takes the joined curve (its UV trims are rebuilt from it) and is split at the joint.
  sCleanCurve.Clear();
  pEdge->SetCurve(pCurve, TRUE);
  pCurve->SetOwner(pEdge);
  pEdge->SetInterval(sJoined);
  SmEdge   *pNewEdge1 = NULL;
  SmEdge   *pNewEdge2 = NULL;
  SmVertex *pNewVertex = NULL;
  SER(m_pOffsetBrep->MakeVertexSplitEdge(pEdge, dJoint, pNewEdge1, pNewEdge2, pNewVertex));
  NER(pNewVertex);
  return SM_SUCCESS;

} // end SmOffsetExecutive::CloseLoopGap

/*******************************************************************//**
PURPOSE: Rebuild every offset face whose loops cross one another.

NOTES: SolidCraft fork. Each loop of an offset face is made of the offsets
   of its edges, each where the face's offset surface meets a neighbour's
   (BuildOffsetOfEdgeuses), and a loop is trimmed only where its own
   consecutive edges meet: nothing trims one loop of a face against another.
   An inner loop whose offset grows past the outer loop's offset crosses it.
   Shelled 1.5 mm inward, a window whose top lies 1 mm below a round gives
   the offset wall the window's offset rising above the wall's top edge,
   where the wall meets the round's offset. Such a face is not valid (SmFace
   test #38: the loops of a face are disjoint) and the merge cannot use it:
   imprinted, its loops share the crossings as vertices (test #39) and the
   merged offset encloses no region; before that, the imprint lost the end
   vertex of an edge running over a crossing (SmEdgeDefinition::Imprint,
   SM_ERR_NULL_POINTER). Each such face is rebuilt from the stretches of its
   loops that bound it (ResolveCrossingLoopsOfFace) before the merge.
***********************************************************************/
SmStatus SmOffsetExecutive::ResolveCrossingLoops()
{
  NER(m_pOffsetBrep);
  SmTArray<SmFace*> sFaces;
  m_pOffsetBrep->GetFaces(sFaces);
  for (ULONG ii=0; ii<sFaces.GetSize(); ii++)
    { SER(ResolveCrossingLoopsOfFace(sFaces[ii])); }
  return SM_SUCCESS;

} // end SmOffsetExecutive::ResolveCrossingLoops

namespace
{
  // A stretch of a loop of an offset face between two of its crossings with other loops or the
  // ends of its edgeuse, in the loop's direction.
  struct sm_LoopStretch
  {
    ULONG            lLoop;   // its loop: 0 = the outer loop
    SmEdgeuse      * pEU;     // the edgeuse (of the face's SM_OT_SAME faceuse) it is a stretch of
    SmBSplineCurve * pUV;     // pEU's UV trim curve
    double           dLo;     // its edge-parameter interval
    double           dHi;
    long             lStart;  // the node it starts at: a crossing (>= 0), or the vertex ending
    long             lEnd;    //   edgeuse k of the face (-1 - k); and the node it ends at
    SmBoolean        bKeep;   // TRUE = it bounds the face
    SmPoint3d        sMid;    // its middle point
  };

  // A crossing of the loops of an offset face that ResolveCrossingLoopsOfFace() does not resolve:
  // reported, and the offset fails with SM_ERR_BAD_INTERSECTIONS rather than merge a face whose
  // loops cross.
  SmStatus sm_UnresolvedCrossing(const TCHAR *cpWhy)
  {
    smos_ErrorMessage(SM_ERR_BAD_INTERSECTIONS, FILE_NAME, LINE_NUMBER, cpWhy, (TCHAR*)0, FUNC_NAME);
    return SM_ERR_BAD_INTERSECTIONS;
  }
}

/*******************************************************************//**
PURPOSE: When the loops of an offset face cross one another, replace the
   face with the face(s) of its surface inside the outer loop and outside
   every inner one, bounded by the stretches of the loops between their
   crossings.

NOTES: SolidCraft fork. SMLib has nothing that resolves crossing loops: the
   face builders and trimming tools require disjoint loops
   (SmBrep::CreateFacesWithCurves: "does not allow for loops to touch or
   cross each other"), and the healer's loop fixes reorder, reorient or split
   off loops, never cut them. So:
   - The crossings are found in the face's UV space, as SmFace test #38
     finds them (GlobalCurveIntersect of the UV trim curves, at the
     surface's tolerance mapped to UV).
   - Every loop is cut at its crossings. A stretch of the outer loop bounds
     the face where it lies outside every inner loop; a stretch of an inner
     loop where it lies inside the outer loop and outside the other inner
     ones (SmLoop::ContainsUVPoint at the stretch's middle).
   - The kept stretches, in their loops' directions, close into loops at
     the crossings: at a crossing exactly one kept stretch arrives and one
     leaves. A new loop holding a stretch of the old outer loop bounds the
     face from outside; the others are holes.
   - The face is rebuilt (MakeFaceWithCurves) from each stretch's edge curve
     and UV trim curve trimmed to it (they share the edge's parameter), on a
     copy of its surface. When it has come apart (several outer loops),
     SmBrep::MakeFaceFromLoopOnFace, as the healer does, gives each outer
     loop a face of its own with the holes inside it.
   - The new faces, edges and vertices take the attributes of the face,
     edges and vertices they come from, and the face's places in the offset
     maps.
   A face whose loops cross in a way this does not resolve is not merged:
   SM_ERR_BAD_INTERSECTIONS, with the reason. Those ways: loops meeting at a
   loop vertex or running along one another, three loops through one point,
   a stretch neither inside nor outside a loop, a crossing where the kept
   stretches do not continue one to one, and a face whose edges bound other
   faces too.
***********************************************************************/
SmStatus SmOffsetExecutive::ResolveCrossingLoopsOfFace(SmFace *pFace)
{
  NER(pFace);
  SmSurface *pSurface = pFace->GetSurface(); NER(pSurface);
  SmFaceuse *pFU = pFace->GetUpwardFaceuse(); NER(pFU);
  if (pFU->GetOrientation() != SM_OT_SAME) { pFU = pFU->GetMate(); NER(pFU); }

  // The face's edge loops (the outer loop first) and the points of its vertex loops.
  SmTArray<SmLoopuse*> sLoopuses;
  pFU->GetLoopuses(sLoopuses);
  std::vector<SmLoopuse*>                vLoops;
  std::vector< std::vector<SmEdgeuse*> > vLoopEUs;
  SmTArray<SmPoint3d>                    sLoopPoints;
  for (ULONG ii=0; ii<sLoopuses.GetSize(); ii++)
    {
      SmLoopuse *pLU = sLoopuses[ii];
      if (pLU->IsVertexLoopuse())
        {
          SmTArray<SmVertexuse*> sVUs;
          pLU->GetVertexuses(sVUs);
          sLoopPoints.Add(sVUs[0]->GetVertex()->GetPoint());
          continue;
        }
      SmTArray<SmEdgeuse*> sEUs;
      pLU->GetEdgeuses(sEUs);
      vLoops.push_back(pLU);
      vLoopEUs.push_back(std::vector<SmEdgeuse*>());
      for (ULONG kk=0; kk<sEUs.GetSize(); kk++) { vLoopEUs.back().push_back(sEUs[kk]); }
    }
  if (vLoops.size() < 2) { return SM_SUCCESS; }

  // Every edgeuse, in loop order, with its UV trim curve.
  std::vector<ULONG>           vLoopOf;
  std::vector<SmEdgeuse*>      vEU;
  std::vector<SmBSplineCurve*> vUV;
  std::vector<size_t>          vFirstOfLoop;
  double dTol3d = m_rGeometryCreation.GetThisApproxTol3d();
  for (ULONG ll=0; ll<vLoops.size(); ll++)
    {
      vFirstOfLoop.push_back(vEU.size());
      for (size_t kk=0; kk<vLoopEUs[ll].size(); kk++)
        {
          SmEdgeuse *pEU = vLoopEUs[ll][kk];
          SmBSplineCurve *pUV = NULL;
          SER(pEU->GetOrCreateUVTrimCurve(pUV)); NER(pUV);
          AERN(pUV->GetNaturalInterval().AreEqual(pEU->GetEdge()->GetInterval()), SM_ERR_ASSERT_FAILURE);
          vLoopOf.push_back(ll);
          vEU.push_back(pEU);
          vUV.push_back(pUV);
          dTol3d = smos_Max(dTol3d, (double)pEU->GetEdge()->GetTolerance());
        }
    }
  vFirstOfLoop.push_back(vEU.size());
  dTol3d *= 10.0;

  // The crossings of edgeuses of different loops, in UV space (as SmFace test #38 finds them).
  SmZoneTol3d sZoneTol3d = SmTol::GetZoneTol3d(pSurface);
  SmTol2d     sUVTol     = SmTol::MapTo2d(sZoneTol3d, pSurface->GetNaturalUVDomain().GetMid(), *pSurface);
  std::vector< std::vector< std::pair<double,long> > > vCuts(vEU.size()); // per edgeuse: (parameter, crossing)
  long lNumCrossings = 0;
  for (size_t aa=0; aa<vEU.size(); aa++)
    {
      for (size_t bb=vFirstOfLoop[vLoopOf[aa]+1]; bb<vEU.size(); bb++)
        {
          SmSolutionArray sSols;
          SER(vUV[aa]->GlobalCurveIntersect(vUV[aa]->GetNaturalInterval(), *vUV[bb], vUV[bb]->GetNaturalInterval(), sUVTol, sSols));
          for (ULONG ss=0; ss<sSols.GetSize(); ss++)
            {
              if (sSols[ss].m_eSolutionType != SM_ST_SINGLE_VALUE)
                { return sm_UnresolvedCrossing(_T("Offset face loops run along one another")); }
              double dTa = sSols[ss].m_vStart[0];
              double dTb = sSols[ss].m_vStart[1];
              SmPoint3d sPa, sPb;
              SER(vEU[aa]->GetEdge()->GetCurve()->EvaluatePoint(dTa, sPa));
              SER(vEU[bb]->GetEdge()->GetCurve()->EvaluatePoint(dTb, sPb));
              if (   sPa.DistanceBetween(vEU[aa]->GetEdge()->GetStartVertex()->GetPoint()) <= dTol3d
                  || sPa.DistanceBetween(vEU[aa]->GetEdge()->GetEndVertex()->GetPoint())   <= dTol3d
                  || sPb.DistanceBetween(vEU[bb]->GetEdge()->GetStartVertex()->GetPoint()) <= dTol3d
                  || sPb.DistanceBetween(vEU[bb]->GetEdge()->GetEndVertex()->GetPoint())   <= dTol3d)
                { return sm_UnresolvedCrossing(_T("Offset face loops meet at a loop vertex")); }
              vCuts[aa].push_back(std::make_pair(dTa, lNumCrossings));
              vCuts[bb].push_back(std::make_pair(dTb, lNumCrossings));
              lNumCrossings++;
            }
        }
    }
  if (lNumCrossings == 0) { return SM_SUCCESS; }

  // Until the merge every offset face is a sheet of its own: its edges bound it alone.
  for (size_t gg=0; gg<vEU.size(); gg++)
    {
      if (!vEU[gg]->GetEdge()->IsLamina())
        { return sm_UnresolvedCrossing(_T("Offset face whose loops cross shares an edge with another face")); }
    }

  // Cut every edgeuse at its crossings into stretches, in the loop's direction.
  std::vector<sm_LoopStretch> vStretches;
  for (size_t gg=0; gg<vEU.size(); gg++)
    {
      SmEdgeuse *pEU   = vEU[gg];
      SmCurve   *pCrv  = pEU->GetEdge()->GetCurve();
      SmBoolean  bSame = (pEU->GetOrientation() == SM_OT_SAME);
      SmExtent1d sIvl  = pEU->GetEdge()->GetInterval();
      std::vector< std::pair<double,long> > &rCuts = vCuts[gg];
      std::sort(rCuts.begin(), rCuts.end());
      if (!bSame) { std::reverse(rCuts.begin(), rCuts.end()); }
      ULONG  lLoop = vLoopOf[gg];
      size_t lPrev = (gg == vFirstOfLoop[lLoop]) ? vFirstOfLoop[lLoop+1] - 1 : gg - 1;
      double dFrom = bSame ? sIvl.GetMin() : sIvl.GetMax();
      long   lFrom = -1 - (long)lPrev;
      for (size_t cc=0; cc<=rCuts.size(); cc++)
        {
          double dTo = (cc < rCuts.size()) ? rCuts[cc].first  : (bSame ? sIvl.GetMax() : sIvl.GetMin());
          long   lTo = (cc < rCuts.size()) ? rCuts[cc].second : -1 - (long)gg;
          sm_LoopStretch sS;
          sS.lLoop  = lLoop;
          sS.pEU    = pEU;
          sS.pUV    = vUV[gg];
          sS.dLo    = smos_Min(dFrom, dTo);
          sS.dHi    = smos_Max(dFrom, dTo);
          sS.lStart = lFrom;
          sS.lEnd   = lTo;
          sS.bKeep  = FALSE;
          SmPoint3d sFromPt, sToPt;
          SER(pCrv->EvaluatePoint(sS.dLo, sFromPt));
          SER(pCrv->EvaluatePoint(sS.dHi, sToPt));
          if (!(sS.dHi > sS.dLo) || sFromPt.DistanceBetween(sToPt) <= dTol3d)
            { return sm_UnresolvedCrossing(_T("Three offset face loops cross at one point")); }
          SER(pCrv->EvaluatePoint(0.5*(sS.dLo + sS.dHi), sS.sMid));
          vStretches.push_back(sS);
          dFrom = dTo;
          lFrom = lTo;
        }
    }

  // Which stretches bound the face: of the outer loop, those outside every inner loop; of an
  // inner loop, those inside the outer loop and outside the other inner ones.
  for (size_t ss=0; ss<vStretches.size(); ss++)
    {
      sm_LoopStretch &rS = vStretches[ss];
      SmPoint3d sUV3;
      SER(rS.pUV->EvaluatePoint(0.5*(rS.dLo + rS.dHi), sUV3));
      SmPoint2d sUV(sUV3.x, sUV3.y);
      rS.bKeep = TRUE;
      for (ULONG ll=0; ll<vLoops.size(); ll++)
        {
          if (ll == rS.lLoop) { continue; }
          SmEdge      *pOnEdge  = NULL;
          double       dOnParam = 0.0;
          SmOrientType eOrient  = SM_OT_UNKNOWN;
          SmPointObjectContainmentType eIn = vLoops[ll]->GetLoop()->ContainsUVPoint(sUV, sZoneTol3d, pOnEdge, dOnParam, eOrient);
          if (eIn != SM_POC_INSIDE && eIn != SM_POC_OUTSIDE)
            { return sm_UnresolvedCrossing(_T("Offset face loop stretch is neither inside nor outside another loop")); }
          if (eIn != ((ll == 0) ? SM_POC_INSIDE : SM_POC_OUTSIDE)) { rS.bKeep = FALSE; }
        }
    }

  // Close the kept stretches into loops: each goes on to the one kept stretch that starts where
  // it ends (at a crossing, the stretch of the other loop that bounds the face).
  std::vector<long> vNext(vStretches.size(), -1);
  std::vector<long> vPrevCount(vStretches.size(), 0);
  for (size_t ss=0; ss<vStretches.size(); ss++)
    {
      if (!vStretches[ss].bKeep) { continue; }
      for (size_t tt=0; tt<vStretches.size(); tt++)
        {
          if (!vStretches[tt].bKeep || vStretches[tt].lStart != vStretches[ss].lEnd) { continue; }
          if (vNext[ss] >= 0) { return sm_UnresolvedCrossing(_T("Offset face loops leave a crossing more than one way")); }
          vNext[ss] = (long)tt;
          vPrevCount[tt]++;
        }
      if (vNext[ss] < 0) { return sm_UnresolvedCrossing(_T("Offset face loops leave a crossing no way")); }
    }
  std::vector< std::vector<size_t> > vNewLoops;
  std::vector<SmBoolean>             vNewLoopIsOuter;
  std::vector<SmBoolean>             vDone(vStretches.size(), FALSE);
  for (size_t ss=0; ss<vStretches.size(); ss++)
    {
      if (!vStretches[ss].bKeep || vDone[ss]) { continue; }
      std::vector<size_t> vLoop;
      SmBoolean bOuter = FALSE;
      for (size_t tt=ss; !vDone[tt]; tt=(size_t)vNext[tt])
        {
          if (vPrevCount[tt] != 1) { return sm_UnresolvedCrossing(_T("Offset face loops arrive at a crossing more than one way")); }
          vDone[tt] = TRUE;
          vLoop.push_back(tt);
          if (vStretches[tt].lLoop == 0) { bOuter = TRUE; }
        }
      vNewLoops.push_back(vLoop);
      vNewLoopIsOuter.push_back(bOuter);
    }

  // Outer loops first (MakeFaceWithCurves takes the first loop for the outer one).
  std::vector<size_t> vOrder;
  for (size_t ll=0; ll<vNewLoops.size(); ll++) { if ( vNewLoopIsOuter[ll]) { vOrder.push_back(ll); } }
  size_t lNumOuter = vOrder.size();
  for (size_t ll=0; ll<vNewLoops.size(); ll++) { if (!vNewLoopIsOuter[ll]) { vOrder.push_back(ll); } }

  // The face(s) the old one becomes; none when no stretch of its outer loop bounds it.
  std::vector<SmFace*> vNewFaces;
  if (lNumOuter > 0)
    {
      SmTArray<ULONG>           sLoopCounts;
      SmTArray<SmCurve*>        s3DCurves;
      SmTArray<SmBSplineCurve*> sUVCurves;
      SmTArray<SmOrientType>    sOrients;
      SmObjsDelete<SmCurve*>        sClean3D(&s3DCurves);
      SmObjsDelete<SmBSplineCurve*> sCleanUV(&sUVCurves);
      for (size_t oo=0; oo<vOrder.size(); oo++)
        {
          const std::vector<size_t> &rLoop = vNewLoops[vOrder[oo]];
          sLoopCounts.Add((ULONG)rLoop.size());
          for (size_t kk=0; kk<rLoop.size(); kk++)
            {
              const sm_LoopStretch &rS = vStretches[rLoop[kk]];
              SmCurve *p3D = NULL;
              SER(rS.pEU->GetEdge()->GetCurve()->Copy(m_crContext, p3D)); NER(p3D);
              s3DCurves.Add(p3D);
              SmExtent1d sIvl3D(rS.dLo, rS.dHi);
              SER(p3D->Trim(sIvl3D));
              SmCurve *pUVCopy = NULL;
              SER(rS.pUV->Copy(m_crContext, pUVCopy)); NER(pUVCopy);
              SmBSplineCurve *pUV = SM_CAST_PTR(SmBSplineCurve, pUVCopy);
              if (pUV == NULL) { delete pUVCopy; SER(SM_ERR_ASSERT_FAILURE); }
              sUVCurves.Add(pUV);
              SmExtent1d sIvlUV(rS.dLo, rS.dHi);
              SER(pUV->Trim(sIvlUV));
              sOrients.Add(rS.pEU->GetOrientation());
            }
        }
      SmSurface *pNewSurface = NULL;
      SER(pSurface->Copy(m_crContext, pNewSurface)); NER(pNewSurface);
      SmRegion *pNewRegion = NULL;
      SmShell  *pNewShell  = NULL;
      SmFace   *pNewFace   = NULL;
      sClean3D.Clear(); // the new face owns its curves and its surface
      sCleanUV.Clear();
      SER(m_pOffsetBrep->MakeFaceWithCurves(m_pOffsetBrep->GetInfiniteRegion(), sLoopCounts, &s3DCurves, &sUVCurves, sOrients,
                                            sLoopPoints, pNewSurface, pFace->GetUVDomain(), SM_OT_SAME,
                                            pNewRegion, pNewShell, pNewFace));
      NER(pNewFace);
      vNewFaces.push_back(pNewFace);

      // A face come apart: each further outer loop a face of its own, with the holes inside it.
      if (lNumOuter > 1)
        {
          SmTArray<SmLoop*> sNewLoops;
          pNewFace->GetLoops(sNewLoops);
          for (ULONG ll=1; ll<sNewLoops.GetSize(); ll++)
            {
              // an outer loop holds a stretch of the old outer loop
              SmTArray<SmEdge*> sLoopEdges;
              sNewLoops[ll]->GetEdges(sLoopEdges);
              SmBoolean bOuter = FALSE;
              for (ULONG ee=0; ee<sLoopEdges.GetSize() && !bOuter; ee++)
                {
                  for (size_t ss=0; ss<vStretches.size() && !bOuter; ss++)
                    {
                      if (!vStretches[ss].bKeep || vStretches[ss].lLoop != 0) { continue; }
                      SmBoolean bFound = FALSE;
                      double dParam = 0.0, dDist = 0.0;
                      SER(sLoopEdges[ee]->ClosestPoint(vStretches[ss].sMid, bFound, dParam, dDist));
                      bOuter = (bFound && dDist <= dTol3d);
                    }
                }
              if (!bOuter) { continue; }
              SmFace *pSplitFace = NULL;
              SER(m_pOffsetBrep->MakeFaceFromLoopOnFace(pNewFace, sNewLoops[ll], pSplitFace));
              NER(pSplitFace);
              vNewFaces.push_back(pSplitFace);
            }
          AERN(vNewFaces.size() == lNumOuter, SM_ERR_ASSERT_FAILURE);
          // MakeFaceFromLoopOnFace leaves the new faces in the old one's shell, which they do not
          // connect: each offset face is a sheet of its own until the merge.
          SER(m_pOffsetBrep->SplitDisconnectedShells());
        }
    }

  // The new faces, edges and vertices take the attributes of the face, edges and vertices they
  // come from, and the old edgeuses' places in the edgeuse map go to the new ones.
  for (size_t ff=0; ff<vNewFaces.size(); ff++)
    {
      SmFace *pNewFace = vNewFaces[ff];
      pFace->Notify(SM_NO_COPY, pNewFace, SM_NO_GET_BREP(pNewFace), SM_NO_GET_BREP(pFace));
      SmTArray<SmEdge*> sNewEdges;
      pNewFace->GetEdges(sNewEdges);
      for (ULONG ee=0; ee<sNewEdges.GetSize(); ee++)
        {
          SmEdge *pNewEdge = sNewEdges[ee];
          // the stretch whose middle lies on it
          const sm_LoopStretch *cpFrom = NULL;
          for (size_t ss=0; ss<vStretches.size() && cpFrom == NULL; ss++)
            {
              if (!vStretches[ss].bKeep) { continue; }
              SmBoolean bFound = FALSE;
              double dParam = 0.0, dDist = 0.0;
              SER(pNewEdge->ClosestPoint(vStretches[ss].sMid, bFound, dParam, dDist));
              if (bFound && dDist <= dTol3d) { cpFrom = &vStretches[ss]; }
            }
          NER(cpFrom);
          SmEdge *pOldEdge = cpFrom->pEU->GetEdge();
          pOldEdge->Notify(SM_NO_COPY, pNewEdge, SM_NO_GET_BREP(pNewEdge), SM_NO_GET_BREP(pOldEdge));
          SmVertex *apNew[2] = { pNewEdge->GetStartVertex(), pNewEdge->GetEndVertex() };
          SmVertex *apOld[2] = { pOldEdge->GetStartVertex(), pOldEdge->GetEndVertex() };
          for (int nn=0; nn<2; nn++)
            for (int oo=0; oo<2; oo++)
              if (apNew[nn]->GetPoint().DistanceBetween(apOld[oo]->GetPoint()) <= dTol3d)
                { apOld[oo]->Notify(SM_NO_COPY, apNew[nn], SM_NO_GET_BREP(apNew[nn]), SM_NO_GET_BREP(apOld[oo])); }
          SmTArray<SmEdgeuse*> sNewEUs;
          SER(pNewEdge->GetEdgeuses(sNewEUs));
          SmEdgeuse *apOldEU[2] = { cpFrom->pEU, cpFrom->pEU->GetMate() };
          for (int oo=0; oo<2; oo++)
            {
              SmTArray<SmEdgeuse*> sOrigEUs;
              m_vEUToEU.GetFirsts(apOldEU[oo], sOrigEUs);
              for (ULONG nn=0; nn<sNewEUs.GetSize(); nn++)
                {
                  if (sNewEUs[nn]->GetFaceuse()->GetOrientation() != apOldEU[oo]->GetFaceuse()->GetOrientation()) { continue; }
                  for (ULONG pp=0; pp<sOrigEUs.GetSize(); pp++) { SER(m_vEUToEU.RelatePair(sOrigEUs[pp], sNewEUs[nn])); }
                }
            }
        }
    }

  // The new faces take the old one's places in the face maps; the old face, its edgeuses and
  // their pairs go.
  SmTArray<SmFaceuse*> sFUs;
  m_vFUToF.GetFirsts(pFace, sFUs);
  for (ULONG ii=0; ii<sFUs.GetSize(); ii++)
    {
      SER(m_vFUToF.RemovePair(sFUs[ii], pFace));
      for (size_t ff=0; ff<vNewFaces.size(); ff++) { SER(m_vFUToF.RelatePair(sFUs[ii], vNewFaces[ff])); }
    }
  SmTArray<SmEdge*> sEs;
  m_vEToF.GetFirsts(pFace, sEs);
  for (ULONG ii=0; ii<sEs.GetSize(); ii++)
    {
      SER(m_vEToF.RemovePair(sEs[ii], pFace));
      for (size_t ff=0; ff<vNewFaces.size(); ff++) { SER(m_vEToF.RelatePair(sEs[ii], vNewFaces[ff])); }
    }
  SmTArray<SmVertex*> sVs;
  m_vVToF.GetFirsts(pFace, sVs);
  for (ULONG ii=0; ii<sVs.GetSize(); ii++)
    {
      SER(m_vVToF.RemovePair(sVs[ii], pFace));
      for (size_t ff=0; ff<vNewFaces.size(); ff++) { SER(m_vVToF.RelatePair(sVs[ii], vNewFaces[ff])); }
    }
  SmTArray<SmFace*> sBlend;
  m_vBlendToFaces.GetSeconds(pFace, sBlend);
  for (ULONG ii=0; ii<sBlend.GetSize(); ii++)
    {
      SER(m_vBlendToFaces.RemovePair(pFace, sBlend[ii]));
      for (size_t ff=0; ff<vNewFaces.size(); ff++) { SER(m_vBlendToFaces.RelatePair(vNewFaces[ff], sBlend[ii])); }
    }
  m_vBlendToFaces.GetFirsts(pFace, sBlend);
  for (ULONG ii=0; ii<sBlend.GetSize(); ii++)
    {
      SER(m_vBlendToFaces.RemovePair(sBlend[ii], pFace));
      for (size_t ff=0; ff<vNewFaces.size(); ff++) { SER(m_vBlendToFaces.RelatePair(sBlend[ii], vNewFaces[ff])); }
    }
  for (size_t gg=0; gg<vEU.size(); gg++)
    {
      SmEdgeuse *apOldEU[2] = { vEU[gg], vEU[gg]->GetMate() };
      for (int oo=0; oo<2; oo++)
        {
          SmTArray<SmEdgeuse*> sOrigEUs;
          m_vEUToEU.GetFirsts(apOldEU[oo], sOrigEUs);
          for (ULONG pp=0; pp<sOrigEUs.GetSize(); pp++) { SER(m_vEUToEU.RemovePair(sOrigEUs[pp], apOldEU[oo])); }
        }
    }
  SER(m_pOffsetBrep->DeleteFace(pFace, TRUE, TRUE));
  return SM_SUCCESS;

} // end SmOffsetExecutive::ResolveCrossingLoopsOfFace

/*******************************************************************//**
PURPOSE: Build the offset of an edgeuse.

NOTES: Make cap faces to fill in gaps between offset faces.
***********************************************************************/
SmStatus SmOffsetExecutive::BuildOffsetOfEdgeuses
  (SmEdgeuse *pEdgeuse,       // in : target edgeuse that has been offset
   SmEdgeuse *pOtherEdgeuse,  // in : radial partner to target edgeuse that has been offset
   SmMarkType eMarkType)      // in : uses without increment eMarkType value
{
  // locals
  SmSSIData sSSIData;
  SmEdgeuse * sEU1Data[16];
  SmEdgeuse * sEU2Data[16];
  SmEdgeuse * sEU3Data[16];
  SmEdgeuse * sEU4Data[16];
  SmTArray<SmEdgeuse*>   sOffsetEUS(16,sEU1Data);
  SmTArray<SmEdgeuse*>   sOtherEUS(16,sEU2Data);
  SmTArray<SmEdgeuse*>   sOffGapEUS(16,sEU3Data);
  SmTArray<SmEdgeuse*>   sOtherOffGapEUS(16,sEU4Data);
  SmTArray<SmCurve*>     s3DTrimCurves;
  SmTArray<SmOrientType> sTrimOrientations;
  SmTArray<SmFace*>      sNewFaces;

  // get offset edge and its two faces (manifold assumption)
  SmEdge *pEdge      = pEdgeuse->GetEdge(); NER(pEdge);
  SmFace *pFace      = pEdgeuse->GetFace();
  SmFace *pOtherFace = pOtherEdgeuse->GetFace();
  
#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe = FALSE;
static ULONG lCount      = 1 ; lCount++ ;
static ULONG lDebugCount = 0 ;
  if (bDebugMe || lCount == lDebugCount) 
    {
      m_pOriginalBrep->Dump() ;  SM_ASSERT_VALID(m_pOriginalBrep) ;
      m_pOffsetBrep->Dump();     SM_ASSERT_VALID(m_pOffsetBrep) ;
      //SmBrep *pTgt2Brep = pEdge->GetBrep() ;

      smgfx_Erase();
      smgfx_SetLook(1,2, 0,0,1) ; m_pOriginalBrep->Draw(TRUE); sm_GraphicsLoop();
      smgfx_SetLook(2,3, 0,1,0) ; m_pOffsetBrep->Draw(TRUE); sm_GraphicsLoop();
      smgfx_SetLook(4,5, 1,0,1) ; if(pEdgeuse) pEdgeuse->Draw() ; sm_GraphicsLoop() ;
      smgfx_SetLook(4,5, 1,1,0) ; if(pOtherEdgeuse) pOtherEdgeuse->Draw() ; sm_GraphicsLoop() ;
      smgfx_SetLook(5,6, 1,0,0) ; pEdge->Draw() ; sm_GraphicsLoop() ;
      smgfx_SetLook(1,2, 0,1,1) ; pFace->Draw(SM_DM_CROSSHATCH); sm_GraphicsLoop() ;
      smgfx_SetLook(1,2, 0,1,1) ; pOtherFace->Draw(SM_DM_CROSSHATCH); sm_GraphicsLoop() ;
      sm_GraphicsLoop();
  }
#endif // SM_DEBUG_CODE

  // If both are shell faces (only being copied not offset) 
  // then don't need to do anything here
  ULONG lIndex;
  if (   m_pShellFaces 
      && m_pShellFaces->FindElement(pFace,lIndex) 
      && m_pShellFaces->FindElement(pOtherFace,lIndex)) 
    {
      // and no extension of convex edges
      if (!this->m_bExtendConvexEdges) 
        {
          return SM_SUCCESS;
        }
    }

  // get mapped offset edgeuses for pEdgeuse and pOtherEdgeuse
  SER(SyncronizeOffsetEdgeuses(pEdgeuse,   pOtherEdgeuse,
                               sOffsetEUS, sOtherEUS, 
                               sOffGapEUS, sOtherOffGapEUS));

  // Note: this doesn't help anything. [cbi B01]
  //// let Lamina edge->pOtherEdgeuse == NULL
  //if ( pEdgeuse->GetEdge()->IsLamina() )
  //  { 
  //    pOtherEdgeuse = NULL;      //RCLxx this causes trouble in EdgeuseOffset below
  //  }

  // for every mapped offset edgeuse pair - see if we need to build a capFace
  ULONG ii, lNumOffEUs = sOffsetEUS.GetSize();
  for (ii=0; ii<lNumOffEUs; ii++) 
    {
      SmEdgeuse *pOffsetEU   = sOffsetEUS[ii];
      SmEdgeuse *pOtherOffEU = sOtherEUS[ii];
//        if (!pOffsetEU || !pOtherOffEU) {
//            if (m_bExtendConvexEdges) {
//                SER(SM_ERR);
//            }
//            continue;
//        }

      // clear lamina other offset edges
      if (pEdgeuse->GetEdge()->IsLamina()) 
        {
          pOtherOffEU = NULL;
        }

      if ( pOffsetEU != NULL )
        { sm_NotifyEdgeAndVertices( pEdgeuse->GetEdge(), pOffsetEU->GetEdge(), SM_NO_COPY ); }

      // build Offset CapSurface to bridge any gap made by offseting the pEdge Faces
      SmSurface *pNewSurface = NULL ;
      SmBoolean  bClosingWall = FALSE ;
      SER(m_rGeometryCreation.EdgeuseOffset(m_crContext,       // in : context for new object construction
                                            this,              // in : offset executive managing this offset
                                            pEdgeuse,          // in : edgeuse from original brep
                                            pOtherEdgeuse,     // in : radial mate of edgeuse from original brep
                                            pOffsetEU,         // in : edgeuse in offset brep mapped to OriginalEdgeuse
                                            pOtherOffEU,       // in : edgeuse in offset brep mapped to OtherOriginalEdgeuse
                                            pNewSurface,       // out: 
                                            s3DTrimCurves,     // out: 
                                            sTrimOrientations, // out: 
                                            sSSIData,          // out: currently not used 
                                            eMarkType,        // in : uses without increment eMarkType Value
                                            &bClosingWall));   // out: TRUE = a step wall, never extended
                                                               // note: uses wihtout increment eMarkType value
      if (!pNewSurface) continue;
      SmFace *pNewFace = NULL;

#ifdef SM_DEBUG_CODE
      if (bDebugMe || lCount == lDebugCount) 
        {
          m_pOriginalBrep->Dump() ;  SM_ASSERT_VALID(m_pOriginalBrep) ;
          m_pOffsetBrep->Dump();     SM_ASSERT_VALID(m_pOffsetBrep) ;
          for(ULONG di=0;di<s3DTrimCurves.GetSize();di++)
            {
              SmCurve *pCurve = s3DTrimCurves[di] ;
              pCurve->Dump() ;
            }

          smgfx_Erase();
          smgfx_SetLook(1,2, 0,0,1) ; m_pOriginalBrep->Draw(TRUE); sm_GraphicsLoop();
          smgfx_SetLook(2,3, 0,1,0) ; m_pOffsetBrep->Draw(TRUE); sm_GraphicsLoop();
          smgfx_SetLook(3,4, 1, 0, 0); if(pEdgeuse)      pEdgeuse->Draw();      sm_GraphicsLoop();
          smgfx_SetLook(4,5, 1, 0, 1); if(pOtherEdgeuse) pOtherEdgeuse->Draw(); sm_GraphicsLoop();
          smgfx_SetLook(5,6, 1,.5, 0); if(pOffsetEU)     pOffsetEU->Draw();     sm_GraphicsLoop();
          smgfx_SetLook(6,7, 1, 0,.5); if(pOtherOffEU)   pOtherOffEU->Draw();   sm_GraphicsLoop();

          smgfx_SetLook(1,2, 0,1,1); if(pEdgeuse)      pEdgeuse->GetFace()->Draw(SM_DM_CROSSHATCH);      sm_GraphicsLoop();
          smgfx_SetLook(1,2, 0,1,1); if(pOtherEdgeuse) pOtherEdgeuse->GetFace()->Draw(SM_DM_CROSSHATCH); sm_GraphicsLoop();
          smgfx_SetLook(1,2, 1,1,0); if(pOffsetEU)     pOffsetEU->GetFace()->Draw(SM_DM_CROSSHATCH);     sm_GraphicsLoop();
          smgfx_SetLook(1,2, 1,1,0); if(pOtherOffEU)   pOtherOffEU->GetFace()->Draw(SM_DM_CROSSHATCH);   sm_GraphicsLoop();

          smgfx_SetLook(1,2, 1,0,0) ; if(pNewSurface) pNewSurface->DrawUV(); sm_GraphicsLoop() ;
          sm_GraphicsLoop();
        }
#endif // SM_DEBUG_CODE

      // Build a face in m_pOffsetBrep for the new CapSurface
      if (s3DTrimCurves.GetSize() == 0) 
        {
          SER(m_pOffsetBrep->CreateFaceFromSurface(pNewSurface,
                                                   pNewSurface->GetNaturalUVDomain(),
                                                   pNewFace));
          // SolidCraft: a step wall spans the gap between two edges exactly, as a lamina side wall
          // does, and is kept out of the face extensions as one: extending the ruled wall along its
          // rulings can carry it through its own axis (a step round a circle), where it folds over
          // itself and its edges drop onto the folded sheet.
          if (bClosingWall)
            { m_sClosingFaces.Add( pNewFace ); }
        }
      else // build a trimmed face 
        {
          SmTArray<ULONG> sCurveLoops;
          sCurveLoops.Add(s3DTrimCurves.GetSize());
          SmTArray<SmPoint3d> sLoopPoints;
          SmRegion *pNewRegion = NULL;
          SmShell *pNewShell = NULL;
          SER(m_pOffsetBrep->MakeFaceWithCurves(m_pOffsetBrep->GetInfiniteRegion(), // in : region to contain new topology objects
                                                sCurveLoops,                        // in : 1 entry per loop, value = loop edge count, 1st entry=outer loop
                                                &s3DTrimCurves,                     // in : opt ordered 3d trimming curves assigned to loops per sLoopEUCounts
                                                NULL,                               // in : opt ordered 2d trimming curves assigned to loops per sLoopEUCounts
                                                sTrimOrientations,                  // in : associated orients for each trimming curve, SM_OT_SAME or SM_OT_OPPOSITE
                                                sLoopPoints,                        // in : Point positions to build SmVertex VertexLoops
                                                pNewSurface,                        // in : new face->Surface
                                                pNewSurface->GetNaturalUVDomain(),  // in : domain of Surface used by face
                                                SM_OT_SAME,                         // in : Surface orient, oneof SM_OT_SAME or SM_OT_OPPOSITE
                                                pNewRegion,                         // out: New region if any. NULL when building trimmed surfaces, may be NotNULL for solids.
                                                pNewShell,                          // out: New shell if any.  Trimmed surfaces always create a new shell.
                                                pNewFace));                         // out: the new face

        } // end build a CapFace for the CapSurface needed by this Edgeuse

#ifdef SM_DEBUG_CODE
  if (bDebugMe || lCount == lDebugCount) 
    {
      m_pOriginalBrep->Dump() ;
      m_pOffsetBrep->Dump();

      smgfx_Erase();
      smgfx_SetLook(1,2, 0,0,1) ; m_pOriginalBrep->Draw(TRUE); sm_GraphicsLoop();
      smgfx_SetLook(2,3, 0,1,0) ; m_pOffsetBrep->Draw(TRUE); sm_GraphicsLoop();
      smgfx_SetLook(3,4, 1, 0, 0); if(pEdgeuse)      pEdgeuse->Draw();      sm_GraphicsLoop();
      smgfx_SetLook(4,5, 1, 0, 1); if(pOtherEdgeuse) pOtherEdgeuse->Draw(); sm_GraphicsLoop();
      smgfx_SetLook(5,6, 1,.5, 0); if(pOffsetEU)     pOffsetEU->Draw();     sm_GraphicsLoop();
      smgfx_SetLook(6,7, 1, 0,.5); if(pOtherOffEU)   pOtherOffEU->Draw();   sm_GraphicsLoop();

      smgfx_SetLook(1,2, 0,1,1); if(pEdgeuse)      pEdgeuse->GetFace()->Draw(SM_DM_CROSSHATCH);      sm_GraphicsLoop();
      smgfx_SetLook(1,2, 0,1,1); if(pOtherEdgeuse) pOtherEdgeuse->GetFace()->Draw(SM_DM_CROSSHATCH); sm_GraphicsLoop();
      smgfx_SetLook(1,2, 1,1,0); if(pOffsetEU)     pOffsetEU->GetFace()->Draw(SM_DM_CROSSHATCH);     sm_GraphicsLoop();
      smgfx_SetLook(1,2, 1,1,0); if(pOtherOffEU)   pOtherOffEU->GetFace()->Draw(SM_DM_CROSSHATCH);   sm_GraphicsLoop();

      smgfx_SetLook(1,2, 1,0,0) ; pNewSurface->DrawUV(); sm_GraphicsLoop() ;
      smgfx_SetLook(1,2, 0,0,0) ; pNewFace->Draw(SM_DM_CROSSHATCH); sm_GraphicsLoop() ;
      sm_GraphicsLoop();

      for (ULONG kkk=0; kkk<sSSIData.m_v3DCurves.GetSize(); kkk++) 
        { smgfx_SetLook(4,5, .5,1,.2) ; sSSIData.m_v3DCurves[kkk]->Draw(); sm_GraphicsLoop(); }
      sm_GraphicsLoop(); 
    }
#endif // SM_DEBUG_CODE

      // When the CapSurface is known to self-intersect
      if (sSSIData.m_v3DCurves.GetSize() > 0) 
        {
//            SmTArray<SmEdge*> sNewEdges;
//            SER(m_pOffsetBrep->MergeSelfIntersectionCurves(*pNewSurface,
//                m_pOffsetBrep->GetTolerance(),s3DCurves,sUVCurves1,sUVCurves2,
//                sNewEdges));
        }

// Remove Composites - moved block up from below
      SM_ASSERT_MSG(pNewSurface->GetFace() != NULL,     _T("SmOffsetExecutive::BuildOffsetOfEdgeuses - Remove Composite single face assumption wrong here - needs debug")) ;
      SM_ASSERT_MSG(pNewSurface->GetFace() == pNewFace, _T("SmOffsetExecutive::BuildOffsetOfEdgeuses - pNewSurface->Face != pNewFace - needs debug")) ;
      // pNewFace = pNewSurface->GetFace() ; - removed redundant line
      // map origEdge to its offset capFaces
      sNewFaces.Add(pNewFace) ;
      SER( m_vEToF.RelatePair(pEdge,pNewFace)  );

// Remove Composites
//    // map origEdge to its offset capFaces
//    SmTArray<SmFace*> sFaces;
//    m_pOffsetBrep->GetFacesOfSurface(pNewSurface,sFaces);
//    sNewFaces.Append(sFaces);
//    for (ULONG kkk=0; kkk<sFaces.GetSize(); kkk++) 
//      {
//        SER(m_vEToF.RelatePair(pEdge,sFaces[kkk]));
//      }
    } // end iter every matched offset edgeuse

  // Now Process any edgeuses that are needed to cover the gaps left
  // because of self intersection of base surfaces.
  // We will simply find the adjacent surfaces and create a 3 sided
  // Coons patch to blend it.
  SmBoolean bProblems = FALSE;
  if (FillGapEdgeuses(pEdge,sOffGapEUS,sNewFaces) != SM_SUCCESS) 
    { bProblems = TRUE; }

  if (FillGapEdgeuses(pEdge,sOtherOffGapEUS,sNewFaces) != SM_SUCCESS) 
    { bProblems = TRUE; }

  // 
  if ( bProblems )
    {
      SmBSplineCurve *pCrv = pEdge->CreateTrimmedNURBSCurve(m_crContext);
      SmObjDelete sCleanCrv(pCrv);
      SmPrimitiveCreation sPC(m_pOffsetBrep->GetInfiniteRegion());
      SmTArray<SmFace*> sStartF, sMidF, sEndF;
      SER(sPC.CreatePipeSweep(m_rGeometryCreation.GetOffsetDistance(),
          pCrv,m_rGeometryCreation.GetThisApproxTol3d(),FALSE,FALSE,
          sStartF,sMidF,sEndF));
    }

  return SM_SUCCESS;

} // end SmOffsetExecutive::BuildOffsetOfEdgeuses

/*******************************************************************//**
PURPOSE: Fill a gap left by a self intersecting offset of an edge.

NOTES: 
***********************************************************************/
SmStatus SmOffsetExecutive::FillGapEdgeuses
 (SmEdge                     * pEdge,
  const SmTArray<SmEdgeuse*> & crGapEdgeuses,
  const SmTArray<SmFace*>    & crCanidateFaces)
{
    if (crGapEdgeuses.GetSize() == 0)
      { return SM_SUCCESS; }

#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe = FALSE;
    if (bDebugMe) {
        smgfx_Erase();
        for (ULONG kkk=0; kkk<crGapEdgeuses.GetSize(); kkk++) {
            crGapEdgeuses[kkk]->GetEdge()->Draw();
            sm_GraphicsLoop();
        }
        for (ULONG mmm=0; mmm<crCanidateFaces.GetSize(); mmm++) {
            crCanidateFaces[mmm]->Draw();
            sm_GraphicsLoop();
        }
        smgfx_Erase();
        m_pOffsetBrep->Draw();
        sm_GraphicsLoop();
    }
#endif // SM_DEBUG_CODE


    SmTArray<SmVertex*> sVertices;
    SmTArray<SmVertex*> sAllVertices;
    for (ULONG ii=0; ii<crCanidateFaces.GetSize(); ii++) {
        crCanidateFaces[ii]->GetVertices(sVertices);
        sAllVertices.Append(sVertices);
    }

    SmTArray<SmEdge*> sCanidateEdges;
    SmTArray<SmEdge*> sEdges;
    ULONG ii, jj, kk, ll;
    for (ii=0; ii<crGapEdgeuses.GetSize(); ii++)
    {
        SmEdgeuse *pGapEU = crGapEdgeuses[ii];
        SmPoint3d sMidPoint;
        SER(pGapEU->NormalizedEvaluate(0.5,FALSE,sMidPoint));  // TRUE = UV Eval, FALSE = 3d Eval

        // Look at each end of the Gap EU
        for (jj=0; jj<2; jj++)
        {
            SmVertex *pV = pGapEU->GetVertexuse()->GetVertex();
            if (jj==1) { pV = pGapEU->GetMate()->GetVertexuse()->GetVertex(); }
            SmPoint3d sVPnt = pV->GetPoint();

            // Look at each vertex of the canidate faces
            for (kk=0; kk<sAllVertices.GetSize(); kk++)
            {
                SmVertex *pTestV = sAllVertices[kk];
                double dTol = pTestV->GetTolerance() + pV->GetTolerance();
                SmPoint3d pTestPnt = pTestV->GetPoint();

                // See if a test vertex is on current end of Gap EU
                if (pTestPnt.DistanceBetween(sVPnt) < dTol)
                {
                    pTestV->GetEdges(sEdges);
                    SmEdge *pEClosest = NULL;
                    double dClosestDist = SM_BIG_DOUBLE;
                    for (ll=0; ll<sEdges.GetSize(); ll++) {
                        SmEdge *pE = sEdges[ll];
                        SmVertex *pOtherV = pE->GetOtherVertex(pTestV);
                        double dDist = sMidPoint.DistanceBetween(pOtherV->GetPoint());
                        if (dDist < dClosestDist) {
                            dClosestDist = dDist;
                            pEClosest = pE;
                        }
                    }
                    sCanidateEdges.Add(pEClosest);
                }
            }
        }

        if (sCanidateEdges.GetSize() != 2) {
            // MSG(_T("Possible Problem In Offset"));
            return SM_SUCCESS;
        }

        SmTArray<SmFace*> sBlendedFaces;
        SmTArray<SmEdgeuse*> sEdgeuses;
        sEdgeuses.Add(pGapEU);
        for (jj=0; jj<sCanidateEdges.GetSize(); jj++) {
            SmEdge *pE = sCanidateEdges[jj];
            sEdgeuses.Add(pE->GetSameOrientedEdgeuse());
            sBlendedFaces.AddUnique(
                sEdgeuses.GetLast()->GetFace());
        }

        // We should now have enough information to make a blend.
        SmTArray<SmSurface*> sBlendSurfaces;
        if (m_pOffsetBrep->CreateCornerBlend(m_crContext,sEdgeuses,NULL,
            m_rGeometryCreation.GetThisApproxTol3d(),
            m_rGeometryCreation.GetTangencyTolRadians(),sBlendSurfaces) != SM_SUCCESS) 
          {
            SER(SM_ERR);
          }

        for (jj=0; jj<sBlendSurfaces.GetSize(); jj++) {
            SmSurface *pNewSurface = sBlendSurfaces[jj];
            SmFace *pNewFace = NULL;
            SER(m_pOffsetBrep->CreateFaceFromSurface(pNewSurface,
                pNewSurface->GetNaturalUVDomain(),pNewFace));
            SER( m_vEToF.RelatePair(pEdge,pNewFace) );
            SER(m_vBlendToFaces.RelateToMany(pNewFace,sBlendedFaces));
        }

    } // end iter every GapEdgeuse

    return SM_SUCCESS;

} // end FillGapEdgeuses

/*******************************************************************//**
PURPOSE: Build the offset surface for a vertex.

NOTES: 
***********************************************************************/
SmStatus SmOffsetExecutive::BuildOffsetOfVertex
  (SmVertex                * pVertex,            // in : Target Vertex
   SmTArray<SmCurve*>      & r3DTrimCurves,      // out: 
   SmTArray<SmOrientType>  & rCurveOrientations, // out: 
   SmSurface              *& rpNewSurface)       // out: 
{
  NER(pVertex);

  // init output
  rpNewSurface = NULL;
  r3DTrimCurves.ReSet();
  rCurveOrientations.ReSet();

#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe = FALSE;
static ULONG lCount      = 1 ; lCount++ ;
static ULONG lDebugCount = 0 ;
  if (bDebugMe || lCount == lDebugCount) 
    {
      smgfx_Erase();
      smgfx_SetLook(1,2, 0,0,1) ; m_pOriginalBrep->Draw(TRUE) ; sm_GraphicsLoop() ;
      smgfx_SetLook(1,2, 0,1,0) ; m_pOffsetBrep->Draw(TRUE) ; sm_GraphicsLoop() ;
      smgfx_SetLook(4,5, 1,0,0) ; pVertex->Draw(); sm_GraphicsLoop();
      sm_GraphicsLoop();
    }
#endif // SM_DEBUG_CODE

  // locals
  SmTArray<SmEdge*> sEdges;
  SmTArray<SmEdge*> sEdgesToBlend;
  SmPoint3d sVPnt         = pVertex->GetPoint();
  double    dTol          = m_pOffsetBrep->GetTolerance();
  double    dOffsetRadius = m_rGeometryCreation.GetOffsetDistance();
  pVertex->GetEdges(sEdges);
  ULONG ii, jj;

  // Pick up the edge(s) from each face which is an arc and
  // has its center at the vertex.

  // when vertex is attached to one edge
  if (sEdges.GetSize() == 1) 
    {
      // get Edgeuse starting at this vertex
      SmEdgeuse *pEU = sEdges[0]->GetPrimaryEdgeuse();
      if (pEU->GetVertexuse()->GetVertex() != pVertex) 
        {
          pEU = pEU->GetRadial();
        }

      // when edgeuse is a strut (this vertex is part of a loop not connected to other edges)
      if (pEU->IsStrut()) 
        {
          // get list of offset faces built from this EU->Faceuse.
          SmFaceuse *pFU = pEU->GetFaceuse();
          SmTArray<SmFace*> sOffFaces;
          m_vFUToF.GetSeconds(pFU,sOffFaces);

          // try faceuse->Mate when no offset faces are found
          if (sOffFaces.GetSize() == 0) 
            {
              m_vFUToF.GetSeconds(pFU->GetMate(),sOffFaces);
            }

          // for every offsetFace built from this Faceuse
          // add circular arc face->Edge->Curves with radii less than the offset distance
          //  to sEdgesToBlend list
          for (ii=0; ii<sOffFaces.GetSize(); ii++) 
            {
              SmFace *pOffF = sOffFaces[ii];

#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe2 = FALSE;
              if (bDebugMe2) 
                {
                  smgfx_Erase();
                  smgfx_SetColor(1,0,0);
                  sVPnt.Draw();
                  sm_GraphicsLoop();
                  smgfx_SetColor(0,0,0);
                  pOffF->Draw();
                  sm_GraphicsLoop();
                }
#endif // SM_DEBUG_CODE
              // for every edge attached to this offsetFace
              pOffF->GetEdges(sEdges);
              for (jj=0; jj<sEdges.GetSize(); jj++) 
                {
                  // get Edge->Curve
                  SmBSplineCurve* pCurve = SM_CAST_PTR(SmBSplineCurve, sEdges[jj]->GetCurve());
                  if (!pCurve) SER_MSG(SM_ERR, _T("Found Edge with no Curve"));

                  // when edgeCurve is an arc
                  SmAxis2Placement sRefFrame;
                  double dRadius, dStartAng, dEndAng;
                  SmPoint3d sLinePoint;
                  SmVector3d sLineVector;
                  if (pCurve->IsArc(10, dTol,
                                    sRefFrame, dRadius,
                                    dStartAng, dEndAng)) 
                    {
                      // when the arc radius is less than the offset distance
                      double dOffDistance = m_rGeometryCreation.GetOffsetDistance();
                      if (   dRadius < dOffDistance
                          && sRefFrame.GetOriginRef().DistanceBetween(sVPnt) < dOffDistance + dTol) 
                        {
                          // add edge to list of edgesToBlend
                          sEdgesToBlend.Add(sEdges[jj]);
                        
                        } // end arcRadius less than offsetDistance check
                    } // end is Arc Check
                } // end iter every edge attached to face
            } // end iter every face offset from this faceuse

          // when there is just 1 edge to blend
          if (sEdgesToBlend.GetSize() == 1) 
            {
              // build the Vertex Offset Surface
              SmStatus eStat = m_rGeometryCreation.VertexOffset(m_crContext,
                                                                pVertex,
                                                                sEdgesToBlend,
                                                                r3DTrimCurves,
                                                                rCurveOrientations,
                                                                rpNewSurface);
              if ( eStat != SM_SUCCESS )
                { SER( eStat ); }
            }

          // all done
          return SM_SUCCESS;
        
        } // end vertex IsStrut check
    } // end vertex attached to just 1 edge check


  // gather offset CapFaces built around edges connected to this vertex
  SmTArray<SmFace*> sOffsetFaces;
  SmTArray<SmFace*> sEFaces;
  for (ii=0; ii<sEdges.GetSize(); ii++) 
    {
      SmEdge *pE = sEdges[ii];
      m_vEToF.GetSeconds(pE,sEFaces);

      // when a vertex->Edge has no offset capFace
      if (sEFaces.GetSize() == 0) 
        { // no need to build a vertex offset capFace
          return SM_SUCCESS;
        }

      // else gather the offset capFaces
      sOffsetFaces.Append(sEFaces);

    } // end iter every vertex->Edge

  // for every offset CapFace built for every vertex->edge
  for (ii=0; ii<sOffsetFaces.GetSize(); ii++) 
    {
      SmFace *pF = sOffsetFaces[ii];

      // for every CapFace->Edge
      pF->GetEdges(sEdges);
      for (jj=0; jj<sEdges.GetSize(); jj++) 
        {
          SmBSplineCurve* pCurve = SM_CAST_PTR(SmBSplineCurve, sEdges[jj]->GetCurve());
          if (!pCurve) { SER(SM_ERR); }

          // locals
          SmAxis2Placement sRefFrame;
          double dRadius, dStartAng, dEndAng;
          SmPoint3d sLinePoint;
          SmVector3d sLineVector;

          // when the edge->Curve is an arc
          if (pCurve->IsArc(10,dTol,
              sRefFrame,dRadius,dStartAng,dEndAng)) 
            {
              // when the arc radius is zero within tolerance
              if (sRefFrame.GetOriginRef().DistanceBetween(sVPnt) < dTol) 
                {
                  // add edge to list of edgesToBlend
                  sEdgesToBlend.Add(sEdges[jj]);
                }
            } // end curve is an arc check

          // else when edge->Curve is a short Line
          else if (   pCurve->IsLine(10,dTol,sLinePoint,sLineVector)
                   && smos_Fabs(sLineVector.Length()-dOffsetRadius) < dTol) 
            {
              double dDist, dParam;
              SmPoint3d sLineEnd = sLinePoint + sLineVector;

              // and vertex is within tol of the line segment
              SER(smgu_SegmentPointDistance(sLinePoint,sLineEnd,sVPnt,dDist,dParam));
              if (dDist < dTol) 
                {
                  // add edge to list of edgesToBlend
                  sEdgesToBlend.Add(sEdges[jj]);
                }
            } // end edge->Curve is a short line check

#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe3 = FALSE;
          if (bDebugMe3) 
            {
              smgfx_Erase();
              smgfx_SetColor(1,0,0);
              sVPnt.Draw();
              sm_GraphicsLoop();
              smgfx_SetColor(0,0,0);
              sEdges[jj]->Draw();
              sm_GraphicsLoop();
            }
#endif // SM_DEBUG_CODE  
        } // end iter every CapFace->Edge
    } // end iter every vertex->edge->CapFace

  // make the vertex offset capSurface 
  SER(m_rGeometryCreation.VertexOffset(m_crContext, pVertex, sEdgesToBlend,
                                       r3DTrimCurves, rCurveOrientations, 
                                       rpNewSurface));

  // all done
  return SM_SUCCESS;

} // end SmOffsetExecutive::BuildOffsetOfVertex

/*******************************************************************//**
PURPOSE: Build the fillet offsets for this faceuse.

NOTES: 

GWC NOTE: This is an incomplete function
***********************************************************************/
SmStatus SmOffsetExecutive::BuildFilletOffsetsOfFaceuse
  (SmFaceuse *pFU)               // NotUsed: in : target faceuse to offset
{
  SM_REF1(pFU) ;
#if 0 // Skip this routine:

  NER(pFU);

#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe3 = FALSE;
static ULONG lCount      = 1 ; lCount++ ;
static ULONG lDebugCount = 0 ;
  if (bDebugMe3 || lCount == lDebugCount) {
      smgfx_Erase();
      pFU->GetFace()->Draw();
      sm_GraphicsLoop();
  }
#endif // SM_DEBUG_CODE

  // locals
  SmFace                   *pF = pFU->GetFace();
  SmBoolean                 bIsFillet ;
  SmBoolean                 bNormalIsOutward;
  ULONG                     lFilletCrossSection ;
  ULONG                     lFilletSolverType;
  SmSurfParamType           eFilletRailDirection;
  SmTArray<SmBSplineCurve*> sUVRails;
  SmTArray<SmFace*>         sRailFaces;
  
  // see if pF is a fillet Surface - get its parameters
  SER(pF->ExtractFillet(*pF->GetContext(),pF->GetTolerance(),2.0,bIsFillet,
      lFilletCrossSection, lFilletSolverType, eFilletRailDirection,
      bNormalIsOutward, sUVRails, sRailFaces));

  // all done
  if (!bIsFillet) return SM_SUCCESS;

#endif // 0 = Skip this routine.

  return SM_SUCCESS;

} // end SmOffsetExecutive::BuildFilletOffsetsOfFaceuse

/*******************************************************************//**
PURPOSE: Build the face that represents the offset of this faceuse
            in m_pOffsetBrep and merge in self-intersection data.

            map pFU to NewFace in m_vFUToF
            and map all pFU->edges to NewFace->Edges in m_vEUToEU

NOTES: 
  If the surface(s) in rSurfaces are SmBSplineSurfaces, and have no
  NURB surface, one will be created for them (MakeNurb()).

  When input Surface(s) have a domain smaller than the pFU->pFace->Surface domain,
  the UVTrimcurves of pFU->Loopuses->Edges are booleaned in 2d to find
  the bit of the original face that maps to the new face.
***********************************************************************/
SmStatus SmOffsetExecutive::BuildOffsetOfFaceuse
(
  SmFaceuse            * pFU,                // in : target faceuse being offset
  SmSurface* &           rSurface,           // in : the offset surfaces of pFU->Surface. might have selfIntersections
  SmSSIData            & rSelfIntersections  // in : summary of surface self-intersections
)
{
  NER(pFU);

#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe = FALSE;
static ULONG lCount      = 1 ; lCount++ ;
static ULONG lDebugCount = 0 ;
  if (bDebugMe || lCount == lDebugCount)
    {
      pFU->GetFace()->GetSurface()->Dump() ;
      rSurface->Dump() ;

      SmBrep *pBrep = pFU->GetBrep() ;

      smgfx_Erase();
      smgfx_SetLook(1,2, 0,0,1) ; if(pBrep) pBrep->Draw(TRUE) ; sm_GraphicsLoop() ;
      smgfx_SetLook(1,2, 0,1,1) ; rSurface->DrawUV(); sm_GraphicsLoop() ;
      smgfx_SetLook(1,2, 1,1,0) ; pFU->GetFace()->GetSurface()->DrawUV(); sm_GraphicsLoop() ;
      smgfx_SetLook(1,2, 0,0,0) ; pFU->GetFace()->Draw(SM_DM_CROSSHATCH); sm_GraphicsLoop() ;
      sm_GraphicsLoop();
    }
#endif // SM_DEBUG_CODE

  // no work - no offset surfaces
  if (rSurface == NULL) 
    { return SM_SUCCESS; }

  // locals
  ULONG kk, ll, kkk, lll, jjj, iii ;
  SmTArray<ULONG>           sCurveLoops;
  SmTArray<SmBSplineCurve*> sUVCurves;
  SmTArray<SmCurve*>        s3DCurves;
  SmTArray<SmOrientType>    sCurveOrients;
  SmTArray<SmPoint3d>       sLoopPoints;
  SmFace* sNewFacesData[32] ;  
  SmTArray<SmFace*> sNewFaces(32, (SmFace**)sNewFacesData) ;

  // for every offset surface of pFU->GetSurface
      SmSurface          *pOffsetSurf = rSurface;
      const SmContext    *cpContext   = m_pOffsetBrep->GetContext();
      SmTArray<SmCurve*> *p3DCurves   = NULL;

      // This will not work on a BSplineSurface if it has no Nurb representation.
      SM_ENSURE_SURFACE_MPNURB(pOffsetSurf) ;
      
      // Special case: face is plane - just translate the 3dCurves, currently only planes have 
      //               mapped NurbDomains where OffSurf(uv) = OrigSurf(uv) + dist*Normal(uv) 
#if 0
      if (pOffsetSurf->IsKindOf(SmPlane_TYPE) && m_eOffsetOperation != SM_OO_LOCAL_OPERATION) 
        {
          // construct face boundary curve and vertexLoop point copies
          SER(pFU->GetFace()->CreateCurvesFromFace(*cpContext, pFU->GetOrientation(), sCurveLoops,
                                                   &s3DCurves,     // out: optional list of 3d Curves
                                                   &sUVCurves,     // out: optional list of 2d Curves
                                                   sCurveOrients,sLoopPoints));

          // get offset vector as vector between corresponding points on pFU and offset Surface
          SmExtent2d sDomain = pOffsetSurf->GetNaturalUVDomain();
          SmPoint3d sOrigPnt;
          SmPoint3d sOffPnt;
          SER(pFU->GetFace()->GetSurface()->EvaluatePoint(sDomain.GetMin(),sOrigPnt));
          SER(pOffsetSurf->EvaluatePoint(sDomain.GetMin(),sOffPnt));
          SmVector3d sTranslation = sOffPnt - sOrigPnt;

          SmAxis2Placement sTransform;
          sTransform.Translate(sTranslation);

          // Translate every pFU->Face->BoundaryCurve
          for (ULONG j=0; j<s3DCurves.GetSize(); j++) 
            {
              SmCurve *pCurve = s3DCurves[j];
              SER(pCurve->Transform(sTransform));
            }

          // set p3DCurves = translated copies of origFace boundary curves
          p3DCurves = &s3DCurves;
        }
#endif // obsolete code if face = plane

      // get offsetSurface and origFace domains
      SmExtent2d sOffsetDomain = pOffsetSurf->GetNaturalUVDomain();
      SmExtent2d sFaceDomain   = pFU->GetFace()->GetUVDomain();
      SmBoolean  bFaceFits     = (   sFaceDomain.IsContainedBy(sOffsetDomain, SM_EFF_ZERO)
                                  || m_eOffsetOperation == SM_OO_LOCAL_OPERATION) ;

      // face is not a plane
      // copy face->UVTrimCurves and vertexLoop points
      // SolidCraft change: a face that does not fit the offset surface's domain is clipped to it
      // by a 2D boolean (below) of planar faces made from its UV loops, and a planar face needs
      // its loops counterclockwise, as the face's SM_OT_SAME faceuse runs them. The loops of an
      // inset's faceuse (SM_OT_OPPOSITE) run clockwise: the 2D face was inside out and its
      // intersection with the domain empty, so the face got no offset at all (the bottom of a vase
      // whose profile meets the axis at a slant, its offset's domain trimmed at the axis,
      // shelled inward: tangency-labs/solidcraft#4567).
      SER(pFU->GetFace()->CreateCurvesFromFace(*cpContext,             // in : context for new object construction
                                                bFaceFits ? pFU->GetOrientation() : SM_OT_SAME,
                                                                       // in : use Faceuse whose orientation is the same as eOrientation
                                                sCurveLoops,           // out: number of edges in each edge loop
                                                NULL,                  // out: optional copies of all 3D loop curves - note: one or both of pOpt3DCurves and
                                               &sUVCurves,             // out: optional copies of all 2D loop curves -         pOptUVCurves must be nonNULL
                                                sCurveOrients,         // out: Curve orientation in relation to its loop for each curve
                                                sLoopPoints));         // out: point position of each of face's vertex loops

      // when sOffsetDomain is larger or equal to sFaceDomain or we are doing a local operation
      if (bFaceFits)
      {
          // make a face in m_pOffset which is the offset of pFU
          SmRegion *pNewRegion = NULL;
          SmShell  *pNewShell = NULL;
          SmFace   *pNewFace = NULL;

          // for plane surfaces p3DCurves are already computed
          // for other surfaces p3DCurves is NULL forcing the call to compute the
          //     offset edges by projecting the UVTrimCurves through the offset surface
          if(m_pOffsetBrep->MakeFaceWithCurves( m_pOffsetBrep->GetInfiniteRegion(),  // in : region to contain new topology objects
              sCurveLoops,                         // in : 1 entry per loop, value = loop edge count, 1st entry=outer loop
              p3DCurves,                           // in : opt ordered 3d trimming curves assigned to loops per sLoopEUCounts
              &sUVCurves,                          // in : opt ordered 2d trimming curves assigned to loops per sLoopEUCounts
              sCurveOrients,                       // in : associated orients for each trimming curve, SM_OT_SAME or SM_OT_OPPOSITE
              sLoopPoints,                         // in : Point positions to build SmVertex VertexLoops
              pOffsetSurf,                         // in : new face->Surface
              sOffsetDomain,                       // in : domain of Surface used by face
              pFU->GetOrientation(),               // in : Surface orient, oneof SM_OT_SAME or SM_OT_OPPOSITE
              pNewRegion,                          // out: New region if any. NULL when building trimmed surfaces, may be NotNULL for solids.
              pNewShell,                           // out: New shell if any.  Trimmed surfaces always create a new shell.
              pNewFace ) == SM_SUCCESS)            // out: the new face    
          {

            NER( pNewFace );

            pFU->GetFace()->Notify( SM_NO_COPY, pNewFace, SM_NO_GET_BREP( pNewFace ), SM_NO_GET_BREP( pFU->GetFace() ) );

#ifdef SM_DEBUG_CODE
            if(bDebugMe || lCount == lDebugCount)
            {
              pNewFace->Dump();
              SmBrep *pBrep = pFU->GetBrep();

              smgfx_Erase();
              smgfx_SetLook( 1, 2, 0, 0, 1 ); if(pBrep) pBrep->Draw( TRUE ); sm_GraphicsLoop();
              smgfx_SetLook( 1, 2, 1, 1, 0 ); pFU->GetFace()->GetSurface()->DrawUV(); sm_GraphicsLoop();
              smgfx_SetLook( 1, 2, 0, 0, 0 ); pFU->GetFace()->Draw( SM_DM_CROSSHATCH ); sm_GraphicsLoop();
              smgfx_SetLook( 1, 2, 0, 1, 1 ); pNewFace->GetSurface()->DrawUV(); sm_GraphicsLoop();
              smgfx_SetLook( 1, 2, 0, 0, 0 ); pNewFace->Draw( SM_DM_CROSSHATCH ); sm_GraphicsLoop();
              sm_GraphicsLoop();
            }
#endif // SM_DEBUG_CODE

            // Merge the self-intersection curves into the m_pOffsetBrep offset face.
            SmSurface *pSurface = pNewFace->GetSurface();
            SmTArray<SmEdge*> sNewEdges;
            if(m_bMergeResults
                && rSelfIntersections.m_v3DCurves.GetSize() > 0)
            {
// Remove Composites
//              SmTemporaryChangeValue<SmBoolean> sChange1( m_pOffsetBrep->m_bMakeComposites, FALSE );
              SER( m_pOffsetBrep->MergeSelfIntersectionCurves( *pSurface,
                                                                m_pOffsetBrep->GetTolerance(), 
                                                                rSelfIntersections.m_v3DCurves, 
                                                                rSelfIntersections.m_vUVCurves1,
                                                                rSelfIntersections.m_vUVCurves2,
                                                                sNewEdges ) );
            } // end need to merge self-intersections check

          // save the new face
            sNewFaces.Add( pNewFace );

            // add in any new faces created by merging self-intersection edges
            SmFace* sEFacesData[32];
            SmTArray<SmFace*> sEFaces( 32, (SmFace**)sEFacesData );
            for(kk = 0; kk < sNewEdges.GetSize(); kk++)
            {
              SmEdge *pE = sNewEdges[kk];
              pE->GetFaces( sEFaces );
              for(kkk = 0; kkk < sEFaces.GetSize(); kkk++)
              {
                SmFace *pF = sEFaces[kkk];
                sNewFaces.AddUnique( pF );
              } // end iter every intersection edge face
            } // end iter every self-intersection edge

          // Now Let's try to remove faces if possible.
            SmTArray<SmEdge*> sFEdges;
            SmTArray<SmFace*> sDelFaces;

            // for every new face
            for(ll = 0; ll < sNewFaces.GetSize(); ll++)
            {
              // If we are extending we can skip this test
              if(m_bExtendConvexEdges) { continue; }

              // get this newFace edges
              SmFace *pF = sNewFaces[ll];
              pF->GetEdges( sFEdges );

              // for every newFace->Edge
              for(lll = 0; lll < sFEdges.GetSize(); lll++)
              {
                SmEdge *pE = sFEdges[lll];

                // for lamina edges
                if(pE->IsLamina())
                {
                  SmEdgeuse *pEU = pE->GetSameOrientedEdgeuse();
                  SmPoint3d sPnt;
                  SER( pEU->NormalizedEvaluate( 0.4567, FALSE, sPnt ) );  // TRUE = UV Eval, FALSE = 3d Eval
                  double dOffDistance = m_rGeometryCreation.GetOffsetDistance();
                  double dTestDist = dOffDistance * 1.1;
                  SmExtent2d sUVDomain = pFU->GetFace()->GetUVDomain();
                  SmSolution sData[16];
                  SmSolutionArray sSolutions( 16, sData );
                  SmSurface *pSrf = pFU->GetFace()->GetSurface();
                  SmSurfaceCache *pSC = smsurf_GetSurfaceCache( pSrf ); NER( pSC );
                  SmCacheCheckOutIn sCheckIO( pSC );
                  { // Makes the change values destructors occur before surface is checked back in
                    // and possibly deleted.

                    // Turn off point testing so GlobalPointSolve() will keep all point solutions
                    //   without classifying the solution point against the trim boundaries.
                    SmTemporaryChangeValue<SmBoolean> sChangeSC( pSC->m_bPointTestEnabled, FALSE );
                    SER( pSrf->GlobalPointSolve( sUVDomain, SM_SO_MINIMIZE, sPnt,
                         pF->GetTolerance(), &dTestDist, SM_SR_SINGLE, sSolutions ) );
                  }
                  double dOffTol = m_rGeometryCreation.GetThisApproxTol3d();
                  if(sSolutions.GetSize() < 1) continue;
                  SmSolution & rSol = sSolutions[0];
                  if(rSol.m_vStart.m_dSolutionValue < dOffDistance - 2.0 * dOffTol)
                  {
                    sDelFaces.AddUnique( pF );
                    break;
                  }
                } // end lamina edge check
              } // end iter every new face edge
            } // end iter every new face

          // remove any faces found to be deletable
            if(sDelFaces.GetSize() > 0)
            {
              SmTopologyIntersector sTI( m_crContext, m_pOffsetBrep, m_pOffsetBrep, SM_EFF_ZERO, SM_EFF_ZERO );
              SER( sTI.DeleteFaces( sDelFaces ) );
              for(jjj = 0; jjj < sDelFaces.GetSize(); jjj++)
              {
                ULONG lIndex;
                if(sNewFaces.FindElement( sDelFaces[jjj], lIndex ))
                {
                  sNewFaces.RemoveAt( lIndex );
                }
              }
            } // end any delete face check
          } // end if MakeFaceWithCurves was successful
        } // end sFaceDomain is smaller than offset Face branch

      else  // Need to do a 2D boolean  
        {
          SmBrep *p2DBrep = new (*cpContext) SmBrep();
          p2DBrep->m_bEditingEnabled = TRUE;
          SmZoneTol3d sZoneTol3d = 1.0e-6 ;
          SM_OLDTOL_LINE p2DBrep->SetTolerance(sZoneTol3d);
          SmFace *pNewFace = NULL;

          // get bounding box containing all pFU->Face->UVTrimCurves
          SmExtent3d sTotalBBox;
          for (ULONG ii=0; ii<sUVCurves.GetSize(); ii++) 
            {
              SER( sUVCurves[ii]->ConvertTo3D() );
              SmExtent3d sCrvBBox;
              SER( sUVCurves[ii]->CalculateBoundingBox(sUVCurves[ii]->GetNaturalInterval(),&sCrvBBox) );
              sTotalBBox.Union(sCrvBBox,sTotalBBox);
            }

          // locals
          SmPoint3d        sMin(sOffsetDomain.GetMin());
          SmAxis2Placement sPlace;

          sPlace.Translate(sMin);
          SmVector2d sSize(sOffsetDomain.GetSize());

          // create a 3d line between every offsetSurface domain corner point
          SM_PTR_ARRAY(sBoxCurves,SmCurve,16);

          // line [0,0] to [1,0]
          SmBSplineCurve *pLine = NULL ;
          SER(SmBSplineCurve::CreateLineSegment(m_crContext,3,sOffsetDomain.GetMin(),sOffsetDomain.Evaluate(1.0,0.0),pLine));
          SmExtent1d sIvl(sOffsetDomain.GetMin().x, sOffsetDomain.GetMax().x);
          SER( pLine->EditParameterization(sIvl) );
          sBoxCurves.Add(pLine); pLine = NULL ;

          // line [0,1] to [1,1]
          SER(SmBSplineCurve::CreateLineSegment(m_crContext,3,sOffsetDomain.Evaluate(0.0,1.0),sOffsetDomain.GetMax(),pLine));
          SER( pLine->EditParameterization(sIvl) );
          sBoxCurves.Add(pLine); pLine = NULL ;

          // line [0,0] to [0,1]
          SER(SmBSplineCurve::CreateLineSegment(m_crContext,3,sOffsetDomain.GetMin(),sOffsetDomain.Evaluate(0,1.0),pLine));
          SER( sIvl.SetMinMax(sOffsetDomain.GetMin().y,sOffsetDomain.GetMax().y) );
          SER( pLine->EditParameterization(sIvl) );
          sBoxCurves.Add(pLine); pLine = NULL ;

          // line [1,0] to [1,1]
          SER(SmBSplineCurve::CreateLineSegment(m_crContext,3,sOffsetDomain.Evaluate(1.0,0.0),sOffsetDomain.GetMax(),pLine));
          SER( pLine->EditParameterization(sIvl) );
          sBoxCurves.Add(pLine);

          // make a new Brep containing 1 planar face bounded by the 3d Curves just made
          // Note, it's a uv-space box even though the curves have dimension 3, for the call.
          SmBrep *pDomBrep = new (m_crContext) SmBrep();
          SmObjDelete sClean( pDomBrep );
          SM_OLDTOL_LINE pDomBrep->SetTolerance(p2DBrep->GetTolerance());
          pDomBrep->m_bEditingEnabled = TRUE;
          SER(pDomBrep->CreatePlanarFaceWith3DCurves(pDomBrep->GetInfiniteRegion(),
                                                     sBoxCurves,
                                                     pDomBrep->GetTolerance(),
                                                     pNewFace));

          // copy the plane just made and get its natural domain
          SmTArray<SmSurface*> sSurfaces;
          pDomBrep->GetSurfaces(sSurfaces);

          SmPlane *pPlaneSurf = SM_CAST_PTR(SmPlane,sSurfaces[0]);

          // Temp workaround to appease Linux gcc compiler
          // Original code with typecast has warning: dereferencing type-punned pointer will break strict-aliasing rules [-Wstrict-aliasing]
          // pPlaneSurf->Copy(*cpContext,(SmSurface*&)pCopyPlane);
          SmPlane* pCopyPlane = NULL;
          SER( pPlaneSurf->Copy(*cpContext,pCopyPlane) );

          SmExtent2d sCopyDomain = pCopyPlane->GetNaturalUVDomain();

          // increase the sCopyDomain size to include the projection of every
          // UVTrimCurve bounding box corner onto this plane we just made.
          SmPoint2d sUV;
          SER( pCopyPlane->ProjectPointToUVDomain(sTotalBBox.GetMin(),sUV) );
          sCopyDomain.AddPoint2d(sUV);
          SER( pCopyPlane->ProjectPointToUVDomain(sTotalBBox.GetMax(),sUV) );
          sCopyDomain.AddPoint2d(sUV);
          SER( pCopyPlane->ProjectPointToUVDomain(sTotalBBox.Evaluate(1,0,0),sUV) );
          sCopyDomain.AddPoint2d(sUV);
          SER( pCopyPlane->ProjectPointToUVDomain(sTotalBBox.Evaluate(0,1,0),sUV) );
          sCopyDomain.AddPoint2d(sUV);
          SER( pCopyPlane->AdjustSTEPUVDomain(sCopyDomain) );

          // make a planar face using the Offset UVTrimCurves
          SmRegion *pNewRegion = NULL;
          SmShell *pNewShell = NULL;
          SER(p2DBrep->MakeFaceWithCurves(p2DBrep->GetInfiniteRegion(),     // in : region to contain new topology objects
                                          sCurveLoops,                      // in : 1 entry per loop, value = loop edge count, 1st entry=outer loop
                                          (SmTArray<SmCurve*>*)&sUVCurves,  // in : opt ordered 3d trimming curves assigned to loops per sLoopEUCounts  
                                          NULL,                             // in : opt ordered 2d trimming curves assigned to loops per sLoopEUCounts  
                                          sCurveOrients,                    // in : associated orients for each trimming curve, SM_OT_SAME or SM_OT_OPPOSITE
                                          sLoopPoints,                      // in : Point positions to build SmVertex VertexLoops
                                          pCopyPlane,                       // in : new face->Surface  
                                          pCopyPlane->GetNaturalUVDomain(), // in : domain of Surface used by face
                                          SM_OT_SAME,                       // in : Surface orient, oneof SM_OT_SAME or SM_OT_OPPOSITE
                                          pNewRegion,                       // out: New region if any. NULL when building trimmed surfaces, may be NotNULL for solids.
                                          pNewShell,                        // out: New shell if any.  Trimmed surfaces always create a new shell.
                                          pNewFace));                       // out: the new face

          // Intersect in 2d the pDomBrep with the p2DBrep
          SmBrep *pResult;
          SER(SmPrimitiveCreation::Boolean2D(pDomBrep,p2DBrep,SM_2D_INTERSECTION,pResult));  // note: increments unlocked mark value

          // pResult now contains Faces that represent the uv domain of the (portion of the)
          // original Face that corresponds to the offset Face.  Its Edge curves represent
          // uv trim curves, which we map into the offset Brep, with pOffsetSurface.
          SM_PTR_ARRAY(s2DFaces, SmFace, 16);  // SmTArray<SmFace *>
          pResult->GetFaces(s2DFaces);
          for (kkk=0; kkk<s2DFaces.GetSize(); kkk++) 
            {
              SmFace *pF = s2DFaces[kkk];
              // the loops in the sense of pFU's orientation, which the offset face is made with
              SER(pF->CreateCurvesFromFace(*cpContext, pFU->GetOrientation(), sCurveLoops,
                                           (SmTArray<SmCurve*>*)&sUVCurves,   // out: optional "3d" Curves, which are actually uv curves
                                           NULL,                              // out: optional UVTrimCurves
                                           sCurveOrients, sLoopPoints));

              // Convert face 3dCurves into 2D curves (all their Z components should be zero)
              for (iii=0; iii<sUVCurves.GetSize(); iii++) 
                {
                  SER( sUVCurves[iii]->ConvertTo2D() );
                }

              // make the offset Face with Curves from intersected Brep
              SmFace *pOffFace = NULL;
              SER(m_pOffsetBrep->MakeFaceWithCurves(m_pOffsetBrep->GetInfiniteRegion(), // in : region to contain new topology objects
                                                    sCurveLoops,                        // in : 1 entry per loop, value = loop edge count, 1st entry=outer loop
                                                    NULL,                               // in : opt ordered 3d trimming curves assigned to loops per sLoopEUCounts 
                                                    &sUVCurves,                         // in : opt ordered 2d trimming curves assigned to loops per sLoopEUCounts 
                                                    sCurveOrients,                      // in : associated orients for each trimming curve, SM_OT_SAME or SM_OT_OPPOSITE
                                                    sLoopPoints,                        // in : Point positions to build SmVertex VertexLoops
                                                    pOffsetSurf,                        // in : new face->Surface
                                                    sOffsetDomain,                      // in : domain of Surface used by face
                                                    pFU->GetOrientation(),              // in : Surface orient, oneof SM_OT_SAME or SM_OT_OPPOSITE
                                                    pNewRegion,                         // out: New region if any. NULL when building trimmed surfaces, may be NotNULL for solids.
                                                    pNewShell,                          // out: New shell if any.  Trimmed surfaces always create a new shell.
                                                    pOffFace ));                         // out: the new face
              NER( pOffFace );
              sNewFaces.Add( pOffFace );

              pFU->GetFace()->Notify(SM_NO_COPY, pOffFace, SM_NO_GET_BREP( pOffFace ), SM_NO_GET_BREP(pFU->GetFace()));

#ifdef SM_DEBUG_CODE
              if (bDebugMe || lCount == lDebugCount) 
                {
                  SmBrep *pBrep = pFU->GetBrep() ;

                  smgfx_Erase();
                  smgfx_SetLook(1,2, 0,0,1) ; if(pBrep) pBrep->Draw(TRUE) ; sm_GraphicsLoop() ;
                  smgfx_SetLook(1,2, 1,1,0) ; pFU->GetFace()->GetSurface()->DrawUV(); sm_GraphicsLoop() ;
                  smgfx_SetLook(1,2, 0,0,0) ; pFU->GetFace()->Draw(SM_DM_CROSSHATCH); sm_GraphicsLoop() ;
                  smgfx_SetLook(1,2, 0,1,1) ; pOffFace->GetSurface()->DrawUV(); sm_GraphicsLoop() ;
                  smgfx_SetLook(1,2, 0,0,0) ; pOffFace->Draw(SM_DM_CROSSHATCH); sm_GraphicsLoop() ;
                  sm_GraphicsLoop();
                }
#endif // SM_DEBUG_CODE
            }
          
        } // end Need to do a 2D boolean branch

#ifdef SM_DEBUG_CODE
      if (bDebugMe || lCount == lDebugCount) 
        {
          SmBrep *pBrep = pFU->GetBrep() ;

          smgfx_Erase();
          smgfx_SetLook(1,2, 0,0,1) ; if(pBrep) pBrep->Draw(TRUE) ; sm_GraphicsLoop() ;
          smgfx_SetLook(1,2, 0,1,0) ; m_pOffsetBrep->Draw(TRUE) ; sm_GraphicsLoop() ;
          smgfx_SetLook(1,2, 0,0,0) ; pFU->GetFace()->Draw(SM_DM_CROSSHATCH); sm_GraphicsLoop() ;
          sm_GraphicsLoop();
        }
#endif // SM_DEBUG_CODE
     

  // Create relations for every newFace - 
  // 1. relate pFU and NewFace in m_vFUToF
  // 2. relate all pFU->edges when NewFace->Edges in m_vEUToEU
  for (ULONG ii=0; ii<sNewFaces.GetSize(); ii++) 
    {
      SmFace *pNewF = sNewFaces[ii];
      SM_ASSERT(pNewF->IsKindOf(SmFace_TYPE)) ;
      SM_ASSERT(pFU->IsKindOf(SmFaceuse_TYPE)) ;
      SER( m_vFUToF.RelatePair(pFU,pNewF) );
      SER(MapEdgeuses(pFU,pNewF));
    }

  // all done
  return SM_SUCCESS;

} // end SmOffsetExecutive::BuildOffsetOfFaceuse

/*******************************************************************//**
PURPOSE: Create a map between edgeuses of original offset faceuse and
            copied offset face edgeuses.

NOTES: 
  Since the offset face is a copy of the original face, all
  the original Face UVTrimCurves should have an offsetFace copy
  with the same UV Shape.  Edgeuses of these two different 
  faces are matched based on their UV Locations.

  1st: match edge->UVTrimCurve midPoints, but if a curve has been trimmed or split
  2nd: match offsetFace->Edge->UVTrimCurve midPoint to closest origFace->Edge->UVTrimCurve
***********************************************************************/
SmStatus SmOffsetExecutive::MapEdgeuses
  (SmFaceuse * pOriginalFU,             // in : original Faceuse bounding infiniteRegion in m_pOriginalBrep
   SmFace    * pOffsetFace)             // in : new Offset face in m_pOffsetBrep
{
#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe = FALSE;
  if (bDebugMe) 
    {
      SmBrep *pBrep = pOriginalFU->GetFace()->GetBrep() ;

      smgfx_Erase();
      smgfx_SetLook(1,2, 0,0,1); if(pBrep) pBrep->Draw(TRUE) ; sm_GraphicsLoop() ;
      smgfx_SetLook(1,2, 0,1,1); pOriginalFU->GetFace()->Draw(SM_DM_CROSSHATCH); sm_GraphicsLoop() ;
      smgfx_SetLook(1,2, 0,1,1); pOriginalFU->Draw(); sm_GraphicsLoop() ;
      smgfx_SetLook(1,2, 0,1,0); pOffsetFace->Draw(SM_DM_CROSSHATCH,11,11); sm_GraphicsLoop() ;
      smgfx_SetLook(1,2, 0,1,0); pOffsetFace->Draw(SM_DM_NORMALS); sm_GraphicsLoop() ;
      sm_GraphicsLoop();
    }
#endif // SM_DEBUG_CODE

  // locals
  SmTArray<SmEdgeuse*> sOrigEUS;
  SmTArray<SmEdgeuse*> sOffEUS;

  // gwc:note matching the originalFU with the OffsetFace->UpwardFaceuse without checking orientation
  //     can match opposite sides of the same face with one another.  Matching the UVTrimCurves
  //     based only on their midpoint value won't generally work when oppositely oriented 
  //     faceuse->edgeuse curves are being matched solely on their mid-point values.
  //
  //  Proof:  Consider a line represented with two different parameterizations as
  //   1. Line_1_Upward   = u*P1 + (1-u)*P2                    note: (u)        + (1-u)        = 1 for all values u
  //   2. Line_2_Upward   = (2*u-u**2)*P1 + (1-2*u+u**2)*P2        : (2*u-u**2) + (1-2*u+u**2) = 1 for all values u
  //  On the downward faceuse the parameterizations of these two lines are reversed as:
  //   1. Line_1_Downward = (1-u)*P1 + (u)*P2
  //   2. Line_2_Downward = (1-2*u+u**2)*P1 + (2*u-u**2)*P2
  //  Now consider the u=.5 parameter value for the lines:
  //   1. Line_1_Upward   = .5*P1 + .5*P2
  //      Line_1_Downward = .5*P1 + .5*P2   These are the same and will match by a simple mid-point check
  //   2. Line_2_Upward   = .75*P1 + .25*P2
  //      Line_2_Downward = .25*P1 + .75*P2  These are different and won't match by a simple mid-point check
  // conclusion - the mid point check works generally when the pOriginalFU orientation is matched to the 
  //              similarly oriented pOffsetFaceuse.  We can do that by checking MidSurface surface normal directions.
  //              and matching OriginalFU with upward or downward pOffsetFace edgeuses as appropriate

  // Due to splitting, the OriginalFace and OffsetFace domains may be a subset of one another.
  SmExtent2d sTgtDomain ;
  SmExtent2d sOriginalDomain = pOriginalFU->GetFace()->GetUVDomain() ;
  SmExtent2d sOffsetDomain   = pOffsetFace           ->GetUVDomain() ;
  SER( sOriginalDomain.Intersect(sOffsetDomain, sTgtDomain)  );
  
  if(sTgtDomain.IsDegenerate())
    {
      SM_ASSERT_MSG(!sTgtDomain.IsDegenerate(),_T("SmOffsetExecutive::MapEdgeuses assumption that Offset Surface domain is a subset the Original is FALSE - needs review")) ;
    } 

  // check OriginalFU and OffsetFace orientations by checking common midSurface Normal directions
  SmPoint3d sOrigPN[2], sCopyPN[2] ;
  SmPoint2d sMidUV       = sTgtDomain.Evaluate(0.5,0.5) ;
  SER( pOriginalFU->GetFace()->GetSurface()->EvaluateNormal(sMidUV, TRUE, TRUE, sOrigPN[1])  );
  SER( pOffsetFace           ->GetSurface()->EvaluateNormal(sMidUV, TRUE, TRUE, sCopyPN[1])  );
  SmBoolean bSameSrfNorms = sOrigPN[1].Dot(sCopyPN[1]) > 0.0 ;

#ifdef SM_DEBUG_CODE
  SmZoneTol3d sZoneTol3d = SmTol::GetZoneTol3d(pOriginalFU) ;
  pOriginalFU->GetFace()->GetSurface()->EvaluatePoint(sMidUV, sOrigPN[0]) ;
  pOffsetFace           ->GetSurface()->EvaluatePoint(sMidUV, sCopyPN[0]) ;
  if(    GetOffsetOperation() != SM_OO_LOCAL_OPERATION    // don't do this check for LocalOperations - surface normals won't generally align
     && !(sOrigPN[1].CloserThan(sZoneTol3d, sCopyPN[1])))
    {
      SM_ASSERT_MSG(   GetOffsetOperation() == SM_OO_LOCAL_OPERATION
                    || sOrigPN[1].CloserThan(sZoneTol3d, sCopyPN[1]),
                    _T("SmOffsetExecutive::MapEdgeuses assumption that orig and copy Face domains are the same is FALSE - needs review")) ;
      if (bDebugMe) 
        {
          SmBrep *pBrep = pOriginalFU->GetFace()->GetBrep() ;

          smgfx_Erase();
          smgfx_SetLook(1,2, 0,0,1); if(pBrep) pBrep->Draw(TRUE) ; sm_GraphicsLoop() ;
          smgfx_SetLook(1,2, 0,1,1); pOriginalFU->GetFace()->Draw(SM_DM_CROSSHATCH); sm_GraphicsLoop() ;
          smgfx_SetLook(1,2, 0,1,1); pOriginalFU->Draw(); sm_GraphicsLoop() ;
          smgfx_SetLook(1,2, 0,1,0); pOffsetFace->Draw(SM_DM_CROSSHATCH,11,11); sm_GraphicsLoop() ;
          smgfx_SetLook(1,2, 0,1,0); pOffsetFace->Draw(SM_DM_NORMALS); sm_GraphicsLoop() ;

          smgfx_SetLook(1,2, 0,1,1); pOriginalFU->GetFace()->GetSurface()->DrawUV(); sm_GraphicsLoop() ;
          smgfx_SetLook(1,2, 0,1,1); pOriginalFU->GetFace()->GetSurface()->DrawParams(); sm_GraphicsLoop() ;
          smgfx_SetLook(1,2, 0,1,0); pOffsetFace->GetSurface()->DrawUV(11,11); sm_GraphicsLoop() ;
          smgfx_SetLook(1,2, 0,1,0); pOffsetFace->GetSurface()->DrawParams(); sm_GraphicsLoop() ;
          smgfx_SetLook(3,4, 0,0,1); sOrigPN[0].Draw() ; sOrigPN[1].Draw(&sOrigPN[0]) ; sm_GraphicsLoop() ;
          smgfx_SetLook(3,4, 0,1,0); sCopyPN[0].Draw() ; sCopyPN[1].Draw(&sCopyPN[0]) ; sm_GraphicsLoop() ;
          sm_GraphicsLoop();
        }
    }
#endif // SM_DEBUG_CODE

  // Select side of pOffsetFace that shares same direction as input pOriginalFU
  pOriginalFU->GetEdgeuses(sOrigEUS);
  if(   ( bSameSrfNorms && pOriginalFU->GetOrientation() == SM_OT_SAME)
     || (!bSameSrfNorms && pOriginalFU->GetOrientation() == SM_OT_OPPOSITE)) 
    { pOffsetFace->GetUpwardFaceuse()->GetEdgeuses(sOffEUS); }
  else             
    { pOffsetFace->GetUpwardFaceuse()->GetMate()->GetEdgeuses(sOffEUS); }
  ULONG ii, jj;

  // for every edgeuse - get OriginalFU and OffsetFace->UpwardFU edgeuses
  for (ii=0; ii<sOffEUS.GetSize(); ii++) 
    {
      // Look for corresponding UV curves
      SmEdgeuse *pOffEU = sOffEUS[ii];

      // get offset EU UV midPoint, gwc:note only .5 can be used because matched edgeuses might be running in opposite directions
      SmPoint3d sOffPnt;
      SER(pOffEU->NormalizedEvaluate(0.5,TRUE,sOffPnt));  // TRUE = UV Eval, FALSE = 3d Eval
      
      // while looking for a matching midPoint
      SmBoolean bFound = FALSE;
      double dScale = 1.0;
      while (!bFound && dScale < 10000.0) 
        {
          // search every original EU
          for (jj=0; jj<sOrigEUS.GetSize(); jj++) 
            {
              SmEdgeuse *pOrigEU = sOrigEUS[jj];
              if (!pOrigEU) continue;

              // get origin EU UV midPoint
              SmPoint3d sOrigPnt;
              SER(pOrigEU->NormalizedEvaluate(0.5,TRUE,sOrigPnt));  // TRUE = UV Eval, FALSE = 3d Eval

              // get midPoint UV Gap
              double dDist  = sOrigPnt.DistanceBetween(sOffPnt);

              // When gap < ScaledZero
              double dUVTol = SM_EFF_ZERO_SQRT * (1.0 + sOffPnt.GetMaxDimension());
              if (dDist < dUVTol * dScale) 
                {
                  // found a match: relate the pair
                  bFound = TRUE;
                  SER( m_vEUToEU.RelatePair(pOrigEU,pOffEU) );

                  sm_NotifyEdgeAndVertices( pOrigEU->GetEdge(), pOffEU->GetEdge(), SM_NO_COPY );

                  sOrigEUS[jj] = NULL;
                  break;
                } // end found a matching pair of edges check

              // We are not doing too well with direct mapping because the curves
              // may have been split.  Try doing a minimization.
              if (dScale > 10.0) 
                {
                  SmCurve *pUVCurve = pOrigEU->GetUVTrimCurve();
#ifdef SM_DEBUG_CODE
                  if (bDebugMe) 
                    {
                      SmBrep *pBrep = pOriginalFU->GetFace()->GetBrep() ;

                      smgfx_Erase();
                      smgfx_SetLook(1,2, 0,0,1); if(pBrep) pBrep->Draw(TRUE) ; sm_GraphicsLoop() ;
                      smgfx_SetLook(1,2, 0,1,1); pOriginalFU->GetFace()->Draw(SM_DM_CROSSHATCH); sm_GraphicsLoop() ;
                      smgfx_SetLook(1,2, 0,1,0); pOffsetFace->Draw(SM_DM_CROSSHATCH); sm_GraphicsLoop() ;
                      smgfx_SetLook(3,6, 1,0,0); sOffPnt.Draw(); sm_GraphicsLoop() ;
                      smgfx_SetLook(3,4, 1,0,1); pUVCurve->Draw(); sm_GraphicsLoop() ;
                      sm_GraphicsLoop();
                    }
#endif // SM_DEBUG_CODE

                  // Drop OffPnt to UVCurve.  Use INTERSECT with max drop distance dUVTol*dScale.
                  SmBoolean bDropFound;
                  double dParameter, dDistance;
                  SER(pUVCurve->DropPoint(pUVCurve->GetNaturalInterval(), // in : target curve allowed domain
                                          sOffPnt,                        // in : Point to drop to curve
                                          NULL,                           // in : For Snapping, NULL=NoSnap, Vec pointing to CrvInside from desired EndPt for snapping.
                                                                          //      when DropPt is within dDistTol in cpVecToCrvInside dir of EndPt, Snap to Ivl EndPt.
                                                                          //      Vec handles ambiguities and makes sure that curves are not just touching at the ends.
                                          dUVTol*dScale,                  // in : Used for SnapDist, and SM_SO_INTERSECT sXSectTol3d distance.
                                                                          //      sXSectTol3d = dDistTol ;      SM_SO_INTERSECT 
                                                                          //      sXSectTol3d = SM_BIG_DOUBLE ; SM_SO_MINIMIZE, SM_SO_NORMALIZE
                                                                          //      sSnapTol3d  = dDistTol ;      SM_SO_MINIMIZE, SM_SO_NORMALIZE, SM_SO_INTERSECT    
                                          NULL,                           // in : If given, call cheap LocalPointSolve() else call expensive GlobalPointSolve()
                                          bDropFound,                     // out: TRUE = found a drop point
                                          dParameter,                     // out: found drop curve param
                                          dDistance,                      // out: found drop distance
                                          SM_SO_INTERSECT )) ;            // in : SM_SO_MINIMIZE = allow nonNormal drops near endPoints
                                                                          //      SM_SO_INTERSECT= like Minimize but point must be within dDistTol
                                                                          //      SM_SO_NORMALIZE= exclude nonNormal drops near endPoints
                                                                          //      default:[SM_SO_MINIMIZE] to preserve original behavior
                  if (bDropFound) 
                    {
                      bFound = TRUE;
                      SER( m_vEUToEU.RelatePair(pOrigEU,pOffEU) );

                      sm_NotifyEdgeAndVertices( pOrigEU->GetEdge(), pOffEU->GetEdge(), SM_NO_COPY );

                      sOrigEUS[jj] = NULL;
                      break;
                    } // end found a matching pair of EUs check
                } // end dScale > 10.0 check
            } // end iter every original Edgeuse

          // failed to find a matching EU mate

          // increase UV tolerance and try again
          dScale = 10.0*dScale;

        } // end while searching for a matching edgeuse
    } // end iter every offset edgeuse

  // all done
  return SM_SUCCESS;

} // end SmOffsetExecutive::MapEdgeuses

//------------------------------------------
// Utility methods
//------------------------------------------

/*******************************************************************//**
PURPOSE: Draw the related geometry

NOTES:
***********************************************************************/
SmDisplayList * SmOffsetExecutive::Draw() const
{
  SmDisplayList *pRtn = nullptr ;

#ifdef SM_DEFINED_HASH_ORDER
  static_assert(false, "SM_DEFINED_HASH_ORDER support is not implemented");
#endif

#ifdef SM_GFX_OUTPUT_CODE
  // start drawlist (unless one is already open)
  smgfx_Open(smgfx_GetRuleColor());

  // draw m_vFUToF - a one to many relationship
  const SmMapPtrToPtrs<SmFaceuse, SmFace> &rFUToF = m_vFUToF.GetFirstToSecond();
  SmTArray<SmFaceuse*> rFUToFKeys;
  rFUToF.GetAllKeys(rFUToFKeys);

  for (ULONG i = 0; i < rFUToFKeys.GetSize(); ++i)
  {
    SmFaceuse* pKey = rFUToFKeys.GetAt(i);
    pKey->GetFace()->Draw(SM_DM_CROSSHATCH);

    for (const SmFace* pValue : rFUToF.At(pKey))
    {
      pValue->Draw(SM_DM_CROSSHATCH);
    }
  }

  // draw m_vEUToEU
  const SmMapPtrToPtrs<SmEdgeuse, SmEdgeuse> &rEUToEU = m_vEUToEU.GetFirstToSecond() ;
  SmTArray<SmEdgeuse*> rEUToEUKeys;
  rEUToEU.GetAllKeys(rEUToEUKeys);

  for (ULONG i = 0; i < rEUToEUKeys.GetSize(); ++i)
  {
    SmEdgeuse* pKey = rEUToEUKeys.GetAt(i);
    pKey->Draw();

    for (const SmEdgeuse* pValue : rEUToEU.At(pKey))
    {
      pValue->Draw();
    }
  }

  // draw m_vEToF
  const SmMapPtrToPtrs<SmEdge, SmFace> &rEToF = m_vEToF.GetFirstToSecond() ;
  SmTArray<SmEdge*> rEToFKeys;
  rEToF.GetAllKeys(rEToFKeys);

  for (ULONG i = 0; i < rEToFKeys.GetSize(); ++i)
  {
    SmEdge* pKey = rEToFKeys.GetAt(i);
    pKey->Draw();

    for (const SmFace* pValue : rEToF.At(pKey))
    {
      pValue->Draw();
    }
  }

  // draw m_vVToF
  const SmMapPtrToPtrs<SmVertex, SmFace> &rVToF = m_vVToF.GetFirstToSecond() ;
  SmTArray<SmVertex*> rVToFKeys;
  rVToF.GetAllKeys(rVToFKeys);

  for (ULONG i = 0; i < rVToFKeys.GetSize(); ++i)
  {
    SmVertex* pKey = rVToFKeys.GetAt(i);
    pKey->Draw();

    for (const SmFace* pValue : rVToF.At(pKey))
    {
      pValue->Draw();
    }
  }

  // draw m_vBlendToFaces
  const SmMapPtrToPtrs<SmFace, SmFace> &rBlendToFaces = m_vBlendToFaces.GetFirstToSecond() ;
  SmTArray<SmFace*> rBlendToFacesKeys;
  rBlendToFaces.GetAllKeys(rBlendToFacesKeys);

  for (ULONG i = 0; i < rBlendToFacesKeys.GetSize(); ++i)
  {
    SmFace* pKey = rBlendToFacesKeys.GetAt(i);
    pKey->Draw(SM_DM_CROSSHATCH);

    for (const SmFace* pValue : rBlendToFaces.At(pKey))
    {
      pValue->Draw(SM_DM_CROSSHATCH);
    }
  }

  // draw m_vSurfExtSurf
  SmTArray<SmSurface*> rSurfExtSurfKeys;
  m_vSurfExtSurf.GetAllKeys(rSurfExtSurfKeys);

  for (ULONG i = 0; i < rSurfExtSurfKeys.GetSize(); ++i)
  {
    SmSurface* pKey = rSurfExtSurfKeys.GetAt(i);
    SmSurface* pValue = m_vSurfExtSurf.At(pKey);
    pKey->DrawUV();
    pValue->DrawUV();
  }

  pRtn = smgfx_Close() ;

#endif // end SM_GFX_CODE
  return pRtn;

} // end SmOffsetExecutive::Draw

/*******************************************************************//**
PURPOSE:  Pretty Print SmTopologyIntersector Summary

NOTES: 
***********************************************************************/
void SmOffsetExecutive::Dump
  (SmBoolean bDumpMapObjects)       // NotUsed: in : TRUE = Dump every map object
                                    //      FALSE= don't
 const
{
  SM_REF1(bDumpMapObjects) ;
  // locals
  TCHAR sBuff[SM_TBLOCK_SIZE], sBuffForFile[SM_TBLOCK_SIZE];

  // brep and Other locals
  SmTArray<SmFace*>   sFaces,    sOFaces;
  SmTArray<SmEdge*>   sEdges,    sOEdges;
  SmTArray<SmVertex*> sVertices, sOVertices;
  if(m_pOriginalBrep) m_pOriginalBrep->GetFaces   (sFaces);     
  if(m_pOriginalBrep) m_pOriginalBrep->GetEdges   (sEdges);     
  if(m_pOriginalBrep) m_pOriginalBrep->GetVertices(sVertices);  

  if(m_pOffsetBrep) m_pOffsetBrep->GetFaces   (sOFaces);   
  if(m_pOffsetBrep) m_pOffsetBrep->GetEdges   (sOEdges);   
  if(m_pOffsetBrep) m_pOffsetBrep->GetVertices(sOVertices);

  smos_sprintf(sBuff,       _T("\n\nBegin Dump SmTopologyIntersector  = 0x%p"), this);
  smos_sprintf(sBuffForFile,_T("\n\nBegin Dump SmTopologyIntersector  = %s"), _T("notNULL"));
             smos_WriteBuffer(sBuff, sBuffForFile);

  // output control parameters
  smos_sprintf(sBuff,_T("\nm_eOffsetOperation = %s"),  
               m_eOffsetOperation == SM_OO_SOLID_OFFSET     ? _T("SOLID_OFFSET")
             : m_eOffsetOperation == SM_OO_NMT_OFFSET       ? _T("NMT_OFFSET")   
             : m_eOffsetOperation == SM_OO_SHELL            ? _T("SHELL")   
             : m_eOffsetOperation == SM_OO_LOCAL_OPERATION  ? _T("LOCAL_OPERATION")   
             : m_eOffsetOperation == SM_OO_MACHINE_OFFSET   ? _T("MACHINE_OFFSET")
             : _T("UNKNONWN_TYPE")) ;
             smos_WriteBuffer(sBuff);
  smos_sprintf(sBuff,_T("\nbInset             = %d,  bExtendConvexEdges = %d"),
             m_bInset, m_bExtendConvexEdges);  
             smos_WriteBuffer(sBuff);
  smos_sprintf(sBuff,_T("\nbMergeResults      = %d,  bCreateOffsetSolid = %d"),  
             m_bMergeResults, m_bCreateOffsetSolid);
             smos_WriteBuffer(sBuff);
  smos_sprintf(sBuff,_T("\nNumber of Faces to Shell (copy not offset)  = %ld"),
             m_pShellFaces ? m_pShellFaces->GetSize() : 0);
             smos_WriteBuffer(sBuff);

  // output m_pOriginalBrep arguments - with some entity counts
  smos_sprintf(sBuff,       _T("\n  Original Brep = 0x%p,  # Faces = %ld, # Edges = %ld, # Vertices = %ld"),
             m_pOriginalBrep,sFaces.GetSize(),sEdges.GetSize(),sVertices.GetSize());
  smos_sprintf(sBuffForFile,_T("\n  Original Brep = %s,  # Faces = %ld, # Edges = %ld, # Vertices = %ld"),
             m_pOriginalBrep ? _T("notNULL") : _T("NULL"),sFaces.GetSize(),sEdges.GetSize(),sVertices.GetSize());
             smos_WriteBuffer(sBuff, sBuffForFile);

  // output m_pOffsetBrep arguments - with some entity counts
  smos_sprintf(sBuff,       _T("\n  Output Brep   = 0x%p,  # Faces = %ld, # Edges = %ld, # Vertices = %ld"),
             m_pOffsetBrep,sOFaces.GetSize(),sOEdges.GetSize(),sOVertices.GetSize());
  smos_sprintf(sBuffForFile,_T("\n  output Brep   = %s,  # Faces = %ld, # Edges = %ld, # Vertices = %ld"),
             m_pOffsetBrep ? _T("notNULL") : _T("NULL"),sOFaces.GetSize(),sOEdges.GetSize(),sOVertices.GetSize());
             smos_WriteBuffer(sBuff, sBuffForFile);

  // output mapped entity lists
  smos_sprintf(sBuff,_T("%s"),_T("\n    Types: Region = 16012"  )); smos_WriteBuffer(sBuff);
  smos_sprintf(sBuff,_T("%s"),_T("\n           Face   = 16002"  )); smos_WriteBuffer(sBuff);
  smos_sprintf(sBuff,_T("%s"),_T("\n           CFace  = 16023"  )); smos_WriteBuffer(sBuff);
  smos_sprintf(sBuff,_T("%s"),_T("\n           Edge   = 16016"  )); smos_WriteBuffer(sBuff);
  smos_sprintf(sBuff,_T("%s"),_T("\n           CEdge  = 16017"  )); smos_WriteBuffer(sBuff);
  smos_sprintf(sBuff,_T("%s"),_T("\n           Vertex = 16021")); smos_WriteBuffer(sBuff);

  // output mapped entity list
  SM_ASSERT(m_vFUToF.GetFirstsCount() == m_vFUToF.GetSecondsCount()) ;
  smos_sprintf(sBuff,_T("\n # m_vFUToF Maps = %ld" ),m_vFUToF.GetFirstsCount());
             smos_WriteBuffer(sBuff);

  // Dump in format "type: [key] -> value"
  const SmMapPtrToPtrs<SmFaceuse, SmFace> &rFUToF = m_vFUToF.GetFirstToSecond() ;
  smos_WriteBuffer(_T("\n m_vFUToF MAP      : "));
  rFUToF.DumpAsObjects() ;

  const SmMapPtrToPtrs<SmEdgeuse, SmEdgeuse> &rEUToEU = m_vEUToEU.GetFirstToSecond() ;
  smos_WriteBuffer(_T(" m_vEUToEU MAP     : "));
  rEUToEU.DumpAsObjects() ;

  const SmMapPtrToPtrs<SmEdge, SmFace> &rEToF = m_vEToF.GetFirstToSecond() ;
  smos_WriteBuffer(_T(" m_vEToF MAP       : "));
  rEToF.DumpAsObjects() ;

  const SmMapPtrToPtrs<SmVertex, SmFace> &rVToF = m_vVToF.GetFirstToSecond() ;
  smos_WriteBuffer(_T(" m_vVToF MAP       : "));
  rVToF.DumpAsObjects() ;

  smos_WriteBuffer(_T(" m_vSurfExtSurf MAP: "));
  m_vSurfExtSurf.DumpAsObjects() ;

  // all done
  smos_WriteBuffer(_T("\nEnd Dump SmOffsetExecutive\n")) ;

} // end SmOffsetExecutive::Dump

