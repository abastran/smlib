// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

/**********************************************************************//**
* FILE NAME --- SmOffsetFaces.cpp
* PURPOSE: SolidCraft change: moving faces of a solid along their normals,
*   the tangent blends on them following (SmOffsetFacesExecutive, see
*   SmOffsetFaces.h): each followed blend is taken off to the sharp edge it
*   replaced, the faces are moved by the local operation, and the sharp
*   edges are blended again at the blends' own sizes, of their own kinds.
**********************************************************************/

#include "StdAfx.h"

#include <SmOffsetFaces.h>

#include <SmAttribute.h>
#include <SmBrep.h>
#include <SmEdge.h>
#include <SmEdgeuse.h>
#include <SmFace.h>
#include <SmFaceuse.h>
#include <SmFilletExecutive.h>
#include <SmFilletRemoval.h>
#include <SmFeatureExecutive.h>
#include <SmPlane.h>
#include <SmOffsetGeometryCreation.h>
#include <SmRegion.h>
#include <SmShell.h>
#include <SmSSIData.h>
#include <SmSurface.h>

#include <cmath>
#include <initializer_list>
#include <string>
#include <utility>
#include <vector>

namespace
{

// The tangency tolerance the removals judge rails by (SmFilletRemoval::CollectFilletTopology).
const double s_dTangencyTolDeg = 2.0;

// What a mark says. A face or edge can carry several (a face both blends of a corner roll on).
const long s_lTarget = -1;   // a surface a moved face is to go onto
const long s_lMoved  = -3;   // a face to move
const long s_lBefore = -2;   // a face there before the blends were put back; an edge there before a group came off
// For followed blend i: (i * s_lStride) + role.
const long s_lStride  = 64;
const long s_lSupport = 0;   // + 0, 1: the face across rail 0, rail 1
const long s_lSharp   = 2;   // the sharp edge the blend became
const long s_lEnd     = 3;   // + k: the face across the blend's k-th other edge (an end face)

/*******************************************************************//**
PURPOSE: Marks on faces, edges and surfaces for the length of one
   OffsetFaces() call: which blend each face rolls on, which edge each
   blend became, which surfaces the moved faces go onto.

NOTES: SMLib hands an attribute on through every split, merge and copy it
   makes (the local operation's copies too) and builds new topology without
   any. A merge keeps what either object stood for.
***********************************************************************/
class SmFollowedBlendMark : public SmAttribute
{
public:
  std::vector<long> m_vCodes;

  explicit SmFollowedBlendMark( const std::vector<long> & crCodes )
    : SmAttribute( SM_AI_FOLLOWED_BLEND, SM_AB_COPY ), m_vCodes( crCodes ) {}

  SmAttribute * MakeCopy( const SmContext & crContext ) const override
    { return new (crContext) SmFollowedBlendMark( m_vCodes ); }
  SM_TYPE   GetType() const override { return SmFollowedBlendMark_TYPE; }
  SmBoolean IsKindOf( SM_TYPE t ) const override
    { return t == SmFollowedBlendMark_TYPE ? TRUE : SmAttribute::IsKindOf( t ); }

  void Merge( SmAObject * pA, SmAObject * pB, SmAObject * pResult ) override
  {
    std::vector<long> sAll;
    for ( SmAObject * pX : { pA, pB } )
      {
        SmFollowedBlendMark * pMark = On( pX );
        if ( pMark != NULL )
          { sAll.insert( sAll.end(), pMark->m_vCodes.begin(), pMark->m_vCodes.end() ); }
      }
    for ( long lCode : sAll )
      { Add( pResult, lCode ); }
  }

  static SmFollowedBlendMark * On( const SmAObject * pX )
  {
    SmAttribute * pA = ( pX == NULL ) ? NULL : pX->FindAttribute( SM_AI_FOLLOWED_BLEND );
    return ( pA != NULL && pA->GetType() == SmFollowedBlendMark_TYPE ) ? static_cast<SmFollowedBlendMark*>( pA ) : NULL;
  }
  static SmBoolean Has( const SmAObject * pX, long lCode )
  {
    SmFollowedBlendMark * pMark = On( pX );
    if ( pMark == NULL )
      { return FALSE; }
    for ( long lHave : pMark->m_vCodes )
      {
        if ( lHave == lCode )
          { return TRUE; }
      }
    return FALSE;
  }
  static void Add( SmAObject * pX, long lCode )
  {
    if ( pX == NULL || Has( pX, lCode ) )
      { return; }
    SmFollowedBlendMark * pMark = On( pX );
    if ( pMark != NULL )
      {
        pMark->m_vCodes.push_back( lCode );
        return;
      }
    pX->AddAttribute( new (*pX->GetContext()) SmFollowedBlendMark( std::vector<long>( 1, lCode ) ) );
  }
  static void Clear( SmAObject * pX )
  {
    SmFollowedBlendMark * pMark = On( pX );
    if ( pMark != NULL )
      { pX->RemoveAttribute( pMark ); }
  }
};

// The pieces joined into rText (SM_TBLOCK_SIZE characters), cut off where they don't fit.
void sm_SetText( TCHAR * rText, std::initializer_list<const TCHAR *> cPieces )
{
  std::basic_string<TCHAR> sAll;
  for ( const TCHAR * cPiece : cPieces )
    {
      if ( cPiece != NULL )
        { sAll += cPiece; }
    }
  size_t lLen = smos_Min( sAll.size(), static_cast<size_t>( SM_TBLOCK_SIZE - 1 ) );
  sAll.copy( rText, lLen );
  rText[lLen] = 0;
}

// The faceuse of pFace on the outside of the solid pBrep (bounding its infinite region).
SmFaceuse * sm_OuterFaceuse( const SmBrep * pBrep, const SmFace * pFace )
{
  SmFaceuse * pFU = pFace->GetUpwardFaceuse();
  if ( pFU == NULL )
    { return NULL; }
  if ( pFU->GetShell() == NULL || pFU->GetShell()->GetRegion() != pBrep->GetInfiniteRegion() )
    { pFU = pFU->GetMate(); }
  return pFU;
}

// The face across pEdge from pFace (a manifold edge), or NULL.
SmFace * sm_FaceAcross( const SmEdge * pEdge, const SmFace * pFace )
{
  SmTArray<SmFace*> sFaces;
  pEdge->GetFaces( sFaces );
  if ( sFaces.GetSize() != 2 )
    { return NULL; }
  return ( sFaces[0] == pFace ) ? sFaces[1] : ( sFaces[1] == pFace ? sFaces[0] : NULL );
}

// What an object stands for, kept while it is gone: copies of its attributes that are copied
// with it (SM_AB_COPY, SM_AB_STANDALONE_COPY: history tags are), these marks left out.
typedef std::vector<SmAttribute*> SmKept;

void sm_Keep( const SmAObject * pFrom, const SmContext & crContext, SmKept & rKept )
{
  SmTArray<SmAttribute*> sAttributes;
  pFrom->GetAttributes( sAttributes );
  for ( ULONG ii = 0; ii < sAttributes.GetSize(); ii++ )
    {
      SmAttribute * pA = sAttributes[ii];
      if ( pA->GetAttributeID() == SM_AI_FOLLOWED_BLEND ||
           !( pA->GetBehavior() == SM_AB_COPY || pA->GetBehavior() == SM_AB_STANDALONE_COPY ) )
        { continue; }
      rKept.push_back( pA->MakeCopy( crContext ) );
    }
}

// The kept attributes put on pTo, replacing any it has of the same kinds.
void sm_Put( const SmKept & crKept, SmAObject * pTo )
{
  if ( pTo == NULL )
    { return; }
  for ( SmAttribute * pA : crKept )
    {
      SmAttribute * pOld = pTo->FindAttribute( pA->GetAttributeID() );
      if ( pOld != NULL )
        { pTo->RemoveAttribute( pOld ); }
      pTo->AddAttribute( pA->MakeCopy( *pTo->GetContext() ));
    }
}

void sm_Free( SmKept & rKept )
{
  for ( SmAttribute * pA : rKept )
    { delete pA; }
  rKept.clear();
}

// One blend that follows a moved face.
struct SmFollowed
{
  SmFace                     * pBlend = NULL;          // the blend face (until it is taken off)
  SmEdge                     * aRails[2] = { NULL, NULL };
  SmFilletSolverType           eSolver = SM_FS_CONST_RADIUS;
  double                       dSize = 0.0;            // the radius, or the chord
  SmKept                       sFaceAttrs;             // the blend face's attributes, kept
  SmKept                       aRailAttrs[2];
  SmKept                       sSeamAttrs;             // its seam's
  std::vector<SmKept>          vEndAttrs;              // its k-th other edge's (with the face s_lEnd + k)
  std::vector<long>            vEndAcross;             // the followed blend across that edge, or -1
  SmBoolean                    bTakenOff = FALSE;
};

} // namespace

/*******************************************************************//**
PURPOSE: SmOffsetFacesExecutive constructor.

NOTES:
***********************************************************************/
SmOffsetFacesExecutive::SmOffsetFacesExecutive( const SmContext & crContext, SmBrep * pBrep )
  : m_crContext( crContext ),
    m_pBrep( pBrep ),
    m_eRefusal( SM_OFR_NONE ),
    m_lFollowed( 0 )
{
  m_cReason[0] = 0;
}

SmOffsetFacesExecutive::~SmOffsetFacesExecutive()
{
}

/*******************************************************************//**
PURPOSE: Note why the operation is refused; returns SM_ERR.

NOTES:
***********************************************************************/
SmStatus SmOffsetFacesExecutive::Refuse( SmOffsetFacesRefusal eRefusal, const TCHAR * cReason, const TCHAR * cMore )
{
  m_eRefusal = eRefusal;
  if ( cMore != NULL && cMore[0] != 0 )
    { sm_SetText( m_cReason, { cReason, _T(": "), cMore } ); }
  else
    { sm_SetText( m_cReason, { cReason } ); }
  return SM_ERR;
}

/*******************************************************************//**
PURPOSE: Move crFaces dDistance along their outward normals, the tangent
   blends on them following.

NOTES: See SmOffsetFaces.h.
***********************************************************************/
SmStatus SmOffsetFacesExecutive::OffsetFaces( const SmTArray<SmFace*> & crFaces, double dDistance )
{
  NER( m_pBrep );
  m_eRefusal  = SM_OFR_NONE;
  m_cReason[0] = 0;
  m_lFollowed = 0;
  if ( crFaces.GetSize() == 0 || !( smos_Fabs( dDistance ) > 0.0 ) )
    { SER( SM_ERR_INVALID_INPUT ); }
  SmTemporaryChangeValue<SmBoolean> sEditing( m_pBrep->m_bEditingEnabled, TRUE );
  const double dTol = m_pBrep->GetTolerance();

  // The kept attributes are freed whatever happens. (The marks are removed at the end of a
  // successful call; a failed one leaves the Brep unspecified, marks and all.)
  std::vector<SmFollowed> sBlends;
  struct FreeKept
  {
    std::vector<SmFollowed> & m_rBlends;
    ~FreeKept()
    {
      for ( SmFollowed & rB : m_rBlends )
        {
          sm_Free( rB.sFaceAttrs );
          sm_Free( rB.aRailAttrs[0] );
          sm_Free( rB.aRailAttrs[1] );
          sm_Free( rB.sSeamAttrs );
          for ( SmKept & rK : rB.vEndAttrs )
            { sm_Free( rK ); }
        }
    }
  } sFreeKept{ sBlends };

  // The faces to move, marked: taking blends off can replace a face by another piece of it.
  for ( ULONG ff = 0; ff < crFaces.GetSize(); ff++ )
    { SmFollowedBlendMark::Add( crFaces[ff], s_lMoved ); }

  // 1. The blends that follow: beside a moved face across a tangent edge, not moved
  //    themselves, recognised as fillets.
  for ( ULONG ff = 0; ff < crFaces.GetSize(); ff++ )
    {
      SmFace * pMoved = crFaces[ff];
      NER( pMoved );
      SmTArray<SmEdge*> sEdges;
      pMoved->GetEdges( sEdges );
      for ( ULONG ee = 0; ee < sEdges.GetSize(); ee++ )
        {
          SmEdge * pEdge = sEdges[ee];
          if ( !pEdge->IsManifold() || !pEdge->IsTangentEdge( s_dTangencyTolDeg ) )
            { continue; }
          SmFace * pBlend = sm_FaceAcross( pEdge, pMoved );
          if ( pBlend == NULL || crFaces.IsIn( pBlend ) )
            { continue; }
          SmBoolean bKnown = FALSE;
          for ( const SmFollowed & crB : sBlends )
            {
              if ( crB.pBlend == pBlend )
                { bKnown = TRUE; }
            }
          if ( bKnown )
            { continue; }

          SmBoolean bIsFillet = FALSE;
          ULONG lSection = 0, lSolver = 0;
          SmSurfParamType eRailDir = SM_SP_U;
          double dRadius = 0.0, dChord = 0.0;
          SmTArray<SmEdge*> sRailsMin, sRailsMax;
          SER( pBlend->IdentifyFillet( dTol, s_dTangencyTolDeg, TRUE, TRUE, bIsFillet, lSection, lSolver, eRailDir,
                                       dRadius, dChord, sRailsMin, sRailsMax ));
          // A tangent neighbour that is no blend (a piece of the same surface, a step) is the
          // local operation's, as before.
          if ( !bIsFillet || !( sRailsMin.IsIn( pEdge ) || sRailsMax.IsIn( pEdge ) ) )
            { continue; }
          if ( sRailsMin.GetSize() != 1 || sRailsMax.GetSize() != 1 )
            { return Refuse( SM_OFR_UNSUPPORTED, _T("a blend on a moved face has its contact line in pieces: it cannot be taken off and put back") ); }

          SmFollowed sB;
          sB.pBlend    = pBlend;
          sB.aRails[0] = sRailsMin[0];
          sB.aRails[1] = sRailsMax[0];
          SmFace * apSupports[2];
          for ( ULONG ii = 0; ii < 2; ii++ )
            {
              apSupports[ii] = sm_FaceAcross( sB.aRails[ii], pBlend );
              NER( apSupports[ii] );
            }
          if ( apSupports[0] == apSupports[1] )
            { return Refuse( SM_OFR_UNSUPPORTED, _T("a blend on a moved face rolls on one face on both sides") ); }

          // The faces it rolls on must meet at an angle: a half round between parallel faces
          // (a slot's end) blends no edge.
          {
            SmPoint3d sOn0, sOn1;
            SER( sB.aRails[0]->GetPrimaryEdgeuse()->NormalizedEvaluate( 0.5, FALSE, sOn0 ));
            SmBoolean bDropped = FALSE;
            double dParam = 0.0, dDist = 0.0;
            SER( sB.aRails[1]->GetCurve()->DropPoint( sB.aRails[1]->GetInterval(), sOn0, NULL, dTol, NULL, bDropped, dParam, dDist ));
            if ( !bDropped )
              { SER( SM_ERR ); }
            SER( sB.aRails[1]->GetCurve()->EvaluatePoint( dParam, sOn1 ));
            SmSurface * pSurf = pBlend->GetSurface();  NER( pSurf );
            SmVector3d asNormals[2];
            const SmPoint3d * apOn[2] = { &sOn0, &sOn1 };
            for ( ULONG ii = 0; ii < 2; ii++ )
              {
                SmBoolean bOn = FALSE, bMulti = FALSE;
                SmPoint2d sUV;
                double dGap = 0.0;
                SER( pSurf->DropPoint( *apOn[ii], pSurf->GetNaturalUVDomain(), NULL, bOn, sUV, dGap, bMulti ));
                if ( !bOn )
                  { SER( SM_ERR ); }
                SER( pSurf->EvaluateNormal( sUV, TRUE, TRUE, asNormals[ii] ));
              }
            double dCos = asNormals[0].Dot( asNormals[1] ) / ( asNormals[0].Length() * asNormals[1].Length() );
            if ( !( dCos > -std::cos( SM_DEG2RAD( s_dTangencyTolDeg ) ) ) || !( dCos < std::cos( SM_DEG2RAD( s_dTangencyTolDeg ) ) ) )
              { continue; }
          }

          // Its kind, as SMLib recovers it exactly: a circular section at constant radius (a
          // round), or at constant chord.
          if ( lSection != 0 )
            {
              return Refuse( SM_OFR_UNSUPPORTED, ( lSection >= 3 && lSection <= 5 )
                             ? _T("a curvature-continuous blend on a moved face: SMLib cannot recover its section exactly to put it back")
                             : _T("a blend on a moved face whose section is no circle: SMLib cannot recover it exactly to put it back") );
            }
          if ( lSolver == 0 && dRadius > 0.0 )
            {
              sB.eSolver = SM_FS_CONST_RADIUS;
              sB.dSize   = dRadius;
            }
          else if ( lSolver == 1 && dChord > 0.0 )
            {
              sB.eSolver = SM_FS_CONST_DIST;
              sB.dSize   = dChord;
            }
          else
            { return Refuse( SM_OFR_UNSUPPORTED, _T("a variable round on a moved face: SMLib cannot recover its radius law exactly to put it back") ); }
          if ( !SmFilletExecutive::IsRadiusSupported( sB.dSize, dTol ) )
            { return Refuse( SM_OFR_UNSUPPORTED, _T("a blend on a moved face is too small for SMLib to put back") ); }
          sBlends.push_back( sB );
        } // end iter every edge of a moved face
    } // end iter every moved face
  m_lFollowed = static_cast<ULONG>( sBlends.size() );

  // 2. What each followed blend stands for, kept, and marks on the faces around it.
  for ( ULONG ii = 0; ii < sBlends.size(); ii++ )
    {
      SmFollowed & rB = sBlends[ii];
      sm_Keep( rB.pBlend, m_crContext, rB.sFaceAttrs );
      for ( ULONG jj = 0; jj < 2; jj++ )
        {
          sm_Keep( rB.aRails[jj], m_crContext, rB.aRailAttrs[jj] );
          SmFollowedBlendMark::Add( sm_FaceAcross( rB.aRails[jj], rB.pBlend ), (long)ii * s_lStride + s_lSupport + (long)jj );
        }
      SmTArray<SmEdge*> sEdges;
      rB.pBlend->GetEdges( sEdges );
      for ( ULONG ee = 0; ee < sEdges.GetSize(); ee++ )
        {
          SmEdge * pEdge = sEdges[ee];
          if ( pEdge == rB.aRails[0] || pEdge == rB.aRails[1] )
            { continue; }
          SmFace * pAcross = sm_FaceAcross( pEdge, rB.pBlend );
          if ( pAcross == NULL )
            {
              // its seam (an edge whose only face it is)
              sm_Free( rB.sSeamAttrs );
              sm_Keep( pEdge, m_crContext, rB.sSeamAttrs );
              continue;
            }
          long lK = static_cast<long>( rB.vEndAttrs.size() );
          if ( lK >= s_lStride - s_lEnd )
            { return Refuse( SM_OFR_UNSUPPORTED, _T("a blend on a moved face has too many edges to put back") ); }
          rB.vEndAttrs.push_back( SmKept() );
          sm_Keep( pEdge, m_crContext, rB.vEndAttrs.back() );
          long lAcrossBlend = -1;
          for ( ULONG jj = 0; jj < sBlends.size(); jj++ )
            {
              if ( sBlends[jj].pBlend == pAcross )
                { lAcrossBlend = (long)jj; }
            }
          rB.vEndAcross.push_back( lAcrossBlend );
          SmFollowedBlendMark::Add( pAcross, (long)ii * s_lStride + s_lEnd + lK );
        }
    }

  // 3. Each followed blend taken off, its faces met at the sharp edge. A blend whose removal
  //    leaves another removable (mitred blends) is taken in a later round.
  SmBoolean bNotExact = FALSE;
  for ( SmBoolean bProgress = TRUE; bProgress; )
    {
      bProgress = FALSE;
      for ( ULONG ii = 0; ii < sBlends.size(); ii++ )
        {
          SmFollowed & rB = sBlends[ii];
          if ( rB.bTakenOff )
            { continue; }
          double dFaceTol = static_cast<double>( rB.pBlend->GetTolerance() );
          SmEdge * pSharp = NULL;
          if ( SmRingFilletRemoval::IsRing( rB.pBlend, s_dTangencyTolDeg ) )
            {
              SmRingFilletRemoval sRemoval( rB.pBlend, dFaceTol );
              if ( sRemoval.RemoveFace() != SM_SUCCESS )
                {
                  bNotExact = bNotExact || sRemoval.NotExact();
                  continue;
                }
              pSharp = sRemoval.GetSharpEdge();
            }
          else
            {
              SmFilletRemoval sRemoval( rB.pBlend, dFaceTol );
              sRemoval.SetExactOnly( TRUE );
              if ( sRemoval.RemoveFace() != SM_SUCCESS )
                {
                  bNotExact = bNotExact || sRemoval.NotExact();
                  continue;
                }
              pSharp = sRemoval.GetSharpEdge();
            }
          NER( pSharp );
          rB.pBlend    = NULL; // deleted
          rB.aRails[0] = rB.aRails[1] = NULL;
          rB.bTakenOff = TRUE;
          SmFollowedBlendMark::Add( pSharp, (long)ii * s_lStride + s_lSharp );
          bProgress = TRUE;
        }
    }
  // Blends that meet each other at their ends (mitred rounds along the edges of a face) cannot
  // come off one by one exactly: each one's end lies on the next one's surface. Each connected
  // group of them comes off at once (SmFeatureExecutive: the faces round the group extended
  // until they meet), where every face round it is a plane, which continues exactly.
  {
    std::vector<ULONG> sLeft;
    for ( ULONG ii = 0; ii < sBlends.size(); ii++ )
      {
        if ( !sBlends[ii].bTakenOff )
          { sLeft.push_back( ii ); }
      }
    std::vector<SmBoolean> sGrouped( sBlends.size(), FALSE );
    for ( ULONG lSeed : sLeft )
      {
        if ( sGrouped[lSeed] )
          { continue; }
        // the group: the blends left that share an edge with it
        std::vector<ULONG> sGroup( 1, lSeed );
        sGrouped[lSeed] = TRUE;
        for ( size_t gg = 0; gg < sGroup.size(); gg++ )
          {
            SmTArray<SmEdge*> sBound;
            sBlends[sGroup[gg]].pBlend->GetEdges( sBound );
            for ( ULONG ee = 0; ee < sBound.GetSize(); ee++ )
              {
                SmFace * pAcross = sm_FaceAcross( sBound[ee], sBlends[sGroup[gg]].pBlend );
                for ( ULONG lOther : sLeft )
                  {
                    if ( !sGrouped[lOther] && sBlends[lOther].pBlend == pAcross )
                      {
                        sGrouped[lOther] = TRUE;
                        sGroup.push_back( lOther );
                      }
                  }
              }
          }
        if ( sGroup.size() < 2 )
          { continue; }
        SmTArray<SmFace*> sFaces;
        for ( ULONG lB : sGroup )
          { sFaces.Add( sBlends[lB].pBlend ); }
        // every face round the group a plane
        SmBoolean bPlanes = TRUE;
        for ( ULONG ff = 0; ff < sFaces.GetSize() && bPlanes; ff++ )
          {
            SmTArray<SmEdge*> sBound;
            sFaces[ff]->GetEdges( sBound );
            for ( ULONG ee = 0; ee < sBound.GetSize() && bPlanes; ee++ )
              {
                SmFace * pAcross = sm_FaceAcross( sBound[ee], sFaces[ff] );
                if ( pAcross == NULL )
                  { bPlanes = FALSE; }
                else if ( !sFaces.IsIn( pAcross ) && !pAcross->GetSurface()->IsKindOf( SmPlane_TYPE ) )
                  { bPlanes = FALSE; }
              }
          }
        if ( !bPlanes )
          {
            bNotExact = TRUE;
            continue;
          }
        // the edges there before, but the group's own, keep their topology (faces split on
        // purpose stay apart; the faces round the group, extended, are joined across its edges)
        {
          SmTArray<SmEdge*> sAllEdges, sGroupEdges;
          m_pBrep->GetEdges( sAllEdges );
          for ( ULONG ff = 0; ff < sFaces.GetSize(); ff++ )
            {
              SmTArray<SmEdge*> sBound;
              sFaces[ff]->GetEdges( sBound );
              sGroupEdges.Append( sBound );
            }
          for ( ULONG ee = 0; ee < sAllEdges.GetSize(); ee++ )
            {
              if ( !sGroupEdges.IsIn( sAllEdges[ee] ) )
                { SmFollowedBlendMark::Add( sAllEdges[ee], s_lBefore ); }
            }
        }
        SmFeatureExecutive sDefeature( m_crContext, *m_pBrep );
        SmFeature * pFeature = sDefeature.AddFeature( sFaces, SM_RB_EXTEND_SURFACE, SM_RB_FROM_OTHER_BREP );
        if ( pFeature == NULL )
          { continue; }
        SmBrep * pBase = NULL;
        SmTArray<SmBrep*> sFeatureBreps;
        SmObjsDelete<SmBrep*> sCleanFeatures( &sFeatureBreps );
        if ( sDefeature.DoDefeature( pBase, &sFeatureBreps ) != SM_SUCCESS || pBase != m_pBrep || !pFeature->GetBaseWasApplied() )
          { return Refuse( SM_OFR_UNSUPPORTED, _T("mitred blends on a moved face could not be taken off together") ); }
        // the extended faces were split where the blends had been: join them again
        SmTArray<SmEdge*> sAllEdges, sKeep;
        m_pBrep->GetEdges( sAllEdges );
        for ( ULONG ee = 0; ee < sAllEdges.GetSize(); ee++ )
          {
            if ( SmFollowedBlendMark::Has( sAllEdges[ee], s_lBefore ) )
              { sKeep.Add( sAllEdges[ee] ); }
          }
        SER( m_pBrep->RemoveTopologicalEdgesAndVertices( &sKeep ));
        // as SmApiRemoveFillet's defeature: the rebuilt body valid, and taken in afresh
        if ( !m_pBrep->AssertValid( NULL, SM_LEVEL_2 ) || !m_pBrep->IsManifoldSolid() )
          { return Refuse( SM_OFR_UNSUPPORTED, _T("mitred blends on a moved face came off leaving no valid solid") ); }
        {
          SmBrep * pWork = new (m_crContext) SmBrep( *m_pBrep, NULL, NULL );
          SmObjDelete sCleanWork( pWork );
          SER( m_pBrep->RemoveAllTopology() );
          SER( m_pBrep->MergeBrep( *pWork ));
        }
        // each blend's sharp edge: where its two faces meet now
        SmTArray<SmEdge*> sNow;
        m_pBrep->GetEdges( sNow );
        for ( ULONG lB : sGroup )
          {
            const long lBase = (long)lB * s_lStride;
            SmBoolean bFound = FALSE;
            for ( ULONG ee = 0; ee < sNow.GetSize(); ee++ )
              {
                SmTArray<SmFace*> sOn;
                sNow[ee]->GetFaces( sOn );
                if ( sOn.GetSize() == 2 &&
                     ( ( SmFollowedBlendMark::Has( sOn[0], lBase + s_lSupport ) && SmFollowedBlendMark::Has( sOn[1], lBase + s_lSupport + 1 ) ) ||
                       ( SmFollowedBlendMark::Has( sOn[1], lBase + s_lSupport ) && SmFollowedBlendMark::Has( sOn[0], lBase + s_lSupport + 1 ) ) ) )
                  {
                    SmFollowedBlendMark::Add( sNow[ee], lBase + s_lSharp );
                    bFound = TRUE;
                  }
              }
            if ( !bFound )
              { return Refuse( SM_OFR_UNSUPPORTED, _T("mitred blends on a moved face came off without the edges they rounded") ); }
            sBlends[lB].pBlend    = NULL; // deleted
            sBlends[lB].aRails[0] = sBlends[lB].aRails[1] = NULL;
            sBlends[lB].bTakenOff = TRUE;
          }
      }
  }
  for ( const SmFollowed & crB : sBlends )
    {
      if ( !crB.bTakenOff )
        {
          return bNotExact
            ? Refuse( SM_OFR_NOT_EXACT, _T("a blend on a moved face could not be taken off: its faces would have to be continued past their surfaces or edges, which SMLib cannot do exactly") )
            : Refuse( SM_OFR_UNSUPPORTED, _T("a blend on a moved face could not be taken off its faces (its ends meet other blends or corners SMLib cannot take apart)") );
        }
    }
  if ( !sBlends.empty() && !m_pBrep->IsManifoldSolid() )
    { return Refuse( SM_OFR_UNSUPPORTED, _T("taking the blends off the moved faces left no closed solid") ); }

  // 4. The faces moved: each onto its surface's offset, along its outward normal (inward from
  //    its mate), the faces around extended or trimmed to meet it.
  {
    SmOffsetGeometryCreation sMake( smos_Fabs( dDistance ), dTol );
    SmTArray<SmFace*>    sMoved;
    {
      SmTArray<SmFace*> sAll;
      m_pBrep->GetFaces( sAll );
      for ( ULONG ff = 0; ff < sAll.GetSize(); ff++ )
        {
          if ( SmFollowedBlendMark::Has( sAll[ff], s_lMoved ) )
            { sMoved.Add( sAll[ff] ); }
        }
    }
    SmTArray<SmFace*>    sFaces;
    SmTArray<SmSurface*> sSurfaces;
    for ( ULONG ff = 0; ff < sMoved.GetSize(); ff++ )
      {
        SmFaceuse * pFU = sm_OuterFaceuse( m_pBrep, sMoved[ff] );
        if ( pFU != NULL && dDistance < 0.0 )
          { pFU = pFU->GetMate(); }
        SmSurface * pSurf = NULL;
        SmSSIData   sCrossings;
        SmStatus eStat = ( pFU == NULL ) ? SM_ERR_INVALID_INPUT : sMake.FaceuseOffset( m_crContext, pFU, TRUE, pSurf, sCrossings );
        if ( eStat != SM_SUCCESS || pSurf == NULL )
          {
            delete pSurf;
            for ( ULONG kk = 0; kk < sSurfaces.GetSize(); kk++ )
              { delete sSurfaces[kk]; }
            SER( eStat != SM_SUCCESS ? eStat : SM_ERR );
          }
        SmFollowedBlendMark::Add( pSurf, s_lTarget );
        sFaces.Add( sMoved[ff] );
        sSurfaces.Add( pSurf );
      }
    // The local operation takes the surfaces over, also when it fails.
    SER( m_pBrep->LocalOperation( sFaces, sSurfaces ));

    // Some face must lie on one of the surfaces, or the operation moved none.
    SmTArray<SmFace*> sAll;
    m_pBrep->GetFaces( sAll );
    SmBoolean bMoved = FALSE;
    for ( ULONG ff = 0; ff < sAll.GetSize(); ff++ )
      {
        if ( SmFollowedBlendMark::Has( sAll[ff]->GetSurface(), s_lTarget ) )
          { bMoved = TRUE; }
        SmFollowedBlendMark::Clear( sAll[ff]->GetSurface() );
      }
    if ( !bMoved )
      { return Refuse( SM_OFR_NO_EFFECT, _T("the local operation reported success but moved no face onto its new surface") ); }
  }
  if ( sBlends.empty() )
    {
      SER( m_pBrep->RemoveAttributeFromTopology( SM_AI_FOLLOWED_BLEND ));
      return SM_SUCCESS;
    }

  // 5. The sharp edges blended again, all in one call, each at its blend's size and of its kind.
  SmTArray<SmEdge*>                    sEdges;
  SmTArray<SmFilletSurfaceGeneratorType> sSections;
  SmTArray<SmFilletSolverType>         sSolvers;
  SmTArray<double>                     sSizes;
  {
    SmTArray<SmEdge*> sAll;
    m_pBrep->GetEdges( sAll );
    for ( ULONG ii = 0; ii < sBlends.size(); ii++ )
      {
        ULONG lFound = 0;
        for ( ULONG ee = 0; ee < sAll.GetSize(); ee++ )
          {
            SmEdge * pEdge = sAll[ee];
            if ( !SmFollowedBlendMark::Has( pEdge, (long)ii * s_lStride + s_lSharp ) )
              { continue; }
            SmTArray<SmFace*> sOn;
            pEdge->GetFaces( sOn );
            // still where the blend's two faces meet
            if ( sOn.GetSize() != 2 || !pEdge->IsManifold() ||
                 !( ( SmFollowedBlendMark::Has( sOn[0], (long)ii * s_lStride + s_lSupport )     && SmFollowedBlendMark::Has( sOn[1], (long)ii * s_lStride + s_lSupport + 1 ) ) ||
                    ( SmFollowedBlendMark::Has( sOn[0], (long)ii * s_lStride + s_lSupport + 1 ) && SmFollowedBlendMark::Has( sOn[1], (long)ii * s_lStride + s_lSupport ) ) ) )
              { continue; }
            sEdges.Add( pEdge );
            sSections.Add( SM_FSG_CIRCULAR );
            sSolvers.Add( sBlends[ii].eSolver );
            sSizes.Add( sBlends[ii].dSize );
            lFound++;
          }
        if ( lFound == 0 )
          { return Refuse( SM_OFR_NO_ROOM, _T("the move took a blend's faces apart: the edge it rolls along is gone, so the blend no longer fits") ); }
      }
  }
  {
    SmTArray<SmFace*> sAll;
    m_pBrep->GetFaces( sAll );
    for ( ULONG ff = 0; ff < sAll.GetSize(); ff++ )
      { SmFollowedBlendMark::Add( sAll[ff], s_lBefore ); }
  }
  {
    SmFilletExecutive sExec( m_crContext, m_pBrep );
    SmStatus eStat = sExec.SetFilletParameters( m_pBrep, sEdges, sSections, sSolvers, sSizes, NULL, -1.0, -1.0, -1.0, -1.0, 1, 1.0 );
    if ( eStat == SM_SUCCESS )
      { eStat = sExec.DoFillet(); }
    if ( eStat != SM_SUCCESS )
      {
        SmFilletErrorInfo * pInfo = sExec.GetFilletErrorInfo();
        const TCHAR * cWhy = ( pInfo != NULL ) ? pInfo->m_cComment : NULL;
        if ( pInfo != NULL && pInfo->GetErrorCode() == SM_FILERR_NO_ROOM )
          { return Refuse( SM_OFR_NO_ROOM, _T("a blend that follows a moved face no longer fits"), cWhy ); }
        if ( pInfo != NULL && pInfo->GetErrorCode() == SM_FILERR_NOT_EXACT )
          { return Refuse( SM_OFR_NOT_EXACT, _T("a blend that follows a moved face overflows where SMLib cannot continue its faces exactly"), cWhy ); }
        if ( cWhy != NULL && cWhy[0] != 0 )
          { sm_SetText( m_cReason, { _T("putting back a blend that follows a moved face: "), cWhy } ); }
        SER( eStat );
      }
  }
  SER( m_pBrep->StitchAndOrient() );

  // 6. The new blend faces stand for the blends they replace, their edges for the old ones: a
  //    rail for the rail on the same face, a seam for the seam, an end edge for the end edge
  //    on the same face (or, between two blends, the edge where those two blends met).
  {
    SmTArray<SmFace*> sAll;
    m_pBrep->GetFaces( sAll );
    std::vector<std::pair<SmFace*, ULONG>> sNewFaces;
    for ( ULONG ff = 0; ff < sAll.GetSize(); ff++ )
      {
        SmFace * pNew = sAll[ff];
        if ( SmFollowedBlendMark::Has( pNew, s_lBefore ) )
          { continue; }
        SmTArray<SmEdge*> sBound;
        pNew->GetEdges( sBound );
        for ( ULONG ii = 0; ii < sBlends.size(); ii++ )
          {
            const long lBase = (long)ii * s_lStride;
            // the blend's: beside both its faces
            SmBoolean abBeside[2] = { FALSE, FALSE };
            for ( ULONG ee = 0; ee < sBound.GetSize(); ee++ )
              {
                SmFace * pAcross = sm_FaceAcross( sBound[ee], pNew );
                for ( ULONG jj = 0; jj < 2; jj++ )
                  {
                    if ( SmFollowedBlendMark::Has( pAcross, lBase + s_lSupport + (long)jj ) )
                      { abBeside[jj] = TRUE; }
                  }
              }
            if ( abBeside[0] && abBeside[1] )
              {
                sNewFaces.push_back( std::make_pair( pNew, ii ));
                break;
              }
          }
      }
    for ( const std::pair<SmFace*, ULONG> & crNew : sNewFaces )
      {
        SmFace * pNew = crNew.first;
        const ULONG ii = crNew.second;
        const long lBase = (long)ii * s_lStride;
        SmFollowed & rB = sBlends[ii];
        sm_Put( rB.sFaceAttrs, pNew );
        SmTArray<SmEdge*> sBound;
        pNew->GetEdges( sBound );
        for ( ULONG ee = 0; ee < sBound.GetSize(); ee++ )
          {
            SmEdge * pEdge = sBound[ee];
            SmFace * pAcross = sm_FaceAcross( pEdge, pNew );
            if ( pAcross == NULL )
              {
                sm_Put( rB.sSeamAttrs, pEdge );
                continue;
              }
            for ( ULONG jj = 0; jj < 2; jj++ )
              {
                if ( SmFollowedBlendMark::Has( pAcross, lBase + s_lSupport + (long)jj ) )
                  { sm_Put( rB.aRailAttrs[jj], pEdge ); }
              }
            // the blend across, if the face across is one of the new blend faces
            long lAcrossBlend = -1;
            for ( const std::pair<SmFace*, ULONG> & crOther : sNewFaces )
              {
                if ( crOther.first == pAcross )
                  { lAcrossBlend = (long)crOther.second; }
              }
            for ( ULONG kk = 0; kk < rB.vEndAttrs.size(); kk++ )
              {
                SmBoolean bSame = ( lAcrossBlend >= 0 ) ? ( rB.vEndAcross[kk] == lAcrossBlend )
                                                        : SmFollowedBlendMark::Has( pAcross, lBase + s_lEnd + (long)kk );
                if ( bSame )
                  { sm_Put( rB.vEndAttrs[kk], pEdge ); }
              }
          }
      }
  }
  SER( m_pBrep->RemoveAttributeFromTopology( SM_AI_FOLLOWED_BLEND ));
  return SM_SUCCESS;

} // end SmOffsetFacesExecutive::OffsetFaces
