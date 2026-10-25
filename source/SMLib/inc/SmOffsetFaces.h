// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

/**********************************************************************//**
* FILE NAME --- SmOffsetFaces.h
* PURPOSE: SolidCraft change: moving faces of a solid along their normals,
*          the tangent blends on them following.
**********************************************************************/

#ifndef __SMOFFSETFACES_H__
#define __SMOFFSETFACES_H__

#ifndef __SMOS_TYPES_H__
#include <SmTypes.h>
#endif

#ifndef __SMMESSAGES_H__
#include <SmMessages.h>
#endif

class SmContext;
class SmBrep;
class SmFace;
template<class TYPE> class SmTArray;

/*******************************************************************//**
PURPOSE: Why SmOffsetFacesExecutive::OffsetFaces() refused.

NOTES:
***********************************************************************/
enum SmOffsetFacesRefusal
{
  SM_OFR_NONE,         // not refused (it succeeded, or failed for another reason: see the status)
  SM_OFR_NO_ROOM,      // a blend that follows a moved face no longer fits the faces it rolls on
  SM_OFR_UNSUPPORTED,  // a blend on a moved face SMLib cannot take off and put back as it was
  SM_OFR_NOT_EXACT,    // following a blend needs a surface or edge continued that SMLib cannot continue exactly
  SM_OFR_NO_EFFECT     // the local operation reported success but moved no face
};

/*******************************************************************//**
PURPOSE: SolidCraft change. Moves faces of a solid dDistance along their
   outward normals (negative: into the body), each onto its surface's
   offset, the faces around them extended or trimmed to meet them
   (SmBrep::LocalOperation), and the tangent blends on the moved faces
   following them.

NOTES:
  A blend follows a moved face when it is a face beside it across a tangent
  edge, not moved itself, that SmFace::IdentifyFillet recognises: two tangent
  rail edges on opposite isoparameter boundaries of its surface, one on the
  moved face, the other on a face it rolls on (moved or not), the two faces
  meeting at an angle (a half round between parallel faces, a slot's end,
  blends no edge and is not one). It is put back where the moved faces meet
  its other face, at its own size and of its own kind, tangent to both:

    the body = blend( move( unblend( body ) ) )

  1. Each followed blend is taken off, its faces met at the sharp edge it
     replaced (SmBrep::DeleteFilletFace's removals: SmRingFilletRemoval for a
     ring, SmFilletRemoval, exact only, for a four-sided one).
  2. The faces are moved (the local operation, as before).
  3. The sharp edges, where they still are, are blended again by the fillet
     executive in one call: a constant-radius round of the radius found
     (circular section, the radius constant), or a constant-chord round of
     the chord found (circular section, the chord constant, the radius not).
     A circular blend whose radius and chord are both constant cannot be
     told from geometry alone; it is put back as the round (the same surface
     wherever the faces still meet at a constant angle).
  4. Each new blend face carries the attributes of the blend it replaces,
     its rails those of the old rails on the same faces, its seam and end
     edges those of the old ones (SmAObject::Notify, SM_NO_COPY): what
     follows identity through SMLib's attributes (history tags) follows it.

  What is refused, and why:
  - A blend whose kind SMLib cannot recover exactly: a variable radius (the
    law along the edge), a curvature-continuous section, a section that is
    no circle: SM_OFR_UNSUPPORTED. So is a blend the removals cannot take
    off (its rails in pieces, its ends at corners of more than three edges,
    blends meeting at its ends in a way neither removal takes).
  - A blend whose removal needs a surface continued past its domain (other
    than a plane) or a curved edge continued past its curve:
    SM_OFR_NOT_EXACT.
  - A blend that no longer fits once the faces have moved: the fillet
    executive's own rules for the blend at that size, as for any blend
    (its rails off the faces they roll on, a concave blend overflowing, its
    rolling ball wider than a face's curve allows, the edge gone):
    SM_OFR_NO_ROOM. A convex blend overflowing a cliff is built by the cliff
    rule (blended as if the faces went on, the material it removes
    subtracted), as the executive does for any blend.
  - Anything the local operation refuses, with its status.

  On failure the Brep is left in an unspecified state (as for
  SmBrep::LocalOperation): work on a copy.
***********************************************************************/
class SM_EXPORT SmOffsetFacesExecutive
{
public:
  SmOffsetFacesExecutive( const SmContext & crContext, SmBrep * pBrep );
  ~SmOffsetFacesExecutive();

  // Move crFaces (faces of the Brep given to the constructor) dDistance along their outward
  // normals; the tangent blends on them follow.
  SmStatus OffsetFaces( const SmTArray<SmFace*> & crFaces, double dDistance );

  SmOffsetFacesRefusal GetRefusal()        const { return m_eRefusal; }
  const TCHAR *        GetReason()         const { return m_cReason; }
  ULONG                GetFollowedBlends() const { return m_lFollowed; }

private:
  SmStatus Refuse( SmOffsetFacesRefusal eRefusal, const TCHAR * cReason, const TCHAR * cMore = NULL );

  const SmContext    & m_crContext;
  SmBrep             * m_pBrep;
  SmOffsetFacesRefusal m_eRefusal;
  TCHAR                m_cReason[SM_TBLOCK_SIZE];
  ULONG                m_lFollowed;

}; // end class SmOffsetFacesExecutive

#endif // !__SMOFFSETFACES_H__
