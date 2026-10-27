// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
// SPDX-License-Identifier: Apache-2.0

/*******************************************************************//**
* FILE NAME --- SmTObjArray.h
* PURPOSE: Header file for the SmTObjArray object.
**********************************************************************/

#ifndef __SMTOBJARRAY_H__
#define __SMTOBJARRAY_H__

#ifndef __SMCORE_TYPES_H__
#include <SmCoreTypes.h>
#endif

#ifndef __SMOS_MATH_H_
#include <SmMath.h>   // SM_BIG_ULONG, SM_UNDEF_ULONG, SM_ASSERT_BREAK (from SmMessages.h)
#endif

#include <vector>

/*******************************************************************//**
PURPOSE: A dynamic array of objects that constructs, copies and destroys
         its elements.

NOTES: SmTArray manages its memory with calloc/memset/memcpy/free and runs no
       constructors or destructors, so it holds only pointers and static types
       (see the WARNING in SmTArray.h). SmTObjArray is for the element types
       SmTArray cannot hold: those with virtual functions, and those with
       members that do (such as an SmTArray member, whose base SmObject has a
       virtual destructor), or that own memory.

       SetSize() default-constructs the elements it adds and destroys the ones
       it removes; ReSet() destroys every element. A copy of the array copies
       every element with its copy constructor.

       Its interface follows the SmTArray methods the callers use, so an
       SmTArray<TYPE> of such a type is replaced by an SmTObjArray<TYPE>
       without changing its callers. Allocation failure throws std::bad_alloc,
       as SMLib's operator new does.
***********************************************************************/
template<class TYPE> class SmTObjArray
{
 public:
    // number of elements in the array
    ULONG GetSize() const { return static_cast<ULONG>(m_vData.size()); }

    // set the number of elements: default-construct added elements, destroy removed ones.
    // RETURNS: old size, or SM_BIG_ULONG (array unchanged) for a size SmTArray also refuses.
    ULONG SetSize(ULONG lNewSize)
    {
      if (lNewSize >= SM_BIG_ULONG || lNewSize == SM_UNDEF_ULONG)
        {
          SM_ASSERT_BREAK(lNewSize < SM_BIG_ULONG && lNewSize != SM_UNDEF_ULONG) ;
          return SM_BIG_ULONG;
        }
      ULONG lOldSize = GetSize();
      m_vData.resize(lNewSize);
      return lOldSize;
    }

    // copy crNewElement to the end of the array
    // RETURNS: index of the new element, or SM_BIG_ULONG (array unchanged) when the array is full.
    ULONG Add(const TYPE & crNewElement)
    {
      ULONG lIndex = GetSize();
      if (lIndex + 1 >= SM_BIG_ULONG)
        {
          SM_ASSERT_BREAK(lIndex + 1 < SM_BIG_ULONG) ;
          return SM_BIG_ULONG;
        }
      m_vData.push_back(crNewElement);
      return lIndex;
    }

    // destroy every element, keep the allocated memory
    void ReSet() { m_vData.clear(); }

    // element access - the lIndex must be less than GetSize()
    TYPE & operator[](ULONG lIndex)
    {
      SM_ASSERT_BREAK(lIndex < GetSize());
      return m_vData[lIndex];
    }
    const TYPE & operator[](ULONG lIndex) const
    {
      SM_ASSERT_BREAK(lIndex < GetSize());
      return m_vData[lIndex];
    }

 private:
    std::vector<TYPE> m_vData;

} ; // end template<class TYPE> class SmTObjArray

#endif // __SMTOBJARRAY_H__
