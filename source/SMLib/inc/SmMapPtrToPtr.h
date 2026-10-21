// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
// SPDX-License-Identifier: Apache-2.0

/*******************************************************************//**
* FILE NAME --- SmMapPtrToPtr.h
* PURPOSE: Header file for SmMapPtrToPtr object.  
**********************************************************************/

/* Array Object definitions */
#ifndef __SMMAPPTRTOPTR_H__
#define __SMMAPPTRTOPTR_H__

#ifndef __SMCORE_TYPES_H__
#include <SmCoreTypes.h>
#endif

#ifndef __SMOBJECT_H__
#include <SmObject.h>
#endif

#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

/*******************************************************************/ /**
 PURPOSE: A set of pointers that lists its members in the order they were
 inserted.

 NOTES: The maps below hand out their keys and values in insertion order, never
 in the order of their addresses: an algorithm that walks them then does the
 same work whatever addresses the heap gives its objects. (std::unordered_set
 iterates in hash order, and the hash of a pointer is its address.)
 ***********************************************************************/
template<class PtrType>
class SmInsertionOrderedPtrSet
{
private:
  std::vector<PtrType>                m_vMembers;  // in insertion order
  std::unordered_map<PtrType, size_t> m_vIndex;    // member -> its place in m_vMembers

public:
  using const_iterator = typename std::vector<PtrType>::const_iterator;

  bool insert(const PtrType& rValue)
  {
    if (m_vIndex.find(rValue) != m_vIndex.end()) { return false; }
    m_vIndex.emplace(rValue, m_vMembers.size());
    m_vMembers.push_back(rValue);
    return true;
  }
  size_t erase(const PtrType& rValue)
  {
    auto it = m_vIndex.find(rValue);
    if (it == m_vIndex.end()) { return 0; }
    size_t lAt = it->second;
    m_vIndex.erase(it);
    m_vMembers.erase(m_vMembers.begin() + lAt);
    for (size_t ii = lAt; ii < m_vMembers.size(); ii++) { m_vIndex[m_vMembers[ii]] = ii; }
    return 1;
  }
  size_t count(const PtrType& rValue) const { return m_vIndex.count(rValue); }
  size_t size()  const { return m_vMembers.size(); }
  bool   empty() const { return m_vMembers.empty(); }
  const_iterator begin() const { return m_vMembers.begin(); }
  const_iterator end()   const { return m_vMembers.end(); }
};

/*******************************************************************/ /**
 PURPOSE: A templated class for maps between pointers of KeyTYPE and ValTYPE.

 NOTES: 

 ***********************************************************************/
template<class KeyTYPE, class ValTYPE>
class SM_EXPORT SmMapPtrToPtr final : public SmObject
{
private:
  static_assert(std::is_pointer<KeyTYPE>::value == false, _T(""));
  static_assert(std::is_pointer<ValTYPE>::value == false, _T(""));
  using KeyType     = KeyTYPE*;
  using ValueType   = ValTYPE*;

  using CRKeyType   = const KeyType&;
  using RValueType  = ValueType&;
  using CRValueType = const ValueType&;
  // Each key's value and the order in which the key was first inserted.
  struct Entry
  {
    ValueType          m_Value;
    unsigned long long m_lSeq;
  };
  // below uses SmAllocator to force memory allocation for maps through same allocators as the rest of the library.
  using MapType = std::unordered_map<KeyType, Entry, std::hash<KeyType>, std::equal_to<KeyType>, SmAllocator<std::pair<const KeyType, Entry>>> ;

  MapType* m_pMap; // not a value to avoid warning C4251
  unsigned long long m_lNextSeq = 0;

  // the map's entries in insertion order
  std::vector<typename MapType::const_iterator> Ordered() const
  {
    std::vector<typename MapType::const_iterator> sItems;
    sItems.reserve(m_pMap->size());
    for (auto it = m_pMap->begin(); it != m_pMap->end(); ++it) { sItems.push_back(it); }
    std::sort(sItems.begin(), sItems.end(),
              [](const typename MapType::const_iterator& a, const typename MapType::const_iterator& b)
              { return a->second.m_lSeq < b->second.m_lSeq; });
    return sItems;
  }
  Entry& Slot(CRKeyType rKey)
  {
    auto sFound = m_pMap->try_emplace(rKey, Entry{ nullptr, m_lNextSeq });
    if (sFound.second) { m_lNextSeq++; }
    return sFound.first->second;
  }

 public:
  SmMapPtrToPtr(const SmContext* cpContext = nullptr);
  SmMapPtrToPtr(const SmMapPtrToPtr&);
  virtual ~SmMapPtrToPtr();
  SmMapPtrToPtr& operator=(const SmMapPtrToPtr&);

  void Insert(CRKeyType rKey, CRValueType rValue) { Slot(rKey).m_Value = rValue; }
  // adds rOther's keys that this map lacks (keeping its own values), in rOther's order
  void Insert(const SmMapPtrToPtr& rOther)
  {
    for (const auto& it : rOther.Ordered())
    {
      if (m_pMap->find(it->first) == m_pMap->end()) { Slot(it->first).m_Value = it->second.m_Value; }
    }
  }
  void Remove(CRKeyType rKey)                     { m_pMap->erase(rKey); }

  ValueType At(CRKeyType rKey) const
  {
    auto it = m_pMap->find(rKey);
    if (it != m_pMap->end()) { return it->second.m_Value; }
    else                     { return nullptr; }
  }
  // inserts new element if rKey is not contained
  RValueType operator[](CRKeyType rKey)       { return Slot(rKey).m_Value; }
  ValueType  operator[](CRKeyType rKey) const { return At(rKey); }
  SmBoolean Contains(CRKeyType rKey) const { return m_pMap->find(rKey) != m_pMap->end(); }

  void RemoveAll(); // does not unset m_cpContext

  ULONG Count() const { return (ULONG)m_pMap->size(); }
  SmBoolean IsEmpty() const { return m_pMap->empty(); }

  // the first key inserted that maps to rVal
  KeyType Get1stKeyFor(CRValueType rVal) const
  {
      for (const auto& it : Ordered())
      {
          if (it->second.m_Value == rVal) { return it->first; }
      }
      return nullptr;
  };

  void GetAllKeys  (SmTArray<KeyType> & rKeys) const ;
  void GetAllValues(SmTArray<ValueType> & rValues) const ;
  void GetAllKeyValuePairs(SmTArray<KeyType> & rKeys, SmTArray<ValueType> & rValues) const;

  SM_COMMON(SmMapPtrToPtr,SmObject,SmMapPtrToPtr_TYPE);

  void Dump(SmBoolean bTerse) const override;
  void DumpAsObjects() const;

} ; // end class SmMapPtrToPtr

/*******************************************************************/ /**
 PURPOSE: A templated class for maps between pointers of KeyTYPE and
 unordered sets of  ValTYPE.

 NOTES:

 ***********************************************************************/

template<class KeyTYPE, class ValTYPE>
class SM_EXPORT SmMapPtrToPtrs final : public SmObject
{
private:
  static_assert(std::is_pointer<KeyTYPE>::value == false, _T(""));
  static_assert(std::is_pointer<ValTYPE>::value == false, _T(""));
  using KeyType = KeyTYPE*;
  using ValueType = ValTYPE*;

  using CRKeyType = const KeyType&;
  using RValueType = ValueType&;
  using CRValueType = const ValueType&;
  using ValueSetType = SmInsertionOrderedPtrSet<ValueType>;
  using RValueSetType = ValueSetType&;
  using CRValueSetType = const ValueSetType&;
  // Each key's values and the order in which the key was first inserted.
  struct Entry
  {
    ValueSetType       m_Values;
    unsigned long long m_lSeq;
  };
  // below uses SmAllocator to force memory allocation for maps through same allocators as the rest of the library.
  using MapType = std::unordered_map<KeyType, Entry, std::hash<KeyType>, std::equal_to<KeyType>, SmAllocator<std::pair<const KeyType, Entry>>> ;

  MapType* m_pMap; // not a value to avoid warning C4251
  unsigned long long m_lNextSeq = 0;

  // the map's entries in insertion order
  std::vector<typename MapType::const_iterator> Ordered() const
  {
    std::vector<typename MapType::const_iterator> sItems;
    sItems.reserve(m_pMap->size());
    for (auto it = m_pMap->begin(); it != m_pMap->end(); ++it) { sItems.push_back(it); }
    std::sort(sItems.begin(), sItems.end(),
              [](const typename MapType::const_iterator& a, const typename MapType::const_iterator& b)
              { return a->second.m_lSeq < b->second.m_lSeq; });
    return sItems;
  }
  Entry& Slot(CRKeyType rKey)
  {
    auto sFound = m_pMap->try_emplace(rKey, Entry{ ValueSetType(), m_lNextSeq });
    if (sFound.second) { m_lNextSeq++; }
    return sFound.first->second;
  }

 public:
  SmMapPtrToPtrs(const SmContext* cpContext = nullptr);
  SmMapPtrToPtrs(const SmMapPtrToPtrs&);
  virtual ~SmMapPtrToPtrs();
  SmMapPtrToPtrs& operator=(const SmMapPtrToPtrs&);

  void Insert(CRKeyType rKey, CRValueType rValue) { Slot(rKey).m_Values.insert(rValue); }
  // adds rOther's keys that this map lacks (keeping its own values), in rOther's order
  void Insert(const SmMapPtrToPtrs& rOther)
  {
    for (const auto& it : rOther.Ordered())
    {
      if (m_pMap->find(it->first) == m_pMap->end()) { Slot(it->first).m_Values = it->second.m_Values; }
    }
  };
  void Remove(CRKeyType rKey)                     { m_pMap->erase(rKey); }
  void Remove(CRKeyType rKey, CRValueType rValue)
  {
    auto it = m_pMap->find(rKey);
    if (it != m_pMap->end())
    {
      it->second.m_Values.erase(rValue);
      if (it->second.m_Values.empty())
      {
        m_pMap->erase(rKey);
      }
    }
  }

  CRValueSetType At(CRKeyType rKey) const
  {
    static const ValueSetType emptyValueSet{};
    auto it = m_pMap->find(rKey);
    if (it != m_pMap->end()) { return it->second.m_Values; }
    else                     { return emptyValueSet; }
  }
  // inserts new element if rKey is not contained
  RValueSetType operator[](CRKeyType rKey) { return Slot(rKey).m_Values; }
  CRValueSetType  operator[](CRKeyType rKey) const { return At(rKey); }
  SmBoolean Contains(CRKeyType rKey) const { return m_pMap->find(rKey) != m_pMap->end(); }

  void RemoveAll(); // does not unset m_cpContext

  ULONG Count() const { return (ULONG)m_pMap->size(); }
  SmBoolean IsEmpty() const { return m_pMap->empty(); }

  void GetAllKeys(SmTArray<KeyType>& rKeys) const;
  void GetAllValues(SmTArray<ValueType>& rValues) const;
  void GetAllKeyValuePairs(SmTArray<KeyType>& rKeys, SmTArray<ValueType>& rValues) const;

  SM_COMMON(SmMapPtrToPtrs,SmObject,SmMapPtrToPtrs_TYPE);

  void Dump(SmBoolean bTerse) const override;
  void DumpAsObjects() const;

}; // end class SmMapPtrToPtrs

template<class KeyTYPE, class ValTYPE>
SmMapPtrToPtr<KeyTYPE, ValTYPE>::SmMapPtrToPtr(const SmContext* cpContext)
{
    m_cpContext = cpContext;
    m_pMap = new MapType;
}

template<class KeyTYPE, class ValTYPE>
SmMapPtrToPtr<KeyTYPE, ValTYPE>::SmMapPtrToPtr(const SmMapPtrToPtr& other)
{
    m_cpContext = other.m_cpContext;
    m_pMap = new MapType(*other.m_pMap);
    m_lNextSeq = other.m_lNextSeq;
}

template<class KeyTYPE, class ValTYPE>
SmMapPtrToPtr<KeyTYPE, ValTYPE>::~SmMapPtrToPtr()
{
    m_cpContext = nullptr;
    delete m_pMap;
    m_pMap = nullptr;
}

template<class KeyTYPE, class ValTYPE>
SmMapPtrToPtr<KeyTYPE, ValTYPE>& SmMapPtrToPtr<KeyTYPE, ValTYPE>::operator=(const SmMapPtrToPtr<KeyTYPE, ValTYPE>& rOther)
{
    m_cpContext = rOther.m_cpContext;
    MapType* pCopy = new MapType(*rOther.m_pMap);
    delete m_pMap;
    m_pMap = pCopy;
    m_lNextSeq = rOther.m_lNextSeq;
    return *this;
}

template<class KeyTYPE, class ValTYPE>
void SmMapPtrToPtr<KeyTYPE, ValTYPE>::RemoveAll()
{
    m_pMap->clear();
}

template<class KeyTYPE, class ValTYPE>
void SmMapPtrToPtr<KeyTYPE, ValTYPE>::GetAllKeys(SmTArray<KeyType>& rKeys) const
{
    rKeys.ReSet();
    rKeys.SetDataSize(Count());
    for (const auto& it : Ordered())
    {
        rKeys.Add(it->first);
    }
}

template<class KeyTYPE, class ValTYPE>
void SmMapPtrToPtr<KeyTYPE, ValTYPE>::GetAllValues(SmTArray<ValueType>& rValues) const
{
    rValues.ReSet();
    rValues.SetDataSize(Count());
    for (const auto& it : Ordered())
    {
        rValues.Add(it->second.m_Value);
    }
}

template<class KeyTYPE, class ValTYPE>
void SmMapPtrToPtr<KeyTYPE, ValTYPE>::GetAllKeyValuePairs(SmTArray<KeyType>& rKeys, SmTArray<ValueType>& rValues) const
{
    rKeys.ReSet();
    rKeys.SetDataSize(Count());
    rValues.ReSet();
    rValues.SetDataSize(Count());
    for (const auto& it : Ordered())
    {
        rKeys.Add(it->first);
        rValues.Add(it->second.m_Value);
    }
}

template<class KeyTYPE, class ValTYPE>
SmBoolean SmMapPtrToPtr<KeyTYPE, ValTYPE>::IsKindOf(SM_TYPE t) const
{
    return ((SmMapPtrToPtr_TYPE == t) ? TRUE : SmObject::IsKindOf((t)));
}

template<class KeyTYPE, class ValTYPE>
void SmMapPtrToPtr<KeyTYPE, ValTYPE>::Dump() const
{
    Dump(TRUE);
}

template<class KeyTYPE, class ValTYPE>
void SmMapPtrToPtr<KeyTYPE, ValTYPE>::Dump
(SmBoolean bTerse)   // in : TRUE = omit pair list
//      FALSE= list pairs
const
{
    TCHAR sBuff[SM_TBLOCK_SIZE];

    smos_sprintf(sBuff, _T("\nSmMapPtrToPtr: NumElems[%zu]\n"), m_pMap->size());
    smos_WriteBuffer(sBuff);
    SmObject::Dump();

    if (bTerse == TRUE)
    {
        return;
    }

    for (const auto& it : Ordered())
    {
        smos_sprintf(sBuff, _T("\t [%p] = %p\n"), it->first, it->second.m_Value);
        smos_WriteBuffer(sBuff);
    }
}

template<class KeyTYPE, class ValTYPE>
void SmMapPtrToPtr<KeyTYPE, ValTYPE>::DumpAsObjects() const
{
    SM_ASSERT(m_pMap != nullptr);

    TCHAR sBuff[SM_TBLOCK_SIZE];

    smos_sprintf(sBuff, _T("\nSmMapPtrToPtr: NumElems[%zu]\n"), m_pMap->size());
    smos_WriteBuffer(sBuff);
    SmObject::Dump();

    for (const auto& it : Ordered())
    {
        SmObject* pKey = (SmObject*)it->first;
        SmObject* pValue = (SmObject*)it->second.m_Value;

        ULONG lKeyType = (ULONG)pKey->GetType();
        ULONG lValType = (ULONG)pValue->GetType();

        smos_sprintf(sBuff, _T("\t %s-%s : [%p] = %p\n"), lKeyType / 1000 == TOPO_BASE_TYPE / 1000 ? SM_TOPO_TYPENAME(lKeyType)
            : lKeyType / 1000 == SURF_BASE_TYPE / 1000 ? SM_SURF_TYPENAME(lKeyType)
            : lKeyType / 1000 == CURV_BASE_TYPE / 1000 ? SM_CURV_TYPENAME(lKeyType)
            : _T("Not A Named TYPE"),
            lValType / 1000 == TOPO_BASE_TYPE / 1000 ? SM_TOPO_TYPENAME(lValType)
            : lValType / 1000 == SURF_BASE_TYPE / 1000 ? SM_SURF_TYPENAME(lValType)
            : lValType / 1000 == CURV_BASE_TYPE / 1000 ? SM_CURV_TYPENAME(lValType)
            : _T("Not A Named TYPE"),
            pKey, pValue);
        smos_WriteBuffer(sBuff);
    }
}

template<class KeyTYPE, class ValTYPE>
SmMapPtrToPtrs<KeyTYPE, ValTYPE>::SmMapPtrToPtrs(const SmContext* cpContext)
{
    m_cpContext = cpContext;
    m_pMap = new MapType;
}

template<class KeyTYPE, class ValTYPE>
SmMapPtrToPtrs<KeyTYPE, ValTYPE>::SmMapPtrToPtrs(const SmMapPtrToPtrs& other)
{
    m_cpContext = other.m_cpContext;
    m_pMap = new MapType(*other.m_pMap);
    m_lNextSeq = other.m_lNextSeq;
}

template<class KeyTYPE, class ValTYPE>
SmMapPtrToPtrs<KeyTYPE, ValTYPE>::~SmMapPtrToPtrs()
{
    m_cpContext = nullptr;
    delete m_pMap;
    m_pMap = nullptr;
}

template<class KeyTYPE, class ValTYPE>
SmMapPtrToPtrs<KeyTYPE, ValTYPE>& SmMapPtrToPtrs<KeyTYPE, ValTYPE>::operator=(const SmMapPtrToPtrs<KeyTYPE, ValTYPE>& rOther)
{
    m_cpContext = rOther.m_cpContext;
    MapType* pCopy = new MapType(*rOther.m_pMap);
    delete m_pMap;
    m_pMap = pCopy;
    m_lNextSeq = rOther.m_lNextSeq;
    return *this;
}

template<class KeyTYPE, class ValTYPE>
void SmMapPtrToPtrs<KeyTYPE, ValTYPE>::RemoveAll()
{
    m_pMap->clear();
}

template<class KeyTYPE, class ValTYPE>
void SmMapPtrToPtrs<KeyTYPE, ValTYPE>::GetAllKeys(SmTArray<KeyType>& rKeys) const
{
    rKeys.ReSet();
    rKeys.SetDataSize(Count());
    for (const auto& it : Ordered())
    {
        rKeys.Add(it->first);
    }
}

template<class KeyTYPE, class ValTYPE>
void SmMapPtrToPtrs<KeyTYPE, ValTYPE>::GetAllValues(SmTArray<ValueType>& rValues) const
{
    rValues.ReSet();
    rValues.SetDataSize(Count());
    for (const auto& it : Ordered())
    {
        for (const auto& rValue : it->second.m_Values)
        {
            rValues.Add(rValue);
        }
    }
}

template<class KeyTYPE, class ValTYPE>
void SmMapPtrToPtrs<KeyTYPE, ValTYPE>::GetAllKeyValuePairs(SmTArray<KeyType>& rKeys, SmTArray<ValueType>& rValues) const
{
    rKeys.ReSet();
    rKeys.SetDataSize(Count());
    rValues.ReSet();
    rValues.SetDataSize(Count());
    for (const auto& it : Ordered())
    {
        rKeys.Add(it->first);
        for (const auto& rValue : it->second.m_Values)
        {
            rValues.Add(rValue);
        }
    }
}

template<class KeyTYPE, class ValTYPE>
SmBoolean SmMapPtrToPtrs<KeyTYPE, ValTYPE>::IsKindOf(SM_TYPE t) const
{
    return ((SmMapPtrToPtrs_TYPE == t) ? TRUE : SmObject::IsKindOf((t)));
}

template<class KeyTYPE, class ValTYPE>
void SmMapPtrToPtrs<KeyTYPE, ValTYPE>::Dump() const
{
    Dump(TRUE);
}

template<class KeyTYPE, class ValTYPE>
void SmMapPtrToPtrs<KeyTYPE, ValTYPE>::Dump
(SmBoolean bTerse)   // in : TRUE = omit pair list
//      FALSE= list pairs
const
{
    TCHAR sBuff[SM_TBLOCK_SIZE];

    smos_sprintf(sBuff, _T("\nSmMapPtrToPtrs: NumElems[%zu]\n"), m_pMap->size());
    smos_WriteBuffer(sBuff);
    SmObject::Dump();

    if (bTerse == TRUE)
    {
        return;
    }

    for (const auto& it : Ordered())
    {
        KeyType rKey = it->first;
        smos_sprintf(sBuff, _T("\t [%p] = {"), rKey);
        for (const auto& rValue : it->second.m_Values)
        {
            smos_sprintf(sBuff, _T("%p, "), rValue);
            smos_WriteBuffer(sBuff);
        }
        smos_sprintf(sBuff, _T("}\n"), rKey);
        smos_WriteBuffer(sBuff);
    }
}

template<class KeyTYPE, class ValTYPE>
void SmMapPtrToPtrs<KeyTYPE, ValTYPE>::DumpAsObjects() const
{
    SM_ASSERT(m_pMap != nullptr);

    TCHAR sBuff[SM_TBLOCK_SIZE];

    smos_sprintf(sBuff, _T("\nSmMapPtrToPtrs: NumElems[%zu]\n"), m_pMap->size());
    smos_WriteBuffer(sBuff);
    SmObject::Dump();

    for (const auto& it : Ordered())
    {
        KeyType rKey = it->first;
        SmObject* pKey = (SmObject*)rKey;
        ULONG lKeyType = (ULONG)pKey->GetType();

        smos_sprintf(sBuff, _T("\t %s : [%p] = {"), lKeyType / 1000 == TOPO_BASE_TYPE / 1000 ? SM_TOPO_TYPENAME(lKeyType)
            : lKeyType / 1000 == SURF_BASE_TYPE / 1000 ? SM_SURF_TYPENAME(lKeyType)
            : lKeyType / 1000 == CURV_BASE_TYPE / 1000 ? SM_CURV_TYPENAME(lKeyType)
            : _T("Not A Named TYPE"),
            pKey);

        for (const auto& rValue : it->second.m_Values)
        {
            SmObject* pValue = (SmObject*)rValue;
            ULONG lValType = (ULONG)pValue->GetType();
            smos_sprintf(sBuff, _T("%s : %p, "), lValType / 1000 == TOPO_BASE_TYPE / 1000 ? SM_TOPO_TYPENAME(lValType)
                : lValType / 1000 == SURF_BASE_TYPE / 1000 ? SM_SURF_TYPENAME(lValType)
                : lValType / 1000 == CURV_BASE_TYPE / 1000 ? SM_CURV_TYPENAME(lValType)
                : _T("Not A Named TYPE"),
                pValue);
            smos_WriteBuffer(sBuff);
        }
        smos_sprintf(sBuff, _T("}\n"), rKey);
        smos_WriteBuffer(sBuff);
    }
}



#endif // no __SMMAPPTRTOPTR_H__

