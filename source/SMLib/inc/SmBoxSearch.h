// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

/*******************************************************************//**
* FILE NAME --- SmBoxSearch.h
* PURPOSE: A static tree of axis-aligned boxes that finds the boxes near a
*          query box, so pairwise tests over n entities need not visit all
*          n^2 pairs.
**********************************************************************/

#ifndef __SMBOXSEARCH_H__
#define __SMBOXSEARCH_H__

#include <SmExtent3d.h>

#include <algorithm>
#include <cmath>
#include <vector>

/*******************************************************************//**
PURPOSE: Boxes of a list of items (vertices, edges), searched for the items
         whose box a query box may come near: the candidates of the pairwise
         tests of the validity checks (SmBrep::ContainsCoincidentTopology,
         SmFace's coincident vertex and edge tests).

NOTES: A query returns, in ascending item order, a superset of the items whose
       box lies within dPad of the query box (it pads by 2 * dPad and a little
       for rounding), and every item whose box is not finite or not known.
       The caller then runs its own test on each, in the order the
       all-pairs loop would have, so the results stay the same; only pairs
       too far apart for that test to pass are no longer visited.
***********************************************************************/
class SmBoxSearch
{
public:
  // Item i has box rBoxes[i]; rKnown[i] == false makes it a candidate of every query.
  void Build(const std::vector<SmExtent3d> & rBoxes, const std::vector<bool> & rKnown)
  {
    const size_t n = rBoxes.size() ;
    m_vLo.assign(3 * n, 0.0) ;
    m_vHi.assign(3 * n, 0.0) ;
    m_vItems.clear() ;
    m_vAlways.clear() ;
    m_vNodes.clear() ;
    m_dScale = 1.0 ;
    for (size_t i = 0; i < n; i++)
      {
        const SmPoint3d sLo = rBoxes[i].GetMin(), sHi = rBoxes[i].GetMax() ;
        const double aLo[3] = { sLo.x, sLo.y, sLo.z }, aHi[3] = { sHi.x, sHi.y, sHi.z } ;
        bool bFinite = i < rKnown.size() && rKnown[i] ;
        for (int k = 0; k < 3 && bFinite; k++)
          { bFinite = std::isfinite(aLo[k]) && std::isfinite(aHi[k]) ; }
        if (!bFinite)
          { m_vAlways.push_back(static_cast<ULONG>(i)) ; continue ; }
        for (int k = 0; k < 3; k++)
          {
            m_vLo[3 * i + k] = aLo[k] ;
            m_vHi[3 * i + k] = aHi[k] ;
            m_dScale = std::max(m_dScale, std::max(std::fabs(aLo[k]), std::fabs(aHi[k]))) ;
          }
        m_vItems.push_back(static_cast<ULONG>(i)) ;
      }
    if (m_vItems.empty())
      { return ; }
    // Median splits along the longest side of the centres' box, leaves of a few items.
    struct Range { size_t begin, end, node ; } ;
    std::vector<Range> sStack ;
    m_vNodes.push_back(Node()) ;
    sStack.push_back(Range{ 0, m_vItems.size(), 0 }) ;
    while (!sStack.empty())
      {
        const Range r = sStack.back() ;
        sStack.pop_back() ;
        Node sNode ;
        sNode.begin = r.begin ;
        sNode.end = r.end ;
        double aCLo[3] = { SM_BIG_DOUBLE, SM_BIG_DOUBLE, SM_BIG_DOUBLE } ;
        double aCHi[3] = { -SM_BIG_DOUBLE, -SM_BIG_DOUBLE, -SM_BIG_DOUBLE } ;
        for (int k = 0; k < 3; k++)
          { sNode.lo[k] = SM_BIG_DOUBLE ; sNode.hi[k] = -SM_BIG_DOUBLE ; }
        for (size_t j = r.begin; j < r.end; j++)
          {
            const ULONG i = m_vItems[j] ;
            for (int k = 0; k < 3; k++)
              {
                sNode.lo[k] = std::min(sNode.lo[k], m_vLo[3 * i + k]) ;
                sNode.hi[k] = std::max(sNode.hi[k], m_vHi[3 * i + k]) ;
                const double c = 0.5 * (m_vLo[3 * i + k] + m_vHi[3 * i + k]) ;
                aCLo[k] = std::min(aCLo[k], c) ;
                aCHi[k] = std::max(aCHi[k], c) ;
              }
          }
        if (r.end - r.begin > kLeafSize)
          {
            int iAxis = 0 ;
            for (int k = 1; k < 3; k++)
              { if (aCHi[k] - aCLo[k] > aCHi[iAxis] - aCLo[iAxis]) iAxis = k ; }
            const size_t mid = r.begin + (r.end - r.begin) / 2 ;
            const std::vector<double> & rLo = m_vLo, & rHi = m_vHi ;
            std::nth_element(m_vItems.begin() + r.begin, m_vItems.begin() + mid, m_vItems.begin() + r.end,
                             [&rLo, &rHi, iAxis](ULONG a, ULONG b)
                               { return rLo[3 * a + iAxis] + rHi[3 * a + iAxis] < rLo[3 * b + iAxis] + rHi[3 * b + iAxis] ; }) ;
            sNode.left = m_vNodes.size() ;
            m_vNodes.push_back(Node()) ;
            sNode.right = m_vNodes.size() ;
            m_vNodes.push_back(Node()) ;
            sStack.push_back(Range{ r.begin, mid, sNode.left }) ;
            sStack.push_back(Range{ mid, r.end, sNode.right }) ;
          }
        m_vNodes[r.node] = sNode ;
      }
  }

  // Sets rOut to the candidates of crBox at distance dPad (see the class notes). Returns
  // false, with rOut empty, when crBox or dPad is not finite: the caller must then test all.
  bool Query(const SmExtent3d & crBox, double dPad, std::vector<ULONG> & rOut) const
  {
    rOut.clear() ;
    const SmPoint3d sLo = crBox.GetMin(), sHi = crBox.GetMax() ;
    const double dGrow = 2.0 * dPad + 1.0e-9 * m_dScale ;
    const double aLo[3] = { sLo.x - dGrow, sLo.y - dGrow, sLo.z - dGrow } ;
    const double aHi[3] = { sHi.x + dGrow, sHi.y + dGrow, sHi.z + dGrow } ;
    for (int k = 0; k < 3; k++)
      {
        if (!std::isfinite(aLo[k]) || !std::isfinite(aHi[k]))
          { return false ; }
      }
    std::vector<size_t> sStack ;
    if (!m_vNodes.empty())
      { sStack.push_back(0) ; }
    while (!sStack.empty())
      {
        const Node & rNode = m_vNodes[sStack.back()] ;
        sStack.pop_back() ;
        if (   rNode.lo[0] > aHi[0] || rNode.lo[1] > aHi[1] || rNode.lo[2] > aHi[2]
            || rNode.hi[0] < aLo[0] || rNode.hi[1] < aLo[1] || rNode.hi[2] < aLo[2])
          { continue ; }
        if (rNode.left != 0)
          {
            sStack.push_back(rNode.left) ;
            sStack.push_back(rNode.right) ;
            continue ;
          }
        for (size_t j = rNode.begin; j < rNode.end; j++)
          {
            const ULONG i = m_vItems[j] ;
            if (   m_vLo[3 * i] > aHi[0] || m_vLo[3 * i + 1] > aHi[1] || m_vLo[3 * i + 2] > aHi[2]
                || m_vHi[3 * i] < aLo[0] || m_vHi[3 * i + 1] < aLo[1] || m_vHi[3 * i + 2] < aLo[2])
              { continue ; }
            rOut.push_back(i) ;
          }
      }
    rOut.insert(rOut.end(), m_vAlways.begin(), m_vAlways.end()) ;
    std::sort(rOut.begin(), rOut.end()) ;
    return true ;
  }

private:
  static const size_t kLeafSize = 8 ;
  struct Node
  {
    double lo[3] = { 0, 0, 0 }, hi[3] = { 0, 0, 0 } ;
    size_t begin = 0, end = 0, left = 0, right = 0 ; // left == 0: a leaf
  } ;
  std::vector<double> m_vLo, m_vHi ;
  std::vector<ULONG>  m_vItems, m_vAlways ;
  std::vector<Node>   m_vNodes ;
  double              m_dScale = 1.0 ;
} ;

// [lFrom, lCount) as a candidate list: what a query gives when it cannot be made.
inline void sm_AllFrom(ULONG lFrom, ULONG lCount, std::vector<ULONG> & rOut)
{
  rOut.clear() ;
  for (ULONG i = lFrom; i < lCount; i++)
    { rOut.push_back(i) ; }
}

class SmVertex ;

/*******************************************************************//**
PURPOSE: A face's vertices, each once in the order SmFace::GetVertices lists
         them (loop by loop), with a box search over their points and their
         largest tolerance: SmFace::PointClassify's vertex test, answered
         without walking the whole face (SmTopologyListWalkScope::FaceVertices).
***********************************************************************/
struct SmFaceVertexIndex
{
  std::vector<SmVertex*> m_vVertices ;
  SmBoxSearch            m_sSearch ;
  double                 m_dMaxVertexTol = 0.0 ;
  SmBoolean              m_bOK = FALSE ;  // FALSE: a point or tolerance was not finite; walk the face
  const void           * m_pFaceuse = NULL ; // the face's upward faceuse when gathered
} ;

#endif // !__SMBOXSEARCH_H__
