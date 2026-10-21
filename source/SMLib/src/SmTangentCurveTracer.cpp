// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
// SPDX-License-Identifier: Apache-2.0

/*******************************************************************//**
* FILE NAME --- SmTangentCurveTracer.cpp
* PURPOSE: SolidCraft change: tracing of tangential surface/surface
*          intersection curves (curves along which the two surfaces touch,
*          n1 x n2 = 0 everywhere on the curve, including high order
*          contact where the sheets cross).
*
* The regular marcher steps along n1 x n2, which vanishes on such a curve,
* and its corrector (a surface/surface/plane Newton solve) is singular there,
* so a tangential curve used to come out as its two end points, or as a
* fragment.  Here a tangential curve is traced as the curve of points where
* the surfaces meet AND their normals are parallel:
*
*   unknowns  x = (u1, v1, u2, v2)
*   residual  S1(u1,v1) - S2(u2,v2)              (3)
*             L (n1 - s n2), s = +-1              (3)
*             plus one or two constraints that pick a single point of the
*             curve: the two surface points on a plane through a predicted
*             point (marching), or one parameter fixed on a domain boundary
*             (the curve's end).
*
* solved by damped Gauss-Newton (Levenberg-Marquardt), which converges on
* this rank deficient system (linearly for odd/high order contact).  Start
* directions at a contact point are found by probing a small circle in the
* tangent plane with the plane-constrained solve: a tangential curve answers
* in two opposite directions, a coincident patch in all, a touching point in
* none.
***********************************************************************/

#include "StdAfx.h"

#include <SmSurfaceIntersector.h>
#include <SmTree.h>
#include <vector>
#include <algorithm>
#include <cmath>

namespace {

struct TcV
{
  double x, y, z;
  TcV() : x(0), y(0), z(0) {}
  TcV(double a, double b, double c) : x(a), y(b), z(c) {}
  TcV(const SmVector3d & v) : x(v.x), y(v.y), z(v.z) {}
  TcV operator+(const TcV & o) const { return TcV(x+o.x, y+o.y, z+o.z); }
  TcV operator-(const TcV & o) const { return TcV(x-o.x, y-o.y, z-o.z); }
  TcV operator*(double d)      const { return TcV(x*d, y*d, z*d); }
  double dot(const TcV & o)    const { return x*o.x + y*o.y + z*o.z; }
  TcV cross(const TcV & o)     const { return TcV(y*o.z - z*o.y, z*o.x - x*o.z, x*o.y - y*o.x); }
  double len()                 const { return std::sqrt(dot(*this)); }
  TcV unit()                   const { double l = len(); return l > 0.0 ? (*this) * (1.0/l) : *this; }
};

// one surface evaluation with what the solver needs
struct TcSrfEval
{
  TcV P, Su, Sv, Suu, Suv, Svv;
  TcV n;          // unit normal
  TcV nu, nv;     // derivatives of the unit normal
  SmBoolean bOk;
};

// a point of a tangential curve
struct TcPt
{
  double x[4];    // u1 v1 u2 v2
  TcV    P;       // 3d point (mean of the two surface points)
  TcV    n;       // unit normal of surface 1
  double dGap;    // |S1 - S2|
};

class TcTracer
{
public:
  const SmSurface * m_apSrf[2];
  SmExtent2d        m_aDom[2];
  double            m_dTol;       // 3d approximation tolerance
  double            m_dSize;      // model size
  double            m_dL;         // weight of the normal residual (a length)
  double            m_dSign;      // n2 orientation relative to n1 on the current curve

  TcTracer(const SmSurface * const apSrf[2], const SmExtent2d aDom[2], double dTol, double dSize)
    : m_dTol(dTol), m_dSize(dSize), m_dL(0.05*dSize), m_dSign(1.0)
  {
    m_apSrf[0] = apSrf[0]; m_apSrf[1] = apSrf[1];
    m_aDom[0]  = aDom[0];  m_aDom[1]  = aDom[1];
  }

  SmBoolean Eval(ULONG lSrf, double u, double v, TcSrfEval & e) const
  {
    SmVector3d aD[9];
    e.bOk = FALSE;
    if (m_apSrf[lSrf]->Evaluate(SmPoint2d(u,v), 2, 2, TRUE, TRUE, TRUE, aD) != SM_SUCCESS)
      return FALSE;
    e.P   = TcV(aD[0]);
    e.Sv  = TcV(aD[1]);
    e.Svv = TcV(aD[2]);
    e.Su  = TcV(aD[3]);
    e.Suv = TcV(aD[4]);
    e.Suu = TcV(aD[6]);
    TcV N  = e.Su.cross(e.Sv);
    double dLen = N.len();
    if (!(dLen > 1e-300)) return FALSE;
    e.n = N * (1.0/dLen);
    TcV Nu = e.Suu.cross(e.Sv) + e.Su.cross(e.Suv);
    TcV Nv = e.Suv.cross(e.Sv) + e.Su.cross(e.Svv);
    // d(N/|N|) = (I - n n^T) dN / |N|
    e.nu = (Nu - e.n * e.n.dot(Nu)) * (1.0/dLen);
    e.nv = (Nv - e.n * e.n.dot(Nv)) * (1.0/dLen);
    e.bOk = TRUE;
    return TRUE;
  }

  void Clamp(double x[4]) const
  {
    for (ULONG s=0; s<2; s++)
      {
        SmPoint2d p = m_aDom[s].ClampPoint2d(SmPoint2d(x[2*s], x[2*s+1]));
        x[2*s] = p.x; x[2*s+1] = p.y;
      }
  }

  // constraint kinds of Solve()
  enum { TC_FREE = 0, TC_PLANE = 1, TC_FIX = 2 };

  // residual and jacobian; returns number of rows
  int Residual(const double x[4], int eMode, const TcV & Q, const TcV & T, int iFix, double dFix,
               double r[9], double J[9][4], TcSrfEval e[2]) const
  {
    if (!Eval(0, x[0], x[1], e[0]) || !Eval(1, x[2], x[3], e[1])) return 0;
    int m = 0;
    TcV d = e[0].P - e[1].P;
    double ad[3] = {d.x, d.y, d.z};
    double a1u[3] = {e[0].Su.x, e[0].Su.y, e[0].Su.z};
    double a1v[3] = {e[0].Sv.x, e[0].Sv.y, e[0].Sv.z};
    double a2u[3] = {e[1].Su.x, e[1].Su.y, e[1].Su.z};
    double a2v[3] = {e[1].Sv.x, e[1].Sv.y, e[1].Sv.z};
    for (int k=0; k<3; k++, m++)
      { r[m] = ad[k]; J[m][0] = a1u[k]; J[m][1] = a1v[k]; J[m][2] = -a2u[k]; J[m][3] = -a2v[k]; }
    TcV dn = e[0].n - e[1].n * m_dSign;
    double an[3]  = {dn.x, dn.y, dn.z};
    double n1u[3] = {e[0].nu.x, e[0].nu.y, e[0].nu.z};
    double n1v[3] = {e[0].nv.x, e[0].nv.y, e[0].nv.z};
    double n2u[3] = {e[1].nu.x, e[1].nu.y, e[1].nu.z};
    double n2v[3] = {e[1].nv.x, e[1].nv.y, e[1].nv.z};
    for (int k=0; k<3; k++, m++)
      { r[m] = m_dL*an[k];
        J[m][0] = m_dL*n1u[k]; J[m][1] = m_dL*n1v[k];
        J[m][2] = -m_dSign*m_dL*n2u[k]; J[m][3] = -m_dSign*m_dL*n2v[k]; }
    if (eMode == TC_PLANE)
      {
        r[m] = (e[0].P - Q).dot(T); J[m][0] = e[0].Su.dot(T); J[m][1] = e[0].Sv.dot(T); J[m][2] = 0; J[m][3] = 0; m++;
        r[m] = (e[1].P - Q).dot(T); J[m][0] = 0; J[m][1] = 0; J[m][2] = e[1].Su.dot(T); J[m][3] = e[1].Sv.dot(T); m++;
      }
    else if (eMode == TC_FIX)
      {
        // weight the parameter residual to a length
        TcV dS = (iFix == 0) ? e[0].Su : (iFix == 1) ? e[0].Sv : (iFix == 2) ? e[1].Su : e[1].Sv;
        double w = smos_Max(dS.len(), 1e-12);
        r[m] = w*(x[iFix] - dFix);
        J[m][0] = J[m][1] = J[m][2] = J[m][3] = 0.0; J[m][iFix] = w; m++;
      }
    return m;
  }

  static SmBoolean Solve4(double A[4][4], double b[4], double z[4])
  {
    double M[4][5];
    for (int i=0;i<4;i++) { for (int j=0;j<4;j++) M[i][j] = A[i][j]; M[i][4] = b[i]; }
    for (int c=0;c<4;c++)
      {
        int p = c;
        for (int i=c+1;i<4;i++) if (std::fabs(M[i][c]) > std::fabs(M[p][c])) p = i;
        if (!(std::fabs(M[p][c]) > 1e-300)) return FALSE;
        if (p != c) for (int j=0;j<5;j++) { double t = M[c][j]; M[c][j] = M[p][j]; M[p][j] = t; }
        for (int i=c+1;i<4;i++)
          { double f = M[i][c]/M[c][c]; for (int j=c;j<5;j++) M[i][j] -= f*M[c][j]; }
      }
    for (int i=3;i>=0;i--)
      { double s = M[i][4]; for (int j=i+1;j<4;j++) s -= M[i][j]*z[j]; z[i] = s/M[i][i]; }
    return TRUE;
  }

  // Levenberg-Marquardt on the residual above. x is in/out.
  SmBoolean Solve(double x[4], int eMode, const TcV & Q, const TcV & T, int iFix, double dFix,
                  TcPt & rOut) const
  {
    Clamp(x);
    double r[9], J[9][4];
    TcSrfEval e[2];
    int m = Residual(x, eMode, Q, T, iFix, dFix, r, J, e);
    if (m == 0) return FALSE;
    double dCost = 0.0; for (int i=0;i<m;i++) dCost += r[i]*r[i];
    double dLambda = 1e-3;
    for (int it=0; it<200 && dCost > 1e-34; it++)
      {
        double A[4][4], g[4];
        for (int a=0;a<4;a++)
          { g[a] = 0.0; for (int i=0;i<m;i++) g[a] += J[i][a]*r[i];
            for (int b=0;b<4;b++) { A[a][b] = 0.0; for (int i=0;i<m;i++) A[a][b] += J[i][a]*J[i][b]; } }
        double dDiagMax = 0.0; for (int a=0;a<4;a++) dDiagMax = smos_Max(dDiagMax, A[a][a]);
        SmBoolean bStepped = FALSE;
        for (int tries=0; tries<12; tries++)
          {
            double B[4][4], nb[4], z[4];
            for (int a=0;a<4;a++)
              { for (int b=0;b<4;b++) B[a][b] = A[a][b];
                B[a][a] += dLambda*(A[a][a] + 1e-12*dDiagMax) + 1e-30; nb[a] = -g[a]; }
            if (!Solve4(B, nb, z)) { dLambda *= 10.0; continue; }
            double xt[4] = { x[0]+z[0], x[1]+z[1], x[2]+z[2], x[3]+z[3] };
            Clamp(xt);
            double rt[9], Jt[9][4];
            TcSrfEval et[2];
            int mt = Residual(xt, eMode, Q, T, iFix, dFix, rt, Jt, et);
            if (mt == 0) { dLambda *= 10.0; continue; }
            double dCostT = 0.0; for (int i=0;i<mt;i++) dCostT += rt[i]*rt[i];
            if (dCostT < dCost)
              {
                double dStep = std::fabs(xt[0]-x[0]) + std::fabs(xt[1]-x[1]) + std::fabs(xt[2]-x[2]) + std::fabs(xt[3]-x[3]);
                for (int a=0;a<4;a++) x[a] = xt[a];
                for (int i=0;i<mt;i++) { r[i] = rt[i]; for (int a=0;a<4;a++) J[i][a] = Jt[i][a]; }
                e[0] = et[0]; e[1] = et[1];
                dCost   = dCostT;
                dLambda = smos_Max(dLambda/4.0, 1e-12);
                bStepped = TRUE;
                if (dStep < 1e-15) it = 1000;
                break;
              }
            dLambda *= 8.0;
          }
        if (!bStepped) break;
      }

    // judge
    double dGap = (e[0].P - e[1].P).len();
    TcV    dn   = e[0].n - e[1].n * m_dSign;
    if (dGap > 0.1*m_dTol)   return FALSE;
    if (dn.len() > 1e-8)     return FALSE;
    if (eMode == TC_PLANE)
      {
        if (std::fabs((e[0].P - Q).dot(T)) > 0.1*m_dTol) return FALSE;
      }
    else if (eMode == TC_FIX)
      {
        if (std::fabs(x[iFix] - dFix) > 1e-9*(1.0 + std::fabs(dFix))) return FALSE;
      }
    for (int a=0;a<4;a++) rOut.x[a] = x[a];
    rOut.P    = (e[0].P + e[1].P) * 0.5;
    rOut.n    = e[0].n;
    rOut.dGap = dGap;
    return TRUE;
  }

  // drop a 3d tangent vector into both parameter spaces
  SmBoolean DropDir(const double x[4], const TcV & T, double duv[4]) const
  {
    for (ULONG s=0; s<2; s++)
      {
        TcSrfEval e;
        if (!Eval(s, x[2*s], x[2*s+1], e)) return FALSE;
        double a = e.Su.dot(e.Su), b = e.Su.dot(e.Sv), c = e.Sv.dot(e.Sv);
        double p = e.Su.dot(T), q = e.Sv.dot(T);
        double det = a*c - b*b;
        if (!(std::fabs(det) > 1e-300)) return FALSE;
        duv[2*s]   = ( c*p - b*q)/det;
        duv[2*s+1] = (-b*p + a*q)/det;
      }
    return TRUE;
  }

  // largest t in [0,1] with x + t*dx inside both domains; sets which coordinate/bound limits it
  double ExitFraction(const double x[4], const double dx[4], int & iFix, double & dFix) const
  {
    double t = 1.0;
    iFix = -1;
    for (ULONG s=0; s<2; s++)
      {
        double lo[2] = { m_aDom[s].GetMin().x, m_aDom[s].GetMin().y };
        double hi[2] = { m_aDom[s].GetMax().x, m_aDom[s].GetMax().y };
        for (int k=0; k<2; k++)
          {
            int i = (int)(2*s + k);
            if (x[i] + dx[i] > hi[k] && dx[i] > 0.0)
              { double f = (hi[k] - x[i])/dx[i]; if (f < t) { t = smos_Max(f,0.0); iFix = i; dFix = hi[k]; } }
            if (x[i] + dx[i] < lo[k] && dx[i] < 0.0)
              { double f = (lo[k] - x[i])/dx[i]; if (f < t) { t = smos_Max(f,0.0); iFix = i; dFix = lo[k]; } }
          }
      }
    return t;
  }

  // one predictor/corrector step attempt from rFrom along T with step h
  // returns 1 = point found, 2 = end point on a domain boundary found, 0 = failed
  int Step(const TcPt & rFrom, const TcV & T, double h, TcPt & rNext) const
  {
    double duv[4];
    if (!DropDir(rFrom.x, T, duv)) return 0;
    double dx[4] = { h*duv[0], h*duv[1], h*duv[2], h*duv[3] };
    int iFix = -1; double dFix = 0.0;
    double f = ExitFraction(rFrom.x, dx, iFix, dFix);
    if (iFix >= 0)
      {
        if (f*h < 1e-3*m_dTol) return 3; // already on the boundary heading out
        double xg[4] = { rFrom.x[0]+f*dx[0], rFrom.x[1]+f*dx[1], rFrom.x[2]+f*dx[2], rFrom.x[3]+f*dx[3] };
        xg[iFix] = dFix;
        TcPt sEnd;
        if (Solve(xg, TC_FIX, TcV(), TcV(), iFix, dFix, sEnd))
          {
            TcV dP = sEnd.P - rFrom.P;
            double dLen = dP.len();
            if (dLen < 1e-3*m_dTol) return 3;
            if (dP.dot(T) > 0.7*dLen && dLen < 2.0*f*h + 10.0*m_dTol)
              { rNext = sEnd; return 2; }
          }
        return 0;
      }
    double xg[4] = { rFrom.x[0]+dx[0], rFrom.x[1]+dx[1], rFrom.x[2]+dx[2], rFrom.x[3]+dx[3] };
    TcV Q = rFrom.P + T*h;
    if (!Solve(xg, TC_PLANE, Q, T, -1, 0.0, rNext)) return 0;
    TcV dP = rNext.P - rFrom.P;
    if ((rNext.P - Q).len() > 0.3*h) return 0;
    if (dP.dot(T) < 0.9*dP.len()) return 0;
    return 1;
  }

  // march from p0 in direction T0; appends points (not p0) to rPts
  // returns TRUE when the curve closed back onto rStart
  SmBoolean March(const TcPt & p0, const TcV & T0, const TcPt & rStart, std::vector<TcPt> & rPts) const
  {
    double h0   = 2e-3*m_dSize;
    double hMax = 0.05*m_dSize;
    double hMin = smos_Max(1e-7*m_dSize, 10.0*m_dTol);
    double h = h0;
    TcPt cur = p0;
    TcV  T   = T0;
    double dTravel = 0.0;
    for (int iter=0; iter<20000 && rPts.size() < 5000; iter++)
      {
        TcPt nxt;
        int eRes = Step(cur, T, h, nxt);
        if (eRes == 3) return FALSE;
        if (eRes == 0)
          {
            h *= 0.5;
            if (h < hMin) return FALSE;
            continue;
          }
        TcV dP = nxt.P - cur.P;
        double dLen = dP.len();
        if (dLen < 1e-3*m_dTol) return FALSE;
        TcV C = dP * (1.0/dLen);
        // the chord direction is close to the mean of the end tangents
        TcV Tn = C*2.0 - T;
        Tn = (Tn - nxt.n * nxt.n.dot(Tn)).unit();
        double dCos = Tn.dot(T);
        if (dCos < std::cos(SM_DEG2RAD(6.0)) && h > hMin)
          { h *= 0.5; continue; }
        // closing back onto the start point
        if (eRes == 1 && rPts.size() >= 3 && dTravel > 4.0*h)
          {
            TcV dS = rStart.P - cur.P;
            double dDist = dS.len();
            if (dDist < 1.5*h && dS.dot(T) > 0.5*dDist)
              {
                rPts.push_back(rStart);
                return TRUE;
              }
          }
        rPts.push_back(nxt);
        dTravel += dLen;
        if (eRes == 2) return FALSE;
        cur = nxt;
        T   = Tn;
        if (dCos > std::cos(SM_DEG2RAD(1.5))) h = smos_Min(h*1.5, hMax);
      }
    return FALSE;
  }

  // tangent plane directions in which the contact continues from p
  // returns number of directions; -1 = contact in (almost) all directions
  int Probe(const TcPt & p, std::vector<TcV> & rDirs) const
  {
    rDirs.clear();
    TcV n  = p.n;
    TcV a  = (std::fabs(n.x) < 0.6) ? TcV(1,0,0) : TcV(0,1,0);
    TcV e1 = (a - n*n.dot(a)).unit();
    TcV e2 = n.cross(e1);
    double r = 2e-3*m_dSize;
    const int N = 24;
    std::vector<TcV> sHits;
    int lNumExact = 0;
    for (int j=0; j<N; j++)
      {
        double th = 2.0*SM_PI*j/N;
        TcV d = e1*std::cos(th) + e2*std::sin(th);
        double duv[4];
        if (!DropDir(p.x, d, duv)) continue;
        double xg[4] = { p.x[0]+r*duv[0], p.x[1]+r*duv[1], p.x[2]+r*duv[2], p.x[3]+r*duv[3] };
        TcPt q;
        TcV Q = p.P + d*r;
        if (!Solve(xg, TC_PLANE, Q, d, -1, 0.0, q)) continue;
        double dOff = (q.P - Q).len();
        if (dOff > 0.35*r) continue;
        // contact right at the probe point: the surfaces coincide there
        // (on a tangential curve this happens only when probing along it)
        if (dOff < 0.01*r) lNumExact++;
        sHits.push_back((q.P - p.P).unit());
      }
    if (lNumExact >= 5 || (int)sHits.size() > N/2) return -1;
    // cluster
    for (size_t i=0; i<sHits.size(); i++)
      {
        SmBoolean bFound = FALSE;
        for (size_t k=0; k<rDirs.size(); k++)
          if (rDirs[k].dot(sHits[i]) > std::cos(SM_DEG2RAD(30.0))) { bFound = TRUE; break; }
        if (!bFound) rDirs.push_back(sHits[i]);
      }
    return (int)rDirs.size();
  }
};

// distance from a point to a polyline
double sm_TcPolyDist(const TcV & X, const std::vector<TcPt> & rPoly, size_t * pOptSeg = NULL, double * pOptT = NULL)
{
  double dBest = SM_BIG_DOUBLE;
  if (rPoly.size() == 1) { if (pOptSeg) *pOptSeg = 0; if (pOptT) *pOptT = 0.0; return (X - rPoly[0].P).len(); }
  for (size_t i=0; i+1<rPoly.size(); i++)
    {
      TcV a = rPoly[i].P, d = rPoly[i+1].P - a;
      double dd = d.dot(d);
      double t = dd > 0.0 ? (X - a).dot(d)/dd : 0.0;
      t = smos_Max(0.0, smos_Min(1.0, t));
      double dist = (X - (a + d*t)).len();
      if (dist < dBest) { dBest = dist; if (pOptSeg) *pOptSeg = i; if (pOptT) *pOptT = t; }
    }
  return dBest;
}

// unit tangents of a polyline by the derivative of the parabola through three neighbours
void sm_TcTangents(const std::vector<TcPt> & rC, SmBoolean bClosed, std::vector<TcV> & rTan)
{
  size_t n = rC.size();
  rTan.assign(n, TcV());
  std::vector<double> s(n, 0.0);
  for (size_t i=1; i<n; i++) s[i] = s[i-1] + (rC[i].P - rC[i-1].P).len();
  for (size_t i=0; i<n; i++)
    {
      TcV A, B, C;
      double tA, tB, tC, t;
      if (bClosed && n > 3 && (i == 0 || i == n-1))
        { A = rC[n-2].P; B = rC[0].P; C = rC[1].P; tA = -(s[n-1]-s[n-2]); tB = 0.0; tC = s[1]; t = 0.0; }
      else
        {
          size_t i0 = (i == 0) ? 0 : (i == n-1) ? n-3 : i-1;
          A = rC[i0].P; B = rC[i0+1].P; C = rC[i0+2].P;
          tA = s[i0]; tB = s[i0+1]; tC = s[i0+2]; t = s[i];
        }
      double w0 = ((t - tB) + (t - tC)) / ((tA - tB)*(tA - tC));
      double w1 = ((t - tA) + (t - tC)) / ((tB - tA)*(tB - tC));
      double w2 = ((t - tA) + (t - tB)) / ((tC - tA)*(tC - tB));
      rTan[i] = (A*w0 + B*w1 + C*w2).unit();
    }
}

} // end anonymous namespace

/*******************************************************************//**
PURPOSE: SolidCraft change: find and trace tangential intersection curves.

NOTES: Seeds are the touching points and the untraceable start points left
  by the regular marcher, the ends of traced curves where the normals are
  parallel, and the interior node pairs FindInteriorCurves() skipped because
  the normals there are (nearly) parallel.  Each seed is refined to a contact
  point; a contact point from which the contact continues in two opposite
  directions is traced both ways to a domain boundary, a dead end, or back to
  itself.  A new tangential curve replaces the curves and points the regular
  marcher produced along it; it is dropped when an existing curve already
  covers it.
***********************************************************************/
SmStatus SmSurfaceIntersector::FindTangentCurves
  (SmTArray<SmCurve*>         & r3DCurves,
   SmTArray<SmCurve*>         & rSurface1UVCurves,
   SmTArray<SmCurve*>         & rSurface2UVCurves,
   SmTArray<SmTsectCurveType> & rCurveTypes,
   SmTArray<double>           & rDeviations)
{
  if (m_apTrees[0] == NULL || m_apTrees[1] == NULL) return SM_SUCCESS;
  double dSize = smos_Min(m_apTrees[0]->GetTopNode()->m_sBBox.GetSize().Length(),
                          m_apTrees[1]->GetTopNode()->m_sBBox.GetSize().Length());
  if (!(dSize > 0.0)) return SM_SUCCESS;

  TcTracer sTr(m_cpSurface, m_vUVDomain, m_dThisApproxTol3d, dSize);

  // --- collect seeds
  std::vector<SmPoint2d> sSeeds;   // pairs: uv1, uv2
  ULONG lNumOrigCurves = r3DCurves.GetSize();
  for (ULONG ii=0; ii<lNumOrigCurves; ii++)
    {
      if (ii >= rSurface1UVCurves.GetSize() || ii >= rSurface2UVCurves.GetSize()) break;
      SmCurve *pUV1 = rSurface1UVCurves[ii], *pUV2 = rSurface2UVCurves[ii];
      if (pUV1 == NULL || pUV2 == NULL) continue;
      SmExtent1d sI1 = pUV1->GetNaturalInterval(), sI2 = pUV2->GetNaturalInterval();
      SmBoolean bDegen = r3DCurves[ii]->IsDegenerate();
      for (int iEnd=0; iEnd<(bDegen ? 1 : 2); iEnd++)
        {
          SmPoint3d sA, sB;
          if (pUV1->EvaluatePoint(iEnd ? sI1.GetMax() : sI1.GetMin(), sA) != SM_SUCCESS) continue;
          if (pUV2->EvaluatePoint(iEnd ? sI2.GetMax() : sI2.GetMin(), sB) != SM_SUCCESS) continue;
          SmPoint2d a(sA.x, sA.y), b(sB.x, sB.y);
          if (!bDegen)
            {
              // only curve ends where the surfaces are (nearly) tangent
              TcSrfEval e1, e2;
              if (!sTr.Eval(0, a.x, a.y, e1) || !sTr.Eval(1, b.x, b.y, e2)) continue;
              if (std::fabs(e1.n.dot(e2.n)) < std::cos(SM_DEG2RAD(1.0))) continue;
            }
          sSeeds.push_back(a); sSeeds.push_back(b);
        }
    }
  {
    SmTArray<SmTsectPnt*> sThrough;
    m_vThroughPoints.GetAllNodes(sThrough);
    for (ULONG ii=0; ii<sThrough.GetSize(); ii++)
      { sSeeds.push_back(sThrough[ii]->UVPos(0)); sSeeds.push_back(sThrough[ii]->UVPos(1)); }
  }
  for (ULONG ii=0; ii+1<m_vTangentSeeds.GetSize(); ii+=2)
    { sSeeds.push_back(m_vTangentSeeds[ii]); sSeeds.push_back(m_vTangentSeeds[ii+1]); }
  m_vTangentSeeds.ReSet();
  if (sSeeds.empty()) return SM_SUCCESS;

  // --- trace
  std::vector< std::vector<TcPt> > sNewCurves;
  std::vector<SmBoolean>           sNewClosed;
  double dNearTol = smos_Max(100.0*m_dThisApproxTol3d, 1e-4*dSize);
  for (size_t is=0; is+1<sSeeds.size(); is+=2)
    {
      TcPt p;
      double x[4] = { sSeeds[is].x, sSeeds[is].y, sSeeds[is+1].x, sSeeds[is+1].y };
      TcSrfEval e1, e2;
      if (!sTr.Eval(0, x[0], x[1], e1) || !sTr.Eval(1, x[2], x[3], e2)) continue;
      if ((e1.P - e2.P).len() > 0.01*dSize) continue;
      if (std::fabs(e1.n.dot(e2.n)) < std::cos(SM_DEG2RAD(10.0))) continue;
      sTr.m_dSign = (e1.n.dot(e2.n) < 0.0) ? -1.0 : 1.0;

      // skip seeds on tangential curves found already
      TcV X = (e1.P + e2.P)*0.5;
      SmBoolean bKnown = FALSE;
      for (size_t k=0; k<sNewCurves.size() && !bKnown; k++)
        if (sm_TcPolyDist(X, sNewCurves[k]) < dNearTol) bKnown = TRUE;
      if (bKnown) continue;

      if (!sTr.Solve(x, TcTracer::TC_FREE, TcV(), TcV(), -1, 0.0, p)) continue;
      for (size_t k=0; k<sNewCurves.size() && !bKnown; k++)
        if (sm_TcPolyDist(p.P, sNewCurves[k]) < dNearTol) bKnown = TRUE;
      if (bKnown) continue;

      std::vector<TcV> sDirs;
      int lNumDirs = sTr.Probe(p, sDirs);
      if (lNumDirs < 0) break;      // coincident surfaces: not ours to handle
      if (lNumDirs == 0) continue;  // isolated touching point

      if (lNumDirs <= 2)
        {
          // one curve through p (or ending at p on a boundary)
          std::vector<TcPt> sFwd, sBwd;
          SmBoolean bClosed = sTr.March(p, sDirs[0], p, sFwd);
          if (!bClosed)
            {
              TcV T1 = (lNumDirs == 2) ? sDirs[1] : sDirs[0]*(-1.0);
              sTr.March(p, T1, p, sBwd);
            }
          std::vector<TcPt> sAll;
          for (size_t k=sBwd.size(); k>0; k--) sAll.push_back(sBwd[k-1]);
          sAll.push_back(p);
          for (size_t k=0; k<sFwd.size(); k++) sAll.push_back(sFwd[k]);
          if (sAll.size() < 3) continue;
          sNewCurves.push_back(sAll);
          sNewClosed.push_back(bClosed);
        }
      else
        {
          // a node of several tangential branches: trace each branch from it
          for (int k=0; k<lNumDirs; k++)
            {
              std::vector<TcPt> sFwd;
              sTr.March(p, sDirs[(size_t)k], p, sFwd);
              if (sFwd.size() < 2) continue;
              std::vector<TcPt> sAll;
              sAll.push_back(p);
              for (size_t j=0; j<sFwd.size(); j++) sAll.push_back(sFwd[j]);
              sNewCurves.push_back(sAll);
              sNewClosed.push_back(FALSE);
            }
        }
    }

  // --- reconcile with what the regular marcher produced and flush
  for (size_t k=0; k<sNewCurves.size(); k++)
    {
      std::vector<TcPt> & rC = sNewCurves[k];
      double dLen = 0.0;
      for (size_t i=0; i+1<rC.size(); i++) dLen += (rC[i+1].P - rC[i].P).len();
      if (dLen < 10.0*m_dThisApproxTol3d) continue;

      // existing curve covering this one: keep the existing one
      SmBoolean bCovered = FALSE;
      for (ULONG ii=0; ii<r3DCurves.GetSize() && !bCovered; ii++)
        {
          SmCurve *pC = r3DCurves[ii];
          if (pC->IsDegenerate()) continue;
          bCovered = TRUE;
          for (size_t i=0; i<rC.size() && bCovered; i += std::max((size_t)1, rC.size()/40))
            {
              SmExtent1d sIvl = pC->GetNaturalInterval();
              double dParam, dDist; SmBoolean bOk;
              if (pC->DropPoint(sIvl, SmPoint3d(rC[i].P.x, rC[i].P.y, rC[i].P.z), NULL, m_dThisApproxTol3d*10.0,
                                NULL, bOk, dParam, dDist, SM_SO_MINIMIZE) != SM_SUCCESS || !bOk || dDist > dNearTol)
                bCovered = FALSE;
            }
        }
      if (bCovered) continue;

      // existing curves (and points) lying along this one: remove them
      for (ULONG ii=0; ii<r3DCurves.GetSize(); ii++)
        {
          SmCurve *pC = r3DCurves[ii];
          SmExtent1d sIvl = pC->GetNaturalInterval();
          SmBoolean bOn = TRUE;
          for (int j=0; j<=32 && bOn; j++)
            {
              SmPoint3d sP;
              if (pC->EvaluatePoint(sIvl.Evaluate(j/32.0), sP) != SM_SUCCESS) { bOn = FALSE; break; }
              if (sm_TcPolyDist(TcV(sP.x, sP.y, sP.z), rC) > 0.02*dSize) bOn = FALSE;
            }
          if (!bOn) continue;
          delete pC;
          r3DCurves.RemoveAt(ii);
          if (ii < rSurface1UVCurves.GetSize()) { delete rSurface1UVCurves[ii]; rSurface1UVCurves.RemoveAt(ii); }
          if (ii < rSurface2UVCurves.GetSize()) { delete rSurface2UVCurves[ii]; rSurface2UVCurves.RemoveAt(ii); }
          if (ii < rCurveTypes.GetSize()) rCurveTypes.RemoveAt(ii);
          if (ii < rDeviations.GetSize()) rDeviations.RemoveAt(ii);
          ii--;
        }

      // refine spans until the cubic Hermite spans FlushCurve() builds are within tolerance
      SmBoolean bClosed = sNewClosed[k];
      std::vector<TcV> sTan;
      for (int iPass=0; iPass<12; iPass++)
        {
          sm_TcTangents(rC, bClosed, sTan);
          std::vector<TcPt> sRef;
          SmBoolean bInserted = FALSE;
          for (size_t i=0; i+1<rC.size(); i++)
            {
              sRef.push_back(rC[i]);
              if (rC.size() + sRef.size() > 6000) continue;
              double h = (rC[i+1].P - rC[i].P).len();
              if (h < 10.0*m_dThisApproxTol3d) continue;
              TcV H = (rC[i].P + rC[i+1].P)*0.5 + (sTan[i] - sTan[i+1])*(h/8.0);
              TcV T = (rC[i+1].P - rC[i].P).unit();
              double xg[4];
              for (int a2=0; a2<4; a2++) xg[a2] = 0.5*(rC[i].x[a2] + rC[i+1].x[a2]);
              TcPt q;
              if (!sTr.Solve(xg, TcTracer::TC_PLANE, H, T, -1, 0.0, q)) continue;
              if ((q.P - H).len() <= 0.5*m_dThisApproxTol3d) continue;
              if ((q.P - H).len() > 0.25*h) continue;
              sRef.push_back(q);
              bInserted = TRUE;
            }
          sRef.push_back(rC.back());
          rC.swap(sRef);
          if (!bInserted) break;
        }
      sm_TcTangents(rC, bClosed, sTan);

      size_t n = rC.size();
      double dParam = 0.0;
      for (size_t i=0; i<n; i++)
        {
          if (i > 0) dParam += (rC[i].P - rC[i-1].P).len();
          const TcV & T = sTan[i];
          SmTsectPnt *pT = (SmTsectPnt*)m_vTSPntMgr.GetNewElement();
          pT->m_ePointType = SM_IP_TANGENT_CURVE;
          pT->UVPos(0) = SmPoint2d(rC[i].x[0], rC[i].x[1]);
          pT->UVPos(1) = SmPoint2d(rC[i].x[2], rC[i].x[3]);
          double duv[4] = {0,0,0,0};
          sTr.DropDir(rC[i].x, T, duv);
          pT->UVDeriv(0) = SmVector2d(duv[0], duv[1]);
          pT->UVDeriv(1) = SmVector2d(duv[2], duv[3]);
          pT->CrvPos()   = SmPoint3d(rC[i].P.x, rC[i].P.y, rC[i].P.z);
          pT->CrvDeriv() = SmVector3d(T.x, T.y, T.z);
          pT->m_dCurveParameter       = dParam;
          pT->m_dDeviation            = rC[i].dGap;
          pT->m_dTangentPlaneAngleRad = 0.0;
          m_vCurvePoints.Append(pT);
        }
      SER(FlushCurve(r3DCurves, rSurface1UVCurves, rSurface2UVCurves, rCurveTypes, rDeviations));
    }

  return SM_SUCCESS;

} // end SmSurfaceIntersector::FindTangentCurves
