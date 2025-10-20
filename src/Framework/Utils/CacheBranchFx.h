//____________________________________________________________________________
/*!

\class    genie::CacheBranchFx

\brief    A simple cache branch storing the cached data in a TNtuple

\author   Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool
          
          Update May 15, 2022 IK: 
          Now type of spline can be:  TSpline3, TSpline5 and 
          ROOT::Math::GSLInterpolator (LINEAR, POLYNOMIAL, CSPLINE, CSPLINE_PERIODIC,
          AKIMA, AKIMA_PERIODIC)
          
\ref      [1] GENIE docdb 297

\author   Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool \n
          Igor Kakorin <kakorin@jinr.ru>
          Joint Institute for Nuclear Research

\created  November 26, 2004

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _CACHE_BRANCH_FUNC_X_H_
#define _CACHE_BRANCH_FUNC_X_H_

#include <iostream>
#include <string>
#include <map>

#include "Framework/Numerical/Spline.h"
#include "Framework/Utils/CacheBranchI.h"

using std::string;
using std::ostream;
using std::map;

namespace genie {

class CacheBranchFx;
ostream & operator << (ostream & stream, const CacheBranchFx & cbntp);

class CacheBranchFx : public CacheBranchI
{
public:
  using TObject::Print; // suppress clang 'hides overloaded virtual function [-Woverloaded-virtual]' warnings

  CacheBranchFx();
  CacheBranchFx(string name);
  ~CacheBranchFx();

  const map<double,double> & Map (void) const { return fFx;     }
  Spline *                   Spl (void) const { return fSpline; }

  void CreateSpline(string type = "TSpline3");
  void AddValues(double x, double y);

  void Reset (void);
  void Print (ostream & stream) const;

  double operator () (double x) const;
  friend ostream & operator << (ostream & stream, const CacheBranchFx & cbntp);

private:
  void Init    (void);
  void CleanUp (void);

  string             fName;   ///< cache branch name
  map<double,double> fFx;     ///< x->y map
  Spline *           fSpline; ///< spline y = f(x)

ClassDef(CacheBranchFx,1)
};

}      // genie namespace
#endif // _CACHE_BRANCH_FUNC_X_H_
