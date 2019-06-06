//____________________________________________________________________________
/*!

\class    genie::CacheBranchFx

\brief    A simple cache branch storing the cached data in a TNtuple

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  November 26, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
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

  void CreateSpline(void);
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

