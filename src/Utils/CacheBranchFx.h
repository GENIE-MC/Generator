//____________________________________________________________________________
/*!

\class    genie::CacheBranchFx

\brief    A simple cache branch storing the cached data in a TNtuple

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  November 26, 2004

\cpright  Copyright (c) 2003-2008, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _CACHE_BRANCH_FUNC_X_H_
#define _CACHE_BRANCH_FUNC_X_H_

#include <iostream>
#include <string>
#include <TNtupleD.h>

#include "Numerical/Spline.h"
#include "Utils/CacheBranchI.h"

using std::string;
using std::ostream;

namespace genie {

class CacheBranchFx : public CacheBranchI
{
public:
  CacheBranchFx();
  CacheBranchFx(string name);
  ~CacheBranchFx();

  inline TNtupleD * Ntuple (void) const { return fNtp;    }
  inline Spline *   Spl    (void) const { return fSpline; }

  void CreateNtuple(string name);
  void CreateSpline(void);
  void AddValues(double x, double y);

  void Reset (void);
  void Print (ostream & stream) const;

  double operator () (double x) const;
  friend ostream & operator << (ostream & stream, const CacheBranchFx & cbntp);

private:
  void Init    (void);
  void CleanUp (void);

  TNtupleD * fNtp;
  Spline   * fSpline;

ClassDef(CacheBranchFx,1)
};

}      // genie namespace
#endif // _CACHE_BRANCH_FUNC_X_H_

