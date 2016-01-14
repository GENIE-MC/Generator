//____________________________________________________________________________
/*!

\class    genie::CacheBranchNtp

\brief    A simple cache branch storing the cached data in a TNtuple

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  November 26, 2004

\cpright  Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _CACHE_BRANCH_NTP_H_
#define _CACHE_BRANCH_NTP_H_

#include <iostream>
#include <string>
#include <TNtupleD.h>

#include "Utils/CacheBranchI.h"

using std::string;
using std::ostream;

namespace genie {

class CacheBranchNtp : public CacheBranchI
{
public:
  CacheBranchNtp();
  CacheBranchNtp(string name, string brdef);
  ~CacheBranchNtp();

  inline TNtupleD * Ntuple (void) const { return fNtp; }

  void CreateNtuple(string name, string branch_def);

  void Reset (void);
  void Print (ostream & stream) const;

  TNtupleD *       operator () (void) const;
  friend ostream & operator << (ostream & stream, const CacheBranchNtp & cbntp);

private:
  void Init    (void);
  void CleanUp (void);

  TNtupleD * fNtp;

ClassDef(CacheBranchNtp,1)
};

}      // genie namespace
#endif // _CACHE_BRANCH_NTP_H_

