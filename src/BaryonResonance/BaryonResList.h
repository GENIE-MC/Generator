//____________________________________________________________________________
/*!

\class    genie::BaryonResList

\brief    Encapsulates a list of baryon resonances.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  May 03, 2004

\cpright  Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _BARYON_RES_LIST_H_
#define _BARYON_RES_LIST_H_

#include <vector>
#include <iostream>

#include "BaryonResonance/BaryonResDataSetI.h"
#include "BaryonResonance/BaryonResonance.h"

using std::vector;
using std::ostream;

namespace genie {

class BaryonResList
{
public:

  BaryonResList();
  BaryonResList(const BaryonResList & rl);
  virtual ~BaryonResList();

  void DecodeFromNameList (string list, string delimiter = ",");

  unsigned int NResonances        (void)              const;
  string       ResonanceName      (unsigned int ires) const;
  Resonance_t  ResonanceId        (unsigned int ires) const;
  int          ResonancePdgCode   (unsigned int ires) const;
  bool         Find               (Resonance_t res)   const;

  void Clear (void);
  void Copy  (const BaryonResList & rl);
  void Print (ostream & stream) const;

  friend ostream & operator << (ostream & stream, const BaryonResList & rl);

private:

  vector<Resonance_t> * fResVec;
};

}       // genie namepace

#endif  // _BARYON_RES_LIST_H_
