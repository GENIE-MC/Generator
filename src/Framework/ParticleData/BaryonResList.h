//____________________________________________________________________________
/*!

\class    genie::BaryonResList

\brief    Encapsulates a list of baryon resonances.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 03, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _BARYON_RES_LIST_H_
#define _BARYON_RES_LIST_H_

#include <vector>
#include <iostream>
#include <string>

#include "Framework/ParticleData/BaryonResonance.h"

using std::vector;
using std::ostream;
using std::string;

namespace genie {

class BaryonResList;
ostream & operator << (ostream & stream, const BaryonResList & rl);

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
