//____________________________________________________________________________
/*!

\class    genie::FKR_MK

\brief    Simple struct-like class holding the extra Feynmann-Kislinger-Ravndall 
          (FKR) baryon excitation model parameters needed for MK model. 

\author  Igor Kakorin <kakorin@jinr.ru> Joint Institute for Nuclear Research
          

\created  Nov 12, 2019

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _FKR_MK_H_
#define _FKR_MK_H_

#include <iostream>
#include "Physics/Resonance/XSection/FKR.h"

using std::ostream;

namespace genie {

class FKR_MK;
ostream & operator<< (ostream & stream, const FKR_MK & parameters);

class FKR_MK : public FKR  {

public:

  friend ostream & operator<< (ostream & stream, const FKR_MK & parameters);

  double Splus;
  double Bplus;
  double Cplus;
  double Sminus;
  double Bminus;
  double Cminus;

  void Reset (void);
  void Print (ostream & stream) const;

  FKR_MK();
  ~FKR_MK();
};

}        // genie namespace

#endif   // _FKR_MK_H_
