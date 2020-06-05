//____________________________________________________________________________
/*!

\class    genie::DummyPXSec

\brief

\author

\created  May 05, 2009

\cpright  Copyright (c) 2003-2020, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org        
*/
//____________________________________________________________________________

#ifndef _DUMMY_PXSEC_H_
#define _DUMMY_PXSEC_H_

#include "Framework/EventGen/XSecAlgorithmI.h"

namespace genie {

class DummyPXSec : public XSecAlgorithmI {

public:
  DummyPXSec();
  DummyPXSec(string config);
 ~DummyPXSec();

  // XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k) const;
  double Integral        (const Interaction * i) const;
  bool   ValidProcess    (const Interaction * i) const;
};

}       // genie namespace
#endif  //
