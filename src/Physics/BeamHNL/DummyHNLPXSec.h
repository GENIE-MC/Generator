//____________________________________________________________________________
/*!

\class    genie::DummyHNLPXSec

\brief

\author   

\created  May 05, 2009

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org        
*/
//____________________________________________________________________________

#ifndef _DUMMYHNL_PXSEC_H_
#define _DUMMYHNL_PXSEC_H_

#include "Framework/EventGen/XSecAlgorithmI.h"

namespace genie {

class DummyHNLPXSec : public XSecAlgorithmI {

public:
  DummyHNLPXSec();
  DummyHNLPXSec(string config);
 ~DummyHNLPXSec();

  // XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k) const;
  double Integral        (const Interaction * i) const;
  bool   ValidProcess    (const Interaction * i) const;
};

}       // genie namespace
#endif  //
