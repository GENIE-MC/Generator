//____________________________________________________________________________
/*!

\class    genie::NNBarOscDummyPXSec

\brief    

\author   Jeremy Hewes, Georgia Karagiorgi
          University of Manchester

\created  November, 2016

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _NNBAR_OSC_DUMMY_PXSEC_H_
#define _NNBAR_OSC_DUMMY_PXSEC_H_

#include "Framework/EventGen/XSecAlgorithmI.h"

namespace genie {

class NNBarOscDummyPXSec : public XSecAlgorithmI {

public:
  NNBarOscDummyPXSec();
  NNBarOscDummyPXSec(string config);
 ~NNBarOscDummyPXSec();

  // XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k) const;
  double Integral        (const Interaction * i) const;
  bool   ValidProcess    (const Interaction * i) const;
};

}       // genie namespace
#endif  //
