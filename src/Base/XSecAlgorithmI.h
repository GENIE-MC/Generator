//____________________________________________________________________________
/*!

\class    genie::XSecAlgorithmI

\brief    Cross Section Calculation Interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _XSEC_ALGORITHM_I_H_
#define _XSEC_ALGORITHM_I_H_

#include "Algorithm/Algorithm.h"
#include "Conventions/KinePhaseSpace.h"
#include "Interaction/Interaction.h"

namespace genie {

class XSecAlgorithmI : public Algorithm {

public:
  virtual ~XSecAlgorithmI();

  //! the XSecAlgorithmI algorithm

  virtual double XSec (const Interaction* i, KinePhaseSpace_t k=kPSfE) const = 0;

  virtual bool ValidProcess    (const Interaction* i) const = 0;
  virtual bool ValidKinematics (const Interaction* i) const = 0;

protected:
  XSecAlgorithmI();
  XSecAlgorithmI(string name);
  XSecAlgorithmI(string name, string config);
};

}       // genie namespace

#endif  // _XSEC_ALGORITHM_I_H_
