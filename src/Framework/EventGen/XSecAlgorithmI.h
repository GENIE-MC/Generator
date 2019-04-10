//____________________________________________________________________________
/*!

\class    genie::XSecAlgorithmI

\brief    Cross Section Calculation Interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 03, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _XSEC_ALGORITHM_I_H_
#define _XSEC_ALGORITHM_I_H_

#include "Framework/Algorithm/Algorithm.h"
#include "Framework/Conventions/KinePhaseSpace.h"
#include "Framework/Interaction/Interaction.h"

namespace genie {

class XSecAlgorithmI : public Algorithm {

public:
  virtual ~XSecAlgorithmI();

  //! Compute the cross section for the input interaction
  virtual double XSec (const Interaction* i, KinePhaseSpace_t k=kPSfE) const = 0;

  //! Integrate the model over the kinematic phase space available to the
  //! input interaction (kinematical cuts can be included)
  virtual double Integral (const Interaction* i) const = 0;

  //! Can this cross section algorithm handle the input process?
  virtual bool ValidProcess    (const Interaction* i) const = 0;

  //! Is the input kinematical point a physically allowed one?
  virtual bool ValidKinematics (const Interaction* i) const;

protected:
  XSecAlgorithmI();
  XSecAlgorithmI(string name);
  XSecAlgorithmI(string name, string config);
};

}       // genie namespace
#endif  // _XSEC_ALGORITHM_I_H_
