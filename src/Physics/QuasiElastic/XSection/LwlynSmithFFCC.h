//____________________________________________________________________________
/*!

\class    genie::LwlynSmithFFCC

\brief    Is a concrete implementation of the QELFormFactorsModelI:
          Form Factors for Quasi Elastic CC vN scattering according to
          Llewellyn-Smith model.

\ref      H.Budd, A.Bodek, J.Arrington, NuINT02 proceedings

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 03, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _LLEWELLYN_SMITH_CC_FORM_FACTOR_MODEL_H_
#define _LLEWELLYN_SMITH_CC_FORM_FACTOR_MODEL_H_

#include "Physics/QuasiElastic/XSection/LwlynSmithFF.h"

namespace genie {

class LwlynSmithFFCC : public LwlynSmithFF {

public:
  LwlynSmithFFCC();
  LwlynSmithFFCC(string config);
  virtual ~LwlynSmithFFCC();

  // QELFormFactorModelI interface implementation
  double F1V    (const Interaction * interaction) const;
  double xiF2V  (const Interaction * interaction) const;
  double FA     (const Interaction * interaction) const;
  double Fp     (const Interaction * interaction) const;
};

}      // genie namespace

#endif

