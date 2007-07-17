//____________________________________________________________________________
/*!

\class    genie::LlewellynSmithModelCC

\brief    Is a concrete implementation of the QELFormFactorsModelI:
          Form Factors for Quasi Elastic CC vN scattering according to
          Llewellyn-Smith model.

\ref      H.Budd, A.Bodek, J.Arrington, NuINT02 proceedings

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  May 03, 2004

\cpright  Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _LLEWELLYN_SMITH_MODEL_CC_H_
#define _LLEWELLYN_SMITH_MODEL_CC_H_

#include "LlewellynSmith/LlewellynSmithModel.h"

namespace genie {

class LlewellynSmithModelCC : public LlewellynSmithModel {

public:

  LlewellynSmithModelCC();
  LlewellynSmithModelCC(string config);
  virtual ~LlewellynSmithModelCC();

  //-- QELFormFactorModelI interface implementation

  double F1V    (const Interaction * interaction) const;
  double xiF2V  (const Interaction * interaction) const;
  double FA     (const Interaction * interaction) const;
  double Fp     (const Interaction * interaction) const;
};

}      // genie namespace

#endif // _LLEWELLYN_SMITH_MODEL_CC_H_

