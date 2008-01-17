//____________________________________________________________________________
/*!

\class    genie::LlewellynSmithModelNC

\brief    Concrete implementation of the QELFormFactorsModelI :
          Form Factors for Quasi Elastic NC vN scattering according to
          Llewellyn-Smith model.

\ref      E.A.Paschos and J.Y.Yu, hep-ph/0107261

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  May 03, 2004

\cpright  Copyright (c) 2003-2008, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _LLEWELLYN_SMITH_MODEL_NC_H_
#define _LLEWELLYN_SMITH_MODEL_NC_H_

#include "LlewellynSmith/LlewellynSmithModel.h"

namespace genie {

class LlewellynSmithModelNC : public LlewellynSmithModel {

public:

  LlewellynSmithModelNC();
  LlewellynSmithModelNC(string config);
  virtual ~LlewellynSmithModelNC();

  //-- QELFormFactorModelI interface implementation

  double F1V     (const Interaction * interaction) const;
  double xiF2V   (const Interaction * interaction) const;
  double FA      (const Interaction * interaction) const;
  double Fp      (const Interaction * interaction) const;
};

}       // genie namespace

#endif  // _LLEWELLYN_SMITH_MODEL_NC_H_

