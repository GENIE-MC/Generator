//____________________________________________________________________________
/*!

\class    genie::LwlynSmithFFNC

\brief    Concrete implementation of the QELFormFactorsModelI :
          Form Factors for Quasi Elastic NC vN scattering according to
          Llewellyn-Smith model.

\ref      E.A.Paschos and J.Y.Yu, hep-ph/0107261

\author   Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  May 03, 2004

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _LLEWELLYN_SMITH_NC_FORM_FACTOR_MODEL_H_
#define _LLEWELLYN_SMITH_NC_FORM_FACTOR_MODEL_H_

#include "Physics/QuasiElastic/XSection/LwlynSmithFF.h"

namespace genie {

class LwlynSmithFFNC : public LwlynSmithFF {

public:
  LwlynSmithFFNC();
  LwlynSmithFFNC(string config);
  virtual ~LwlynSmithFFNC();

  // QELFormFactorModelI interface implementation
  double F1V     (const Interaction * interaction) const;
  double xiF2V   (const Interaction * interaction) const;
  double FA      (const Interaction * interaction) const;
  double Fp      (const Interaction * interaction) const;
};

}       // genie namespace

#endif
