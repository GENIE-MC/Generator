//____________________________________________________________________________
/*!

\class    genie::LwlynSmithFFDeltaS.h

\brief    Is a concrete implementation of the QELFormFactorsModelI:
          Form Factors for Quasi Elastic CC vN Delta S=1 scattering.

\ref      Equations for the strange form factors are taken from Cabibbo 
          and Chilton, Phys.Rev. 137 (1965) B1628-B1634 and 
          Alam et al., J.Phys. G42 (2015) no.5, 055107.

\author   Hugh Gallagher <Hugh.Gallagher \at tufts.edu>
          Tufts University 

\created  April 10, 2016

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _LLEWELLYN_SMITH_DELTAS_FORM_FACTOR_MODEL_H_
#define _LLEWELLYN_SMITH_DELTAS_FORM_FACTOR_MODEL_H_

#include "Physics/QuasiElastic/XSection/LwlynSmithFF.h"

namespace genie {

class LwlynSmithFFDeltaS : public LwlynSmithFF {

public:
  LwlynSmithFFDeltaS();
  LwlynSmithFFDeltaS(string config);
  virtual ~LwlynSmithFFDeltaS();

  // QELFormFactorModelI interface implementation
  double F1V    (const Interaction * interaction) const;
  double xiF2V  (const Interaction * interaction) const;
  double FA     (const Interaction * interaction) const;
  double Fp     (const Interaction * interaction) const;
};

}      // genie namespace

#endif

