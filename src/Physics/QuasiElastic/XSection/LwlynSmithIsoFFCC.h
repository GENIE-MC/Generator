//____________________________________________________________________________
/*!

\class    genie::LwlynSmithIsoFFCC

\brief    Is a concrete implementation of the QELFormFactorsModelI:
          Form Factors for Quasi Elastic CC vN scattering according to
          Llewellyn-Smith model for isoscalar nucleon.


\author   Igor Kakorin <kakorin@jinr.ru>, Joint Institute for Nuclear Research \n
          based on code of
          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  Nov 12, 2019

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _LLEWELLYN_SMITH_ISO_CC_FORM_FACTOR_MODEL_H_
#define _LLEWELLYN_SMITH_ISO_CC_FORM_FACTOR_MODEL_H_

#include "Physics/QuasiElastic/XSection/LwlynSmithFF.h"

namespace genie {

class LwlynSmithIsoFFCC : public LwlynSmithFF {

public:
  LwlynSmithIsoFFCC();
  LwlynSmithIsoFFCC(string config);
  virtual ~LwlynSmithIsoFFCC();

  // QELFormFactorModelI interface implementation
  double F1V    (const Interaction * interaction) const;
  double xiF2V  (const Interaction * interaction) const;
  double FA     (const Interaction * interaction) const;
  double Fp     (const Interaction * interaction) const;
  double tau    (const Interaction * interaction) const;
  
protected:
  virtual void LoadConfig (void);
  bool fIsCC;
  
  
};

}      // genie namespace

#endif

