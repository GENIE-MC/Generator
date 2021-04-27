//____________________________________________________________________________
/*!

\class    genie::MKFFCC

\brief    Is a concrete implementation of the QELFormFactorsModelI:
          Form Factors for MK SPP model.


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

#ifndef _MK_CC_FORM_FACTOR_MODEL_H_
#define _MK_CC_FORM_FACTOR_MODEL_H_

#include "Physics/QuasiElastic/XSection/LwlynSmithFF.h"

namespace genie {

class MKFFCC : public LwlynSmithFF {

public:
  MKFFCC();
  MKFFCC(string config);
  virtual ~MKFFCC();

  // QELFormFactorModelI interface implementation
  double F1V    (const Interaction * interaction) const;
  double xiF2V  (const Interaction * interaction) const;
  double FA     (const Interaction * interaction) const;
  double Fp     (const Interaction * interaction) const;
  double tau    (const Interaction * interaction) const;
  
  
};

}      // genie namespace

#endif

