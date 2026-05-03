//____________________________________________________________________________
/*!

\class    genie::LwlynSmithFFEM

\brief    Is a concrete implementation of the QELFormFactorsModelI:
          Form Factors for Quasi Elastic EM eN scattering according to
          Llewellyn-Smith model.

\author   Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  March 25, 2019

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _LLEWELLYN_SMITH_EM_FORM_FACTOR_MODEL_H_
#define _LLEWELLYN_SMITH_EM_FORM_FACTOR_MODEL_H_

#include "Physics/QuasiElastic/XSection/LwlynSmithFF.h"

namespace genie {

class LwlynSmithFFEM : public LwlynSmithFF {

public:
  LwlynSmithFFEM();
  LwlynSmithFFEM(string config);
  virtual ~LwlynSmithFFEM();

  // QELFormFactorModelI interface implementation
  virtual double F1V    (const Interaction * interaction) const /*override*/;
  virtual double xiF2V  (const Interaction * interaction) const /*override*/;
  //virtual double F3V  (const Interaction * interaction) const /*override*/;

  inline double FA(const Interaction* /*interaction*/) const /*override*/
    { return 0.; }
  inline double Fp(const Interaction* /*interaction*/) const /*override*/
    { return 0.; }
  //inline double F3A(const Interaction * interaction) const /*override*/
  //{ return 0.; }
};

}      // genie namespace

#endif

