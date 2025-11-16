//____________________________________________________________________________
/*!

\class    genie::DipoleELFormFactorsModel

\brief    Concrete implementation of the ELFormFactorsModelI interface.
          Computes dipole elastic form factors.

\author   Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  Oct 19, 2005

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org        
*/
//____________________________________________________________________________

#ifndef _DIPOLE_EL_FORM_FACTORS_MODEL_H_
#define _DIPOLE_EL_FORM_FACTORS_MODEL_H_

#include "Physics/QuasiElastic/XSection/ELFormFactorsModelI.h"

namespace genie {

class DipoleELFormFactorsModel : public ELFormFactorsModelI {

public:
  DipoleELFormFactorsModel();
  DipoleELFormFactorsModel(string config);
  virtual ~DipoleELFormFactorsModel();

  // implement the ELFormFactorsModelI interface
  double Gep (const Interaction * interaction) const;
  double Gmp (const Interaction * interaction) const;
  double Gen (const Interaction * interaction) const;
  double Gmn (const Interaction * interaction) const;

  // overload Algorithm's Configure()
  void   Configure  (const Registry & config);
  void   Configure  (string param_set);

private:

  void LoadConfig(void);

  double fMv;
  double fMv2;
  double fMuP;
  double fMuN;
};

}         // genie namespace

#endif    // _DIPOLE_EL_FORM_FACTORS_MODEL_H_
