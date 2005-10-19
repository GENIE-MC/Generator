//____________________________________________________________________________
/*!

\class    genie::DipoleELFormFactorsModel

\brief    Concrete implementation of the ELFormFactorsModelI interface.
          Computes dipole elastic form factors.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  October 19, 2005

*/
//____________________________________________________________________________

#ifndef _DIPOLE_EL_FORM_FACTORS_MODEL_H_
#define _DIPOLE_EL_FORM_FACTORS_MODEL_H_

#include "Base/ELFormFactorsModelI.h"

namespace genie {

class DipoleELFormFactorsModel : public ELFormFactorsModelI {

public:
  DipoleELFormFactorsModel();
  DipoleELFormFactorsModel(const char * param_set);
  virtual ~DipoleELFormFactorsModel();

  double Gep (double q2) const;
  double Gmp (double q2) const;
  double Gen (double q2) const;
  double Gmn (double q2) const;
};

}         // genie namespace

#endif    // _DIPOLE_EL_FORM_FACTORS_MODEL_H_
