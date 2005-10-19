//____________________________________________________________________________
/*!

\class    genie::BBA03ELFormFactorsModel

\brief    Concrete implementation of the ELFormFactorsModelI interface.
          Computes elastic form factors using the BBA2003 parameterization.

\ref

\author

\created  October 19, 2005

*/
//____________________________________________________________________________

#ifndef _BBA2003_EL_FORM_FACTORS_MODEL_H_
#define _BBA2003_EL_FORM_FACTORS_MODEL_H_

#include "Base/ELFormFactorsModelI.h"

namespace genie {

class BBA03ELFormFactorsModel : public ELFormFactorsModelI {

public:
  BBA03ELFormFactorsModel();
  BBA03ELFormFactorsModel(const char * param_set);
  virtual ~BBA03ELFormFactorsModel();

  double Gep (double q2) const;
  double Gmp (double q2) const;
  double Gen (double q2) const;
  double Gmn (double q2) const;
};

}         // genie namespace

#endif    // _BBA2003_EL_FORM_FACTORS_MODEL_H_
