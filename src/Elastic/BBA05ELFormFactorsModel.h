//____________________________________________________________________________
/*!

\class    genie::BBA05ELFormFactorsModel

\brief    Concrete implementation of the ELFormFactorsModelI interface.
          Computes elastic form factors using the BBA2005 parameterization.

\ref

\author

\created  October 19, 2005

*/
//____________________________________________________________________________

#ifndef _BBA2005_EL_FORM_FACTORS_MODEL_H_
#define _BBA2005_EL_FORM_FACTORS_MODEL_H_

#include "Base/ELFormFactorsModelI.h"

namespace genie {

class BBA05ELFormFactorsModel : public ELFormFactorsModelI {

public:
  BBA05ELFormFactorsModel();
  BBA05ELFormFactorsModel(string config);
  virtual ~BBA05ELFormFactorsModel();

  //! implement the ELFormFactorsModelI interface
  double Gep (const Interaction * interaction) const;
  double Gmp (const Interaction * interaction) const;
  double Gen (const Interaction * interaction) const;
  double Gmn (const Interaction * interaction) const;
};

}         // genie namespace

#endif    // _BBA2005_EL_FORM_FACTORS_MODEL_H_
