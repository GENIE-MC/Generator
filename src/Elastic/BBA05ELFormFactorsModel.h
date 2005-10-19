//____________________________________________________________________________
/*!

\class    genie::BBA05ELFormFactorsModel

\brief    Concrete implementation of the ELFormFactorsModelI interface.
          Computes elastic form factors using the BBA2005 parameterization.

\ref

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

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
  BBA05ELFormFactorsModel(const char * param_set);
  virtual ~BBA05ELFormFactorsModel();

  double Ge (const Interaction * interaction) const;
  double Gm (const Interaction * interaction) const;
};

}         // genie namespace

#endif    // _BBA2005_EL_FORM_FACTORS_MODEL_H_
