//____________________________________________________________________________
/*!

\class    genie::ELFormFactorsModelI

\brief    Pure abstract base class. Defines the ELFormFactorsModelI interface
          to be implemented by any algorithmic class computing Elastic Form
          Factors.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004

*/
//____________________________________________________________________________

#ifndef _EL_FORM_FACTORS_MODEL_I_H_
#define _EL_FORM_FACTORS_MODEL_I_H_

#include "Algorithm/Algorithm.h"

namespace genie {

class ELFormFactorsModelI : public Algorithm {

public:

  virtual ~ELFormFactorsModelI();

  virtual double Gep (double q2) const = 0;
  virtual double Gmp (double q2) const = 0;
  virtual double Gen (double q2) const = 0;
  virtual double Gmn (double q2) const = 0;

protected:

  ELFormFactorsModelI();
  ELFormFactorsModelI(const char * param_set);
};

}         // genie namespace

#endif    // _EL_FORM_FACTORS_MODEL_I_H_
