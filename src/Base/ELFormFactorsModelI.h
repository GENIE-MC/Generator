//____________________________________________________________________________
/*!

\class    genie::ELFormFactorsModelI

\brief    Pure abstract base class. Defines the ELFormFactorsModelI interface
          to be implemented by any algorithmic class computing Elastic Form
          Factors.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _EL_FORM_FACTORS_MODEL_I_H_
#define _EL_FORM_FACTORS_MODEL_I_H_

#include "Algorithm/Algorithm.h"

namespace genie {

class Interaction;

class ELFormFactorsModelI : public Algorithm {

public:
  virtual ~ELFormFactorsModelI();

  //! Compute the elastic form factor G_{ep} for the input interaction
  virtual double Gep (const Interaction * interaction) const = 0;

  //! Compute the elastic form factor G_{mp} for the input interaction
  virtual double Gmp (const Interaction * interaction) const = 0;

  //! Compute the elastic form factor G_{en} for the input interaction
  virtual double Gen (const Interaction * interaction) const = 0;

  //! Compute the elastic form factor G_{mn} for the input interaction
  virtual double Gmn (const Interaction * interaction) const = 0;

protected:
  ELFormFactorsModelI();
  ELFormFactorsModelI(string name);
  ELFormFactorsModelI(string name, string config);
};

}         // genie namespace
#endif    // _EL_FORM_FACTORS_MODEL_I_H_
