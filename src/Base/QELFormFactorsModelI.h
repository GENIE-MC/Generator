//____________________________________________________________________________
/*!

\class    genie::QELFormFactorsModelI

\brief    Pure abstract base class. Defines the QELFormFactorsModelI interface
          to be implemented by any algorithmic class computing Quasi-Elastic
          Form Factors.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _QEL_FORM_FACTORS_MODEL_I_H_
#define _QEL_FORM_FACTORS_MODEL_I_H_

#include "Algorithm/Algorithm.h"
#include "Interaction/Interaction.h"

namespace genie {

class QELFormFactorsModelI : public Algorithm {

public:

  virtual ~QELFormFactorsModelI();

  virtual double F1V   (const Interaction * interaction) const = 0;
  virtual double xiF2V (const Interaction * interaction) const = 0;
  virtual double FA    (const Interaction * interaction) const = 0;
  virtual double Fp    (const Interaction * interaction) const = 0;

protected:

  QELFormFactorsModelI();
  QELFormFactorsModelI(string name);
  QELFormFactorsModelI(string name, string config);
};

}         // genie namespace 

#endif    // _QEL_FORM_FACTORS_MODEL_I_H_
