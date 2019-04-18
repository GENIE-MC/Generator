//____________________________________________________________________________
/*!

\class    genie::QELFormFactorsModelI

\brief    Pure abstract base class. Defines the QELFormFactorsModelI interface
          to be implemented by any algorithmic class computing Quasi-Elastic
          Form Factors.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 03, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _QEL_FORM_FACTORS_MODEL_I_H_
#define _QEL_FORM_FACTORS_MODEL_I_H_

#include "Framework/Algorithm/Algorithm.h"
#include "Framework/Interaction/Interaction.h"

namespace genie {

class QELFormFactorsModelI : public Algorithm {

public:
  virtual ~QELFormFactorsModelI();

  //! Compute the form factor F1V for the input interaction
  virtual double F1V   (const Interaction * interaction) const = 0;

  //! Compute the form factor xi*F2V for the input interaction
  virtual double xiF2V (const Interaction * interaction) const = 0;

  //! Compute the form factor FA for the input interaction
  virtual double FA    (const Interaction * interaction) const = 0;

  //! Compute the form factor Fp for the input interaction
  virtual double Fp    (const Interaction * interaction) const = 0;

protected:
  QELFormFactorsModelI();
  QELFormFactorsModelI(string name);
  QELFormFactorsModelI(string name, string config);
};

}         // genie namespace 
#endif    // _QEL_FORM_FACTORS_MODEL_I_H_
