//____________________________________________________________________________
/*!

\class    genie::QELFormFactorsModelI

\brief    Pure abstract base class. Defines the QELFormFactorsModelI interface
          to be implemented by any algorithmic class computing Quasi-Elastic
          Form Factors.

\author   Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
          University of Liverpool & STFC Rutherford Appleton Laboratory

\created  May 03, 2004

\cpright  Copyright (c) 2003-2020, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
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

  //! \name Second-class currents
  //! @{
  //! Compute the form factor F3V for the input interaction
  inline virtual double F3V(const Interaction* /*interaction*/)
    const { return 0.; }

  //! Compute the form factor F3A for the input interaction
  inline virtual double F3A(const Interaction* /*interaction*/)
    const { return 0.; }
  //! @}

protected:
  QELFormFactorsModelI();
  QELFormFactorsModelI(string name);
  QELFormFactorsModelI(string name, string config);
};

}         // genie namespace
#endif    // _QEL_FORM_FACTORS_MODEL_I_H_
