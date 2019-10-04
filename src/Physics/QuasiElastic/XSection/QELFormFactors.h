//____________________________________________________________________________
/*!

\class    genie::QELFormFactors

\brief    A class holding Quasi Elastic (QEL) Form Factors.

          This class is using the \b Strategy Pattern. \n
          It can accept requests to calculate itself, for a given interaction,
          that it then delegates to the algorithmic object, implementing the
          QELFormFactorsModelI interface, that it finds attached to itself.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  April 20, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _QEL_FORM_FACTORS_H_
#define _QEL_FORM_FACTORS_H_

#include <iostream>

#include "Physics/QuasiElastic/XSection/QELFormFactorsModelI.h"
#include "Framework/Interaction/Interaction.h"

using std::ostream;

namespace genie {

class QELFormFactors;
ostream & operator << (ostream & stream, const QELFormFactors & ff);

class QELFormFactors {

public:

  QELFormFactors();
  QELFormFactors(const QELFormFactors & form_factors);
  virtual ~QELFormFactors() { }

  //! Attach an algorithm
  void   SetModel  (const QELFormFactorsModelI * model);

  //! Compute the form factors for the input interaction using the attached model
  void   Calculate (const Interaction * interaction);

  //! Get the computed form factor F1V
  double F1V    (void) const { return fF1V;   }

  //! Get the computed form factor xi*F2V
  double xiF2V  (void) const { return fxiF2V; }

  //! Get the computed form factor FA
  double FA     (void) const { return fFA;    }

  //! Get the computed form factor Fp
  double Fp     (void) const { return fFp;    }

  //! Get the attached model
  const QELFormFactorsModelI * Model (void) const {return fModel;}

  void   Reset    (Option_t * opt="");
  void   Copy     (const QELFormFactors & ff);
  bool   Compare  (const QELFormFactors & ff) const;
  void   Print    (ostream & stream) const;

  bool             operator == (const QELFormFactors & ff) const;
  QELFormFactors & operator =  (const QELFormFactors & ff);
  friend ostream & operator << (ostream & stream, const QELFormFactors & ff);

private:

  double fF1V;
  double fxiF2V;
  double fFA;
  double fFp;

  const QELFormFactorsModelI * fModel;
};

}        // genie namespace

#endif   // _QEL_FORM_FACTORS_H_
