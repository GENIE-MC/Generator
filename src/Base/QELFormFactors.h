//____________________________________________________________________________
/*!

\class    genie::QELFormFactors

\brief    A class holding Quasi Elastic (QEL) Form Factors.

          This class is using the \b Strategy Pattern. \n
          It can accept requests to calculate itself, for a given interaction,
          that it then delegates to the algorithmic object, implementing the
          QELFormFactorsModelI interface, that it finds attached to itself.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  April 20, 2004

*/
//____________________________________________________________________________

#ifndef _QEL_FORM_FACTORS_H_
#define _QEL_FORM_FACTORS_H_

#include <iostream>

#include "Base/QELFormFactorsModelI.h"
#include "Interaction/Interaction.h"

using std::ostream;

namespace genie {

class QELFormFactors {

public:

  QELFormFactors();
  QELFormFactors(const QELFormFactors & form_factors);
  virtual ~QELFormFactors() { }

  void   SetModel  (const QELFormFactorsModelI * model);
  void   Calculate (const Interaction * interaction);

  double F1V    (void) const { return fF1V;   }
  double xiF2V  (void) const { return fxiF2V; }
  double FA     (void) const { return fFA;    }
  double Fp     (void) const { return fFp;    }

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
