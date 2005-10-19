//____________________________________________________________________________
/*!

\class    genie::ELFormFactors

\brief    A class holding the Elastic Form Factors Ge,Gm

          This class is using the \b Strategy Pattern. \n

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  April 20, 2004

*/
//____________________________________________________________________________

#ifndef _EL_FORM_FACTORS_H_
#define _EL_FORM_FACTORS_H_

#include <iostream>

#include "Base/ELFormFactorsModelI.h"
#include "Interaction/Interaction.h"

using std::ostream;

namespace genie {

class ELFormFactors {

public:

  ELFormFactors();
  ELFormFactors(const ELFormFactors & form_factors);
  virtual ~ELFormFactors() { }

  void   SetModel  (const ELFormFactorsModelI * model);
  void   Calculate (const Interaction * interaction);

  double Ge  (void) const { return fGe; }
  double Gm  (void) const { return fGm; }

  friend ostream & operator << (ostream & stream, const ELFormFactors & ff);

  void Print(ostream & stream) const;

private:

  void   InitFormFactors(void);

  double fGe;
  double fGm;

  const ELFormFactorsModelI * fModel;
};

}        // genie namespace

#endif   // _QEL_FORM_FACTORS_H_
