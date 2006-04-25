//____________________________________________________________________________
/*!

\class    genie::ELFormFactors

\brief    A class holding the Elastic Form Factors Ge,Gm

          This class is using the \b Strategy Pattern. \n

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  April 20, 2004

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _EL_FORM_FACTORS_H_
#define _EL_FORM_FACTORS_H_

#include <iostream>

#include "Base/ELFormFactorsModelI.h"

using std::ostream;

namespace genie {

class Interaction;

class ELFormFactors {

public:

  ELFormFactors();
  ELFormFactors(const ELFormFactors & form_factors);
  virtual ~ELFormFactors() { }

  void   SetModel  (const ELFormFactorsModelI * model);
  void   Calculate (const Interaction * interaction);

  double Gep (void) const { return fGep; }
  double Gmp (void) const { return fGmp; }
  double Gen (void) const { return fGen; }
  double Gmn (void) const { return fGmn; }

  const ELFormFactorsModelI * Model (void) const {return fModel;}

  void   Reset    (Option_t * opt="");
  void   Copy     (const ELFormFactors & ff);
  bool   Compare  (const ELFormFactors & ff) const;
  void   Print    (ostream & stream) const;

  bool             operator == (const ELFormFactors & ff) const;
  ELFormFactors &  operator =  (const ELFormFactors & ff);
  friend ostream & operator << (ostream & stream, const ELFormFactors & ff);

private:

  double fGep;
  double fGmp;
  double fGen;
  double fGmn;

  const ELFormFactorsModelI * fModel;
};

}        // genie namespace

#endif   // _QEL_FORM_FACTORS_H_
