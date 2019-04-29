//____________________________________________________________________________
/*!

\class    genie::AxialFormFactor

\brief    A class holding the Axial Form Factor

          This class is using the \b Strategy Pattern. \n

\author   Aaron Meyer <asmeyer2012 \at uchicago.edu>
          based off AxialFormFactorModelI by
          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  August 19, 2013

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _AXIAL_FORM_FACTOR_H_
#define _AXIAL_FORM_FACTOR_H_

#include <iostream>

#include "Physics/QuasiElastic/XSection/AxialFormFactorModelI.h"

using std::ostream;

namespace genie {

class Interaction;
class AxialFormFactor;
ostream & operator << (ostream & stream, const AxialFormFactor & ff);

class AxialFormFactor {

public:
  AxialFormFactor();
  AxialFormFactor(const AxialFormFactor & form_factors);
  virtual ~AxialFormFactor() { }

  //! Attach an algorithm
  void   SetModel  (const AxialFormFactorModelI * model);

  //! Calculate the form factors for the input interaction using the attached algorithm
  void   Calculate (const Interaction * interaction);

  //! Get the computed axial form factor
  double FA (void) const { return fFA; }

  //! Get the attached model
  const AxialFormFactorModelI * Model (void) const {return fModel;}

  void   Reset    (Option_t * opt="");
  void   Copy     (const AxialFormFactor & ff);
  bool   Compare  (const AxialFormFactor & ff) const;
  void   Print    (ostream & stream) const;

  bool              operator == (const AxialFormFactor & ff) const;
  AxialFormFactor & operator =  (const AxialFormFactor & ff);
  friend ostream  & operator << (ostream & stream, const AxialFormFactor & ff);

private:

  double fFA;

  const AxialFormFactorModelI * fModel;
};

}        // genie namespace

#endif   // _AXIAL_FORM_FACTOR_H_
