//____________________________________________________________________________
/*!

\class    genie::AxialFormFactorModelI

\brief    Pure abstract base class. Defines the AxialFormFactorModelI interface
          to be implemented by LlewellynSmith Algorithm for calculating the 
          Axial Form Factor.

\author   Aaron Meyer <asmeyer2012 \at uchicago.edu>
          based off ELFormFactorsModelI by
          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  August 16, 2013

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _AXIAL_FORM_FACTOR_MODEL_I_H_
#define _AXIAL_FORM_FACTOR_MODEL_I_H_

#include "Framework/Algorithm/Algorithm.h"

namespace genie {

class Interaction;

class AxialFormFactorModelI : public Algorithm {

public:
  virtual ~AxialFormFactorModelI();

  //! Compute the axial form factor
  virtual double FA (const Interaction * interaction) const = 0;

protected:
  AxialFormFactorModelI();
  AxialFormFactorModelI(string name);
  AxialFormFactorModelI(string name, string config);
};

}         // genie namespace
#endif    // _AXIAL_FORM_FACTOR_MODEL_I_H_
