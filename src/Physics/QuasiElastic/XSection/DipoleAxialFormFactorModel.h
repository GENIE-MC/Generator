//____________________________________________________________________________
/*!

\class    genie::DipoleAxialFormFactorModel

\brief    Concrete implementation of the AxialFormFactorModelI interface.
          Computes the axial form factor using the dipole form factor
          approximation.

\author   Aaron Meyer <asmeyer2012 \at uchicago.edu>
          based off DipoleELFormFactorsModel by
          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  August 16, 2013

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _DIPOLE_AXIAL_FORM_FACTOR_MODEL_H_
#define _DIPOLE_AXIAL_FORM_FACTOR_MODEL_H_

#include "Physics/QuasiElastic/XSection/AxialFormFactorModelI.h"

namespace genie {

class DipoleAxialFormFactorModel : public AxialFormFactorModelI {

public:
  DipoleAxialFormFactorModel();
  DipoleAxialFormFactorModel(string config);
  virtual ~DipoleAxialFormFactorModel();

  // implement the AxialFormFactorModelI interface
  double FA (const Interaction * interaction) const;

  // overload Algorithm's Configure() 
  void   Configure  (const Registry & config);
  void   Configure  (string param_set);

private:

  void LoadConfig(void);

  double fMa;  ///< axial mass
  double fMa2;
  double fFA0; ///< FA(q2=0)
};

}         // genie namespace

#endif    // _DIPOLE_AXIAL_FORM_FACTOR_MODEL_H_
