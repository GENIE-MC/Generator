//____________________________________________________________________________
/*!

\class    genie::ZExpAxialFormFactorModel

\brief    Concrete implementation of the AxialFormFactorModelI interface.
          Computes the axial form factor using the model-independent
          z-expansion technique

\ref      Hill et al.
          arXiv:1008.4619
          DOI: 10.1103/PhysRevD.82.113005

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

#ifndef _Z_EXPANSION_AXIAL_FORM_FACTOR_MODEL_H_
#define _Z_EXPANSION_AXIAL_FORM_FACTOR_MODEL_H_

#include "Physics/QuasiElastic/XSection/AxialFormFactorModelI.h"

namespace genie {

class ZExpAxialFormFactorModel : public AxialFormFactorModelI {

public:
  ZExpAxialFormFactorModel();
  ZExpAxialFormFactorModel(string config);
  virtual ~ZExpAxialFormFactorModel();

  // implement the AxialFormFactorModelI interface
  double FA (const Interaction * interaction) const;

  // overload Algorithm's Configure() 
  void   Configure  (const Registry & config);
  void   Configure  (string param_set);

private:

  // calculate z parameter used in expansion
  double CalculateZ(double q2) const;
  void FixCoeffs (void);
  void FixA0     (void);
  void FixQ4Limit(void);
  void LoadConfig(void);

  bool   fQ4limit;
  int    fKmax;
  double fT0;
  double fTcut;
  double fFA0;
  //double fZ_An[11];
  double* fZ_An;
};

}         // genie namespace

#endif    // _Z_EXPANSION_AXIAL_FORM_FACTOR_MODEL_H_
