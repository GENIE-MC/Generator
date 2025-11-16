//____________________________________________________________________________
/*!

\class    genie::ZExpELFormFactorModel

\brief    Concrete implementation of the ELFormFactorModelI interface.
          Computes the EL form factor using the model-independent
          z-expansion technique

\ref      Hill et al.
          arXiv:1008.4619
          DOI: 10.1103/PhysRevD.82.113005

\author   Kaushik Borah <kaushik.borah99 \at uky.edu>
          based off DipoleELFormFactorsModel by
          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  August 16, 2013

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _Z_EXPANSION_EL_FORM_FACTOR_MODEL_H_
#define _Z_EXPANSION_EL_FORM_FACTOR_MODEL_H_

#include "Physics/QuasiElastic/XSection/ELFormFactorsModelI.h"

namespace genie {

class ZExpELFormFactorModel : public ELFormFactorsModelI {

public:
  ZExpELFormFactorModel();
  ZExpELFormFactorModel(string config);
  virtual ~ZExpELFormFactorModel();

  // implement the ELFormFactorModelI interface
  double Gep (const Interaction * interaction) const override;
  double Gen (const Interaction * interaction) const override;
  double Gmp (const Interaction * interaction) const override;
  double Gmn (const Interaction * interaction) const override;

  // overload Algorithm's Configure()
  void   Configure  (const Registry & config) override;
  void   Configure  (string param_set) override;

private:

  // calculate z parameter used in expansion
  double CalculateZ(double q2) const;
  void FixCoeffs (void);
  void FixEL0     (void);
  void FixQ4Limit(void);
  void LoadConfig(void);

  bool   fQ4limit;
  int    fKmax;
  double fT0;
  double fTcut;
  double fGep0;
  double fGmp0;
  double fGen0;
  double fGmn0;
  //double fZ_An[11];
  std::vector<double> fZ_APn;
  std::vector<double> fZ_BPn;
  std::vector<double> fZ_ANn;
  std::vector<double> fZ_BNn;
};

}         // genie namespace

#endif    // _Z_EXPANSION_EL_FORM_FACTOR_MODEL_H_

