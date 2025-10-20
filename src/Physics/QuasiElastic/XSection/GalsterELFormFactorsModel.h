//____________________________________________________________________________
/*!

\class    genie::GalsterELFormFactorsModel

\brief    Concrete implementation of the ELFormFactorsModelI interface.
          Computes elastic form factors using the Galster parameterization.

\ref      S.Galster et al., Nuclear Physics B32 (1971) 221-237

\author   Igor Kakorin <kakorin@jinr.ru>, Joint Institute for Nuclear Research

\created  May 19, 2019

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _GALSTER_EL_FORM_FACTORS_MODEL_H_
#define _GALSTER_EL_FORM_FACTORS_MODEL_H_

#include "Physics/QuasiElastic/XSection/ELFormFactorsModelI.h"

namespace genie {


class GalsterELFormFactorsModel : public ELFormFactorsModelI {

public:
  GalsterELFormFactorsModel();
  GalsterELFormFactorsModel(string config);
  virtual ~GalsterELFormFactorsModel();

  // implement the ELFormFactorsModelI interface
  double Gep (const Interaction * interaction) const;
  double Gmp (const Interaction * interaction) const;
  double Gen (const Interaction * interaction) const;
  double Gmn (const Interaction * interaction) const;

  // overload Algorithm's Configure() to load the BBA2003Fit_t
  // structs from the configuration Registry
  void   Configure  (const Registry & config);
  void   Configure  (string param_set);

private:

  // fill data members from the configuration Registry
  void LoadConfig(void);

  // model parameters.
  double       fGenp;  ///< parameter for Gen
  double       fMv;    ///< Elactic vector mass 
  double       fMv2;   ///< Elactic vector mass 
  double       fMuP;   ///< Anomalous proton magnetic moment
  double       fMuN;   ///< Anomalous neutron magnetic moment
  
  bool fIsIsoscalarNucleon;  ///< Is assuming isoscalar nucleon?
};

}         // genie namespace

#endif    // _GALSTER_EL_FORM_FACTORS_MODEL_H_
