//____________________________________________________________________________
/*!

\class    genie::BBA07ELFormFactorsModel

\brief    Computes elastic form factors using the BBA2007 parameterization.
          Concrete implementation of the ELFormFactorsModelI interface.

\ref      A.Bodek, R.Bradford, H.Budd and S.Avvakumov, 
          Euro.Phys.J.C53 (2008)

	  Adapted by code provided by authors at:
          http://www.pas.rochester.edu/~bodek/FF/

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 31, 2008

\cpright  Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _BBA2007_EL_FORM_FACTORS_MODEL_H_
#define _BBA2007_EL_FORM_FACTORS_MODEL_H_

#include "ElFF/ELFormFactorsModelI.h"

namespace genie {

class BBA07ELFormFactorsModel : public ELFormFactorsModelI {

public:
  BBA07ELFormFactorsModel();
  BBA07ELFormFactorsModel(string config);
  virtual ~BBA07ELFormFactorsModel();

  // implement the ELFormFactorsModelI interface
  double Gep (const Interaction * interaction) const;
  double Gmp (const Interaction * interaction) const;
  double Gen (const Interaction * interaction) const;
  double Gmn (const Interaction * interaction) const;

  // overload Algorithm's Configure() to load the BBA2007 parameters
  // structs from the configuration Registry
  void   Configure  (const Registry & config);
  void   Configure  (string param_set);

private:

  // fill data members from the configuration Registry
  void LoadConfig(void);

  // various parametrizations
  double Lagrange      (const Interaction * interaction, double* par) const;
  double Kelly         (const Interaction * interaction, double* par) const;
  double GalsterFactor (const Interaction * interaction, double* par) const;

  // various kinematical params 
  double Tau (const Interaction * interaction) const;
  double Xi  (const Interaction * interaction) const;

  // model parameters.
  double fMuP;   ///< Anomalous proton  magnetic moment
  double fMuN;   ///< Anomalous neutron magnetic moment

  static double fsParGalsterFactor  [3]; ///<
  static double fsParGepKelly       [5]; ///<
  static double fsParGmpKelly       [5]; ///<
  static double fsParGepLagrange    [8]; ///<
  static double fsParGmpLagrange    [8]; ///<
  static double fsParGmnLagrange_25 [8]; ///<
  static double fsParGmnLagrange_43 [8]; ///<
  static double fsParGenLagrange_25 [8]; ///<
  static double fsParGenLagrange_43 [8]; ///<
};

}         // genie namespace
#endif    // _BBA2007_EL_FORM_FACTORS_MODEL_H_
