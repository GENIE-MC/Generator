//____________________________________________________________________________
/*!

\class    genie::BBA05ELFormFactorsModel

\brief    Concrete implementation of the ELFormFactorsModelI interface.
          Computes elastic form factors using the BBA2005 parameterization.

\ref      R.Bradford, A.Bodek, H.Budd and J.Arrington, hep-ex/0602017

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  Oct 19, 2005

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _BBA2005_EL_FORM_FACTORS_MODEL_H_
#define _BBA2005_EL_FORM_FACTORS_MODEL_H_

#include "Physics/QuasiElastic/XSection/ELFormFactorsModelI.h"

namespace genie {

typedef struct SBBA2005Fit
{
  double a0, a1, a2, b1, b2, b3, b4;
}
BBA2005Fit_t;

class BBA05ELFormFactorsModel : public ELFormFactorsModelI {

public:
  BBA05ELFormFactorsModel();
  BBA05ELFormFactorsModel(string config);
  virtual ~BBA05ELFormFactorsModel();

  // implement the ELFormFactorsModelI interface
  double Gep (const Interaction * interaction) const;
  double Gmp (const Interaction * interaction) const;
  double Gen (const Interaction * interaction) const;
  double Gmn (const Interaction * interaction) const;

  // overload Algorithm's Configure() to load the BBA2005Fit_t
  // structs from the configuration Registry
  void   Configure  (const Registry & config);
  void   Configure  (string param_set);

private:

  // fill data members from the configuration Registry
  void LoadConfig(void);

  // the actual BBA2005 fit function
  double BBA05Fit (double tau, const BBA2005Fit_t & fp) const;
  double tau      (const Interaction * interaction) const;

  // model parameters.
  BBA2005Fit_t fGep;   ///< BBA2005 fit coefficients for Gep
  BBA2005Fit_t fGen;   ///< BBA2005 fit coefficients for Gen
  BBA2005Fit_t fGmp;   ///< BBA2005 fit coefficients for Gmp
  BBA2005Fit_t fGmn;   ///< BBA2005 fit coefficients for Gmn
  double       fMuP;   ///< Anomalous proton magnetic moment
  double       fMuN;   ///< Anomalous neutron magnetic moment
};

}         // genie namespace

#endif    // _BBA2005_EL_FORM_FACTORS_MODEL_H_
