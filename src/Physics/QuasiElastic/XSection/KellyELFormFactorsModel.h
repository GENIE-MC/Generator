//____________________________________________________________________________
/*!

\class    genie::KellyELFormFactorsModel

\brief    Concrete implementation of the ELFormFactorsModelI interface.
          Computes elastic form factors using the Kelly parameterization.
 	  Based on J.J. Kelly, Phys.Rev.C 70 (2004) 068202

\author   Noah Steinberg <nsteinbe \at fnal.gov>
	  Fermi National Accelerator Laboratory 	

\created  Sept 26, 2023

\cpright  Copyright (c) 2003-2023, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org        
*/
//____________________________________________________________________________

#ifndef _KELLY_EL_FORM_FACTORS_MODEL_H_
#define _KELLY_EL_FORM_FACTORS_MODEL_H_

#include "Physics/QuasiElastic/XSection/ELFormFactorsModelI.h"

namespace genie {

typedef struct SKellyFit
{
  double a0, a1, a2, b1, b2, b3, b4;
}
KellyFit_t;

class KellyELFormFactorsModel : public ELFormFactorsModelI {

public:
  KellyELFormFactorsModel();
  KellyELFormFactorsModel(string config);
  virtual ~KellyELFormFactorsModel();

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

  // the actual Kelly fit function
  double KellyFit (double tau, const KellyFit_t & fp) const;
  double tau      (const Interaction * interaction) const;
  double KellyVectorDipole (const Interaction * interaction) const;

  // model parameters.
  KellyFit_t fGep;     ///< Kelly fit coefficients for Gep
  KellyFit_t fGen;     ///< Kelly fit coefficients for Gen
  KellyFit_t fGmp;     ///< Kelly fit coefficients for Gmp
  KellyFit_t fGmn;     ///< Kelly fit coefficients for Gmn
  double       fMuP;   ///< Anomalous proton magnetic moment
  double       fMuN;   ///< Anomalous neutron magnetic moment
  double       fMv;    ///< Vector dipole mass
  double       fMv2;   ///< Vector dipole mass squared
};

}         // genie namespace

#endif    // _KELLY_EL_FORM_FACTORS_MODEL_H_
