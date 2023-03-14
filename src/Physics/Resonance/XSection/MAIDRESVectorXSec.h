//____________________________________________________________________________
/*!

\class    genie::MAIDRESVectorXSec

\brief    Computes the double differential cross section for resonance
          electro- or neutrino-production according to the Rein-Sehgal model
	  using MAID Form Factors.

          The computed cross section is the d^2 xsec/ dQ^2 dW \n

          where \n
            \li \c Q^2 : momentum transfer ^ 2
            \li \c W   : invariant mass of the final state hadronic system

          Is a concrete implementation of the XSecAlgorithmI interface.

\ref      

\author   Julia Tena Vidal <jtenavidal@tauex.tau.ac.il>
          Tel Aviv Unversity

\created  Feb 28, 2023

\cpright  Copyright (c) 2023-2033, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _MAID_RES_VEC_XSECMODEL_H_
#define _MAID_RES_VEC_XSECMODEL_H_

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/ParticleData/BaryonResonance.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Physics/Resonance/XSection/RESVectFormFactorsI.h"

namespace genie {

class MAIDRESVectorXSec : public XSecAlgorithmI {

public:
  MAIDRESVectorXSec();
  MAIDRESVectorXSec(string config);
  virtual ~MAIDRESVectorXSec();

  // implement the XSecAlgorithmI interface
  double XSec         (const Interaction * i, KinePhaseSpace_t k) const;
  double Integral     (const Interaction * i) const;
  bool   ValidProcess (const Interaction * i) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void LoadConfig (void);

  RESVectFormFactorsI * fVFFEMp ;
  RESVectFormFactorsI * fVFFEMn ;

  bool     fUsingDisResJoin;   ///< use a DIS/RES joining scheme?
  double   fWcut;              ///< apply DIS/RES joining scheme < Wcut
  double   fN2ResMaxNWidths;   ///< limits allowed phase space for n=2 res
  double   fN0ResMaxNWidths;   ///< limits allowed phase space for n=0 res
  double   fGnResMaxNWidths;   ///< limits allowed phase space for other res
  string   fKFTable;             ///< table of Fermi momentum (kF) constants for various nuclei
  bool     fUseRFGParametrization; ///< use parametrization for fermi momentum insted of table?
  bool     fUsePauliBlocking;      ///< account for Pauli blocking?
  double   fRESScaling ; //< Overall scaling 

  const XSecIntegratorI * fXSecIntegrator;

};

}       // genie namespace

#endif  // _MAID_RES_VEC_XSECMODEL_H_
