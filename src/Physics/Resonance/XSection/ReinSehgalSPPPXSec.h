//____________________________________________________________________________
/*!

\class    genie::ReinSehgalSPPPXSec

\brief    Computes the differential cross section for an exclusive 1-pion 
          reaction through resonance neutrinoproduction according to the 
          Rein-Sehgal model.

          The cross section is computed for an input list of resonances
          as the sum of the Rein-Sehgal single resonance cross sections
          weighted:

          \li With the value of their Breit-Wigner distributions at the given
              W,Q^2 (The code for BW weighting is included in the single
              resonance cross section algorithm. The user needs to make sure
              that he does not run the single resonance cross section code with
              a configuration that inhibits weighting).

          \li With the isospin Glebsch-Gordon coefficient determining the
              contribution of each resonance to the exclusive final state.

          \li With the BR for the produced resonance to decay into the given
              exclusive final state.

          In this algorithm we follow the non-coherent approach: we sum
          the weighted resonance production cross sections rather than the
          resonance production amplitudes.

          Is a concrete implementation of the XSecAlgorithmI interface.

\ref      D.Rein and L.M.Sehgal, Neutrino Excitation of Baryon Resonances
          and Single Pion Production, Ann.Phys.133, 79 (1981)

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  November 22, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _REIN_SEHGAL_EXCLUSIVE_SPP_PXSEC_H_
#define _REIN_SEHGAL_EXCLUSIVE_SPP_PXSEC_H_

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/ParticleData/BaryonResList.h"

namespace genie {

class XSecIntegratorI;

class ReinSehgalSPPPXSec : public XSecAlgorithmI {

public:
  ReinSehgalSPPPXSec();
  ReinSehgalSPPPXSec(string config);
  virtual ~ReinSehgalSPPPXSec();

  //-- XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k) const;
  double Integral        (const Interaction * i) const;
  bool   ValidProcess    (const Interaction * i) const;
	
  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  //-- load algorithm configuration when Algorithm::Configure() 
  void LoadConfig (void);

  double XSecNRES(const Interaction * i, KinePhaseSpace_t k) const;
  double XSec1RES(const Interaction * i, KinePhaseSpace_t k) const;

  //-- private data members
  BaryonResList           fResList;
  const XSecAlgorithmI *  fSingleResXSecModel;
  const XSecIntegratorI * fXSecIntegrator;
};

}       // genie namespace
#endif  // _REIN_SEHGAL_EXCLUSIVE_SPP_PXSEC_H_
