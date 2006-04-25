//____________________________________________________________________________
/*!

\class    genie::RSExclusiveRESPXSec

\brief    Computes the differential cross section for an exclusive 1-pion 
          reaction through resonance neutrinoproduction according to the 
          Rein-Seghal model.

          The computed xsec is the double differential d^2 xsec/ dQ^2 dW \n
          where \n
          \li Q^2 : momentum transfer ^ 2
          \li W   : invariant mass of the final state hadronic system

          The cross section is computed for an input list of resonances
          as the sum of the Rein-Seghal single resonance cross sections
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

          Is a concrete implementation of the XSecAlgorithmI initerface.

\ref      D.Rein and L.M.Seghal, Neutrino Excitation of Baryon Resonances
          and Single Pion Production, Ann.Phys.133, 79 (1981)

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  November 22, 2004

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _REIN_SEGHAL_EXCLUSIVE_SPP_PARTIAL_XSEC_H_
#define _REIN_SEGHAL_EXCLUSIVE_SPP_PARTIAL_XSEC_H_

#include "Base/XSecAlgorithmI.h"
#include "BaryonResonance/BaryonResList.h"

namespace genie {

class ReinSeghalSPPPXSec : public XSecAlgorithmI {

public:

  ReinSeghalSPPPXSec();
  ReinSeghalSPPPXSec(string config);
  virtual ~ReinSeghalSPPPXSec();

  //-- XSecAlgorithmI interface implementation
  double XSec            (const Interaction * interaction) const;
  bool   ValidProcess    (const Interaction * interaction) const;
  bool   ValidKinematics (const Interaction * interaction) const;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  double XSecNRES(const Interaction * interaction) const;
  double XSec1RES(const Interaction * interaction) const;

  void LoadSubAlg       (void);
  void GetResonanceList (void);

  BaryonResList          fResList;
  const XSecAlgorithmI * fSingleResXSecModel;
};

}       // genie namespace
#endif  // _REIN_SEGHAL_EXCLUSIVE_SPP_PARTIAL_XSEC_H_
