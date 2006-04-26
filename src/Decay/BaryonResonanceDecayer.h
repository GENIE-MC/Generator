//____________________________________________________________________________
/*!

\class    genie::BaryonResonanceDecayer

\brief    Baryon resonance decayer.

          A simple decayer based on resonance's branching fractions (BRs) and
          an N-body phase space generator. Since the resonance can be produced
          off-shell, decay channels with total-mass > W are suppressed. \n

          Unlike PythiaDecayer, which relies on PYTHIA, this algorithm is not
          based on any rich underlying model. Should be used for decaying
          baryon resonances only. \n

          The baryon resonance PDG codes follow the MINOS extensions to
          PDG tables. \n

          The BaryonResonanceDecayer is a concrete implementation of the
          DecayModelI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  November 27, 2004

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _BARYON_RESONANCE_DECAYER_H_
#define _BARYON_RESONANCE_DECAYER_H_

#include <TGenPhaseSpace.h>

#include "Decay/DecayModelI.h"

namespace genie {

class BaryonResonanceDecayer : public DecayModelI {

public:

  BaryonResonanceDecayer();
  BaryonResonanceDecayer(string config);
  virtual ~BaryonResonanceDecayer();

  //-- implement the DecayModelI interface

  bool           IsHandled  (int pdgc)                    const;
  void           Initialize (void)                        const;
  TClonesArray * Decay      (const DecayerInputs_t & inp) const;

private:

  TClonesArray * DecayExclusive (int pdgc, TLorentzVector & p, TDecayChannel * ch) const;
  double         FinalStateMass (TDecayChannel * channel) const;

  mutable TGenPhaseSpace fPhaseSpaceGenerator;
};

}         // genie namespace

#endif    // _BARYON_RESONANCE_DECAYER_H_
