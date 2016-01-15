//____________________________________________________________________________
/*!

\class    genie::BaryonResonanceDecayer

\brief    Baryon resonance decayer.

          A simple decayer based on resonance's branching fractions (BRs) and
          an N-body phase space generator. Since the resonance can be produced
          off-shell, decay channels with total-mass > W are suppressed. \n

          Is a concrete implementation of the DecayModelI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  November 27, 2004

\cpright  Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _BARYON_RESONANCE_DECAYER_H_
#define _BARYON_RESONANCE_DECAYER_H_

#include <TGenPhaseSpace.h>
#include <TLorentzVector.h>

#include "Decay/DecayModelI.h"

namespace genie {

class BaryonResonanceDecayer : public DecayModelI {

public:
  BaryonResonanceDecayer();
  BaryonResonanceDecayer(string config);
  virtual ~BaryonResonanceDecayer();

  // implement the DecayModelI interface
  bool           IsHandled      (int pdgc)                      const;
  void           Initialize     (void)                          const;
  TClonesArray * Decay          (const DecayerInputs_t & inp)   const;
  double         Weight         (void)                          const;
  void           InhibitDecay   (int pdg, TDecayChannel * dc=0) const;
  void           UnInhibitDecay (int pdg, TDecayChannel * dc=0) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void           LoadConfig     (void);
  TClonesArray * DecayExclusive (int pdgc, TLorentzVector & p, TDecayChannel * ch) const;
  double         FinalStateMass (TDecayChannel * channel) const;

  mutable TGenPhaseSpace fPhaseSpaceGenerator;
  mutable double         fWeight;

  bool fGenerateWeighted;
};

}         // genie namespace

#endif    // _BARYON_RESONANCE_DECAYER_H_
