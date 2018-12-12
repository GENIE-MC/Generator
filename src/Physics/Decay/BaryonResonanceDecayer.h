//____________________________________________________________________________
/*!

\class    genie::BaryonResonanceDecayer

\brief    Baryon resonance decayer module.

          A simple resonance decay simulation built using resonance branching
          fraction data and an N-body phase space generator.
          Since the resonance can be produced off-the-mass-shell, decay
          channels with total-mass > W are suppressed.

          Is a concerete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  November 27, 2004

\cpright  Copyright (c) 2003-2018, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _BARYON_RESONANCE_DECAYER_H_
#define _BARYON_RESONANCE_DECAYER_H_

#include <TGenPhaseSpace.h>
#include <TLorentzVector.h>

#include "Physics/Decay/Decayer.h"

namespace genie {

class GHepParticle;
class BaryonResonanceDecayer : protected Decayer {

public:
  BaryonResonanceDecayer();
  BaryonResonanceDecayer(string config);
  virtual ~BaryonResonanceDecayer();

  // Implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event) const;

  virtual void LoadConfig    (void);

private:

  void           Initialize        (void) const;
  bool           IsHandled         (int pdgc) const;
  void           InhibitDecay      (int pdgc, TDecayChannel * ch=0) const;
  void           UnInhibitDecay    (int pdgc, TDecayChannel * ch=0) const;
  double         Weight            (void) const;
  bool           Decay             (int dec_part_id, GHepRecord * event) const;
  TDecayChannel* SelectDecayChannel(int dec_part_id, GHepRecord * event) const;
  void           DecayExclusive    (int dec_part_id, GHepRecord * event, TDecayChannel * ch) const;

  TObjArray *    EvolveDeltaBR        (int dec_part_pdgc, TObjArray * decay_list, double W) const;
  double         EvolveDeltaDecayWidth(int dec_part_pdgc, TDecayChannel * ch, double W) const;

  double         FinalStateMass    ( TDecayChannel * ch ) const;
  bool           IsPiNDecayChannel ( TDecayChannel * ch ) const;

  mutable TGenPhaseSpace fPhaseSpaceGenerator;
  mutable double         fWeight;

  double fProb32, fProb12 ; // parameters for pion angular distribution for Delta -> N + pi

  double fFFScaling ;  // Scaling factor of the form factor of the Delta wrt to Q2

};

}         // genie namespace
#endif    // _BARYON_RESONANCE_DECAYER_H_
