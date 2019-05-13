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

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
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
  TDecayChannel* SelectDecayChannel(int dec_part_id, GHepRecord * event, bool & to_be_deleted ) const;
  // the flag to_be_deleted is referred to the returned decay channel 
  bool           DecayExclusive    (int dec_part_id, GHepRecord * event, TDecayChannel * ch) const;

  // Methods specific for Delta decay
  TObjArray *    EvolveDeltaBR        (int dec_part_pdgc, TObjArray * decay_list, double W) const;
  double         EvolveDeltaDecayWidth(int dec_part_pdgc, TDecayChannel * ch, double W) const;
  bool           AcceptPionDecay( TLorentzVector lab_pion, int dec_part_id, const GHepRecord * event ) const ;

  double         FinalStateMass    ( TDecayChannel * ch ) const;
  bool           IsPiNDecayChannel ( TDecayChannel * ch ) const;

  static bool IsDelta( int dec_part_pdgc ) ; 
  static bool HasEvolvedBRs( int dec_part_pdgc ) ; 

  mutable TGenPhaseSpace fPhaseSpaceGenerator;
  mutable double         fWeight;

  bool   fDeltaThetaOnly ;

  std::vector<double> fR33, fR31, fR3m1 ;
  std::vector<double> fQ2Thresholds ;

  std::vector<double> fW_max ;

  double fFFScaling ;  // Scaling factor of the form factor of the Delta wrt to Q2

};

}         // genie namespace
#endif    // _BARYON_RESONANCE_DECAYER_H_
