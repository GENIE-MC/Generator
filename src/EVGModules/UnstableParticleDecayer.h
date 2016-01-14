//____________________________________________________________________________
/*!

\class    genie::UnstableParticleDecayer

\brief    Decays unstable particles found in the generated event record.
          After the interaction vertex generation it visits the event record
          and it decays the unstable particles using an externally specified
          particle decay model. The decay products are added to the event
          record and the status of parent particle is toggled. \n
          Is a concerete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  November 17, 2004

\cpright  Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _UNSTABLE_PARTICLE_DECAYER_H_
#define _UNSTABLE_PARTICLE_DECAYER_H_

#include <vector>

#include "EVGCore/EventRecordVisitorI.h"
#include "PDG/PDGCodeList.h"

using std::vector;

namespace genie {

class GHepParticle;
class DecayModelI;

class UnstableParticleDecayer : public EventRecordVisitorI {

public :

  UnstableParticleDecayer();
  UnstableParticleDecayer(string config);
  ~UnstableParticleDecayer();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void  LoadConfig        (void);
  bool  ToBeDecayed       (GHepParticle * particle) const;
  bool  IsUnstable        (GHepParticle * particle) const;
  void  CopyToEventRecord (TClonesArray * dp, GHepRecord * ev, GHepParticle * p,
                           int mother_pos, bool in_nucleus) const;

  bool                           fRunBefHadroTransp;   ///< is invoked before or after hadron transport?
  PDGCodeList                    fParticlesToDecay;    ///< list of particles to be decayed
  PDGCodeList                    fParticlesNotToDecay; ///< list of particles for which decay is inhibited
  vector <const DecayModelI *> * fDecayers;            ///< list of all specified decayers
  mutable const DecayModelI *    fCurrDecayer;         ///< current selected decayer

  //double fMaxLifetime; ///< define "unstable" particle
};

}      // genie namespace

#endif // _UNSTABLE_PARTICLE_DECAYER_H_
