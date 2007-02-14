//____________________________________________________________________________
/*!

\class    genie::UnstableParticleDecayer

\brief    Decays unstable particles found in the generated event record.
          After the interaction vertex generation it visits the event record
          and it decays the unstable particles using an externally specified
          particle decay model. The decay products are added to the event
          record and the status of parent particle is toggled. \n
          Is a concerete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  November 17, 2004

\cpright  Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _UNSTABLE_PARTICLE_DECAYER_H_
#define _UNSTABLE_PARTICLE_DECAYER_H_

#include <vector>

#include "EVGCore/EventRecordVisitorI.h"

using std::vector;

namespace genie {

class GHepParticle;
class DecayModelI;

class UnstableParticleDecayer : public EventRecordVisitorI {

public :

  UnstableParticleDecayer();
  UnstableParticleDecayer(string config);
  ~UnstableParticleDecayer();

  //-- implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void  LoadConfig        (void);
  bool  ToBeDecayed       (GHepParticle * particle) const;
  bool  IsUnstable        (GHepParticle * particle) const;
  void  CopyToEventRecord (TClonesArray * dp, GHepRecord * ev, 
                           int mother_pos, bool in_nucleus) const;

  double fMaxLifetime;       ///< define "unstable" particle?
  bool   fRunBefHadroTransp; ///< is being run before or after hadron transport?

  vector<const DecayModelI *> * fDecayers; ///< list of all specified decayers

  mutable const DecayModelI * fCurrDecayer; ///< current selected decayer

};

}      // genie namespace

#endif // _UNSTABLE_PARTICLE_DECAYER_H_
