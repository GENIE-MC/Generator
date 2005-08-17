//____________________________________________________________________________
/*!

\class   genie::UnstableParticleDecayer

\brief   Decays unstable particles found in the generated event record.

         After the interaction vertex generation it visits the event record
         and it decays the unstable particles using an externally specified
         particle decay model. The decay products are added to the event
         record and the status of parent particle is toggled. \n

         Is a concerete implementation of the EventRecordVisitorI interface.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created November 17, 2004

*/
//____________________________________________________________________________

#ifndef _UNSTABLE_PARTICLE_DECAYER_H_
#define _UNSTABLE_PARTICLE_DECAYER_H_

#include "EVGCore/EventRecordVisitorI.h"

namespace genie {

class GHepParticle;

class UnstableParticleDecayer : public EventRecordVisitorI {

public :

  UnstableParticleDecayer();
  UnstableParticleDecayer(const char * param_set);
  ~UnstableParticleDecayer();

  //-- implement the EventRecordVisitorI interface

  void ProcessEventRecord(GHepRecord * event_rec) const;

private:

  bool  ToBeDecayed (GHepParticle * particle) const;
  bool  IsUnstable  (GHepParticle * particle) const;

  void  CopyToEventRecord(TClonesArray * dp,
                            GHepRecord * event_rec, int mother_pos) const;
};

}      // genie namespace

#endif // _UNSTABLE_PARTICLE_DECAYER_H_
