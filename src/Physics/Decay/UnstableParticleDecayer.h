//____________________________________________________________________________
/*!

\class    genie::UnstableParticleDecayer

\brief    A hook for concrete particle decayers in the chain of event
          processing modules.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  November 17, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _UNSTABLE_PARTICLE_DECAYER_H_
#define _UNSTABLE_PARTICLE_DECAYER_H_

#include <vector>

#include "Framework/EventGen/EventRecordVisitorI.h"

using std::vector;

namespace genie {

class GHepParticle;

class UnstableParticleDecayer : public EventRecordVisitorI {

public :

  UnstableParticleDecayer();
  UnstableParticleDecayer(string config);
  ~UnstableParticleDecayer();

  // Implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event) const;

  // Overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void  LoadConfig (void);
  vector <const EventRecordVisitorI *> fDecayers;///< list of all specified decayers
};

}      // genie namespace
#endif // _UNSTABLE_PARTICLE_DECAYER_H_
