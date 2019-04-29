//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab
*/
//____________________________________________________________________________

#include <sstream>

#include "Framework/GHEP/GHepRecord.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Registry/Registry.h"
#include "Physics/Decay/UnstableParticleDecayer.h"

using std::ostringstream;

using namespace genie;
//___________________________________________________________________________
UnstableParticleDecayer::UnstableParticleDecayer() :
EventRecordVisitorI("genie::UnstableParticleDecayer")
{

}
//___________________________________________________________________________
UnstableParticleDecayer::UnstableParticleDecayer(string config) :
EventRecordVisitorI("genie::UnstableParticleDecayer", config)
{

}
//___________________________________________________________________________
UnstableParticleDecayer::~UnstableParticleDecayer()
{
  fDecayers.clear();
}
//___________________________________________________________________________
void UnstableParticleDecayer::ProcessEventRecord(GHepRecord * event) const
{
  vector<const EventRecordVisitorI *>::const_iterator it = fDecayers.begin();
  for( ; it != fDecayers.end(); ++it)
  {
    const EventRecordVisitorI * decayer = *it;
    decayer->ProcessEventRecord(event);
  }
}
//___________________________________________________________________________
void UnstableParticleDecayer::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();

  fAllowReconfig = false;
}
//___________________________________________________________________________
void UnstableParticleDecayer::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();

  fAllowReconfig = false;
}
//___________________________________________________________________________
void UnstableParticleDecayer::LoadConfig(void)
{
  fDecayers.clear();

  // Load particle decayers
  // Order is important if both decayers can handle a specific particle
  // as only the first would get the chance to decay it
  int ndec = 0 ;
  this->GetParam("NDecayers", ndec);
  assert(ndec>0);

  for(int idec = 0; idec < ndec; idec++) {
     ostringstream alg_key;
     alg_key     << "Decayer-" << idec;
     const EventRecordVisitorI * decayer =
        dynamic_cast<const EventRecordVisitorI *>
            (this->SubAlg(alg_key.str()));
     fDecayers.push_back(decayer);
  }
}
//___________________________________________________________________________
