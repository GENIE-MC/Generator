#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_INCL_ENABLED__

#ifndef G4INCLGENIEDISChannel_HH_
#define G4INCLGENIEDISChannel_HH_ 1

#include "G4INCLParticle.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"
#include "Physics/HadronTransport/G4INCLGENIEParticleRecord.h"


namespace G4INCL {
  class GENIEDISChannel : public IChannel {

  public:
    GENIEDISChannel(Particle *p, Nucleus *n, std::vector<GENIEParticleRecord> *eventRecord);
    virtual ~GENIEDISChannel();

    void fillFinalState(FinalState *fs);

  private:
    Particle *hitParticle;
    Nucleus * theNucleus;
    std::vector<GENIEParticleRecord> *genie_evtrec; // GENIE event record

    INCL_DECLARE_ALLOCATION_POOL(GENIEDISChannel)
  };

}

#endif

#endif // __GENIE_INCL_ENABLED__
