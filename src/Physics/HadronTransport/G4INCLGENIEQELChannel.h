#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_INCL_ENABLED__

#ifndef G4INCLGENIEQELChannel_HH_
#define G4INCLGENIEQELChannel_HH_ 1

#include "G4INCLParticle.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"
#include "Physics/HadronTransport/G4INCLGENIEParticleRecord.h"

namespace G4INCL {
  class GENIEQELChannel : public IChannel {

    public:
      GENIEQELChannel(Particle *p, Nucleus *n, std::vector<GENIEParticleRecord> *eventRecord);
      virtual ~GENIEQELChannel();

      void fillFinalState(FinalState *fs);

    private:
      Particle *hitParticle;
      Nucleus * theNucleus;
      std::vector<GENIEParticleRecord> *genie_evtrec; // GENIE event record

      INCL_DECLARE_ALLOCATION_POOL(GENIEQELChannel)
  };

}

#endif

#endif // __GENIE_INCL_ENABLED__
