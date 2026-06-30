#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_INCL_ENABLED__

#include "G4INCLGENIEDISChannel.h"
#include "G4INCLRandom.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLParticleTable.hh"
#include "G4INCLCrossSections.hh"
#include "G4INCLGlobals.hh"
#include "Framework/ParticleData/PDGCodes.h"

namespace G4INCL {

  GENIEDISChannel::GENIEDISChannel(Particle *p, Nucleus *n, std::vector<GENIEParticleRecord> *eventRecord)
    :hitParticle(p), theNucleus(n),
    genie_evtrec(eventRecord)
  {
  }

  GENIEDISChannel::~GENIEDISChannel()
  {
  }

  void GENIEDISChannel::fillFinalState(FinalState *fs)
  {
    std::vector<GENIEParticleRecord>::iterator ip;

    for(ip = genie_evtrec->begin(); ip != genie_evtrec->end(); ip++){
      if(ip->Status() == genie::kIStHadronInTheNucleus){
        if(ip->Pdg() == genie::kPdgNeutron || ip->Pdg() == genie::kPdgProton){
          ip->setID(int(hitParticle->getID()));
          hitParticle->setType(ip->Type());
          hitParticle->setMomentum(ip->P3());
          hitParticle->setPosition(ip->X3());
          hitParticle->adjustEnergyFromMomentum();
          fs->addModifiedParticle(hitParticle);
        }
        else{
          Particle *ihadron = new Particle(ip->Type(), ip->P3(), ip->X3());
          ip->setID(int(ihadron->getID()));
          ihadron->setMass(ip->Mass());
          ihadron->adjustEnergyFromMomentum();
          fs->addCreatedParticle(ihadron);
        }
      }
    }
  }
}

#endif // __GENIE_INCL_ENABLED__
