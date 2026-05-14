#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_INCL_ENABLED__

#include "G4INCLGENIERESChannel.h"
#include "G4INCLRandom.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLParticleTable.hh"
#include "G4INCLCrossSections.hh"
#include "G4INCLGlobals.hh"

namespace G4INCL {

  GENIERESChannel::GENIERESChannel(Particle *p, Nucleus *n, std::vector<GENIEParticleRecord> *eventRecord)
    : hitParticle(p), theNucleus(n),
    genie_evtrec(eventRecord)
  {
  }

  GENIERESChannel::~GENIERESChannel()
  {
  }

  void GENIERESChannel::fillFinalState(FinalState *fs)
  {


    std::vector<GENIEParticleRecord>::iterator ip;
    for(ip = genie_evtrec->begin(); ip != genie_evtrec->end(); ip++){
      if(ip->Status() == genie::kIStHadronInTheNucleus || ip->Status() == genie::kIStPreDecayResonantState){
        if(abs(ip->Pdg()) > 1000){ // is not meson
          ip->setID(int(hitParticle->getID()));
          hitParticle->setType(ip->Type());
          hitParticle->setMass(ip->Mass());
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
