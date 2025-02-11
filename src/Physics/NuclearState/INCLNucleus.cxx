//____________________________________________________________________________
/*!

  \class    genie::INCLNucleus

  \brief    INCLXX nuclear model. Implements the NuclearModelI 
  interface.

  \ref      

  \author   Liang Liu, (liangliu@fnal.gov)

  \created  Oct. 2024

  \cpright  Copyright (c) 2003-2024, The GENIE Collaboration
  For the full text of the license visit http://copyright.genie-mc.org

*/
//____________________________________________________________________________

#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_INCL_ENABLED__

#include <cassert>
#include <iostream>

#include <TSystem.h>
#include <TNtupleD.h>
#include <TTree.h>

#include "Framework/Messenger/Messenger.h"
#include "Physics/NuclearState/INCLNucleus.h"
#include "Framework/Numerical/Spline.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"


#include "G4INCLCascade.hh"
#include "G4INCLRandom.hh"
#include "G4INCLStandardPropagationModel.hh"
#include "G4INCLParticleTable.hh"
#include "G4INCLParticle.hh"
#include "G4INCLNuclearMassTable.hh"
#include "G4INCLGlobalInfo.hh"
#include "G4INCLNucleus.hh"

#include "G4INCLPauliBlocking.hh"

#include "G4INCLCrossSections.hh"

#include "G4INCLPhaseSpaceGenerator.hh"

#include "G4INCLLogger.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLNuclearDensityFactory.hh"

#include "G4INCLINuclearPotential.hh"

#include "G4INCLCoulombDistortion.hh"

#include "G4INCLClustering.hh"

#include "G4INCLIntersection.hh"

#include "G4INCLBinaryCollisionAvatar.hh"

#include "G4INCLCascadeAction.hh"
#include "G4INCLAvatarDumpAction.hh"

#include <cstring> 
#include <cstdlib>
#include <numeric>

#include "G4INCLPbarAtrestEntryChannel.hh"


#include "G4INCLGeant4Compat.hh"
#include "G4INCLCascade.hh"

#include "G4INCLClustering.hh"
#include "G4INCLParticle.hh"
#include "G4INCLIAvatar.hh"
#include "G4INCLIPropagationModel.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLStandardPropagationModel.hh"
#include "G4INCLRandom.hh"
#include "G4INCLRanecu.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLParticleTable.hh"
#include "G4INCLKinematicsUtils.hh"
#include "g4inclpauliblocking.hh"
#include "g4inclphasespacegenerator.hh"
#include "g4inclcoulombdistortion.hh"
#include "g4inclbinarycollisionavatar.hh"

// signal handler (for Linux and GCC)
#include "G4INCLSignalHandling.hh"

// For I/O
#include "IWriter.hh"
#include "ASCIIWriter.hh"
#include "ProtobufWriter.hh"
#include "INCLTree.hh"
#include "ROOTWriter.hh"
#include "HDF5Writer.hh"

// For configuration
#include "G4INCLConfig.hh"

// For logging
#include "G4INCLLogger.hh"

// Generic de-excitation interface
#include "G4INCLIDeExcitation.hh"

// ABLA v3p de-excitation
#ifdef INCL_DEEXCITATION_ABLAXX
#include "G4INCLAblaInterface.hh"
#endif

// ABLACXX de-excitation
#ifdef INCL_DEEXCITATION_ABLACXX
#include "G4INCLAblaxxInterface.hh"
#endif

// ABLA07 de-excitation
#ifdef INCL_DEEXCITATION_ABLA07
#include "G4INCLAbla07Interface.hh"
#endif

// SMM de-excitation
#ifdef INCL_DEEXCITATION_SMM
#include "G4INCLSMMInterface.hh"
#endif

// GEMINIXX de-excitation
#ifdef INCL_DEEXCITATION_GEMINIXX
#include "G4INCLGEMINIXXInterface.hh"
#endif



// INCL++
#include "G4INCLConfig.hh"
#include "G4INCLCascade.hh"
#include "G4INCLConfigEnums.hh"
#include "G4INCLParticle.hh"
// signal handler (for Linux and GCC)
#include "G4INCLSignalHandling.hh"

// Generic de-excitation interface
#include "G4INCLIDeExcitation.hh"

// ABLA v3p de-excitation
#ifdef INCL_DEEXCITATION_ABLAXX
#include "G4INCLAblaInterface.hh"
#endif

// ABLA07 de-excitation
#ifdef INCL_DEEXCITATION_ABLA07
#include "G4INCLAbla07Interface.hh"
#endif

// SMM de-excitation
#ifdef INCL_DEEXCITATION_SMM
#include "G4INCLSMMInterface.hh"
#endif

// GEMINIXX de-excitation
#ifdef INCL_DEEXCITATION_GEMINIXX
#include "G4INCLGEMINIXXInterface.hh"
#endif



using std::cout;
using std::endl;

using namespace genie;

//____________________________________________________________________________
INCLNucleus * INCLNucleus::fInstance = 0;

//____________________________________________________________________________
INCLNucleus::INCLNucleus():propagationModel_(0)
{
  //  this->Load();
  fInstance = 0;
  nucleon_index_ = -1;
  nucleus_ = nullptr;
  theConfig_ = nullptr;
}
//____________________________________________________________________________
INCLNucleus::~INCLNucleus()
{
  //  if(!gAbortingInErr) {
  //    cout << "INCLNucleus singleton dtor: Deleting inputs... " << endl;
  //  }
  //  delete fNuclSupprD2;
}
//____________________________________________________________________________
INCLNucleus * INCLNucleus::Instance()
{
  if(fInstance == 0) {
    LOG("NuclData", pINFO) << "INCLNucleus late initialization";
    //    static INCLNucleus::Cleaner cleaner;
    //    cleaner.DummyMethodAndSilentCompiler();
    fInstance = new INCLNucleus;
    fInstance->theConfig_ = new G4INCL::Config();
    fInstance->nucleus_ = nullptr;
    fInstance->init();
  }
  return fInstance;
}


void INCLNucleus::init(){
  LOG("NuclData", pINFO) << "init()";
  theConfig_->init();
  theConfig_->setINCLXXDataFilePath("/root/inclxx/inclxx-v6.33.1-e5857a1/data"); // FIXME:: using config to set path
  // initialize INCL model
  G4INCL::Random::initialize(theConfig_);
  // Select the Pauli and CDPP blocking algorithms
  G4INCL::Pauli::initialize(theConfig_);
  // Set the phase-space generator
  G4INCL::PhaseSpaceGenerator::initialize(theConfig_);
  // Select the Coulomb-distortion algorithm:
  G4INCL::CoulombDistortion::initialize(theConfig_);
  // Select the clustering algorithm:
  G4INCL::Clustering::initialize(theConfig_);
  // Initialize the INCL particle table:
  G4INCL::ParticleTable::initialize(theConfig_);
  // Initialize the value of cutNN in BinaryCollisionAvatar
  G4INCL::BinaryCollisionAvatar::setCutNN(theConfig_->getCutNN());
  // Initialize the value of strange cross section bias
  G4INCL::BinaryCollisionAvatar::setBias(theConfig_->getBias());
  // Set the cross-section set
  G4INCL::CrossSections::initialize(theConfig_);
  //theConfig_->setLocalEnergyBBType(G4INCL::NeverLocalEnergy);
  //theConfig_->setLocalEnergyPiType(G4INCL::NeverLocalEnergy);


  // Propagation model is responsible for finding avatars and
  // transporting the particles. In principle this step is "hidden"
  // behind an abstract interface and the rest of the system does not
  // care how the transportation and avatar finding is done. This
  // should allow us to "easily" experiment with different avatar
  // finding schemes and even to support things like curved
  // trajectories in the future.
  propagationModel_ = new G4INCL::StandardPropagationModel(theConfig_->getLocalEnergyBBType(),theConfig_->getLocalEnergyPiType(),theConfig_->getHadronizationTime());
  if(theConfig_->getCascadeActionType() == G4INCL::AvatarDumpActionType)
    cascadeAction_ = new G4INCL::AvatarDumpAction();
  else
    cascadeAction_ = new G4INCL::CascadeAction();
//  cascadeAction_->beforeRunAction(theConfig_);



}


void INCLNucleus::initialize(const GHepRecord * evrec){

  LOG("NuclData", pINFO) << "initialize()";

  // initialize according process Event in INCL

  GHepParticle * nucleus = evrec->TargetNucleus();
  G4INCL::ParticleSpecies targetSpecies = G4INCL::ParticleSpecies(nucleus->Name());
  theConfig_->setTargetA(targetSpecies.theA);
  theConfig_->setTargetZ(targetSpecies.theZ);
  theConfig_->setTargetS(targetSpecies.theS);
  LOG("NuclData", pINFO) << "initialize()";
  // define Nucleus and initialize it

  if(nucleus_){
    if(!nucleus_->getStore()->getParticles().empty()){
      if(nucleus_->getA() == evrec->TargetNucleus()->A() && nucleus_->getZ() == evrec->TargetNucleus()->Z()){
//  LOG("NuclData", pINFO)  << nucleus_->print();
	return;
      }
      else{
	nucleus_->deleteParticles();
	nucleus_->getStore()->clear();
	delete nucleus_;
      }
    }
  }

  // ReInitialize the bias vector
  G4INCL::Particle::INCLBiasVector.clear();
  //Particle::INCLBiasVector.Clear();
  G4INCL::Particle::nextBiasedCollisionID = 0;

  // Set the target and the projectile 
  // implement prepare reaction
  
  // Reset the forced-transparent flag
  // forceTransparent = false; FIXME
  //
  // Initialise the maximum universe radius
  // INCL initialize universe radius according to particle species,
  // kenetic energy, and nucleus type. void INCL::initUniverseRadius(ParticleSpecies const &p, const double kineticEnergy, const int A, const int Z)
  // FIXME: I'm not sure if we need the interaction distance for neutrino scattering
  double rMax = 0.0;
  const double pMaximumRadius = G4INCL::ParticleTable::getMaximumNuclearRadius(G4INCL::Proton, targetSpecies.theA, targetSpecies.theZ);
  const double nMaximumRadius = G4INCL::ParticleTable::getMaximumNuclearRadius(G4INCL::Neutron, targetSpecies.theA, targetSpecies.theZ);
  const double maximumRadius = std::max(pMaximumRadius, nMaximumRadius);
  rMax = std::max(maximumRadius, rMax);
  maxUniverseRadius_ = rMax;

  // FIXME the last two parameters need to be configed
  // theConfig_, G4INCL::ParticleTable::getMaximumNuclearRadius(G4INCL::Proton, targetSpecies.theA, targetSpecies.theZ)
  // G4INCL::NType
  nucleus_ = new G4INCL::Nucleus(targetSpecies.theA, targetSpecies.theZ, targetSpecies.theS, 
  theConfig_, maxUniverseRadius_, G4INCL::Def);
  nucleus_->getStore()->getBook().reset();
  nucleus_->initializeParticles();
  propagationModel_->setNucleus(nucleus_);

  // initialize max interaction distance
  // FIXME: in INCL, composite has non-zero max interaction distance.
  maxInteractionDistance_ = 0;

  // set the min Remnant size to be 4
  // the min remnant is alpha particle
  // FIXME: it is only works for nuclei with large A
  minRemnantSize_ = 4;

  // cascade action is not related to simulation
  // it is just output the casade to file FIXME
  //cascadeAction_->beforeCascadeAction(propagationModel_);
  //
  // INCL need to decide whether the cascade can be ran or not
  // For genie, we need to run casecade for every events
  // const bool canRunCascade = preCascade(projectileSpecies, kineticEnergy);

  LOG("INCLNucleus", pDEBUG) << nucleus_->getStore()->getParticles().at(2)->getPotentialEnergy() ;
  LOG("INCLNucleus", pDEBUG) << nucleus_->getStore()->getParticles().at(2)->getEnergy() -  nucleus_->getStore()->getParticles().at(2)->getPotentialEnergy()  ;
  LOG("INCLNucleus", pDEBUG) << nucleus_->getStore()->getParticles().at(2)->getMomentum().print() ;
  LOG("INCLNucleus", pDEBUG) << nucleus_->getStore()->getParticles().at(2)->getPosition().print() ;

  GHepParticle * nucleon = evrec->HitNucleon();
  LOG("INCLNucleus", pDEBUG) << "hit nucleon pdg : " << nucleon->Pdg();

  RandomGen * rnd = RandomGen::Instance();
  if(pdg::IsProton(nucleon->Pdg()))
    nucleon_index_ = rnd->RndGen().Integer(targetSpecies.theZ);
  else if(pdg::IsNeutron(nucleon->Pdg()))
    nucleon_index_ = rnd->RndGen().Integer(targetSpecies.theA - targetSpecies.theZ) + targetSpecies.theZ;
  else
    exit(1);
  hitNucleon_ = nucleus_->getStore()->getParticles().at(nucleon_index_);
}

void INCLNucleus::reset(const GHepRecord * evrec){
  // nucleus must exsit!
  if(!nucleus_)  LOG("INCLNucleus", pFATAL) << "nucleus doesn't exsit!";
  // can't reset a nucleus with different type
  if(!(nucleus_->getA() == evrec->TargetNucleus()->A() && nucleus_->getZ() == evrec->TargetNucleus()->Z())) 
    LOG("INCLNucleus", pFATAL) << "you are try to reset a nucleus with different type!";
  // reset the nucleus
  if(nucleus_){
    nucleus_->deleteParticles();
    nucleus_->getStore()->clear();
    nucleus_->initializeParticles();
    nucleus_->getStore()->getBook().reset();

  GHepParticle * nucleon = evrec->HitNucleon();
  LOG("INCLNucleus", pDEBUG) << "hit nucleon pdg : " << nucleon->Pdg();
  RandomGen * rnd = RandomGen::Instance();
  if(pdg::IsProton(nucleon->Pdg()))
    nucleon_index_ = rnd->RndGen().Integer(evrec->TargetNucleus()->Z());
  else if(pdg::IsNeutron(nucleon->Pdg()))
    nucleon_index_ = rnd->RndGen().Integer(evrec->TargetNucleus()->A() - evrec->TargetNucleus()->Z()) + evrec->TargetNucleus()->Z();
  else
    exit(1);
  hitNucleon_ = nucleus_->getStore()->getParticles().at(nucleon_index_);


  delete propagationModel_;
  propagationModel_ = new G4INCL::StandardPropagationModel(theConfig_->getLocalEnergyBBType(),theConfig_->getLocalEnergyPiType(),theConfig_->getHadronizationTime());

  }
}
TVector3 INCLNucleus::getHitNucleonPosition(){
  if(nucleus_ && nucleon_index_ != -1){
    TVector3 v3_(999999.,999999.,999999.);
    v3_.SetXYZ(nucleus_->getStore()->getParticles().at(nucleon_index_)->getPosition().getX(),
	nucleus_->getStore()->getParticles().at(nucleon_index_)->getPosition().getY(),
	nucleus_->getStore()->getParticles().at(nucleon_index_)->getPosition().getZ());
    return v3_;
  }
  else 
    exit(1);
}
TVector3 INCLNucleus::getHitNucleonMomentum(){
  if(nucleus_ && nucleon_index_ != -1){
    TVector3 v3_(999999.,999999.,999999.);
    v3_.SetXYZ(nucleus_->getStore()->getParticles().at(nucleon_index_)->getMomentum().getX(),
	nucleus_->getStore()->getParticles().at(nucleon_index_)->getMomentum().getY(),
	nucleus_->getStore()->getParticles().at(nucleon_index_)->getMomentum().getZ());
    return v3_;
  }
  else 
    exit(1);

}
double INCLNucleus::getHitNucleonEnergy(){
  if(nucleus_ && nucleon_index_ != -1){
    return nucleus_->getStore()->getParticles().at(nucleon_index_)->getEnergy();
  }
  else 
    exit(1);
}

double INCLNucleus::getHitNucleonMass(){
  if(nucleus_ && nucleon_index_ != -1){
    return nucleus_->getStore()->getParticles().at(nucleon_index_)->getMass();
  }
  else 
    exit(1);
}

double INCLNucleus::getMass(){
  if(nucleus_ && nucleon_index_ != -1){
    return nucleus_->getMass();
  }
  else 
    exit(1);
}

G4INCL::Nucleus * INCLNucleus::getNuclues(){
  if(nucleus_ && nucleon_index_ != -1)
    return nucleus_;
  else 
    exit(1);
}

G4INCL::Particle * INCLNucleus::getHitParticle(){
    return hitNucleon_;
}

G4INCL::StandardPropagationModel * INCLNucleus::getPropagationModel(){
  return propagationModel_;
}

double INCLNucleus::getRemovalEnergy(){
  if(nucleus_ && nucleon_index_ != -1){
    double removal_energy = 0;
    double nucleon_mass = nucleus_->getStore()->getParticles().at(nucleon_index_)->getRealMass();
    double mag = nucleus_->getStore()->getParticles().at(nucleon_index_)->getMomentum().mag();
    removal_energy = TMath::Sqrt(mag*mag + nucleon_mass*nucleon_mass) - nucleus_->getStore()->getParticles().at(nucleon_index_)->getEnergy();
    return removal_energy;
  }
  else 
    exit(1);
}

#endif // __GENIE_INCL_ENABLED__

