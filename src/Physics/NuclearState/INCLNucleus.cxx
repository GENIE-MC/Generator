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
// GENIE headers
#include "Framework/Messenger/Messenger.h"
#include "Physics/NuclearState/INCLNucleus.h"
#include "Framework/Numerical/Spline.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"

// INCL headers
#include "G4INCLNuclearMassTable.hh"
#include "G4INCLGlobalInfo.hh"
#include "G4INCLPauliBlocking.hh"
#include "G4INCLCrossSections.hh"
#include "G4INCLPhaseSpaceGenerator.hh"
#include "G4INCLLogger.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLNuclearDensityFactory.hh"
#include "G4INCLINuclearPotential.hh"
#include "G4INCLCoulombDistortion.hh"
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
#include "G4INCLConfigEnums.hh"


using std::cout;
using std::endl;

using namespace genie;

//____________________________________________________________________________
INCLNucleus * INCLNucleus::fInstance = 0;

//____________________________________________________________________________
// Private constructor: direct instantiation is disallowed.
// Use the static factory method(s) above to create instances.
// This enforces controlled construction
// and prevents clients from bypassing it.
INCLNucleus::INCLNucleus():propagationModel_(0)
{
  //  this->Load();
  fInstance = 0;
  nucleon_index_ = -1;
  nucleus_ = nullptr;
  theConfig_ = nullptr;
  hitNucleon_ = nullptr;
  model_type_ = kNucmINCL; // default value, will be override in NucleusGenHybridStruck
}
//____________________________________________________________________________
INCLNucleus::~INCLNucleus()
{
  // ...
}
//____________________________________________________________________________
INCLNucleus * INCLNucleus::Instance()
{
  if(fInstance == 0) {
    LOG("INCLNucleus", pINFO) << "INCLNucleus initialization later";
    fInstance = new INCLNucleus;
    fInstance->theConfig_ = new G4INCL::Config();
  }
  return fInstance;
}

void INCLNucleus::configure(){

  theConfig_->init();
  theConfig_->setINCLXXDataFilePath(INCLXXDataFilePath_);
  theConfig_->setABLAXXDataFilePath(ablaxxDataFilePath_);
  theConfig_->setcABLA07DataFilePath(abla07DataFilePath_);
  theConfig_->setGEMINIXXDataFilePath(geminixxDataFilePath_);
  theConfig_->setDeExcitationType(deExcitationType_);
  theConfig_->setPotentialType(potentialType_);
  theConfig_->setPauliType(pauliType_);
  theConfig_->setPauliString(pauliString_);
  theConfig_->setLocalEnergyBBType(localEnergyTypeBB_);
  theConfig_->setLocalEnergyPiType(localEnergyTypePi_);
  theConfig_->setHadronizationTime(hadronizationTime_);
  theConfig_->setClusterAlgorithm(clusterAlgorithmType_);
  theConfig_->setClusterAlgorithmString(clusterAlgorithmString_);
  //theConfig_->setsrcPairConfig(true);
  std::cout << "DEBUG: " << __FILE__ << ":" << __LINE__ << " clasuter algorithm: " <<  clusterAlgorithmString_ << "  " << clusterAlgorithmType_ << std::endl;

  // TODO:  finish the configuration options
  //theConfig_->setRPCorrelationCoefficient(1.0); // Using r-p correlation without fuzzy
  //theConfig_->setLocalEnergyBBType(G4INCL::NeverLocalEnergy);
  //theConfig_->setLocalEnergyPiType(G4INCL::NeverLocalEnergy);

  // initialize INCL model
  G4INCL::Random::initialize(theConfig_);
  // Select the Pauli and CDPP blocking algorithms
  G4INCL::Pauli::initialize(theConfig_);
  // Set the cross-section set
  G4INCL::CrossSections::initialize(theConfig_);
  // Set the phase-space generator
  G4INCL::PhaseSpaceGenerator::initialize(theConfig_);
  // Initialize the INCL particle table:
  G4INCL::ParticleTable::initialize(theConfig_);
  // Select the Coulomb-distortion algorithm:
  G4INCL::CoulombDistortion::initialize(theConfig_);
  // Select the clustering algorithm:
  G4INCL::Clustering::initialize(theConfig_);
  // Initialize the value of cutNN in BinaryCollisionAvatar
  G4INCL::BinaryCollisionAvatar::setCutNN(theConfig_->getCutNN());
  // Initialize the value of strange cross section bias
  G4INCL::BinaryCollisionAvatar::setBias(theConfig_->getBias());


}

void INCLNucleus::initialize(const Target * tgt){

  // initialize the nucleus
  // If we need to re-initialize the nucleus, and nucleus_ != nullptr, delete it.
  if(nucleus_){
    delete nucleus_;
    nucleus_ = nullptr;
  }

  // initialize according process Event in INCL
  G4INCL::ParticleSpecies targetSpecies = G4INCL::ParticleSpecies(tgt->A(), tgt->Z());
  theConfig_->setTargetA(targetSpecies.theA);
  theConfig_->setTargetZ(targetSpecies.theZ);
  theConfig_->setTargetS(targetSpecies.theS);
  // define Nucleus and initialize it

  // ReInitialize the bias vector
  G4INCL::Particle::INCLBiasVector.clear();
  //Particle::INCLBiasVector.Clear();
  G4INCL::Particle::nextBiasedCollisionID = 0;


  // Initialise the maximum universe radius
  // INCL initialize universe radius according to particle species,
  // kenetic energy, and nucleus type. void INCL::initUniverseRadius(ParticleSpecies const &p, const double kineticEnergy, const int A, const int Z)
  // FIXME: I'm not sure if we need the interaction distance for neutrino scattering
  this->initUniverseRadius(targetSpecies.theA, targetSpecies.theZ);

  // FIXME the last two parameters need to be configed
  // theConfig_, G4INCL::ParticleTable::getMaximumNuclearRadius(G4INCL::Proton, targetSpecies.theA, targetSpecies.theZ)
  // G4INCL::NType
  nucleus_ = new G4INCL::Nucleus(targetSpecies.theA, targetSpecies.theZ, targetSpecies.theS, 
      theConfig_, maxUniverseRadius_, G4INCL::Def);
  nucleus_->getStore()->getBook().reset();
  nucleus_->initializeParticles();
  theDensity = nucleus_->getDensity();
  thePotential = nucleus_->getPotential();

  // maxInteractionDistance_ = 0;

  // set the min Remnant size to be 4
  // the min remnant is alpha particle
  // FIXME: it only works for nuclei with large A
  //
  // minRemnantSize_ = 4;

}

void INCLNucleus::reset(const Target * tgt){
  // initialize the target
  // reset the nucleus
  if(propagationModel_){
    delete propagationModel_;
  }
  this->initialize(tgt);

  // sample the index of nucleon hitted by lepton
  int nucleon_pdg = tgt->HitNucPdg();

  if(nucleon_pdg == kPdgClusterNN){
    // randomly pick a neutron and then use the closest neutron to form the cluster NN
    clusterNN_ = getNNCluster(kPdgNeutron, kPdgNeutron);
  }
  else if(nucleon_pdg == kPdgClusterNP){
    // randomly pick a neutron and then use the closest proton to form the cluster NP
    clusterNN_ = getNNCluster(kPdgNeutron, kPdgProton);
  }
  else if(nucleon_pdg == kPdgClusterPP){
    // randomly pick a proton and then use the closest proton to form the cluster PP
    clusterNN_ = getNNCluster(kPdgProton, kPdgProton);
  }
  else if(pdg::IsProton(nucleon_pdg) || pdg::IsNeutron(nucleon_pdg)){
    hitNucleon_ = this->getNucleon(nucleon_pdg);
  }
  else{
    // TODO:  might be 3p3h and NpNh channels
    LOG("INCLNucleus", pFATAL) << "Can't get a valid nucleon! " << nucleon_pdg;
    exit(1);
  }

  // Propagation model is responsible for finding avatars and
  // transporting the particles. In principle this step is "hidden"
  // behind an abstract interface and the rest of the system does not
  // care how the transportation and avatar finding is done. This
  // should allow us to "easily" experiment with different avatar
  // finding schemes and even to support things like curved
  // trajectories in the future.

  propagationModel_ = new G4INCL::StandardPropagationModel(
      theConfig_->getLocalEnergyBBType(),
      theConfig_->getLocalEnergyPiType(),
      theConfig_->getHadronizationTime());
  propagationModel_->setNucleus(nucleus_);
  // TODO: Using cascade action to do book keeping 
  //if(theConfig_->getCascadeActionType() == G4INCL::AvatarDumpActionType)
  //  cascadeAction_ = new G4INCL::AvatarDumpAction();
  //else
  //  cascadeAction_ = new G4INCL::CascadeAction();
}

TVector3 INCLNucleus::getHitNucleonPosition(){
  TVector3 v3(999999.,999999.,999999.);
  v3.SetXYZ(hitNucleon_->getPosition().getX(),
      hitNucleon_->getPosition().getY(),
      hitNucleon_->getPosition().getZ());
  return v3;
}

TVector3 INCLNucleus::getHitNucleonMomentum(){
  // INCL initial state;
  // we need to subtract the local energy from INCL nucleon before interaction
  double localEnergy = G4INCL::KinematicsUtils::getLocalEnergy(nucleus_, hitNucleon_);
  double oldEnergy = hitNucleon_->getEnergy();
  // subtract the local energy
  hitNucleon_->setEnergy(oldEnergy - localEnergy);
  hitNucleon_->adjustMomentumFromEnergy();
  TVector3 p3(999999.,999999.,999999.);
  p3.SetXYZ(hitNucleon_->getMomentum().getX(),
      hitNucleon_->getMomentum().getY(),
      hitNucleon_->getMomentum().getZ());
  // put it back to old energy
  hitNucleon_->setEnergy(oldEnergy);
  hitNucleon_->adjustMomentumFromEnergy();
  return p3;
}
double INCLNucleus::getHitNucleonEnergy(){
  double localEnergy = G4INCL::KinematicsUtils::getLocalEnergy(nucleus_, hitNucleon_);
  double oldEnergy = hitNucleon_->getEnergy();
  return (oldEnergy - localEnergy);
}

double INCLNucleus::getHitNucleonMass(){
  return hitNucleon_->getMass();
}

double INCLNucleus::getMass(){
  return nucleus_->getMass();
}

G4INCL::Nucleus * INCLNucleus::getNuclues(){
  return nucleus_;
}

G4INCL::Particle * INCLNucleus::getHitParticle(){
  return hitNucleon_;
}
std::shared_ptr<G4INCL::Cluster> INCLNucleus::getHitNNCluster(){
  return clusterNN_;
}


G4INCL::StandardPropagationModel * INCLNucleus::getPropagationModel(){
  return propagationModel_;
}

double INCLNucleus::getRemovalEnergy(){
  // FIXME: need to find the correct way for removal energy
  //   double removal_energy = 0;
  //   double nucleon_mass = hitNucleon_->getRealMass();
  //   double mag = hitNucleon_->getMomentum().mag();
  //   removal_energy = TMath::Sqrt(mag*mag + nucleon_mass*nucleon_mass) - hitNucleon_->getEnergy();
  // return removal_energy;
  //return hitNucleon_->getPotentialEnergy();
  //return nucleus_->getPotential()->getSeparationEnergy(hitNucleon_->getType());
  // FIXME: removal = potential - qvalue
  double removal = hitNucleon_->getPotentialEnergy() - hitNucleon_->getEmissionQValueCorrection(nucleus_->getA(),nucleus_->getZ(),nucleus_->getS());
  return removal;

}

void INCLNucleus::initUniverseRadius(const int A, const int Z){
  // This function is analogy to function in incl_physics/src/G4INCLCascade.cc
  // void INCL::initUniverseRadius(ParticleSpecies const &p, 
  //                const double kineticEnergy, const int A, 
  //                const int Z)
  double rMax = 0.0;
  // A should be large than 0
  // FIXME: 
  // 1. do we need to consider the isotopes?
  // 2. do we need to consider the extra-impact parameter for neutrino?
  // 	The xsec for neutrino-nucleus is ~10 fb
  // 	the xsec for hadron-nucleus is ~800 mb
  if(!(A > 0)) 
    throw std::runtime_error("Mass number A is not real!");
  const double pMaximumRadius = G4INCL::ParticleTable::getMaximumNuclearRadius(G4INCL::Proton,  A, Z);
  const double nMaximumRadius = G4INCL::ParticleTable::getMaximumNuclearRadius(G4INCL::Neutron, A, Z);
  const double maximumRadius = std::max(pMaximumRadius, nMaximumRadius);
  rMax = std::max(maximumRadius, rMax);
  maxUniverseRadius_ = rMax;
  //  LOG("INCLNucleus", pINFO) << "max Universe Radius : " << maxUniverseRadius_; 
}

std::shared_ptr<G4INCL::Cluster> INCLNucleus::getNNCluster(const int pdg1, const int pdg2){

  LOG("INCLNucleus", pINFO) << "get cluster";
  cluster_index1_ = -1;
  cluster_index2_ = -1;
  RandomGen * rnd = RandomGen::Instance();
  G4INCL::Particle *cluster_N1;
  G4INCL::Particle *cluster_N2;
  if(pdg::IsProton(pdg1))
    cluster_index1_ = rnd->RndGen().Integer(nucleus_->getZ());
  else if(pdg::IsNeutron(pdg1))
    cluster_index1_ = rnd->RndGen().Integer(nucleus_->getA() - nucleus_->getZ()) + nucleus_->getZ();
  else {
    LOG("INCLNucleus", pFATAL) << "Can't get a valid nucleon! " << pdg1;
    exit(1);
  }

  cluster_N1 = nucleus_->getStore()->getParticles().at(cluster_index1_);
  G4INCL::ParticleList const &particles = nucleus_->getStore()->getParticles();
  if(pdg::IsProton(pdg2)){
    double size = 10e16;
    for(G4INCL::ParticleIter i=particles.begin(), e=particles.end(); i!=e; ++i) {
      if((*i)->getID() == cluster_N1->getID()) continue;
      if((*i)->getType() != G4INCL::Proton) continue;
      double space    = ((*i)->getPosition() - cluster_N1->getPosition()).mag2();
      //double momentum = ((*i)->getMomentum() - cluster_N1->getMomentum()).mag2();
      double temp_size     = space;
      if(temp_size < size){ // TODO: maybe need to find a reasonable way to get the cluster
        size =  temp_size;
        cluster_N2 = (*i);
      }
    }
  }
  else if(pdg::IsNeutron(pdg2)){
    double size = 10e16;
    for(G4INCL::ParticleIter i=particles.begin(), e=particles.end(); i!=e; ++i) {
      if((*i)->getID() == cluster_N1->getID()) continue;
      if((*i)->getType() != G4INCL::Neutron) continue;
      double space    = ((*i)->getPosition() - cluster_N1->getPosition()).mag2();
      //double momentum = ((*i)->getMomentum() - cluster_N1->getMomentum()).mag2();
      double temp_size     = space;
      if(temp_size < size){
        size =  temp_size;
        cluster_N2 = (*i);
      }
    }
  }
  else{
    LOG("INCLNucleus", pFATAL) << "Can't get a valid nucleon! " << pdg2;
    exit(1);
  }
  G4INCL::ParticleList selectedParticles;
  selectedParticles.push_back(cluster_N1);
  selectedParticles.push_back(cluster_N2);
  return std::make_shared<G4INCL::Cluster>(selectedParticles.begin(), selectedParticles.end());

}

G4INCL::Particle * INCLNucleus::getNucleon(const int pdg){

  RandomGen * rnd = RandomGen::Instance();
  nucleon_index_ = -1;
  if(pdg::IsProton(pdg))
    nucleon_index_ = rnd->RndGen().Integer(nucleus_->getZ());
  else if(pdg::IsNeutron(pdg))
    nucleon_index_ = rnd->RndGen().Integer(nucleus_->getA() - nucleus_->getZ()) + nucleus_->getZ();
  else{
    LOG("INCLNucleus", pFATAL) << "Can't get a valid nucleon!";
    exit(1);
  }
  // set the index of nucleon hitted by lepton before initialize a nucleus
  //    nucleus_->setLeptonScatteringDensity(theDensityForLepton);
  //    nucleus_->setLeptonHitNucleonIndex(nucleon_index_);
  // reset the hit nucleon
  return nucleus_->getStore()->getParticles().at(nucleon_index_);
}

bool INCLNucleus::isRPValid(double r, double p){
  (void)r;
  double theFermiMomentum = thePotential->getFermiMomentum(hitNucleon_->getType());
  // local energy
  double locE = G4INCL::KinematicsUtils::getLocalEnergy(nucleus_, hitNucleon_);

  double theFermiEnergy = std::sqrt(theFermiMomentum*theFermiMomentum + hitNucleon_->getMass()*hitNucleon_->getMass());
  double MaxMomAtR = std::sqrt((theFermiEnergy - locE) *  (theFermiEnergy - locE) - hitNucleon_->getMass()*hitNucleon_->getMass());

  return (std::abs(p) < MaxMomAtR);

}

TVector3 INCLNucleus::ResamplingVertex(const int pdg){
  G4INCL::ParticleType t = G4INCL::Proton;
  if(pdg::IsProton(pdg)){
    t = G4INCL::Proton;
  }else if(pdg::IsNeutron(pdg)){
    t = G4INCL::Neutron;
  }
  double rpCorrelationCoefficient[G4INCL::UnknownParticle];
  std::fill(rpCorrelationCoefficient, rpCorrelationCoefficient + G4INCL::UnknownParticle, 1.);
  rpCorrelationCoefficient[G4INCL::Proton] = G4INCL::ParticleTable::getRPCorrelationCoefficient(G4INCL::Proton);
  rpCorrelationCoefficient[G4INCL::Neutron] = G4INCL::ParticleTable::getRPCorrelationCoefficient(G4INCL::Neutron);
  rpCorrelationCoefficient[G4INCL::Lambda] = G4INCL::ParticleTable::getRPCorrelationCoefficient(G4INCL::Lambda);
  assert(theDensity && thePotential);
  std::pair<double,double> ranNumbers = G4INCL::Random::correlatedUniform(rpCorrelationCoefficient[t]);
  const double x = G4INCL::Math::pow13(ranNumbers.first);
  //const double y = G4INCL::Math::pow13(ranNumbers.second);
  //const double theFermiMomentum = thePotential->getFermiMomentum(t);
  //const G4INCL::ThreeVector momentumVector = G4INCL::Random::normVector(y*theFermiMomentum);
  const double reflectionRadius = theDensity->getMaxRFromP(t, x);
  const G4INCL::ThreeVector positionVector = G4INCL::Random::sphereVector(reflectionRadius);
  const G4INCL::ThreeVector nucleus_position = nucleus_->getPosition();
  const G4INCL::ThreeVector new_position = positionVector + nucleus_position;
  return TVector3(new_position.getX(), new_position.getY(), new_position.getZ());
}

void INCLNucleus::ResamplingHitNucleon(){
  // local energy
  double locE = G4INCL::KinematicsUtils::getLocalEnergy(nucleus_, hitNucleon_);

  int iteration_count = 0;
  const double theFermiMomentum = thePotential->getFermiMomentum(hitNucleon_->getType());
  while(true){
    const G4INCL::ThreeVector momentumVector = G4INCL::Random::sphereVector(theFermiMomentum);
    const double momentumAbs = momentumVector.mag();
    hitNucleon_->setMomentum(momentumVector);
    hitNucleon_->setUncorrelatedMomentum(momentumAbs);
    hitNucleon_->adjustEnergyFromMomentum();
    double KE = hitNucleon_->getEnergy() - hitNucleon_->getMass();
    iteration_count++;
    if(KE > locE){
      break;
    }
    if(iteration_count > 10000){
      LOG("INCLNucleus", pFATAL) << "Resamping the momentum of struck nucleon more than 10000 times!";
      exit(1);
    }
  }
}
/*
void INCLNucleus::setHitParticle(const int pdg, TVector3 &posi){
  LOG("INCLNucleus", pINFO) << "set hit nucleon according to position";
  nucleon_index_ = -1;
  G4INCL::ThreeVector hitposi(posi.X(), posi.Y(), posi.Z());

  G4INCL::ParticleList const &particles = nucleus_->getStore()->getParticles();
  if(pdg::IsProton(pdg)){
    double size = 1e16;
    for(G4INCL::ParticleIter i=particles.begin(), e=particles.end(); i!=e; ++i) {
      if((*i)->getType() != G4INCL::Proton) continue;
      double space    = ((*i)->getPosition() - hitposi).mag2();
      // we don't consider the momentum now. 
      // double momentum = ...... ; 
      double temp_size     = space;
      if(temp_size < size){ // TODO: maybe need to find a reasonable way to get the cluster
        size =  temp_size;
        hitNucleon_ = (*i);
      }
    }
  }
  else if(pdg::IsNeutron(pdg)){
    double size = 1e16;
    for(G4INCL::ParticleIter i=particles.begin(), e=particles.end(); i!=e; ++i) {
      if((*i)->getType() != G4INCL::Neutron) continue;
      double space    = ((*i)->getPosition() - hitposi).mag2();
      // we don't consider the momentum now. 
      // double momentum = ...... ; 
      double temp_size     = space;
      if(temp_size < size){
        size =  temp_size;
        hitNucleon_ = (*i);
      }
    }
  }
  else{
    LOG("INCLNucleus", pFATAL) << "Can't get a valid nucleon! " << pdg;
    exit(1);
  }
  propagationModel_->setNucleus(nucleus_);

}
*/
void INCLNucleus::setHitNNCluster(const int pdg1, const int pdg2, TVector3 &posi){
  // FIXME: find the closest nucleon pair might basing the nuclear density. 
  // We will randomly drop two nucleon according to pdg number.
  LOG("INCLNucleus", pINFO) << "get cluster";
  //RandomGen * rnd = RandomGen::Instance();
  G4INCL::ThreeVector hitposi(posi.X(), posi.Y(), posi.Z());

  int cluster_pdg[2];
  cluster_pdg[0] = pdg1;
  cluster_pdg[1] = pdg2;
  G4INCL::Particle *cluster_N[2];
  for(int idx = 0; idx < 2; idx++){
    int pdg_ = cluster_pdg[idx];
    G4INCL::ParticleList const &particles = nucleus_->getStore()->getParticles();
    if(pdg::IsProton(pdg_)){
      double size = 1e16;
      for(G4INCL::ParticleIter i=particles.begin(), e=particles.end(); i!=e; ++i) {
        if(idx == 1){ //  the second nucleon and the first nucleon should have different ID.
          if((*i)->getID() == cluster_N[0]->getID()) continue;
        }
        if((*i)->getType() != G4INCL::Proton) continue;
        double space    = ((*i)->getPosition() - hitposi).mag2();
        // double momentum = ((*i)->getMomentum() - cluster_N1->getMomentum()).mag2();
        double temp_size     = space;
        if(temp_size < size){ // TODO: maybe need to find a reasonable way to get the cluster
          size =  temp_size;
          cluster_N[idx] = (*i);
        }
      }
    }
    else if(pdg::IsNeutron(pdg_)){
      double size = 1e16;
      for(G4INCL::ParticleIter i=particles.begin(), e=particles.end(); i!=e; ++i) {
        if(idx == 1){ //  the second nucleon and the first nucleon should have different ID.
          if((*i)->getID() == cluster_N[0]->getID()) continue;
        }
        if((*i)->getType() != G4INCL::Neutron) continue;
        double space    = ((*i)->getPosition() - hitposi).mag2();
        //double momentum = ((*i)->getMomentum() - cluster_N1->getMomentum()).mag2();
        double temp_size     = space;
        if(temp_size < size){
          size =  temp_size;
          cluster_N[idx] = (*i);
        }
      }
    }
    else{
      LOG("INCLNucleus", pFATAL) << "Can't get a valid nucleon! " << pdg2;
      exit(1);
    }
  }

  G4INCL::ParticleList selectedParticles;
  selectedParticles.push_back(cluster_N[0]);
  selectedParticles.push_back(cluster_N[1]);
  clusterNN_ = std::make_shared<G4INCL::Cluster>(selectedParticles.begin(), selectedParticles.end());
}

#endif // __GENIE_INCL_ENABLED__

