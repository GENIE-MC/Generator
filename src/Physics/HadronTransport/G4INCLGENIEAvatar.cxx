#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_INCL_ENABLED__

#include "G4INCLGENIEAvatar.h"
#include "G4INCLGENIEQELChannel.h"
#include "G4INCLGENIEMECChannel.h"
#include "G4INCLGENIEDISChannel.h"
#include "G4INCLGENIERESChannel.h"
#include "G4INCLPauliBlocking.hh"
#include <sstream>
#include <string>
#include <cassert>


namespace G4INCL {

  GENIEAvatar::GENIEAvatar(double time, Particle *p, Nucleus *n, std::vector<GENIEParticleRecord> *eventRecord, bool ishybrid)
    : InteractionAvatar(time, n, p),
    particle1(p), theNucleus(n),
    genie_evtrec(eventRecord), fHybridModel(ishybrid)
  {
    delta_Z = 0;
  }


  GENIEAvatar::GENIEAvatar(double time, Cluster *p, Nucleus *n, std::vector<GENIEParticleRecord> *eventRecord, bool ishybrid)
    : InteractionAvatar(time, n, nullptr),
    theNucleus(n), cluster(p),
    genie_evtrec(eventRecord), fHybridModel(ishybrid)
  {
    delta_Z = 0;
    fHybridModel = false;
  }

  GENIEAvatar::~GENIEAvatar() {

  }

  G4INCL::IChannel* GENIEAvatar::getChannel() {

    // using genie event record to fill final states of INCL
    // the initial states is hit nucleons, (for MEC channel, initial states
    // is NN cluster)
    // the final states should be hadron in nucleus
    switch((*genie_evtrec)[0].ScatteringType()){
      case 1: return new GENIEQELChannel(particle1, theNucleus, genie_evtrec);
      case 4: return new GENIERESChannel(particle1, theNucleus, genie_evtrec);
      case 3: return new GENIEDISChannel(particle1, theNucleus, genie_evtrec);
      case 10: return new GENIEMECChannel(cluster, theNucleus, genie_evtrec);
      default: return NULL;
    }

  }


  void GENIEAvatar::preInteraction() {


    if((*genie_evtrec)[0].ScatteringType() != genie::kScMEC){
      int index = 0;
      double lepton_initial_energy = 0;
      ThreeVector leptonInitialMom;
      std::vector<GENIEParticleRecord>::iterator ip;
      for(ip = genie_evtrec->begin(); ip != genie_evtrec->end(); ip++){
        if(ip->RecordCode() == kProbe){
          lepton_initial_energy = std::sqrt(ip->P3().mag2() + ip->Mass()*ip->Mass());
          leptonInitialMom = ip->P3();
          delta_Z += ip->Charge();
        }
        else if(ip->RecordCode() == kFinalStateLepton){
          leptonE = std::sqrt(ip->P3().mag2() + ip->Mass()*ip->Mass());
          leptonMom = ip->P3();
          delta_Z -= ip->Charge();
        }
        else if(ip->RecordCode() == kHitNucleon){
          if(!fHybridModel){
            ThreeVector n_mom = particle1->getMomentum();
            ip->setMomentum(n_mom);
            ip->setMass(particle1->getMass());
          }
        }
        index++;
      }
      oldTotalEnergy = lepton_initial_energy + particle1->getEnergy() - particle1->getPotentialEnergy();

      // transfrom the target nucleon to local energy frame
      KinematicsUtils::transformToLocalEnergyFrame(theNucleus, particle1);

      // make boost vector
      boostVector = (leptonInitialMom + particle1->getMomentum())/(lepton_initial_energy + particle1->getEnergy());
    }
    else{
      int index = 0;
      double lepton_initial_energy = 0;
      ThreeVector leptonInitialMom;
      std::vector<GENIEParticleRecord>::iterator ip;
      for(ip = genie_evtrec->begin(); ip != genie_evtrec->end(); ip++){
        if(ip->RecordCode() == kProbe){
          lepton_initial_energy = std::sqrt(ip->P3().mag2() + ip->Mass()*ip->Mass());
          leptonInitialMom = ip->P3();
          delta_Z += ip->Charge();
        }
        else if(ip->RecordCode() == kFinalStateLepton){
          leptonE = std::sqrt(ip->P3().mag2() + ip->Mass()*ip->Mass());
          leptonMom = ip->P3();
          delta_Z -= ip->Charge();
        }
        else if(ip->RecordCode() == kHitNucleon){
        }
        index++;
      }

      G4INCL::ParticleList particles = cluster->getParticleList();
      oldTotalEnergy = lepton_initial_energy;

      ThreeVector local_mom = leptonInitialMom;
      double local_energy = lepton_initial_energy;
      for(G4INCL::ParticleIter i=particles.begin(), e=particles.end(); i!=e; ++i) {
        oldTotalEnergy += (*i)->getEnergy() - (*i)->getPotentialEnergy();
        KinematicsUtils::transformToLocalEnergyFrame(theNucleus, (*i));
        local_mom += (*i)->getMomentum();
        local_energy += (*i)->getEnergy();
      }
      boostVector = local_mom / local_energy;
    }
  }

  void GENIEAvatar::postInteractionHybridModel(FinalState *fs){
    // Make sure we have at least two particles in the final state
    // Removed because of neutral kaon decay

    // InteractionAvatar::postInteraction(fs);
    // TODO: modify the InteractionAvatar::postInteraction(FinalState *fs)
    // in here to handle the final state from primary neutrino interaction
    // QEL RES MEC DIS
    modified = fs->getModifiedParticles();
    created = fs->getCreatedParticles();
    Destroyed = fs->getDestroyedParticles();
    modifiedAndCreated = modified;     
    modifiedAndCreated.insert(modifiedAndCreated.end(), created.begin(), created.end());

    // If there is no Nucleus, just return
    if(!theNucleus) return;
    // put hadrons into potential and make them off-shell
    // FIXME: need to add the Q-value for hadron from primary vertex.

    for(ParticleIter i=modifiedAndCreated.begin(), e=modifiedAndCreated.end(); i!=e; ++i ){
      double Qval = (*i)->getEmissionQValueCorrection(theNucleus->getA(),theNucleus->getZ(),theNucleus->getS());
      bool success = this->putIntoPotential(Qval, (*i));
      (*i)->rpCorrelate();
      if(!success){
        fs->reset();
        fs->makeNoEnergyConservation();
        fs->setTotalEnergyBeforeInteraction(0.0);
        return; // Interaction is blocked. Return an empty final state.
      }
    }


    /// update the event record after call enforceEnergyConservation;
    int index = 0;
    std::vector<GENIEParticleRecord>::iterator ip;
    ParticleIter imc = modifiedAndCreated.begin();
    for(ip = genie_evtrec->begin(); ip != genie_evtrec->end(); ip++){
      if(ip->Status() == genie::kIStHadronInTheNucleus || ip->Status() == genie::kIStPreDecayResonantState){
        ThreeVector p_mom = (*imc)->getMomentum();
        ip->setMomentum(p_mom);
        ip->setMass((*imc)->getMass());
        imc++;
      }
      index++;
    }

    // Mark pions and kaons that have been created outside their well (we will force them
    // to be emitted later).
    for(ParticleIter i=created.begin(), e=created.end(); i!=e; ++i ){
      if(((*i)->isPion() || (*i)->isKaon() || (*i)->isAntiKaon()) && (*i)->getPosition().mag() > theNucleus->getSurfaceRadius(*i)) {
        (*i)->makeParticipant();
        (*i)->setOutOfWell();
        fs->addOutgoingParticle(*i);
        INCL_DEBUG("Pion was created outside its potential well." << '\n'
            << (*i)->print());
      }
    }


    // Test Pauli blocking
    bool isBlocked = Pauli::isBlocked(modifiedAndCreated, theNucleus);

    if(isBlocked) {
      // Restore the state of the initial particles
      //restoreParticles();

      // Delete newly created particles
      for(ParticleIter i=created.begin(), e=created.end(); i!=e; ++i )
        delete *i;

      fs->reset();
      fs->makePauliBlocked();
      fs->setTotalEnergyBeforeInteraction(0.0);

      return; // Interaction is blocked. Return an empty final state.
    }



    // Test CDPP blocking
    bool isCDPPBlocked = Pauli::isCDPPBlocked(created, theNucleus);
    if(isCDPPBlocked) {

      // Restore the state of the initial particles
      //restoreParticles();

      // Delete newly created particles
      for(ParticleIter i=created.begin(), e=created.end(); i!=e; ++i ){
        delete *i;
      }

      fs->reset();
      fs->makePauliBlocked();
      fs->setTotalEnergyBeforeInteraction(0.0);

      return; // Interaction is blocked. Return an empty final state.
    }


    // If all went well, try to bring particles inside the nucleus...
    for(ParticleIter i=modifiedAndCreated.begin(), e=modifiedAndCreated.end(); i!=e; ++i ){
      // ...except for pions beyond their surface radius.
      if((*i)->isOutOfWell()) continue;

      //const bool successBringParticlesInside = InteractionAvatar::bringParticleInside(*i);
    }

    // Collision accepted!
    // Biasing of the final state
    std::vector<int> newBiasCollisionVector;
    newBiasCollisionVector = ModifiedAndDestroyed.getParticleListBiasVector();
    if(std::fabs(weight-1.) > 1E-6){
      newBiasCollisionVector.push_back(Particle::nextBiasedCollisionID);
      Particle::FillINCLBiasVector(1./weight);
      weight = 1.; // useless?
    }

    for(ParticleIter i=modifiedAndCreated.begin(), e=modifiedAndCreated.end(); i!=e; ++i ) {
      (*i)->setBiasCollisionVector(newBiasCollisionVector);
      if(!(*i)->isOutOfWell()) {
        // Decide if the particle should be made into a spectator
        // (Back to spectator)
        bool goesBackToSpectator = false;
        if((*i)->isNucleon() && theNucleus->getStore()->getConfig()->getBackToSpectator()) {
          double threshold = (*i)->getPotentialEnergy();
          if((*i)->getType()==Proton)
            threshold += Math::twoThirds*theNucleus->getTransmissionBarrier(*i);
          if((*i)->getKineticEnergy() < threshold)
            goesBackToSpectator = true;
        }
        // Thaw the particle propagation
        (*i)->thawPropagation();

        // Increment or decrement the participant counters
        if(goesBackToSpectator) {
          INCL_DEBUG("The following particle goes back to spectator:" << '\n'
              << (*i)->print() << '\n');
          if(!(*i)->isTargetSpectator()) {
            theNucleus->getStore()->getBook().decrementCascading();
          }
          (*i)->makeTargetSpectator();
        } else {
          if((*i)->isTargetSpectator()) {
            theNucleus->getStore()->getBook().incrementCascading();
          }
          (*i)->makeParticipant();
        }
      }
    }
    ParticleList destroyed = fs->getDestroyedParticles();
    for(ParticleIter i=destroyed.begin(), e=destroyed.end(); i!=e; ++i )
      if(!(*i)->isTargetSpectator())
        theNucleus->getStore()->getBook().decrementCascading();

    theNucleus->setZ(theNucleus->getZ() + delta_Z);
    theNucleus->getStore()->getBook().incrementAcceptedCollisions();

    return;

  }

  void GENIEAvatar::postInteraction(FinalState *fs) {
    // Make sure we have at least two particles in the final state
    // Removed because of neutral kaon decay

    // InteractionAvatar::postInteraction(fs);
    // TODO: modify the InteractionAvatar::postInteraction(FinalState *fs)
    // in here to handle the final state from primary neutrino interaction
    // QEL RES MEC DIS
    modified = fs->getModifiedParticles();
    created = fs->getCreatedParticles();
    Destroyed = fs->getDestroyedParticles();
    modifiedAndCreated = modified;     
    modifiedAndCreated.insert(modifiedAndCreated.end(), created.begin(), created.end());

    // If there is no Nucleus, just return
    if(!theNucleus) return;

    if(fHybridModel){
      return this->postInteractionHybridModel(fs);
    }


    // Mark pions and kaons that have been created outside their well (we will force them
    // to be emitted later).
    for(ParticleIter i=created.begin(), e=created.end(); i!=e; ++i ){
      if(((*i)->isPion() || (*i)->isKaon() || (*i)->isAntiKaon()) && (*i)->getPosition().mag() > theNucleus->getSurfaceRadius(*i)) {
        (*i)->makeParticipant();
        (*i)->setOutOfWell();
        fs->addOutgoingParticle(*i);
        INCL_DEBUG("Pion was created outside its potential well." << '\n'
            << (*i)->print());
      }
    }


    // using INCL nuclear model for primary vertex
    bool success = enforceEnergyConservation(fs);
    if(!success){
      std::cout << "enforceEnergyConservation fs wrong!" << std::endl;
      fs->reset();
      fs->makeNoEnergyConservation();
      fs->setTotalEnergyBeforeInteraction(0.0);
      return; // Interaction is blocked. Return an empty final state.
    }

    // Test Pauli blocking
    bool isBlocked = Pauli::isBlocked(modifiedAndCreated, theNucleus);

    if(isBlocked) {
      // Restore the state of the initial particles
      //restoreParticles();

      // Delete newly created particles
      for(ParticleIter i=created.begin(), e=created.end(); i!=e; ++i )
        delete *i;

      fs->reset();
      fs->makePauliBlocked();
      fs->setTotalEnergyBeforeInteraction(0.0);

      return; // Interaction is blocked. Return an empty final state.
    }



    // Test CDPP blocking
    bool isCDPPBlocked = Pauli::isCDPPBlocked(created, theNucleus);
    if(isCDPPBlocked) {

      // Restore the state of the initial particles
      //restoreParticles();

      // Delete newly created particles
      for(ParticleIter i=created.begin(), e=created.end(); i!=e; ++i ){
        delete *i;
      }

      fs->reset();
      fs->makePauliBlocked();
      fs->setTotalEnergyBeforeInteraction(0.0);

      return; // Interaction is blocked. Return an empty final state.
    }


    // If all went well, try to bring particles inside the nucleus...
    for(ParticleIter i=modifiedAndCreated.begin(), e=modifiedAndCreated.end(); i!=e; ++i ){
      // ...except for pions beyond their surface radius.
      if((*i)->isOutOfWell()) continue;

      //const bool successBringParticlesInside = InteractionAvatar::bringParticleInside(*i);
    }


    /// update the event record after call enforceEnergyConservation;
    int index = 0;
    std::vector<GENIEParticleRecord>::iterator ip;
    ParticleIter imc = modifiedAndCreated.begin();
    for(ip = genie_evtrec->begin(); ip != genie_evtrec->end(); ip++){
      if(ip->RecordCode() == kFinalStateLepton){
        ip->setMomentum(leptonMom);
        ip->setMass(std::sqrt(std::max(leptonE*leptonE - leptonMom.mag2(), 0.)));
      }
      else if(ip->Status() == genie::kIStHadronInTheNucleus || ip->Status() == genie::kIStPreDecayResonantState){
        ThreeVector p_mom = (*imc)->getMomentum();
        ip->setMomentum(p_mom);
        ip->setMass((*imc)->getMass());
        imc++;
      }
      index++;
    }

    // Collision accepted!
    // Biasing of the final state
    std::vector<int> newBiasCollisionVector;
    newBiasCollisionVector = ModifiedAndDestroyed.getParticleListBiasVector();
    if(std::fabs(weight-1.) > 1E-6){
      newBiasCollisionVector.push_back(Particle::nextBiasedCollisionID);
      Particle::FillINCLBiasVector(1./weight);
      weight = 1.; // useless?
    }
    for(ParticleIter i=modifiedAndCreated.begin(), e=modifiedAndCreated.end(); i!=e; ++i ) {
      (*i)->setBiasCollisionVector(newBiasCollisionVector);
      if(!(*i)->isOutOfWell()) {
        // Decide if the particle should be made into a spectator
        // (Back to spectator)
        bool goesBackToSpectator = false;
        if((*i)->isNucleon() && theNucleus->getStore()->getConfig()->getBackToSpectator()) {
          double threshold = (*i)->getPotentialEnergy();
          if((*i)->getType()==Proton)
            threshold += Math::twoThirds*theNucleus->getTransmissionBarrier(*i);
          if((*i)->getKineticEnergy() < threshold)
            goesBackToSpectator = true;
        }
        // Thaw the particle propagation
        (*i)->thawPropagation();

        // Increment or decrement the participant counters
        if(goesBackToSpectator) {
          INCL_DEBUG("The following particle goes back to spectator:" << '\n'
              << (*i)->print() << '\n');
          if(!(*i)->isTargetSpectator()) {
            theNucleus->getStore()->getBook().decrementCascading();
          }
          (*i)->makeTargetSpectator();
        } else {
          if((*i)->isTargetSpectator()) {
            theNucleus->getStore()->getBook().incrementCascading();
          }
          (*i)->makeParticipant();
        }
      }
    }
    ParticleList destroyed = fs->getDestroyedParticles();
    for(ParticleIter i=destroyed.begin(), e=destroyed.end(); i!=e; ++i )
      if(!(*i)->isTargetSpectator())
        theNucleus->getStore()->getBook().decrementCascading();

    theNucleus->setZ(theNucleus->getZ() + delta_Z);
    theNucleus->getStore()->getBook().incrementAcceptedCollisions();
    return;
  }

  bool GENIEAvatar::enforceEnergyConservation(FinalState * const fs){
    (void) fs;
    // Set up the violationE calculation
    violationEFunctor = new ViolationLeptonEMomentumFunctor(theNucleus, modifiedAndCreated, leptonE, leptonMom, boostVector, oldTotalEnergy, true);
    const RootFinder::Solution theSolution = RootFinder::solve(violationEFunctor, 1.0);
    if(theSolution.success) { // Apply the solution
      (*violationEFunctor)(theSolution.x);
    } else if(theNucleus){
      theNucleus->getStore()->getBook().incrementEnergyViolationInteraction();
    }
    delete violationEFunctor;
    violationEFunctor = NULL;
    return theSolution.success;

  }

  /* ***                                                      ***
   * *** InteractionAvatar::ViolationEMomentumFunctor methods ***
   * ***                                                      ***/

  GENIEAvatar::ViolationLeptonEMomentumFunctor::ViolationLeptonEMomentumFunctor(Nucleus * const nucleus, ParticleList const &modAndCre,
      double &leptonE, ThreeVector &leptonMom, ThreeVector const &boost,
      const double totalEnergyBeforeInteraction, const bool localE) :
    RootFunctor(0., 1E6),
    finalParticles(modAndCre),
    leptonEnergy(leptonE),
    leptonMomentum(leptonMom),
    initialEnergy(totalEnergyBeforeInteraction),
    theNucleus(nucleus),
    boostVector(boost),
    shouldUseLocalEnergy(localE)
  {
    // Store the particle momenta (necessary for the calls to
    // scaleParticleMomenta() to work)
    lepton = new Lepton(leptonMomentum, leptonEnergy);
    lepton->boost(boostVector);
    leptonMomentumCM = lepton->getMomentum();
    for(ParticleIter i=finalParticles.begin(), e=finalParticles.end(); i!=e; ++i) {
      (*i)->boost(boostVector);
      particleMomenta.push_back((*i)->getMomentum());
    }
  }

  GENIEAvatar::ViolationLeptonEMomentumFunctor::~ViolationLeptonEMomentumFunctor() {
    particleMomenta.clear();
    delete lepton;
  }


  double GENIEAvatar::ViolationLeptonEMomentumFunctor::operator()(const double alpha) const {
    //LOG("INCLCascadeIntranuke", pNOTICE) << "alpha : " << alpha;
    scaleParticleMomenta(alpha);

    double deltaE = 0.0;
    for(ParticleIter i=finalParticles.begin(), e=finalParticles.end(); i!=e; ++i){
      deltaE += (*i)->getEnergy() - (*i)->getPotentialEnergy();
    }
    deltaE += leptonEnergy;
    deltaE -= initialEnergy;
    return deltaE;
  }

  void GENIEAvatar::ViolationLeptonEMomentumFunctor::scaleParticleMomenta(const double alpha) const {

    lepton->setMomentum(leptonMomentumCM*alpha);
    lepton->adjustEnergyFromMomentum();
    lepton->boost(-boostVector);
    leptonEnergy = lepton->getEnergy();
    leptonMomentum = lepton->getMomentum();

    std::vector<ThreeVector>::const_iterator iP = particleMomenta.begin();
    for(ParticleIter i=finalParticles.begin(), e=finalParticles.end(); i!=e; ++i, ++iP) {
      (*i)->setMomentum((*iP)*alpha);
      (*i)->adjustEnergyFromMomentum();
      (*i)->rpCorrelate();
      (*i)->boost(-boostVector);
      if(theNucleus){
        theNucleus->updatePotentialEnergy(*i);
      } else {
        (*i)->setPotentialEnergy(0.);
      }

      //jcd      if(shouldUseLocalEnergy && !(*i)->isPion())  // This translates AECSVT's loops 1, 3 and 4
      if(shouldUseLocalEnergy && !(*i)->isPion() && !(*i)->isEta() && !(*i)->isOmega() &&
          !(*i)->isKaon() && !(*i)->isAntiKaon()  && !(*i)->isSigma() && !(*i)->isPhoton() && !(*i)->isLambda() && !(*i)->isAntiNucleon()) { // This translates AECSVT's loops 1, 3 and 4
        assert(theNucleus); // Local energy without a nucleus doesn't make sense
        const double energy = (*i)->getEnergy(); // Store the energy of the particle
        double locE = KinematicsUtils::getLocalEnergy(theNucleus, *i); // Initial value of local energy
        double locEOld;
        double deltaLocE = InteractionAvatar::locEAccuracy + 1E3;
        for(int iterLocE=0;
            deltaLocE>InteractionAvatar::locEAccuracy && iterLocE<InteractionAvatar::maxIterLocE;
            ++iterLocE) {
          locEOld = locE;
          (*i)->setEnergy(energy + locE); // Update the energy of the particle...
          (*i)->adjustMomentumFromEnergy();
          theNucleus->updatePotentialEnergy(*i); // ...update its potential energy...
          locE = KinematicsUtils::getLocalEnergy(theNucleus, *i); // ...and recompute locE.
          deltaLocE = std::abs(locE-locEOld);
        }
      }

      //jlrs  For lambdas and nuclei with masses higher than 19 also local energy
      if(shouldUseLocalEnergy && (*i)->isLambda() && theNucleus->getA()>19) {
        assert(theNucleus); // Local energy without a nucleus doesn't make sense
        const double energy = (*i)->getEnergy(); // Store the energy of the particle
        double locE = KinematicsUtils::getLocalEnergy(theNucleus, *i); // Initial value of local energy
        double locEOld;
        double deltaLocE = InteractionAvatar::locEAccuracy + 1E3;
        for(int iterLocE=0;
            deltaLocE>InteractionAvatar::locEAccuracy && iterLocE<InteractionAvatar::maxIterLocE;
            ++iterLocE) {
          locEOld = locE;
          (*i)->setEnergy(energy + locE); // Update the energy of the particle...
          (*i)->adjustMomentumFromEnergy();
          theNucleus->updatePotentialEnergy(*i); // ...update its potential energy...
          locE = KinematicsUtils::getLocalEnergy(theNucleus, *i); // ...and recompute locE.
          deltaLocE = std::abs(locE-locEOld);
        }
      }
    }
  }

  void GENIEAvatar::ViolationLeptonEMomentumFunctor::cleanUp(const bool success) const {
    if(!success)
      scaleParticleMomenta(1.);
  }


  std::string GENIEAvatar::dump() const {
    std::stringstream ss;
    ss << "(avatar " << theTime << " 'decay" << '\n'
      << "(list " << '\n'
      << particle1->dump()
      << "))" << '\n';
    return ss.str();
  }

  bool GENIEAvatar::putIntoPotential(const double theQValueCorrection, Particle * const theParticle){

    // TODO {this is the place to add refraction}
    theParticle->setINCLMass(); // Will automatically put the particle on shell

    // Add the nuclear potential to the kinetic energy when entering the
    // nucleus

    //double vv = 0.0;
    //double energyInside = std::max(theParticle->getMass(), theParticle->getEnergy() + vv - theQValueCorrection);
    //theParticle->setEnergy(energyInside);
    //theParticle->setPotentialEnergy(vv);
    //theParticle->adjustMomentumFromEnergy();
    //return true;


    class IncomingEFunctor : public RootFunctor {
      public:
        IncomingEFunctor(Particle * const p, Nucleus const * const n, const double correction) :
          RootFunctor(0., 1E6),
          theParticle(p),
          thePotential(n->getPotential()),
          theEnergy(theParticle->getEnergy()),
          theMass(theParticle->getMass()),
          theQValueCorrection(correction),
          theMomentumDirection(theParticle->getMomentum())
      {

      }
        ~IncomingEFunctor() {}
        double operator()(const double v) const {
          double energyInside = std::max(theMass, theEnergy + v - theQValueCorrection);
          theParticle->setEnergy(energyInside);
          theParticle->setPotentialEnergy(v);
          theParticle->setMomentum(theMomentumDirection); // keep the same direction
          theParticle->adjustMomentumFromEnergy();
          return v - thePotential->computePotentialEnergy(theParticle);
        }
        void cleanUp(const bool /*success*/) const {}
      private:
        Particle *theParticle;
        NuclearPotential::INuclearPotential const *thePotential;
        const double theEnergy;
        const double theMass;
        const double theQValueCorrection;
        const ThreeVector theMomentumDirection;
    } theIncomingEFunctor(theParticle,theNucleus,theQValueCorrection);

    double v = theNucleus->getPotential()->computePotentialEnergy(theParticle);
    if(theParticle->getKineticEnergy()+v-theQValueCorrection<0.) { // Particle entering below 0. Die gracefully
      INCL_DEBUG("Particle " << theParticle->getID() << " is trying to enter below 0" << '\n');
      return false;
    }

    const RootFinder::Solution theSolution = RootFinder::solve(&theIncomingEFunctor, v);
    if(theSolution.success) { // Apply the solution
      theIncomingEFunctor(theSolution.x);
      INCL_DEBUG("Particle successfully entered:\n" << theParticle->print() << '\n');
    } else {
      INCL_WARN("Couldn't compute the potential for incoming particle, root-finding algorithm failed." << '\n');
    }
    return theSolution.success;

  }

}

#endif // __GENIE_INCL_ENABLED__
