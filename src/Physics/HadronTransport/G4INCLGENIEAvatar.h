#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_INCL_ENABLED__

#ifndef G4INCLGENIEAVATAR_HH_
#define G4INCLGENIEAVATAR_HH_

#include "G4INCLInteractionAvatar.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLParticle.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLAllocationPool.hh"
#include "Physics/HadronTransport/G4INCLGENIEParticleRecord.h"

namespace G4INCL {

  /**
   * Decay avatar
   *
   * The reflection avatar is created when a particle reaches the boundary of the nucleus.
   * At this point it can either be reflected from the boundary or exit the nucleus.
   */
  class GENIEAvatar: public InteractionAvatar {
    public:
      GENIEAvatar(double time, Particle *p, Nucleus *n, std::vector<GENIEParticleRecord> *eventRecord, bool ishybrid);
      GENIEAvatar(double time, Cluster *p, Nucleus *n, std::vector<GENIEParticleRecord> *eventRecord, bool ishybrid);
      virtual ~GENIEAvatar();

      IChannel* getChannel();
      void fillFinalState(FinalState *fs);

      virtual void preInteraction();
      virtual void postInteraction(FinalState *fs);

      ParticleList getParticles() const {
        ParticleList theParticleList;
        theParticleList.push_back(particle1);
        return theParticleList;
      }

      std::vector<GENIEParticleRecord> *getEventRecord() {
        return genie_evtrec;
      }

      std::string dump() const;
    private:
      Particle *particle1;
      Nucleus  *theNucleus;
      Cluster  *cluster;
      std::vector<GENIEParticleRecord> *genie_evtrec;
      ThreeVector boostVector;
      ThreeVector leptonMom;
      double leptonE;
      double oldTotalEnergy;

      ParticleList modified, created, modifiedAndCreated, Destroyed;

      bool enforceEnergyConservation(FinalState * const fs);

      bool putIntoPotential(const double theQValueCorrection, Particle * const theParticle);

      int delta_Z;

      bool fHybridModel;

      void postInteractionHybridModel(FinalState *fs);


    private:

      /// \brief An object that hold the energy and momentum of lepton
      class Lepton {

        public:
          Lepton(ThreeVector mom, double e):
            theMomentum(mom), theEnergy(e){
              theMass = std::sqrt(std::max( theEnergy*theEnergy - theMomentum.mag2(), 0.));
            }
          ~Lepton(){
          }
          void boost(const ThreeVector &aBoostVector){

            const double beta2 = aBoostVector.mag2();
            const double gamma = 1.0 / std::sqrt(1.0 - beta2);
            const double bp = theMomentum.dot(aBoostVector);
            const double alpha = (gamma*gamma)/(1.0 + gamma);

            theMomentum = theMomentum + aBoostVector * (alpha * bp - gamma * theEnergy);
            theEnergy = gamma * (theEnergy - bp);
          }

          ThreeVector getMomentum(){
            return theMomentum;
          }
          double getEnergy(){
            return theEnergy;
          }

          void setMomentum(const ThreeVector &momentum){
            theMomentum = momentum;
          }
          double adjustEnergyFromMomentum() {
            theEnergy = std::sqrt(theMomentum.mag2() + theMass*theMass);
            return theEnergy;
          }

        private:
          ThreeVector theMomentum;
          double      theEnergy;
          double      theMass;

      };


      ///  \brief RootFunctor-derived object for enforcing energy conservation in lepton-nucleon
      class ViolationLeptonEMomentumFunctor : public RootFunctor {
        public:
          /** \brief Prepare for calling the () operator and scaleParticleMomenta
           *
           * The constructor sets the private class members.
           */

          ViolationLeptonEMomentumFunctor(Nucleus * const nucleus, ParticleList const &modAndCre, 
              double &leptonE, ThreeVector &leptonMom, ThreeVector const &boost,
              const double totalEnergyBeforeInteraction, const bool localE);
          virtual ~ViolationLeptonEMomentumFunctor();

          /** \brief Compute the energy-conservation violation.
           *
           * \param x scale factor for the particle momenta
           * \return the energy-conservation violation
           */
          double operator()(const double x) const;

          /// \brief Clean up after root finding
          void cleanUp(const bool success) const;

        private:
          /// \brief List of final-state particles.
          ParticleList finalParticles;
          /// \brief energy of final state lepton
          double &leptonEnergy;
          /// \brief momentum of final state lepton
          ThreeVector &leptonMomentum;
          /// \brief Total energy before the interaction.
          double initialEnergy;
          /// \brief Pointer to the nucleus
          Nucleus *theNucleus;
          /// \brief Pointer to the boost vector
          ThreeVector const &boostVector;
          /// \brief True if we should use local energy
          const bool shouldUseLocalEnergy;
          /// \brief CM particle momenta, as determined by the channel.
          std::vector<ThreeVector> particleMomenta;



          /// \brief CM momentum of final state lepton
          ThreeVector leptonMomentumCM;


          Lepton *lepton;

          /** \brief Scale the momenta of the modified and created particles.
           *
           * Set the momenta of the modified and created particles to alpha times
           * their original momenta (stored in particleMomenta). You must call
           * init() before using this method.
           *
           * \param alpha scale factor
           */
          void scaleParticleMomenta(const double alpha) const;

      };

      RootFunctor *violationEFunctor;



      INCL_DECLARE_ALLOCATION_POOL(GENIEAvatar)
  };

}

#endif /* G4INCLDECAYAVATAR_HH_ */
#endif // __GENIE_INCL_ENABLED__
