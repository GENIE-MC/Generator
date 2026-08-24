/*
 * =====================================================================================
 *
 *       Filename:  G4INCLGENIECascadeAction.h
 *
 *    Description: : 
 *
 *        Version:  1.0
 *        Created:  04/22/2026 04:11:48 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Liang Liu (L. Liu), liangliu@fnal.gov
 *		    Fermi National Accelerator Laboratory
 *  Collaboration:  GENIE
 *
 * =====================================================================================
 */

#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_INCL_ENABLED__

#ifndef G4INCLGENIECASCADEACTION_HH
#define G4INCLGENIECASCADEACTION_HH 1

#include "Physics/NuclearState/INCLNucleus.h"
#include "G4INCLCascadeAction.hh"
#include <fstream>

#include "Framework/GHEP/GHepRecord.h"
#include "G4INCLGENIEParticleRecord.h"
#include "Framework/Conventions/Units.h"
#include "Framework/ParticleData/PDGCodes.h"

namespace G4INCL {

  class GENIECascadeAction : public CascadeAction {

    // class CascadeAction make class INCL a friend because it needs to call private methods
    // but for GENIE interaction, we don't want to change the code inside G4INCL::INCL yet
    // Althrough GENIECascadeAction is inherited from CascadeAction, I only using the virtual 
    // functions.

    public:
      GENIECascadeAction();
      virtual ~GENIECascadeAction();

      virtual void beforeRunUserAction(Config const *);
      virtual void beforeCascadeUserAction(IPropagationModel *);
      virtual void beforePropagationUserAction(IPropagationModel *);
      virtual void beforeAvatarUserAction(IAvatar *, Nucleus *);
      virtual void afterAvatarUserAction(IAvatar *, Nucleus *, FinalState *);
      virtual void afterNPVAvatarUserAction(IAvatar *, Nucleus *, FinalState *);
      virtual void afterPropagationUserAction(IPropagationModel *, IAvatar *);
      virtual void afterCascadeUserAction(Nucleus *);
      virtual void afterRunUserAction();

      void setGHepRecord(genie::GHepRecord *evr){
        evrec = evr;
      }

      std::vector<G4INCL::GENIEParticleRecord> *getEventRecord(){
        return &eventRecord;
      }
      void fillStep(Particle *par, std::vector<INCLRecord> &stepList, G4INCLFinalStateType type, double time);
      void fillEventRecord(FinalState *fs, ParticleList mother_list, double time, G4INCL::AvatarType avaType);
      /// Snapshot the TARGET REMNANT (the G4INCL::Nucleus itself; decayMe() only) BEFORE
      /// ClusterDecay::decay() mutates it in place into its last nucleon. The remnant is the only
      /// decaying object that is never in the GHEP record; outgoing clusters are already recorded
      /// when emitted and must NOT use this. Consumed (and cleared) by fillEventRecord() on the
      /// decayMe DecayAvatarType step, matched by INCL ID, to emit the kIStPreDeExNuclearRemnant entry.
      void snapshotTargetRemnant(Particle *remnant, double time);

    private:
      //std::ofstream *oFile;
      long eventCounter;
      long stepCounter;
      ParticleList mlist;

      genie::GHepRecord * evrec;
      static constexpr double MeV = genie::units::MeV;
      static constexpr double GeV = genie::units::GeV;

      std::vector<G4INCL::GENIEParticleRecord> eventRecord;
      std::vector<INCLRecord> tempFinalState;
      std::map<int, std::vector<INCLRecord>> stepFinalState;
      std::vector<INCLRecord> backup_mother; // this cantainer is used to store the initial momentum and position of mother particles in binary collision



  };

}
#endif // G4INCLGENIECASCADEACTION_HH

#endif // __GENIE_INCL_ENABLED__
