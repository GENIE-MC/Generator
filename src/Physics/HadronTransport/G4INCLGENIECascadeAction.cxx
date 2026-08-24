/*
 * =====================================================================================
 *
 *       Filename:  G4INCLGENIECascadeAction.cxx
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  04/22/2026 04:11:54 PM
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

#include "G4INCLGENIECascadeAction.h"
#include "G4INCLGENIEParticleRecord.h"
#include "Physics/NuclearState/INCLNucleus.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Interaction/ProcessInfo.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include <sstream>
#include <string>

namespace G4INCL {
  using namespace genie;

  GENIECascadeAction::GENIECascadeAction():
    eventCounter(0), stepCounter(0)
  {
  }

  GENIECascadeAction::~GENIECascadeAction() {}

  void GENIECascadeAction::beforeRunUserAction(Config const *){
  }

  void GENIECascadeAction::beforeCascadeUserAction(IPropagationModel * /*pm*/) {

    // reset and initial
    eventRecord.clear();
    tempFinalState.clear();
    stepFinalState.clear();
    backup_mother.clear();   // no target-remnant snapshot may cross an event boundary

    const ProcessInfo & proc_info = evrec->Summary()->ProcInfo();
    // convert ghep event record to INCL Style.
    // G4INCL::GENIEParticleRecord is the bridge
    TObjArrayIter piter(evrec);
    GHepParticle * p = nullptr;
    eventRecord.clear();
    tempFinalState.clear();
    while ( (p = (GHepParticle *) piter.Next() ) ) {
      // the code of the particles in primary neutrino interaction
      G4INCL::GENIERecordCode recordCode;
      if     (eventRecord.size() == static_cast<size_t>(evrec->ProbePosition()))                         { recordCode = G4INCL::kProbe; }
      else if(eventRecord.size() == static_cast<size_t>(evrec->TargetNucleusPosition()))            { recordCode = G4INCL::kTarget;}
      else if(eventRecord.size() == static_cast<size_t>(evrec->HitNucleonPosition()))               { recordCode = G4INCL::kHitNucleon;}
      else if(eventRecord.size() == static_cast<size_t>(evrec->RemnantNucleusPosition()))           { recordCode = G4INCL::kRemnant;}
      else if(eventRecord.size() == static_cast<size_t>(evrec->FinalStatePrimaryLeptonPosition()))  { recordCode = G4INCL::kFinalStateLepton;}
      else { recordCode = G4INCL::kUnknown;}
      eventRecord.emplace_back(p, int(proc_info.ScatteringTypeId()), recordCode);
    }
  }

  void GENIECascadeAction::beforePropagationUserAction(IPropagationModel *){
  }
  void GENIECascadeAction::beforeAvatarUserAction(IAvatar *avatar, Nucleus *) { 
    mlist.clear();
    mlist = avatar->getParticles();
    backup_mother.clear();
    for(G4INCL::ParticleIter imom =  mlist.begin(); imom != mlist.end(); imom++){
      this->fillStep(*imom, backup_mother, kUnknownType, -1);
    }
  }

  void GENIECascadeAction::afterAvatarUserAction(IAvatar *avatar, Nucleus *nucleus, FinalState *finalState) {
    (void) nucleus;
    this->fillEventRecord(finalState, mlist, avatar->getTime(), avatar->getType());
  }
  void GENIECascadeAction::afterNPVAvatarUserAction(IAvatar *avatar, Nucleus *nucleus, FinalState *finalState) {
    (void) avatar;
    (void) nucleus;


    // update the event record after INCL postInteraction
    // INCL might rescale the four momentum of final states
    // particles, we update p4 in GHep event record

    TObjArrayIter piter(evrec);
    piter.Reset();   // rewind
    GHepParticle * p = nullptr;
    auto evr = eventRecord.begin();
    int idx =0;
    while ( (p = (GHepParticle *) piter.Next() ) ) {
      TLorentzVector *p4 = p->P4();
      p4->SetPx(evr->P3().getX() * MeV / GeV);
      p4->SetPy(evr->P3().getY() * MeV / GeV);
      p4->SetPz(evr->P3().getZ() * MeV / GeV);
      p4->SetE(std::sqrt(evr->P3().mag2() + evr->Mass()*evr->Mass()) * MeV / GeV);
      tempFinalState.emplace_back(evr->ID(), evr->Pdg(), evr->FirstMother(), idx++, kUnknownType);
      evr++;
    }

    // std::vector<GENIEParticleRecord> *eventRecord = avatar->getEventRecord();
    ParticleList outgoing = finalState->getOutgoingParticles();

    if(!outgoing.empty()){

      std::unique_ptr<G4INCL::FinalState> tmpfs(new FinalState);
      G4INCL::ParticleList modified  = finalState->getModifiedParticles();
      G4INCL::ParticleList created   = finalState->getCreatedParticles();
      G4INCL::ParticleList Destroyed = finalState->getDestroyedParticles();

      for(ParticleIter i = modified.begin(), e = modified.end(); i!=e; ++i){
        tmpfs->addModifiedParticle((*i));
      }
      for(ParticleIter i = Destroyed.begin(), e = Destroyed.end(); i!=e; ++i){
        tmpfs->addDestroyedParticle((*i));
      }
      for(ParticleIter i = created.begin(), e = created.end(); i!=e; ++i){
        if(!((*i)->isOutOfWell())){
          tmpfs->addCreatedParticle((*i));
        }
      }
      for(ParticleIter i = outgoing.begin(), e = outgoing.end(); i!=e; ++i){
        tmpfs->addOutgoingParticle((*i));
      }

      finalState->reset();
      modified  = tmpfs->getModifiedParticles();
      created   = tmpfs->getCreatedParticles();
      Destroyed = tmpfs->getDestroyedParticles();
      outgoing =  tmpfs->getOutgoingParticles();

      for(ParticleIter i = modified.begin(), e = modified.end(); i!=e; ++i){
        finalState->addModifiedParticle((*i));
      }
      for(ParticleIter i = Destroyed.begin(), e = Destroyed.end(); i!=e; ++i){
        finalState->addDestroyedParticle((*i));
      }
      for(ParticleIter i = created.begin(), e = created.end(); i!=e; ++i){
        finalState->addCreatedParticle((*i));
      }
      for(ParticleIter i = outgoing.begin(), e = outgoing.end(); i!=e; ++i){
        finalState->addOutgoingParticle((*i));
      }
    }


    for(ParticleIter iter=outgoing.begin(); iter!=outgoing.end(); ++iter){
      int outp_mother_idx = -1;
      int tmp_idx_ = -1;
      int pdg = 0;
      for(auto er = eventRecord.begin(); er != eventRecord.end(); er++){
        tmp_idx_++;
        if((*iter)->getID() == er->ID()){
          outp_mother_idx = tmp_idx_;
          pdg = er->Pdg();
        }
      }
      GHepParticle fsp(pdg, kIStStableFinalState, outp_mother_idx, -1, -1, -1, 
          TLorentzVector((*iter)->getMomentum().getX()  * MeV / GeV,
            (*iter)->getMomentum().getY()  * MeV / GeV,
            (*iter)->getMomentum().getZ()  * MeV / GeV,
            (*iter)->getEnergy()  * MeV / GeV),
          TLorentzVector((*iter)->getPosition().getX(),
            (*iter)->getPosition().getY(),
            (*iter)->getPosition().getZ(),
            0)
          );
      evrec->AddParticle(fsp);
      tempFinalState.emplace_back((*iter)->getID(), pdg, outp_mother_idx, idx++, kUnknownType);
    }

  }
  void GENIECascadeAction::afterPropagationUserAction(IPropagationModel *, IAvatar *) {
  }

  void GENIECascadeAction::afterCascadeUserAction(Nucleus * /*nucleus*/) {
    std::cout << "DEBUG: " << __FILE__ << ":" << __LINE__ << std::endl;
  }
  void GENIECascadeAction::afterRunUserAction() {
    std::cout << "DEBUG: " << __FILE__ << ":" << __LINE__ << std::endl;
  }

  void GENIECascadeAction::fillStep(Particle *par, std::vector<INCLRecord> &stepList, G4INCLFinalStateType type, double time){

    ParticleSpecies ptype = par->getSpecies();
    TLorentzVector p4mom;
    p4mom.SetPx(par->getMomentum().getX()  * MeV / GeV);
    p4mom.SetPy(par->getMomentum().getY()  * MeV / GeV);
    p4mom.SetPz(par->getMomentum().getZ()  * MeV / GeV);
    p4mom.SetE(par->getEnergy()  * MeV / GeV);
    TLorentzVector p4posi;
    p4posi.SetX(par->getPosition().getX());
    p4posi.SetY(par->getPosition().getY());
    p4posi.SetZ(par->getPosition().getZ());
    p4posi.SetT(time);
    stepList.emplace_back(par->getID(), ptype.getPDGCode(), -2, 0, type,  p4mom, p4posi, ptype.theType);

  }

  void GENIECascadeAction::snapshotTargetRemnant(Particle *remnant, double time){
    // Cluster::getSpecies() reports theType==Composite with the current (A,Z,S), so the record
    // carries the INCL composite code (A + 1000*Z - 1e6*S) and the pre-decay four-momentum.
    backup_mother.clear();
    this->fillStep(remnant, backup_mother, kComposite, time);
    LOG("INCLCascadeIntranuke", pINFO) << "snapshotTargetRemnant: INCL ID " << remnant->getID()
      << " A=" << remnant->getA() << " Z=" << remnant->getZ() << " S=" << remnant->getS()
      << " composite code " << backup_mother.back().pdgid
      << " E=" << backup_mother.back().p4mom.E() << " GeV";
  }

  void GENIECascadeAction::fillEventRecord(FinalState *fs, ParticleList mother_list, double time, G4INCL::AvatarType avaType){


    std::vector<INCLRecord> stepParticleList;
    stepParticleList.clear();
    stepCounter++;
    LOG("INCLCascadeIntranuke", pDEBUG) << "stepCounter : " << stepCounter;
    for(ParticleIter iter=mother_list.begin(); iter!=mother_list.end(); ++iter){
      this->fillStep(*iter, stepParticleList, kMother, time);
      ParticleSpecies ptype = (*iter)->getSpecies();
      stepFinalState[stepCounter].emplace_back((*iter)->getID(), ptype.getPDGCode(), -2, 0, kMother);
      LOG("INCLCascadeIntranuke", pDEBUG) << "mother list ID : " << (*iter)->getID() << " pdg : " << ptype.getPDGCode(); 
    }
    ParticleList modified = fs->getModifiedParticles();
    for(ParticleIter iter=modified.begin(); iter!=modified.end(); ++iter){
      this->fillStep(*iter, stepParticleList, kModified, time);
      ParticleSpecies ptype = (*iter)->getSpecies();
      stepFinalState[stepCounter].emplace_back((*iter)->getID(), ptype.getPDGCode(), -2, 0, kModified);
      LOG("INCLCascadeIntranuke", pDEBUG) << "Modified ID : " << (*iter)->getID() << " pdg : " << ptype.getPDGCode(); 
    }
    ParticleList outgoing = fs->getOutgoingParticles();
    for(ParticleIter iter=outgoing.begin(); iter!=outgoing.end(); ++iter){
      this->fillStep(*iter, stepParticleList, kOutgoing, time);
      ParticleSpecies ptype = (*iter)->getSpecies();
      stepFinalState[stepCounter].emplace_back((*iter)->getID(), ptype.getPDGCode(), -2, 0, kOutgoing);
      LOG("INCLCascadeIntranuke", pDEBUG) << "Outgoing ID : " << (*iter)->getID() << " pdg : " << ptype.getPDGCode(); 
    }
    ParticleList destroyed = fs->getDestroyedParticles();
    for(ParticleIter iter=destroyed.begin();  iter!=destroyed.end(); ++iter){
      this->fillStep(*iter, stepParticleList, kDestroyed, time);
      ParticleSpecies ptype = (*iter)->getSpecies();
      stepFinalState[stepCounter].emplace_back((*iter)->getID(), ptype.getPDGCode(), -2, 0, kDestroyed);
      LOG("INCLCascadeIntranuke", pDEBUG) << "Destroyed ID : " << (*iter)->getID() << " pdg : " << ptype.getPDGCode(); 
    }
    ParticleList created = fs->getCreatedParticles();
    for(ParticleIter iter=created.begin(); iter!=created.end(); ++iter){
      this->fillStep(*iter, stepParticleList, kCreated, time);
      ParticleSpecies ptype = (*iter)->getSpecies();
      stepFinalState[stepCounter].emplace_back((*iter)->getID(), ptype.getPDGCode(), -2, 0, kCreated);
      LOG("INCLCascadeIntranuke", pDEBUG) << "Created ID : " << (*iter)->getID() << " pdg : " << ptype.getPDGCode(); 
    }

    int num_partiles = tempFinalState.size();
    if(created.size() == 0 && outgoing.size() == 0 && modified.size() == 0) return;
    int index_ = num_partiles;
    if(mother_list.size() == 1){
      // FIXME: channel with mother size = 1 could be reflection, outgoing and decay
      // Try to match the mother in step to the event record history
      // if the mother of this step is not in the event record history and it is a reflection
      // avatar, we will not put this step into event record history
      int mother_position = -1;
      for(auto ip = stepParticleList.begin(); ip != stepParticleList.end(); ++ip){
        if(ip->fsType != kMother) continue;
        for(auto it = tempFinalState.begin(); it != tempFinalState.end(); ++it){
          if(it->global_index == ip->global_index){
            mother_position = it->local_index;
          }
        }
      }

      if(mother_position != -1){

        if(avaType == G4INCL::SurfaceAvatarType){
          //
          for(auto ip = stepParticleList.begin(); ip != stepParticleList.end(); ++ip){
            if(ip->fsType == kMother || ip->fsType == kDestroyed) continue;
            LOG("INCLCascadeIntranuke", pNOTICE) << "created: " << created.size();
            LOG("INCLCascadeIntranuke", pNOTICE) << "outgoing: " << outgoing.size();
            LOG("INCLCascadeIntranuke", pNOTICE) << "modified: " << modified.size();
            LOG("INCLCascadeIntranuke", pNOTICE) << "destroyed: " << destroyed.size();
            // put this particle in event record
            tempFinalState.emplace_back(ip->global_index, ip->pdgid, mother_position, index_++, ip->fsType, ip->p4mom, ip->p4posi);
            evrec->Particle(mother_position)->SetRescatterCode(int(avaType));

            // get the pdg in genie style
            int pdg = ip->pdgid;
            LOG("INCLCascadeIntranuke", pINFO) << "PDG : " << ip->pdgid;
            if(ip->theType == G4INCL::Composite){
              if(pdg != genie::kPdgProton && pdg != genie::kPdgNeutron && pdg != genie::kPdgLambda){
                int S = pdg / 1000000;
                int pdg_no_s = pdg % 1000000;
                int A = pdg_no_s % 1000;
                int Z = pdg_no_s / 1000;
                pdg = genie::pdg::IonPdgCode( A , Z, S, 0);
                pdg = GENIEINCLUtil::INCLPDG_to_GHEPPDG(pdg, A, Z, S);
              }
            }
            LOG("INCLCascadeIntranuke", pINFO) << "PDG : " << pdg;
            // get the type of the particle
            EGHepStatus ptype;

            if(ip->fsType == kOutgoing){
              //ptype = kIStPreDecayResonantState;
              ptype = kIStStableFinalState;
            }
            else{
              ptype = kIStHadronInTheNucleus;
            }

            GHepParticle p(pdg, ptype, mother_position, -1, -1, -1, ip->p4mom, ip->p4posi);
            evrec->AddParticle(p);
          }
        }
        else if(avaType == G4INCL::DecayAvatarType){

          /* If the decay was Pauli-blocked, make sure the propagation model
           * generates a new decay avatar on the next call to propagate().
           *
           * \bug{Note that we don't generate new decay avatars for deltas that
           * could not satisfy energy conservation. This is in keeping with
           * INCL4.6, but doesn't seem to make much sense to me (DM), as energy
           * conservation can be impossible to satisfy due to weird local-energy
           * conditions, for example, that evolve with time.}
           */

          /* decay channel without daughters, skip it
           * a decay must have more than 1 daughters
           */

          if((created.size() + outgoing.size() + modified.size()) != 1){
            evrec->Particle(mother_position)->SetRescatterCode(int(avaType));
            evrec->Particle(mother_position)->SetStatus(kIStDecayedState);
            for(auto ip = stepParticleList.begin(); ip != stepParticleList.end(); ++ip){
              if(ip->fsType == kMother || ip->fsType == kDestroyed) continue;
              LOG("INCLCascadeIntranuke", pNOTICE) << "created: " << created.size();
              LOG("INCLCascadeIntranuke", pNOTICE) << "outgoing: " << outgoing.size();
              LOG("INCLCascadeIntranuke", pNOTICE) << "modified: " << modified.size();
              LOG("INCLCascadeIntranuke", pNOTICE) << "destroyed: " << destroyed.size();
              // put this particle in event record
              tempFinalState.emplace_back(ip->global_index, ip->pdgid, mother_position, index_++, ip->fsType, ip->p4mom, ip->p4posi);

              // get the pdg in genie style
              int pdg = ip->pdgid;
              LOG("INCLCascadeIntranuke", pINFO) << "PDG : " << ip->pdgid;
              if(ip->theType == G4INCL::Composite){
                if(pdg != genie::kPdgProton && pdg != genie::kPdgNeutron && pdg != genie::kPdgLambda){
                  int S = pdg / 1000000;
                  int pdg_no_s = pdg % 1000000;
                  int A = pdg_no_s % 1000;
                  int Z = pdg_no_s / 1000;
                  pdg = genie::pdg::IonPdgCode( A , Z, S, 0);
                }
              }
              LOG("INCLCascadeIntranuke", pINFO) << "PDG : " << pdg;
              // get the type of the particle
              EGHepStatus ptype;

              if(ip->fsType == kOutgoing){
                //ptype = kIStPreDecayResonantState;
                ptype = kIStStableFinalState;
              }
              else{
                ptype = kIStHadronInTheNucleus;
              }

              GHepParticle p(pdg, ptype, mother_position, -1, -1, -1, ip->p4mom, ip->p4posi);
              evrec->AddParticle(p);
            }
          }
        }
      }
      else if(avaType == G4INCL::DecayAvatarType){
        // The mother is NOT in the event-record history. The legitimate case is decayMe(): the
        // target remnant itself (never a cascade particle, hence never in tempFinalState) was
        // phase-space decayed, and decayMe() snapshotted it via snapshotTargetRemnant() BEFORE
        // ClusterDecay::decay() turned it into its last nucleon. The snapshot is applied only to
        // that remnant (matched by INCL ID). Any other decay whose mother was never registered
        // (e.g. an outgoing cluster whose emission was not recorded) must not kill the job: its
        // products are attached to the intermediate remnant nucleus (GHEP index 3, the convention
        // used for the pre-de-excitation remnant) and enough is logged to find the root cause.
        const int motherID = mother_list.empty() ? -1 : (*mother_list.begin())->getID();
        const INCLRecord * snap = nullptr;
        for(auto imom = backup_mother.begin(); imom != backup_mother.end(); ++imom){
          if(imom->fsType == kComposite && imom->theType == G4INCL::Composite && imom->global_index == motherID){
            snap = &(*imom);
            break;
          }
        }
        if(snap){
          int A = snap->pdgid%1000;
          int Z = (snap->pdgid/1000)%1000;
          int S = (snap->pdgid/1000000)%1000;
          int pdg = genie::pdg::IonPdgCode(A,Z,S,0);
          pdg = GENIEINCLUtil::INCLPDG_to_GHEPPDG(pdg, A, Z, S);
          LOG("INCLCascadeIntranuke", pINFO) << "decayMe: pre-decay target remnant A=" << A << " Z=" << Z
            << " S=" << S << " -> PDG " << pdg;
          GHepParticle p(pdg, kIStPreDeExNuclearRemnant, 3, -1, -1, -1, snap->p4mom, snap->p4posi);
          tempFinalState.emplace_back(snap->global_index, snap->pdgid, 3, index_++, snap->fsType, snap->p4mom, snap->p4posi, snap->theType);
          evrec->AddParticle(p);
          mother_position = index_ - 1;
        }
        else{
          mother_position = 3;
          LOG("INCLCascadeIntranuke", pERROR)
            << "fillEventRecord: the mother (INCL ID " << motherID << ") of a decay step is not in the"
            << " event-record history and no target-remnant snapshot matches it; attaching "
            << (created.size() + outgoing.size() + modified.size()) << " product(s) to the intermediate"
            << " remnant nucleus (GHEP idx " << mother_position << "). mothers=" << mother_list.size()
            << " created=" << created.size() << " outgoing=" << outgoing.size()
            << " modified=" << modified.size() << " destroyed=" << destroyed.size();
          for(auto ip = stepParticleList.begin(); ip != stepParticleList.end(); ++ip){
            LOG("INCLCascadeIntranuke", pERROR) << "   step particle: INCL ID " << ip->global_index
              << " INCL pdg " << ip->pdgid << " fsType " << ip->fsType << " INCL type " << ip->theType
              << " E=" << ip->p4mom.E() << " GeV";
          }
        }
        backup_mother.clear();   // the snapshot is single-use
        LOG("INCLCascadeIntranuke", pINFO) << "the index of the decay mother : " << mother_position;
        for(auto ip = stepParticleList.begin(); ip != stepParticleList.end(); ++ip){
          if(ip->fsType == kMother || ip->fsType == kDestroyed) continue;
          LOG("INCLCascadeIntranuke", pNOTICE) << "created: " << created.size();
          LOG("INCLCascadeIntranuke", pNOTICE) << "outgoing: " << outgoing.size();
          LOG("INCLCascadeIntranuke", pNOTICE) << "modified: " << modified.size();
          LOG("INCLCascadeIntranuke", pNOTICE) << "destroyed: " << destroyed.size();
          // put this particle in event record
          tempFinalState.emplace_back(ip->global_index, ip->pdgid, mother_position, index_++, ip->fsType, ip->p4mom, ip->p4posi);
          evrec->Particle(mother_position)->SetRescatterCode(int(avaType));

          // get the pdg in genie style
          int pdg = ip->pdgid;
          LOG("INCLCascadeIntranuke", pINFO) << "PDG : " << ip->pdgid;
          if(ip->theType == G4INCL::Composite){
            if(pdg != genie::kPdgProton && pdg != genie::kPdgNeutron && pdg != genie::kPdgLambda){
              int S = pdg / 1000000;
              int pdg_no_s = pdg % 1000000;
              int A = pdg_no_s % 1000;
              int Z = pdg_no_s / 1000;
              pdg = genie::pdg::IonPdgCode( A , Z, S, 0);
              pdg = GENIEINCLUtil::INCLPDG_to_GHEPPDG(pdg, A, Z, S);
            }
          }
          LOG("INCLCascadeIntranuke", pINFO) << "PDG : " << pdg;
          // get the type of the particle
          EGHepStatus ptype;

          if(ip->fsType == kOutgoing){
            //ptype = kIStPreDecayResonantState;
            ptype = kIStStableFinalState;
          }
          else{
            ptype = kIStHadronInTheNucleus;
          }

          GHepParticle p(pdg, ptype, mother_position, -1, -1, -1, ip->p4mom, ip->p4posi);
          evrec->AddParticle(p);
        }
      }
    }
    else if(mother_list.size() == 2){
      int first_mother_position = -1;
      int second_mother_position = -1;
      // find parents for nn-collision
      for(auto ip = stepParticleList.begin(); ip != stepParticleList.end(); ++ip){
        if(ip->fsType != kMother) continue;
        int mother_position = -1;
        bool is_spectator = true;
        for(auto it = tempFinalState.begin(); it != tempFinalState.end(); ++it){
          if(it->global_index == ip->global_index){
            mother_position = it->local_index;
            is_spectator = false;
          }
        }
        if(is_spectator){
          for(auto imom = backup_mother.begin(); imom != backup_mother.end(); ++imom){
            if(imom->global_index == ip->global_index){
              tempFinalState.emplace_back(ip->global_index, imom->pdgid, -1, index_++, imom->fsType, imom->p4mom, imom->p4posi);
              GHepParticle p(imom->pdgid, kIStSpectator, -1, -1, -1, -1, imom->p4mom, imom->p4posi);
              evrec->AddParticle(p);
            }
          }
          mother_position = index_ - 1;
        }
        // assign the mother index
        if(first_mother_position == -1)
          first_mother_position = mother_position;
        else
          second_mother_position = mother_position;
      }

      if(first_mother_position > second_mother_position){
        int temp_position = first_mother_position;
        first_mother_position = second_mother_position;
        second_mother_position = temp_position;
      }
      LOG("INCLCascadeIntranuke", pNOTICE) << first_mother_position;
      LOG("INCLCascadeIntranuke", pNOTICE) << second_mother_position;
      evrec->Particle(first_mother_position)->SetRescatterCode(int(avaType));
      evrec->Particle(second_mother_position)->SetRescatterCode(int(avaType));

      for(auto ip = stepParticleList.begin(); ip != stepParticleList.end(); ++ip){
        if(ip->fsType == kMother || ip->fsType == kDestroyed) continue;
        LOG("INCLCascadeIntranuke", pNOTICE) << "created: " << created.size();
        LOG("INCLCascadeIntranuke", pNOTICE) << "outgoing: " << outgoing.size();
        LOG("INCLCascadeIntranuke", pNOTICE) << "modified: " << modified.size();
        LOG("INCLCascadeIntranuke", pNOTICE) << "destroyed: " << destroyed.size();
        // put this particle in event record
        tempFinalState.emplace_back(ip->global_index, ip->pdgid, first_mother_position, index_++, ip->fsType, ip->p4mom, ip->p4posi);

        // get the pdg in genie style
        int pdg = ip->pdgid;
        LOG("INCLCascadeIntranuke", pINFO) << "PDG : " << ip->pdgid;
        if(ip->theType == G4INCL::Composite){
          if(pdg != genie::kPdgProton && pdg != genie::kPdgNeutron && pdg != genie::kPdgLambda){
            int S = pdg / 1000000;
            int pdg_no_s = pdg % 1000000;
            int A = pdg_no_s % 1000;
            int Z = pdg_no_s / 1000;
            pdg = genie::pdg::IonPdgCode( A , Z, S, 0);
          }
        }
        LOG("INCLCascadeIntranuke", pINFO) << "PDG : " << pdg;
        // get the type of the particle
        EGHepStatus ptype;
        ptype = kIStHadronInTheNucleus;
        GHepParticle p(pdg, ptype, first_mother_position, second_mother_position, -1, -1, ip->p4mom, ip->p4posi);
        evrec->AddParticle(p);
      }
      // update the daughter's indices for spectator
      if(evrec->Particle(first_mother_position)->FirstDaughter() == -1){
        evrec->Particle(first_mother_position)->SetFirstDaughter(evrec->Particle(second_mother_position)->FirstDaughter());
        evrec->Particle(first_mother_position)->SetLastDaughter(evrec->Particle(second_mother_position)->LastDaughter());
      }
      else if(evrec->Particle(second_mother_position)->FirstDaughter() == -1){
        evrec->Particle(second_mother_position)->SetFirstDaughter(evrec->Particle(first_mother_position)->FirstDaughter());
        evrec->Particle(second_mother_position)->SetLastDaughter(evrec->Particle( first_mother_position)->LastDaughter());
      }

    }
    else{
      LOG("INCLCascadeIntranuke", pNOTICE) << "wrong mother size : " << mother_list.size();
      LOG("INCLCascadeIntranuke", pNOTICE) << "created: " << created.size();
      LOG("INCLCascadeIntranuke", pNOTICE) << "outgoing: " << outgoing.size();
      LOG("INCLCascadeIntranuke", pNOTICE) << "modified: " << modified.size();
      LOG("INCLCascadeIntranuke", pNOTICE) << "destroyed: " << destroyed.size();
      exit(1);
    }



  }

}
#endif // __GENIE_INCL_ENABLED__
