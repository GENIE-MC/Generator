//____________________________________________________________________________
/*
 Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Jim Dobson <j.dobson07 \at imperial.ac.uk>
         Imperial College London

 For documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Mar 04, 2009 - JD
   Was first added in v2.5.1.
 @ Mar 05, 2009 - CA
   Modified ReconstructHadronFateHA() to work with hadron+A event files in
   addition to neutrino event files.

*/
//____________________________________________________________________________

#include "Algorithm/AlgConfigPool.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepParticle.h"
#include "HadronTransport/INukeUtils.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGCodes.h"
#include "Registry/Registry.h"
#include "Utils/NuclearUtils.h"

using namespace genie;

//____________________________________________________________________________
INukeFateHA_t genie::utils::intranuke::ReconstructHadronFateHA(
   GHepRecord * event, int i, bool hA_mode)
{
// Reconstruct the INTRANUKE/hA model fate for the hadron at position i.
// Returns one of
// - kIHAFtUndefined  -> undefined (no re-interaction, or i not a hadron)   
// - kIHAFtCEx        -> charge exchange            
// - kIHAFtElas       -> elastic
// - kIHAFtInelas     -> inelastic
// - kIHAFtAbsNP      -> absorption followed by emission of np
// - kIHAFtAbsPP      -> absorption followed by emission of pp
// - kIHAFtAbsNPP      -> absorption followed by emission of npp
// - kIHAFtAbsNNP     -> absorption followed by emission of nnp
// - kIHAFtAbs2N2P    -> absorption followed by emission of 2n2p
// - kIHAFtAbs2N3P    -> absorption followed by emission of 2n3p
// - kIHAFtNPip       -> pi production : n pi+
// - kIHAFtNPipPi0    -> pi production : n pi+ pi0
//


  // get the particle from the event record
  GHepParticle * p = dynamic_cast<GHepParticle *>(event->Particle(i)) ;
  if(!p) return kIHAFtUndefined;

  //
  //
  if(!hA_mode) {
    int ist   = p->Status();
    if(ist != kIStHadronInTheNucleus) return kIHAFtUndefined;
  }

  int pdg_i = p->Pdg();
  bool is_handled = pdg_i == kPdgProton  || 
                    pdg_i == kPdgNeutron ||
                    pdg_i == kPdgPiP     || 
                    pdg_i == kPdgPi0     || 
                    pdg_i == kPdgPiM;
  if(!is_handled) return kIHAFtUndefined;

  // determine the binding energy for given nuclear target

  double binding_energy = 0.0;
  
  AlgConfigPool * conf_pool = AlgConfigPool::Instance();
  Registry * gpl = conf_pool->GlobalParameterList();
  
  // get target pdgc code
  int target_index = 1; // target index is always 1
  int target_pdgc  = event->Particle(target_index)->Pdg();
  
  // make string to pass to fUserPhysicsConfig->Get(...)
  string get_string = "RFG-NucRemovalE@Pdg=";
  char buffer[20];
  sprintf (buffer, "%i", target_pdgc);  
  get_string.append(buffer);
  
  if(gpl->Exists(get_string)) {
    binding_energy = gpl->GetDouble(get_string);
    LOG("INukeUtils", pDEBUG)
      << "Target pdgc code: " << target_pdgc 
      << " has binding energy = " << binding_energy;
  } else {
    Target * target_nucl = new Target(target_pdgc);  
    binding_energy = utils::nuclear::BindEnergy(*target_nucl);
  }

  //
  // tolerances have been set so as to catch events with zero energy change
  // and those that have had binding energy subtracted but to be narrow 
  // enough that almost no other events will be included due to machine 
  // floating point precision differences
  //
  double epsilon_zero_energy_max    = 2.0E-7;
  double epsilon_binding_energy_min = binding_energy - 0.00000001;
  double epsilon_binding_energy_max = binding_energy + 0.00000001;
  double epsilon_delta_P_squared    = 1.0E-12;

  //
  // Determine the fate
  //

  INukeFateHA_t hadron_fate = kIHAFtUndefined;
  
  int fd = p->FirstDaughter();
  int ld = p->LastDaughter();
  int nd = ld - fd + 1;
  if(nd==1) {
    int j = fd;
    if(event->Particle(j)->Status() == 3) {
      int fd_j = event->Particle(j)->FirstDaughter();
      int ld_j = event->Particle(j)->LastDaughter(); 
      int fd_j_pdg = event->Particle(fd_j)->Pdg();
      int ld_j_pdg = event->Particle(ld_j)->Pdg();
      bool is_pi0_gg = (event->Particle(j)->Pdg() == kPdgPi0 ) && 
        ( fd_j_pdg == kPdgGamma
          || fd_j_pdg == kPdgElectron
          || fd_j_pdg == kPdgPositron)
        && (ld_j_pdg == kPdgGamma
            || ld_j_pdg == kPdgElectron
            || ld_j_pdg == kPdgPositron);
      
      if(!is_pi0_gg) {
        fd = fd_j;
        ld = ld_j;
        nd = ld - fd + 1;
      }
    }                                                                                     
  }

   LOG("INukeUtils", pDEBUG) 
      << "Has " << nd << " daughters at [" << fd << ", " << ld << "]";
  
    int np   = 0;
    int nn   = 0;
    int npi0 = 0;
    int npip = 0;
    int npim = 0;
    for(int id = fd; id <= ld; id++) {
      int pdg_id = event->Particle(id)->Pdg();
       if (pdg_id == kPdgProton  ) np++;
       if (pdg_id == kPdgNeutron ) nn++;
       if (pdg_id == kPdgPi0     ) npi0++;
       if (pdg_id == kPdgPiP     ) npip++;
       if (pdg_id == kPdgPiM     ) npim++;
    }

    LOG("INukeUtils", pDEBUG) 
      << "p: " << np << ", n: " << nn 
      << ", pi+: " << npip << ", pi-" << npim << ", pi0: " << npi0;

    //
    // check whether it is absorption followed by nuclear breakup
    // hadron dissapears & many nucleons (>=2) appear in the final state
    // possible modes: pi + A -> A' + np/pp/npp/nnp/nnpp
    //                 N  + A -> A' + np/pp/npp/nnp/nnppp
    //
    if      (npip+npim+npi0==0 && nn==1 && np==1) hadron_fate = kIHAFtAbsNP;
    else if (npip+npim+npi0==0 && nn==0 && np==2) hadron_fate = kIHAFtAbsPP;
    else if (npip+npim+npi0==0 && nn==1 && np==2) hadron_fate = kIHAFtAbsNPP;
    else if (npip+npim+npi0==0 && nn==2 && np==1) hadron_fate = kIHAFtAbsNNP;
    else if (npip+npim+npi0==0 && nn==2 && np==2) hadron_fate = kIHAFtAbs2N2P;
    else if (npip+npim+npi0==0 && nn==2 && np==3) hadron_fate = kIHAFtAbs2N3P;

    //
    // check whether it is pion production
    // possible modes: pi + A -> A' + n + pi+ + pi0
    //                 N  + A -> A' + n + pi+
    //                 N  + A -> A' + n + pi+ + pi0
    //
    else if (nn==1 && np==0 && npi0==0 && npip==1 && npim==0) hadron_fate = kIHAFtNPip;
    else if (nn==1 && np==0 && npi0==1 && npip==1 && npim==0) hadron_fate = kIHAFtNPipPi0;
    
    //
    // check whether it is charge exchange
    // possible modes: p   + A -> A' + n
    //                 n   + A -> A' + p
    //                 pi- + A -> A' + pi0
    //                 pi+ + A -> A' + pi0
    //                 pi0 + A -> A' + pi+
    //                 pi0 + A -> A' + pi-
    //
    else if((nd==1)
             &&
             ((pdg_i == kPdgProton  && nn   == 1) ||
              (pdg_i == kPdgNeutron && np   == 1) ||
              (pdg_i == kPdgPiM     && npi0 == 1) ||
              (pdg_i == kPdgPiP     && npi0 == 1) ||
              (pdg_i == kPdgPi0     && npip == 1) ||
              (pdg_i == kPdgPi0     && npim == 1))  )
      {
        hadron_fate = kIHAFtCEx;
      }
    
    // if it is the same particle before and after rescattering
    //
    else if ((nd==1) && (pdg_i == event->Particle(fd)->Pdg())){
      double delta_E = TMath::Abs(event->Particle(i)->P4()->Energy() - event->Particle(fd)->P4()->Energy()); 
      double delta_P_squared = TMath::Power(event->Particle(i)->P4()->Px()-event->Particle(fd)->P4()->Px(), 2.) +   
                               TMath::Power(event->Particle(i)->P4()->Py()-event->Particle(fd)->P4()->Py(), 2.) + 
                               TMath::Power(event->Particle(i)->P4()->Pz()-event->Particle(fd)->P4()->Pz(), 2.);
      double mod_Pfd_squared = TMath::Power(event->Particle(fd)->P4()->Px(), 2.) +   
                               TMath::Power(event->Particle(fd)->P4()->Py(), 2.) + 
                               TMath::Power(event->Particle(fd)->P4()->Pz(), 2.);

      // and if there is no energy loss
      //
      if(delta_E < epsilon_zero_energy_max){
        // then it must be either an elastic
        //
        if(TMath::Abs(delta_P_squared) > epsilon_delta_P_squared){
          hadron_fate = kIHAFtElas;
        }
        // or no interaction at all
        //
        else {  
          hadron_fate = kIHAFtUndefined;
        }
      }
      
      // if deltaE > binding energy then must be inelastic
      else if( delta_E > epsilon_binding_energy_max){
        hadron_fate = kIHAFtInelas;
      }
      // if deltaE == binding energy 
      else if(delta_E > epsilon_binding_energy_min && delta_E < epsilon_binding_energy_max){
        // at this point know that only enery difference is nuclear binding 
        // energy (except in the unlikely case of inelastic scattering with 
        // coincidental energy loss == nucl binding energy of target.
        
        // if only had a nuclear binding energy subtraction and no elastic scattering
        // then all scale_px == scale_py == scale_pz 
        double scale_px = -1.0; //initial momentum in x / final momentum in x
        double scale_py = -2.0; //     "    "         y        "   "        y 
        double scale_pz = -3.0; //     "    "         y        "   "        y 
        
        scale_px = (double) (event->Particle(i)->P4()->Px() / event->Particle(fd)->P4()->Px());
        scale_py = (double) (event->Particle(i)->P4()->Py() / event->Particle(fd)->P4()->Py());
        scale_pz = (double) (event->Particle(i)->P4()->Pz() / event->Particle(fd)->P4()->Pz());
        
        double epsilon_parrallel = 1E-6;
        
        if((TMath::Abs(scale_px-scale_py) < epsilon_parrallel) &&
           (TMath::Abs(scale_py-scale_pz) < epsilon_parrallel) &&
           (TMath::Abs(scale_pz-scale_px) < epsilon_parrallel) ){

          hadron_fate = kIHAFtUndefined;
        }
        
        // at this point we know that there was no energy change except for the binding energy 
        else {hadron_fate = kIHAFtElas;
        }
        
      }

      // energy change not consistent with binding energy
      // three possibilities are that 1) it is inelastic event. 2) It is elastic but had energy 
      // so low that not possible to subtract all binding energy and now has E = m --> P = 0
      // or 3) same as 2 but was inelastic. 2) and 3) are indistinguishable with current information
      else if (delta_E >= epsilon_zero_energy_max){
        // there is ambiguity
        if(mod_Pfd_squared < epsilon_delta_P_squared)  {
           LOG("INukeUtils", pWARN) << "In the ambiguous case!";
          hadron_fate = kIHAFtElas; // this could be inelastic but setting to elastic due to greater abundance
        }
        // it is inelastic
        else { hadron_fate = kIHAFtInelas;}
      }
    } // end f single particle if
    //    // otherwise its inelastic   
    else {
      hadron_fate = kIHAFtElas;
    }
      
  return hadron_fate;
}
//____________________________________________________________________________

