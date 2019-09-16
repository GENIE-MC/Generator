//____________________________________________________________________________
/*
 Copyright (c) 2003-2017, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - Nov 30, 2008

 For documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Nov 30, 2008 - CA
   Added in version 2.5.1

*/
//____________________________________________________________________________


#include "Framework/Messenger/Messenger.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepUtils.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"

//____________________________________________________________________________
int genie::utils::ghep::NeutReactionCode(const GHepRecord * event)
{
// Ryan Terri, Yoshinari Hayato, Costas Andreopoulos
//
// A description of NEUT event types can be seen here: 
// http://t2k.phy.duke.edu/bin/view/Main/NeutModes
// Any extension used here has been agreed with SK (Hayato et al)
//
// Updated for NEUT 5.4.0 by Christophe Bronner

  
  if(!event) {
    LOG("GHepUtils", pWARN) << "Null event!";
    return 0;
  }

  int evtype = 0;

  Interaction * interaction = event->Summary();  
          
  const ProcessInfo &  proc = interaction->ProcInfo();
  const InitialState & init = interaction->InitState();
  const XclsTag &      xcls = interaction->ExclTag();
  const Kinematics &   kine = interaction->Kine();
  const Target &       tgt  = init.Tgt();
        
  bool is_cc    = proc.IsWeakCC();
  bool is_nc    = proc.IsWeakNC();
  bool is_charm = xcls.IsCharmEvent();
  bool is_qel   = proc.IsQuasiElastic();
  bool is_dis   = proc.IsDeepInelastic();
  bool is_res   = proc.IsResonant();
  bool is_coh_pr = proc.IsCoherentProduction();
  bool is_ve    = proc.IsNuElectronElastic();
  bool is_mec   = proc.IsMEC();
  bool is_imd   = proc.IsInverseMuDecay();
  bool is_ask   = proc.IsSingleKaon();
  bool is_p     = tgt.HitNucIsSet() ? tgt.HitNucPdg()==kPdgProton  : false;
  bool is_n     = tgt.HitNucIsSet() ? tgt.HitNucPdg()==kPdgNeutron : false;
  bool is_nu    = pdg::IsNeutrino    (init.ProbePdg());
  bool is_nubar = pdg::IsAntiNeutrino(init.ProbePdg());
  bool W_gt_2   = kine.KVSet(kKVW) ?  (kine.W() > 2.0) : false;
        
  // (quasi-)elastic, nc+cc, nu+nubar
  //
  if      (is_qel && !is_charm && is_cc && is_nu           ) evtype =   1;
  else if (is_qel && !is_charm && is_nc && is_nu && is_p   ) evtype =  51;   
  else if (is_qel && !is_charm && is_nc && is_nu && is_n   ) evtype =  52;
  else if (is_qel && !is_charm && is_cc && is_nubar        ) evtype =  -1;
  else if (is_qel && !is_charm && is_nc && is_nubar && is_p) evtype = -51;
  else if (is_qel && !is_charm && is_nc && is_nubar && is_n) evtype = -52;
           
  // MEC - only CC implemented in NEUT
  else if      (is_mec && !is_charm && is_cc && is_nu           ) evtype =   2;
  else if      (is_mec && !is_charm && is_cc && is_nubar        ) evtype =  -2;

  // quasi-elastic charm production
  // Part of the DIS W>2GeV mode in NEUT - CB
  else if (is_qel && is_charm && is_cc && is_nu    ) evtype =   26;
  else if (is_qel && is_charm && is_cc && is_nubar ) evtype =  -26;

  // inverse mu- (tau-) decay and ve- elastic
  //Those modes don't actually exist in NEUT, 9 and 59 used as place holders
  else if ( is_imd ) evtype =  9;
  else if ( is_ve  ) evtype = 59;
               
  // coherent pi, nc+cc, nu+nubar
  //
  else if (is_coh_pr && is_cc && is_nu   ) evtype =  16;
  else if (is_coh_pr && is_cc && is_nubar) evtype = -16;
  else if (is_coh_pr && is_nc && is_nu   ) evtype =  36;
  else if (is_coh_pr && is_nc && is_nubar) evtype = -36;
                     
  // dis, W>2, nc+cc, nu+nubar
  // (charm DIS not simulated by NEUT, will bundle GENIE charm DIS into this category)
  //
  else if (is_dis && W_gt_2 && is_cc && is_nu   ) evtype =  26;
  else if (is_dis && W_gt_2 && is_nc && is_nu   ) evtype =  46;
  else if (is_dis && W_gt_2 && is_cc && is_nubar) evtype = -26; 
  else if (is_dis && W_gt_2 && is_nc && is_nubar) evtype = -46; 

  // resonance or dis with W < 2 GeV or single kaon
  //
  else if ( is_res || (is_dis && !W_gt_2) || is_ask ) {
        
     LOG("GHepUtils", pNOTICE) << "Current event is RES or DIS with W<2";
        
     // check the number of pions and nucleons in the primary hadronic system
     // (_before_ intranuclear rescattering)
     //
     int nn=0, np=0, npi0=0, npip=0, npim=0, nKp=0, nKm=0, nK0=0, neta=0, nlambda=0, ngamma=0;
     bool nuclear_target = init.Tgt().IsNucleus();

     TIter event_iter(event);
     GHepParticle * p = 0;

     while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) )
     {
         GHepStatus_t ghep_ist = (GHepStatus_t) p->Status();
         int ghep_pdgc    = p->Pdg();
         int ghep_fm      = p->FirstMother();
         int ghep_fmpdgc  = (ghep_fm==-1) ? 0 : event->Particle(ghep_fm)->Pdg();
        
         // For nuclear targets use hadrons marked as 'hadron in the nucleus'
         // which are the ones passed in the intranuclear rescattering
         // For free nucleon targets use particles marked as 'final state'
         // but make an exception for decayed pi0's,eta's (count them and not their daughters)

         bool decayed         = (ghep_ist==kIStDecayedState && (ghep_pdgc==kPdgPi0 || ghep_pdgc==kPdgEta));
         bool parent_included = (ghep_fmpdgc==kPdgPi0 || ghep_fmpdgc==kPdgEta);

         bool count_it =
               ( nuclear_target && ghep_ist==kIStHadronInTheNucleus) ||
               (!nuclear_target && decayed) ||
               (!nuclear_target && ghep_ist==kIStStableFinalState && !parent_included);

         if(!count_it) continue;
                
         if(ghep_pdgc == kPdgProton )    np++;            // p
         if(ghep_pdgc == kPdgNeutron)    nn++;            // n
         if(ghep_pdgc == kPdgPiP)        npip++;          // pi+
         if(ghep_pdgc == kPdgPiM)        npim++;          // pi-
         if(ghep_pdgc == kPdgPi0)        npi0++;          // pi0
         if(ghep_pdgc == kPdgEta)        neta++;          // eta0
         if(ghep_pdgc == kPdgKP)         nKp++;           // K+
         if(ghep_pdgc == kPdgKM)         nKm++;           // K-
         if(ghep_pdgc == kPdgK0)         nK0++;           // K0
         if(ghep_pdgc == kPdgAntiK0)     nK0++;           // K0
         if(ghep_pdgc == kPdgLambda)     nlambda++;       // Lamda
         if(ghep_pdgc == kPdgAntiLambda) nlambda++;       // Lamda
         if(ghep_pdgc == kPdgGamma)      ngamma++;        // photon
     }
     LOG("GHepUtils", pNOTICE)
           << "Num of primary particles: \n p = " << np << ", n = " << nn
              << ", pi+ = " << npip << ", pi- = " << npim << ", pi0 = " << npi0 
              << ", eta = " << neta 
              << ", K+ = " << nKp << ", K- = " << nKm << ", K0 = " << nK0 
              << ", Lambda's = " << nlambda
              << ", gamma's = " << ngamma;
              
     int nnuc = np + nn;
     int npi  = npi0 + npip + npim;
     int nK   = nK0 + nKp + nKm;
     int neKL = neta + nK + nlambda;
              
     bool is_radiative_dec = (nnuc==1) && (npi==0) && (ngamma==1);
              
     //
     // single gamma from resonances
     //
              
     if      (is_res && is_nu    && is_cc && is_n && is_radiative_dec) evtype =  17;
     else if (is_res && is_nu    && is_nc && is_n && is_radiative_dec) evtype =  38;
     else if (is_res && is_nu    && is_nc && is_p && is_radiative_dec) evtype =  39;
               
     else if (is_res && is_nubar && is_cc && is_p && is_radiative_dec) evtype = -17;
     else if (is_res && is_nubar && is_nc && is_n && is_radiative_dec) evtype = -38;
     else if (is_res && is_nubar && is_nc && is_p && is_radiative_dec) evtype = -39;
               
     //
     // single pi (res + non-res bkg)
     //
            
     // nu CC
     else if (is_nu    && is_cc && is_p && np==1 && nn==0 && npip==1 && npim==0 && npi0==0 && neKL==0) evtype =  11;
     else if (is_nu    && is_cc && is_n && np==1 && nn==0 && npip==0 && npim==0 && npi0==1 && neKL==0) evtype =  12;
     else if (is_nu    && is_cc && is_n && np==0 && nn==1 && npip==1 && npim==0 && npi0==0 && neKL==0) evtype =  13;
           
     // nu NC
     else if (is_nu    && is_nc && is_n && np==0 && nn==1 && npip==0 && npim==0 && npi0==1 && neKL==0) evtype =  31;
     else if (is_nu    && is_nc && is_p && np==1 && nn==0 && npip==0 && npim==0 && npi0==1 && neKL==0) evtype =  32;
     else if (is_nu    && is_nc && is_n && np==1 && nn==0 && npip==0 && npim==1 && npi0==0 && neKL==0) evtype =  33;
     else if (is_nu    && is_nc && is_p && np==0 && nn==1 && npip==1 && npim==0 && npi0==0 && neKL==0) evtype =  34;
           
     //nubar CC
     else if (is_nubar && is_cc && is_n && np==0 && nn==1 && npip==0 && npim==1 && npi0==0 && neKL==0) evtype = -11;
     else if (is_nubar && is_cc && is_p && np==0 && nn==1 && npip==0 && npim==0 && npi0==1 && neKL==0) evtype = -12;
     else if (is_nubar && is_cc && is_p && np==1 && nn==0 && npip==0 && npim==1 && npi0==0 && neKL==0) evtype = -13;
                     
     //nubar NC
     else if (is_nubar && is_nc && is_n && np==0 && nn==1 && npip==0 && npim==0 && npi0==1 && neKL==0) evtype = -31;
     else if (is_nubar && is_nc && is_p && np==1 && nn==0 && npip==0 && npim==0 && npi0==1 && neKL==0) evtype = -32;
     else if (is_nubar && is_nc && is_n && np==1 && nn==0 && npip==0 && npim==1 && npi0==0 && neKL==0) evtype = -33;
     else if (is_nubar && is_nc && is_p && np==0 && nn==1 && npip==1 && npim==0 && npi0==0 && neKL==0) evtype = -34;
              
     //
     // single eta from res
     //
              
     else if (is_res &&  is_nu    && is_cc && is_n && np==1 && nn==0 && npi==0 && nK==0 && nlambda==0 && neta==1) evtype =  22;
     else if (is_res &&  is_nu    && is_nc && is_n && np==0 && nn==1 && npi==0 && nK==0 && nlambda==0 && neta==1) evtype =  42;
     else if (is_res &&  is_nu    && is_nc && is_p && np==1 && nn==0 && npi==0 && nK==0 && nlambda==0 && neta==1) evtype =  43;
              
     else if (is_res &&  is_nubar && is_cc && is_p && np==0 && nn==1 && npi==0 && nK==0 && nlambda==0 && neta==1) evtype = -22;
     else if (is_res &&  is_nubar && is_nc && is_n && np==0 && nn==1 && npi==0 && nK==0 && nlambda==0 && neta==1) evtype = -42;
     else if (is_res &&  is_nubar && is_nc && is_p && np==1 && nn==0 && npi==0 && nK==0 && nlambda==0 && neta==1) evtype = -43;
              
     //
     // single K from res (dS=0)
     //
              
     else if (is_res &&  is_nu    && is_cc && is_n && nnuc==0 && npi==0 && nK==1 && nlambda==1 && neta==0) evtype =  23;
     else if (is_res &&  is_nu    && is_nc && is_n && nnuc==0 && npi==0 && nK==1 && nlambda==1 && neta==0) evtype =  44;
     else if (is_res &&  is_nu    && is_nc && is_p && nnuc==0 && npi==0 && nK==1 && nlambda==1 && neta==0) evtype =  45;
              
     else if (is_res &&  is_nubar && is_cc && is_p && nnuc==0 && npi==0 && nK==1 && nlambda==1 && neta==0) evtype = -23;
     else if (is_res &&  is_nubar && is_nc && is_n && nnuc==0 && npi==0 && nK==1 && nlambda==1 && neta==0) evtype = -44;
     else if (is_res &&  is_nubar && is_nc && is_p && nnuc==0 && npi==0 && nK==1 && nlambda==1 && neta==0) evtype = -45;
     
     //
     // single K from AtharSingleKaon (dS=1)
     //
     //Those modes are assigned but not used (xsec=0) in NEUT
     else if (is_ask &&  is_nu && is_cc && is_n && nn==1 && np==0 && nKp==1 && neKL==1) evtype =  18;
     else if (is_ask &&  is_nu && is_cc && is_n && nn==0 && np==1 && nK0==1 && neKL==1) evtype =  19;
     else if (is_ask &&  is_nu && is_cc && is_p && nn==0 && np==1 && nKp==1 && neKL==1) evtype =  20;
     

     // antineutrino modes not yet implemented
     //else if (is_ask &&  is_nubar && is_cc && is_n && nn==1 && np==0 && nKp==1 && neKL==1) evtype = -18;
     //else if (is_ask &&  is_nubar && is_cc && is_n && nn==0 && np==1 && nK0==1 && neKL==1) evtype = -19;
     //else if (is_ask &&  is_nubar && is_cc && is_p && nn==0 && np==1 && nKp==1 && neKL==1) evtype = -20;

     //
     // multi-pi (res or dis (W<2GeV)
     //
              
     else if (is_nu    && is_cc && npi>1) evtype =  21;
     else if (is_nu    && is_nc && npi>1) evtype =  41;
     else if (is_nubar && is_cc && npi>1) evtype = -21;
     else if (is_nubar && is_nc && npi>1) evtype = -41;

     //
     // rare final state for RES or low-W (<2GeV) DIS events
     // (eg K0\bar{K0} final states, N0(1720) -> Sigma- K+ res decays, etc)
     // bundled-in with multi-pi
     //
     else {              
       LOG("GHepUtils", pWARN)
         << "Rare RES/low-W DIS final state: Bundled-in with multi-pi events";

             if (is_nu    && is_cc) evtype =  21;
        else if (is_nu    && is_nc) evtype =  41;
        else if (is_nubar && is_cc) evtype = -21;
        else if (is_nubar && is_nc) evtype = -41;
     }
  }

  return evtype;
}
//____________________________________________________________________________
int genie::utils::ghep::NuanceReactionCode(const GHepRecord * event)
{
// Josh Spitz, Costas Andreopoulos
//  
  if(!event) {
    LOG("GHepUtils", pWARN) << "Null event!";
    return 0;
  }

  int evtype = 0;

  Interaction * interaction = event->Summary();  

  const ProcessInfo &  proc = interaction->ProcInfo();
  const InitialState & init = interaction->InitState();
  if      (proc.IsQuasiElastic()   && proc.IsWeakCC()) evtype =  1;
  else if (proc.IsQuasiElastic()   && proc.IsWeakNC()) evtype =  2;
  else if (proc.IsDeepInelastic()  && proc.IsWeakCC()) evtype = 91;
  else if (proc.IsDeepInelastic()  && proc.IsWeakNC()) evtype = 92;
  else if (proc.IsCoherentProduction() && proc.IsWeakNC()) evtype = 96;
  else if (proc.IsCoherentProduction() && proc.IsWeakCC()) evtype = 97;
  else if (proc.IsNuElectronElastic())                 evtype = 98;
  else if (proc.IsInverseMuDecay())                    evtype = 99;
  else if (proc.IsResonant()) {
     int nn=0, np=0, npi0=0, npip=0, npim=0; 
     bool nuclear_target = init.Tgt().IsNucleus();
     GHepStatus_t matched_ist = (nuclear_target) ?
                kIStHadronInTheNucleus : kIStStableFinalState;

     TIter event_iter(event);
     GHepParticle * p = 0;
  
     while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) )
     {
         GHepStatus_t ghep_ist = (GHepStatus_t) p->Status();
         if(ghep_ist != matched_ist) continue;

         int ghep_pdgc = p->Pdg();
         if(ghep_pdgc == kPdgProton ) np++;
         if(ghep_pdgc == kPdgNeutron) nn++;
         if(ghep_pdgc == kPdgPi0)     npi0++;
         if(ghep_pdgc == kPdgPiP)     npip++;
         if(ghep_pdgc == kPdgPiM)     npim++;
     }
     if(proc.IsWeakCC() && init.IsNuP()) {
         // v p -> l- p pi+
         if(np==1 && nn==0 && npip==1 && npi0==0 && npim==0) evtype = 3;
     }
     if(proc.IsWeakCC() && init.IsNuN()) {
         // v n -> l- p pi0  
         if(np==1 && nn==0 && npip==0 && npi0==1 && npim==0) evtype = 4;
         // v n -> l- n pi+
         if(np==0 && nn==1 && npip==1 && npi0==0 && npim==0) evtype = 5;
     }
     if(proc.IsWeakNC() && init.IsNuP()) {
         // v p -> v p pi0
         if(np==1 && nn==0 && npip==0 && npi0==1 && npim==0) evtype = 6;
         // v p -> v n pi+
         if(np==0 && nn==1 && npip==1 && npi0==0 && npim==0) evtype = 7;
     }
     if(proc.IsWeakNC() && init.IsNuN()) {
         // v n -> v n pi0
         if(np==0 && nn==1 && npip==0 && npi0==1 && npim==0) evtype = 8;
         // v n -> v p pi-
         if(np==1 && nn==0 && npip==0 && npi0==0 && npim==1) evtype = 9;
     }
     if(proc.IsWeakCC() && init.IsNuBarN()) {
         // \bar{v} n -> l+ n pi-
         if(np==1 && nn==0 && npip==1 && npi0==0 && npim==0) evtype = 10;
     }
     if(proc.IsWeakCC() && init.IsNuBarP()) {
         // \bar{v} p -> l+ n pi0
         if(np==1 && nn==0 && npip==0 && npi0==1 && npim==0) evtype = 11;
         // \bar{v} p -> l+ p pi-
         if(np==0 && nn==1 && npip==1 && npi0==0 && npim==0) evtype = 12;
     }
     if(proc.IsWeakNC() && init.IsNuBarP()) {
         // \bar{v} p -> \bar{v} p pi0
         if(np==1 && nn==0 && npip==0 && npi0==1 && npim==0) evtype = 13;
         // \bar{v} p -> \bar{v} n pi+
         if(np==0 && nn==1 && npip==1 && npi0==0 && npim==0) evtype = 14;
      }
      if(proc.IsWeakNC() && init.IsNuBarN()) {
         // \bar{v} n -> \bar{v} n pi0
         if(np==0 && nn==1 && npip==0 && npi0==1 && npim==0) evtype = 15;
         // \bar{v} n -> \bar{v} p pi-
         if(np==1 && nn==0 && npip==0 && npi0==0 && npim==1) evtype = 16;
      }
  }

  return evtype;
}
//____________________________________________________________________________
