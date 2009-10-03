//____________________________________________________________________________
/*
 Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Jim Dobson <j.dobson07 \at imperial.ac.uk>
         Imperial College London

         Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

 For documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Mar 04, 2009 - JD
   Was first added in v2.5.1. Adapted from the T2K GENIE reweighting tool.
 @ Mar 05, 2009 - CA
   Modified ReconstructHadronFateHA() to work with hadron+A event files in
   addition to neutrino event files.
 @ Sep 10, 2009 - CA
   Added MeanFreePath(), Dist2Exit(), Dist2ExitMFP()
 @ Sep 30, 2009 - CA
   Added StepParticle() from Intranuke.cxx
 @ Oct 02, 2009 - CA
   Added test MeanFreePath_Delta().
*/
//____________________________________________________________________________

#include <TLorentzVector.h>
#include <TMath.h>

#include "Algorithm/AlgConfigPool.h"
#include "Conventions/Units.h"
#include "Conventions/GBuild.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepParticle.h"
#include "HadronTransport/INukeUtils.h"
#include "HadronTransport/INukeHadroData.h"
#include "Messenger/Messenger.h"
#include "Numerical/Spline.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGCodes.h"
#include "Registry/Registry.h"
#include "Utils/NuclearUtils.h"

using namespace genie;

//____________________________________________________________________________
double genie::utils::intranuke::MeanFreePath(
   int pdgc, const TLorentzVector & x4, const TLorentzVector & p4, 
   double A, double nRpi, double nRnuc)
{
// Calculate the mean free path (in fm) for a pions and nucleons in a nucleus
//
// Inputs
//  pdgc : Hadron PDG code
//  x4   : Hadron 4-position in the nucleus coordinate system (units: fm)
//  p4   : Hadron 4-momentum (units: GeV)
//  A    : Nucleus atomic mass number
//  nRpi : Controls the pion ring size in terms of de-Broglie wavelengths
//  nRnuc: Controls the nuclepn ring size in terms of de-Broglie wavelengths
//
  bool is_pion    = pdgc == kPdgPiP || pdgc == kPdgPi0 || pdgc == kPdgPiM;
  bool is_nucleon = pdgc == kPdgProton || pdgc == kPdgNeutron;

  if(!is_pion && !is_nucleon) return 0.;
        
  // before getting the nuclear density at the current position
  // check whether the nucleus has to become larger by const times the
  // de Broglie wavelength -- that is somewhat empirical, but this
  // is what is needed to get piA total cross sections right.
  // The ring size is different for light nuclei (using gaus density) /
  // heavy nuclei (using woods-saxon density).
  // The ring size is different for pions / nucleons.
  //
  double momentum = p4.Vect().Mag(); // hadron momentum in GeV
  double ring = (momentum>0) ? 1.240/momentum : 0; // de-Broglie wavelength
 
  if(A<=20) { ring /= 2.; }
 
  if      (is_pion   ) { ring *= nRpi;  }
  else if (is_nucleon) { ring *= nRnuc; }

  // get the nuclear density at the current position
  double rnow = x4.Vect().Mag(); 
  double rho  = A * utils::nuclear::Density(rnow,A,ring);

  // the hadron+nucleon cross section will be evaluated within the range
  // of the input spline and assumed to be const outside that range
  //
  double ke = (p4.Energy() - p4.M()) / units::MeV;  // kinetic energy in MeV
  ke = TMath::Max(INukeHadroData::fMinKinEnergy,   ke);
  ke = TMath::Min(INukeHadroData::fMaxKinEnergyHN, ke);

  // get total xsection for the incident hadron at its current
  // kinetic energy
  double sigtot = 0;

  INukeHadroData * fHadroData = INukeHadroData::Instance();

  if (pdgc == kPdgPiP)
      sigtot = fHadroData -> XSecPipN_Tot() -> Evaluate(ke);
  else if (pdgc == kPdgPi0)
      sigtot = fHadroData -> XSecPi0N_Tot() -> Evaluate(ke);
  else if (pdgc == kPdgPiM)
      sigtot = fHadroData -> XSecPimN_Tot() -> Evaluate(ke);
  else if (pdgc == kPdgProton)
      sigtot = fHadroData -> XSecPN_Tot()   -> Evaluate(ke);
  else if (pdgc == kPdgNeutron)
      sigtot = fHadroData -> XSecNN_Tot()   -> Evaluate(ke);
  else {
     return 0;
  }

  // the xsection splines in INukeHadroData return the hadron x-section in
  // mb -> convert to fm^2
  sigtot *= (units::mb / units::fm2);

  // compute the mean free path
  double lamda = 1. / (rho * sigtot);

  // exits if lamda is InF (if cross section is 0)
  if( ! TMath::Finite(lamda) ) {
     return -1;
  }

  return lamda;
}
//____________________________________________________________________________
double genie::utils::intranuke::MeanFreePath_Delta(
   int pdgc, const TLorentzVector & x4, const TLorentzVector & p4, double A)
{
//
// **test**
//

// Calculate the mean free path (in fm) for Delta's in a nucleus
//
// Inputs
//  pdgc : Hadron PDG code
//  x4   : Hadron 4-position in the nucleus coordinate system (units: fm)
//  p4   : Hadron 4-momentum (units: GeV)
//  A    : Nucleus atomic mass number
//
  bool is_deltapp = (pdgc==kPdgP33m1232_DeltaPP);
  if(!is_deltapp) return 0.;
        
  // get the nuclear density at the current position
  double rnow = x4.Vect().Mag(); 
  double rho  = A * utils::nuclear::Density(rnow,A);

  // the Delta+N->N+N cross section will be evaluated within the range
  // of the input spline and assumed to be const outside that range
  double ke = (p4.Energy() - p4.M()) / units::MeV;  // kinetic energy in MeV
  ke = TMath::Max(1500., ke);
  ke = TMath::Min(   0., ke);

  // get the Delta+N->N+N cross section
  double sig = 0;
  if      (ke< 500) sig=20;
  else if (ke<1000) sig=40;
  else              sig=30;

  // value is in mb -> convert to fm^2
  sig *= (units::mb / units::fm2);

  // compute the mean free path
  double lamda = 1. / (rho * sig);

  // exits if lamda is InF (if cross section is 0)
  if( ! TMath::Finite(lamda) ) {
     return -1;
  }

  return lamda;
}
//____________________________________________________________________________
double genie::utils::intranuke::ProbSurvival(
  int pdgc, const TLorentzVector & x4, const TLorentzVector & p4, double A,
  double mfp_scale_factor, double nRpi, double nRnuc, double NR, double R0)
{
// Calculate the survival probability for a hadron inside a nucleus
//
// Inputs
//  pdgc : Hadron PDG code
//  x4 : Hadron 4-position in the nucleus coordinate system (units: fm)
//  p4 : Hadron 4-momentum (units: GeV)
//  A : Nucleus atomic mass number
//  mfp_scale_factor: Tweaks the mean free path (mfp -> mfp*scale). Def: 1.0
//  nRpi: Controls the pion ring size in terms of de-Broglie wavelengths
//  nRnuc: Controls the nuclepn ring size in terms of de-Broglie wavelengths
//  NR: How far away to track the hadron, in terms of the corresponding 
//      nuclear radius. Def: 3
//  R0: R0 in R=R0*A^1/3 (units:fm). Def. 1.4

   double prob = 1.0;

   double step = 0.05; // fermi
   double R    = NR * R0 * TMath::Power(A, 1./3.);

   TVector3 dr3 = p4.Vect().Unit();  // unit vector along its direction
   TLorentzVector dr4(dr3,0);
     
   TLorentzVector x4_curr(x4); // current position

   while(1) {
     double rnow  = x4_curr.Vect().Mag();
     if (rnow > (R+step)) break;

     x4_curr += (step*dr4);
     rnow = x4_curr.Vect().Mag();
     double mfp = 
       genie::utils::intranuke::MeanFreePath(pdgc,x4_curr,p4,A,nRpi,nRnuc);
     mfp *= mfp_scale_factor;

     double dprob = (mfp>0) ? TMath::Exp(-step/mfp) : 0.;
     prob*=dprob;
   }

   return prob;
}
//____________________________________________________________________________
INukeFateHA_t genie::utils::intranuke::ReconstructHadronFateHA(
   GHepRecord * event, int i, bool hA_mode)
{
// Reconstruct the INTRANUKE/hA model fate for the hadron at position i.
//
// Returns one of:
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
double genie::utils::intranuke::Dist2Exit(
   const TLorentzVector & x4, const TLorentzVector & p4, 
   double A, double NR, double R0)
{
// Calculate distance within a nucleus (units: fm) before we stop tracking 
// the hadron.
// See previous functions for a description of inputs.
//
   double R    = NR * R0 * TMath::Power(A, 1./3.);
   double step = 0.05; // fermi
   
   TVector3 dr3 = p4.Vect().Unit();  // unit vector along its direction
   TLorentzVector dr4(dr3,0);

   TLorentzVector x4_curr(x4); // current position
  
   double d=0;
   while(1) {
        double rnow  = x4_curr.Vect().Mag();
        x4_curr += (step*dr4);
        d+=step;
        rnow = x4_curr.Vect().Mag();
        if (rnow > R) break;
   }
   return d;
}
//____________________________________________________________________________
double genie::utils::intranuke::Dist2ExitMFP(
   int pdgc, const TLorentzVector & x4, const TLorentzVector & p4, 
   double A, double NR, double R0)
{
// Calculate distance within a nucleus (expressed in terms of 'mean free 
// paths') before we stop tracking the hadron.
// See previous functions for a description of inputs.
//

// distance before exiting in mean free path lengths
// 
   double R    = NR * R0 * TMath::Power(A, 1./3.);
   double step = 0.05; // fermi

   TVector3 dr3 = p4.Vect().Unit();  // unit vector along its direction
   TLorentzVector dr4(dr3,0);

   TLorentzVector x4_curr(x4); // current position

   double d=0;
   double d_mfp=0;
   while(1) {
        double rnow  = x4_curr.Vect().Mag();
        x4_curr += (step*dr4);
        d+=step;
        rnow = x4_curr.Vect().Mag();

        double lambda = genie::utils::intranuke::MeanFreePath(pdgc,x4_curr,p4,A);
        d_mfp += (step/lambda);

        if (rnow > R) break;
   }
   return d_mfp;
}
//____________________________________________________________________________
void genie::utils::intranuke::StepParticle(
           GHepParticle * p, double step, double nuclear_radius)
{
// Steps a particle starting from its current position (in fm) and moving
// along the direction of its current momentum by the input step (in fm).
// The particle is stepped in a straight line.
// If a nuclear radius is set then the following check is performed:
// If the step is too large and takes the the particle far away from the
// nucleus then its position is scaled back so that the escaped particles are
// always within a ~1fm from the "outer nucleus surface"

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("INukeUtils", pDEBUG)
      << "Stepping particle [" << p->Name() << "] by dr = " << step << " fm";
#endif

  // Step particle
  TVector3 dr = p->P4()->Vect().Unit();      // unit vector along its direction
  dr.SetMag(step);                           // spatial step size
  double dt = 0;                             // temporal step:
  TLorentzVector dx4(dr,dt);                 // 4-vector step
  TLorentzVector x4new = *(p->X4()) + dx4;   // new position
  
  if(nuclear_radius > 0.) {
    // Check position against nuclear boundary. If the particle was stepped
    // too far away outside the nuclear boundary bring it back to within 
    // 1fm from that boundary
    double epsilon = 1; // fm  
    double r       = x4new.Vect().Mag(); // fm
    double rmax    = nuclear_radius+epsilon;
    if(r > rmax) {
       LOG("INukeUtils", pINFO)
         << "Particle was stepped too far away (r = " << r << " fm)";
       LOG("INukeUtils", pINFO)
         << "Placing it " << epsilon 
         << " fm outside the nucleus (r' = " << rmax << " fm)";
       double scale = rmax/r;
       x4new *= scale;
    }//r>rmax
  }//nucl radius set

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("INukeUtils", pDEBUG)
      << "\n Init direction = " << print::Vec3AsString(&dr)
      << "\n Init position (in fm,nsec) = " << print::X4AsString(p->X4())
      << "\n Fin  position (in fm,nsec) = " << print::X4AsString(&x4new);
#endif

  p->SetPosition(x4new);
}
//___________________________________________________________________________
