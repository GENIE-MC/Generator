//____________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory
*/
//____________________________________________________________________________
#include <cmath>

#include <cmath>

#include <TClonesArray.h>
#include <TDecayChannel.h>
#include <TMath.h>

#include "Framework/GHEP/GHepFlags.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
// #include "Framework/ParticleData/BaryonResUtils.h" // TODO: do I need to create a similar file to this?
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Utils/PrintUtils.h"
#include "Framework/Utils/StringUtils.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/Decay/DarkSectorDecayer.h"
#include "Math/GSLMinimizer.h"
#include <Math/WrappedParamFunction.h>

using namespace genie;
using namespace genie::controls;
using namespace genie::constants;
//____________________________________________________________________________
DarkSectorDecayer::DarkSectorDecayer() :
Decayer("genie::DarkSectorDecayer")
{
  this->Initialize();
}
//____________________________________________________________________________
DarkSectorDecayer::DarkSectorDecayer(string config) :
Decayer("genie::DarkSectorDecayer", config)
{
  this->Initialize();
}
//____________________________________________________________________________
DarkSectorDecayer::~DarkSectorDecayer()
{

  for ( unsigned int i = 0; i < fRParams.size() ; ++i ) {
    delete fRParams[i] ;
  }

}
//____________________________________________________________________________
void DarkSectorDecayer::ProcessEventRecord(GHepRecord * event) const
{
  LOG("DarkSectorDecay", pINFO)
    << "Running dark sector decayer ";
    // << ((fRunBefHadroTransp) ? "*before*" : "*after*") << " FSI";

  // Loop over particles, find unstable ones and decay them
  TObjArrayIter piter(event);
  GHepParticle * p = 0;
  int ipos = -1;

  while( (p = (GHepParticle *) piter.Next()) ) {

    ipos++;
    LOG("DarkSectorDecay", pDEBUG) << "Checking: " << p->Name();

    int pdg_code = p->Pdg();
    GHepStatus_t status_code = p->Status();

    //    std::cout << "Decaing particle " << ipos << " with PDG " << pdg_code << std::endl ; 

    if(!this->IsHandled  (pdg_code)) continue;
    if(!this->ToBeDecayed(pdg_code, status_code)) continue;

    LOG("DarkSectorDecay", pNOTICE)
          << "Decaying unstable particle: " << p->Name();


    // TODO: here check pdg of particle
    // if it's dark mediator call DecayDarkMediator
    // if it's dark neutrino call DecayDarkNeutrino
    
    //TODO: this block below
    // bool decayed = this->Decay(ipos, event);
    // if ( ! decayed ) {
    //   LOG("DarkSectorDecay", pWARN) << "Dark stuff not decayed!" ;
    //   LOG("DarkSectorDecay", pWARN) << "Quitting the current event generation thread" ;

    //   event -> EventFlags() -> SetBitNumber(kHadroSysGenErr, true);

    //   genie::exceptions::EVGThreadException exception;
    //   exception.SetReason("Dark stuff not decayed"); // TODO
    //   exception.SwitchOnFastForward();
    //   throw exception;

    //   return ;
    // }

  } // loop over particles

  LOG("DarkSectorDecay", pNOTICE)
     << "Done finding & decaying dark sector particles";
}
//____________________________________________________________________________
// bool DarkSectorDecayer::Decay(
//   int decay_particle_id, GHepRecord * event) const
// {
//   // Reset previous decay weight
//   fWeight = 1.;

//   // Get particle to be decayed
//   GHepParticle * decay_particle = event->Particle(decay_particle_id);
//   if( ! decay_particle) {
//     LOG("DarkSectorDecay", pERROR)
//       << "Particle to be decayed not in the event record. Particle ud: " << decay_particle_id ; 
//     return false;
//   }

//   bool to_be_deleted ;

//   // Select a decay channel
//   TDecayChannel * selected_decay_channel =
//     this->SelectDecayChannel(decay_particle_id, event, to_be_deleted ) ;

//   if(!selected_decay_channel) {
//     LOG("DarkSectorDecay", pERROR)
//       << "No decay channel for particle " << decay_particle_id ;
//     LOG("DarkSectorDecay", pERROR)
//       << *event ;

//     return false;
//   }

//   // Decay the exclusive state and copy daughters in the event record
//   bool decayed = this->DecayExclusive(decay_particle_id, event, selected_decay_channel);

//   if ( to_be_deleted ) 
//     delete selected_decay_channel ; 

//   if ( ! decayed ) return false ;

//   // Update the event weight for each weighted particle decay
//   double weight = event->Weight() * fWeight;
//   event->SetWeight(weight);

//   // Mark input particle as a 'decayed state' & add its daughter links
//   decay_particle->SetStatus(kIStDecayedState);

//   return true;
// }
//____________________________________________________________________________
// double DarkSectorDecayer::Weight(void) const
// {
//   return fWeight;
// }
//____________________________________________________________________________
// double DarkSectorDecayer::FinalStateMass( TDecayChannel * ch ) const
// {
// // Computes the total mass of the final state system

//   double mass = 0;
//   unsigned int nd = ch->NDaughters();

//   for(unsigned int iparticle = 0; iparticle < nd; iparticle++) {

//      int daughter_code = ch->DaughterPdgCode(iparticle);
//      TParticlePDG * daughter = PDGLibrary::Instance()->Find(daughter_code);
//      assert(daughter);

//      double md = daughter->Mass();

//      // hack to switch off channels giving rare  occurences of |1114| that has
//      // no decay channels in the pdg table (08/2007)
//      if(TMath::Abs(daughter_code) == 1114) {
//          LOG("DarkSectorDecay", pNOTICE)
//             << "Disabling decay channel containing resonance 1114";;
//          md = 999999999;
//      }
//      mass += md;
//   }
//   return mass;
// }
//____________________________________________________________________________
void DarkSectorDecayer::Initialize(void) const
{

}
//____________________________________________________________________________
bool DarkSectorDecayer::IsHandled(int pdg_code) const
{
  // bool is_handled = utils::res::IsBaryonResonance(pdg_code); //TODO: add utils::darks::isdark
  bool is_handled = true;

  LOG("DarkSectorDecay", pDEBUG)
      << "Can decay particle with PDG code = " << pdg_code
      << "? " << ((is_handled)? "Yes" : "No");

  return is_handled ;
}
//____________________________________________________________________________
void DarkSectorDecayer::InhibitDecay(int pdgc, TDecayChannel * dc) const
{
  if(! this->IsHandled(pdgc)) return;
  if(!dc) return;

  //
  // Not implemented
  //
}
//____________________________________________________________________________
void DarkSectorDecayer::UnInhibitDecay(int pdgc, TDecayChannel * dc) const
{
  if(! this->IsHandled(pdgc)) return;
  if(!dc) return;

  //
  // Not implemented
  //
}
//____________________________________________________________________________
// TDecayChannel * DarkSectorDecayer::SelectDecayChannel( int decay_particle_id,
//                                                        GHepRecord * event,
//                                                        bool & to_be_deleted ) const
// {
//   // Get particle to be decayed
//   GHepParticle * decay_particle = event->Particle(decay_particle_id);
//   if(!decay_particle) return 0;

//   // Get the particle 4-momentum and PDG code
//   TLorentzVector decay_particle_p4 = *(decay_particle->P4());
//   int decay_particle_pdg_code = decay_particle->Pdg();

//   // Find the particle in the PDG library & quit if it does not exist
//   TParticlePDG * mother =
//      PDGLibrary::Instance()->Find(decay_particle_pdg_code);
//   if(!mother) {
//      LOG("DarkSectorDecay", pERROR)
//         << "\n *** The particle with PDG code = " << decay_particle_pdg_code
//          << " was not found in PDGLibrary";
//      return 0;
//   }
//   LOG("DarkSectorDecay", pINFO)
//     << "Decaying a " << mother->GetName()
//     << " with P4 = " << utils::print::P4AsString(&decay_particle_p4);

//   // Get the resonance mass W (generally different from the mass associated
//   // with the input PDG code, since the it is produced off the mass shell)
//   double W = decay_particle_p4.M();
//   LOG("DarkSectorDecay", pINFO) << "Available mass W = " << W;

//   // Get all decay channels
//   TObjArray * original_decay_list = mother->DecayList();

//   unsigned int nch = original_decay_list -> GetEntries();
//   LOG("DarkSectorDecay", pINFO)
//     << mother->GetName() << " has: " << nch << " decay channels";

//   // Loop over the decay channels (dc) and write down the branching
//   // ratios to be used for selecting a decay channel.
//   // Since a baryon resonance can be created at W < Mres, explicitly
//   // check and inhibit decay channels for which W > final-state-mass

//   bool has_evolved_brs = DarkSectorDecayer::HasEvolvedBRs( decay_particle_pdg_code ) ; 
  
//   TObjArray * actual_decay_list = nullptr ;

//   if ( has_evolved_brs ) {
//     actual_decay_list = EvolveDeltaBR( decay_particle_pdg_code, original_decay_list, W ) ;
//     if ( ! actual_decay_list ) return nullptr ;
//     nch = actual_decay_list -> GetEntries() ;
//     to_be_deleted = true ;
//   }
//   else {
//     actual_decay_list = original_decay_list ;
//     to_be_deleted = false ;
//   }

//   double BR[nch], tot_BR = 0;

//   for(unsigned int ich = 0; ich < nch; ich++) {

//     TDecayChannel * ch = (TDecayChannel *) actual_decay_list -> At(ich);

//     double fsmass = this->FinalStateMass(ch) ;
//     if ( fsmass < W ) {

//       SLOG("DarkSectorDecay", pDEBUG)
//                 << "Using channel: " << ich
//                 << " with final state mass = " << fsmass << " GeV";

//       tot_BR += ch->BranchingRatio();

//     } else {
//       SLOG("DarkSectorDecay", pINFO)
//                 << "Suppresing channel: " << ich
//                 << " with final state mass = " << fsmass << " GeV";
//     } // final state mass

//     BR[ich] = tot_BR;
//   }//channel loop

//   if( tot_BR <= 0. ) {
//     SLOG("DarkSectorDecay", pWARN)
//           << "None of the " << nch << " decay channels is available @ W = " << W;
//     return 0;
//   }

//   // Select a resonance based on the branching ratios
//   unsigned int ich = 0, sel_ich; // id of selected decay channel
//   RandomGen * rnd = RandomGen::Instance();
//   double x = tot_BR * rnd->RndDec().Rndm();
//   do {
//     sel_ich = ich;
//   } while (x > BR[ich++]);

//   TDecayChannel * sel_ch = (TDecayChannel *) actual_decay_list -> At(sel_ich);

//   LOG("DarkSectorDecay", pINFO)
//     << "Selected " << sel_ch->NDaughters() << "-particle decay channel ("
//     << sel_ich << ") has BR = " << sel_ch->BranchingRatio();

//   if ( has_evolved_brs ) {

//     for ( unsigned int i = 0; i < nch; ++i) {
//       if ( sel_ich != i ) delete actual_decay_list -> At(i);
//     }

//     delete actual_decay_list ;
//   }

//   return sel_ch;
// }
//____________________________________________________________________________
// bool DarkSectorDecayer::DecayExclusive(
//   int decay_particle_id, GHepRecord * event, TDecayChannel * ch) const
// {
//   // Find the particle to be decayed in the event record
//   GHepParticle * decay_particle = event->Particle(decay_particle_id);
//   if(!decay_particle) return false ;

//   // Get the decayed particle 4-momentum, 4-position and PDG code
//   TLorentzVector decay_particle_p4 = *(decay_particle->P4());
//   TLorentzVector decay_particle_x4 = *(decay_particle->X4());
//   int decay_particle_pdg_code = decay_particle->Pdg();

//   // Get the final state mass spectrum and the particle codes
//   // for the selected decay channel
//   unsigned int nd = ch->NDaughters();

//   int    pdgc[nd];
//   double mass[nd];

//   for(unsigned int iparticle = 0; iparticle < nd; iparticle++) {

//      int daughter_code = ch->DaughterPdgCode(iparticle);
//      TParticlePDG * daughter = PDGLibrary::Instance()->Find(daughter_code);
//      assert(daughter);

//      pdgc[iparticle] = daughter_code;
//      mass[iparticle] = daughter->Mass();

//      SLOG("DarkSectorDecay", pINFO)
//          << "+ daughter[" << iparticle << "]: "
//          << daughter->GetName() << " (pdg-code = "
//          << pdgc[iparticle] << ", mass = " << mass[iparticle] << ")";
//   }

//   // Check whether the expected channel is Delta->pion+nucleon
//   bool is_delta = (decay_particle_pdg_code == kPdgP33m1232_DeltaPP ||
//                    decay_particle_pdg_code == kPdgP33m1232_DeltaP  ||
//                    decay_particle_pdg_code == kPdgP33m1232_Delta0);

//   bool is_delta_N_Pi_decay = is_delta && this->IsPiNDecayChannel(ch);

//   // Decay the resonance using an N-body phase space generator
//   // The particle will be decayed in its rest frame and then the daughters
//   // will be boosted back to the original frame.

//   bool is_permitted = fPhaseSpaceGenerator.SetDecay(decay_particle_p4, nd, mass);
//   if ( ! is_permitted ) return false ;

//   // Find the maximum phase space decay weight
//   // double wmax = fPhaseSpaceGenerator.GetWtMax();
//   double wmax = -1;
//   for(int i=0; i<50; i++) {
//      double w = fPhaseSpaceGenerator.Generate();
//      wmax = TMath::Max(wmax,w);
//   }
//   assert(wmax>0);
//   LOG("DarkSectorDecay", pINFO)
//     << "Max phase space gen. weight for current decay: " << wmax;

//   if(fGenerateWeighted)
//   {
//     // Generating weighted decays
//     // Do a single draw of momentum 4-vectors and then stop,
//     // taking into account the weight for this particular draw
//     double w = fPhaseSpaceGenerator.Generate();
//     fWeight *= TMath::Max(w/wmax, 1.);
//   }
//   else
//   {
//     // Generating un-weighted decays
//     RandomGen * rnd = RandomGen::Instance();
//     wmax *= 2;
//     bool accept_decay=false;
//     unsigned int itry=0;

//     while(!accept_decay)
//     {
//       itry++;
//       assert(itry<kMaxUnweightDecayIterations);

//       double w  = fPhaseSpaceGenerator.Generate();
//       double gw = wmax * rnd->RndDec().Rndm();

//       if(w>wmax) {
//          LOG("DarkSectorDecay", pWARN)
//             << "Current decay weight = " << w << " > wmax = " << wmax;
//       }
//       LOG("DarkSectorDecay", pINFO)
//         << "Current decay weight = " << w << " / R = " << gw;

//       accept_decay = (gw<=w);

//       // Extra logic that applies only for Delta -> N + pi
//       if( accept_decay && is_delta_N_Pi_decay ) {

//         // We don't want the decay Delta -> Pi + N to be isotropic in the Delta referemce frame
//         // as generated by the simple phase space generator.
//         // In order to sample pion distribution according to W(Theta, phi) in the Delta resonance decay,
//         // we make use of the following.
//     	// Note that Theta and Phi are defined in a reference frame which depends on the whole event
//         // For each event generated from a Delta -> N + Pi event with Pi emitted at
//         // at angles Theta and Phi (in the Delta rest frame), attach a random number to it.
//         // then we calculate W(Theta, Phi).
//         // Each possible final state is used to evaluate (Theta, Phi),
//         // then a random number is thrown, if the the random number is higher than W(Theta, Phi) drop this event and go
//         // back to re-generate an event and repeat the procedure.
//         // Otherwise keep this event to the record.
//     	// For efficiency reasons the maxium of the function is Q2 dependent

//         // Locate the pion in the decay products
//         // at this point we already know that the pion is unique so the first pion we find is our pion
//         unsigned int pi_id = 0 ;

//         for(unsigned int iparticle = 0; iparticle < nd; iparticle++) {

//           if ( genie::pdg::IsPion( ch->DaughterPdgCode(iparticle) ) ) {
//             pi_id = iparticle ;
//             break ;
//           }
//         }//iparticle

//         TLorentzVector * lab_pion = fPhaseSpaceGenerator.GetDecay(pi_id);

//         accept_decay = AcceptPionDecay( *lab_pion, decay_particle_id, event) ;

//       } //if it is a Delta -> N + Pi

//     }//accept_decay

//   }//fGenerateWeighted

//   // A decay was generated - Copy to the event record

//   // Check whether the interaction is off a nuclear target or free nucleon
//   // Depending on whether this module is run before or after the hadron
//   // transport module it would affect the daughters status code
//   GHepParticle * target_nucleus = event->TargetNucleus();
//   bool in_nucleus = (target_nucleus!=0);

//   // Loop over daughter list and add corresponding GHepParticles
//   for(unsigned int id = 0; id < nd; id++) {

//      int daughter_pdg_code = pdgc[id];

//      TLorentzVector * daughter_p4 = fPhaseSpaceGenerator.GetDecay(id);
//      LOG("DarkSectorDecay", pDEBUG)
//         << "Adding daughter particle with PDG code = " << pdgc[id]
//         << " and mass = " << mass[id] << " GeV";

//      bool is_hadron = pdg::IsHadron(daughter_pdg_code);
//      bool hadron_in_nuc = (in_nucleus && is_hadron && fRunBefHadroTransp);

//      GHepStatus_t daughter_status_code = (hadron_in_nuc) ?
//           kIStHadronInTheNucleus : kIStStableFinalState;

//      event->AddParticle(
//        daughter_pdg_code, daughter_status_code, decay_particle_id,-1,-1,-1,
//        *daughter_p4, decay_particle_x4);
//   }

//   return true ;
// }
//____________________________________________________________________________
void DarkSectorDecayer::LoadConfig(void) {

  Decayer::LoadConfig() ;

  this -> GetParam( "FFScaling", fFFScaling ) ;
  this -> GetParam( "DarkMediatorMass", fDarkMediatorMass ) ;

  // this -> GetParamDef( "Delta-ThetaOnly", fDeltaThetaOnly, true ) ;

  // this -> GetParamDef( "DeltaDecayMaximumTolerance", fMaxTolerance, 0.0005 ) ;

  // bool invalid_configuration = false ;

  // // load R33 parameters
  // this -> GetParamVect( "Delta-R33", fR33 ) ; 

  // // load Q2 thresholds if necessary
  // if ( fR33.size() > 1 ) {
  //   this -> GetParamVect("Delta-Q2", fQ2Thresholds ) ;
  // }
  // else {
  //   fQ2Thresholds.clear() ;
  // }

  // // check if the number of Q2 matches the number of R33
  // if ( fQ2Thresholds.size() != fR33.size() -1 ) {
  //   invalid_configuration = true ;
  //   LOG("DarkSectorDecayer", pFATAL) << "Delta-Q2 and Delta-R33 have wrong sizes" ;
  //   LOG("DarkSectorDecayer", pFATAL) << "Delta-Q2  -> " << fQ2Thresholds.size() ;
  //   LOG("DarkSectorDecayer", pFATAL) << "Delta-R33 -> " << fR33.size() ;
  // }

  // if ( fDeltaThetaOnly ) {

  //   // check the parameters validity
  //   for ( unsigned int i = 0 ; i < fR33.size(); ++i ) {
  //     if ( (fR33[i] < -0.5) ||  (fR33[i] > 1.) ) {
  //   	invalid_configuration = true ;
  //   	LOG("DarkSectorDecayer", pFATAL) << "Delta-R33[" << i << "] out of valid range [-0.5 ; 1 ]" ;
  //   	LOG("DarkSectorDecayer", pFATAL) << "Delta-R33[" << i << "] = " << fR33[i] ;
  //   	break ;
  //     }
  //   }
    
  //   // set appropriate maxima
  //   fW_max.resize( fR33.size(), 0. ) ;
  //   for ( unsigned int i = 0 ; i < fR33.size(); ++i ) {
  //     fW_max[i] = ( fR33[i] < 0.5 ? 2. * ( 1. - fR33[i] ) : fR33[i] + 0.5 ) + fMaxTolerance ;
  //   }
    
  // } // Delta Theta Only
  
  // else {

  //   // load R31 and R3m1 parameters
  //   this -> GetParamVect( "Delta-R31", fR31 ) ;

  //   this -> GetParamVect( "Delta-R3m1", fR3m1 ) ;

  //   // check if they match the numbers of R33
  //   if ( (fR31.size() != fR33.size()) || (fR3m1.size() != fR33.size()) ) {
  //     LOG("DarkSectorDecayer", pFATAL) << "Delta-R31 or Delta-R3m1 sizes don't match Delta-R33" ;
  //     LOG("DarkSectorDecayer", pFATAL) << "R31: " << fR31.size()
  //                                      << ", R3m1: " << fR31.size()
  //                                      << " while R33: " << fR33.size() ;
  //     invalid_configuration = true ;
  //   }

  //   for ( unsigned int i = 0; i < fRParams.size() ; ++i ) {
  //     delete [] fRParams[i] ;
  //   }
  //   fRParams.clear() ; 
  //   // fill the container by Q2 bin instead of the parmaeter bin
  //   for ( unsigned int i = 0 ; i < fR33.size(); ++i ) {
  //     fRParams.push_back( new double[3]{ fR33[i], fR31[i], fR3m1[i] } ) ;
  //   }


  //   // check if they are physical
  //   fW_max.resize( fR33.size(), 0. ) ;
  //   for ( unsigned int i = 0 ; i < fR33.size(); ++i ) {
  //     double temp_min = FindDistributionExtrema( i, false ) ;
  //     if ( temp_min < 0. ) {
  //       LOG("DarkSectorDecayer", pFATAL) << "pion angular distribution minimum is negative for Q2 bin " << i ;
  //       invalid_configuration = true ;
  //       break ;
  //     }
      
  //     double temp_max = FindDistributionExtrema( i, true ) ;
  //     if ( temp_max <= 0. ) {
  //       LOG("DarkSectorDecayer", pFATAL) << "pion angular distribution maximum is non positive for Q2 bin " << i ;
  //       invalid_configuration = true ;
  //       break ;
  //     }
      
  //     fW_max[i] = temp_max + fMaxTolerance ; 
      
  //   }
    
  // }

  // if ( invalid_configuration ) {

  //   LOG("DarkSectorDecayer", pFATAL)
  //     << "Invalid configuration: Exiting" ;

  //   // From the FreeBSD Library Functions Manual
  //   //
  //   // EX_CONFIG (78)   Something was found in an unconfigured or miscon-
  //   //                  figured state.

  //   exit( 78 ) ;

  // }
  
}
