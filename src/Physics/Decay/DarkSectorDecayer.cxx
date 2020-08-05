//____________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Author: Iker de Icaza <i.de-icaza-astiz \at sussex.ac.uk>
 University of Sussex

 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory
*/
//____________________________________________________________________________

#include <cmath>
#include <numeric>

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
  EventRecordVisitorI("genie::DarkSectorDecayer")
{

}
//____________________________________________________________________________
DarkSectorDecayer::DarkSectorDecayer(string config) :
  EventRecordVisitorI("genie::DarkSectorDecayer", config)
{

}
//____________________________________________________________________________
DarkSectorDecayer::~DarkSectorDecayer()
{

}
//____________________________________________________________________________
void DarkSectorDecayer::ProcessEventRecord(GHepRecord * event) const
{
  LOG("DarkSectorDecayer", pINFO)
    << "Running dark sector decayer ";

  // Loop over particles, find unstable ones and decay them
  TObjArrayIter piter(event);
  GHepParticle * p = 0;
  int ipos = -1;

  while( (p = (GHepParticle *) piter.Next()) ) {
    ipos++;
    LOG("DarkSectorDecayer", pDEBUG) << "Checking: " << p->Name();

    if(!this->ToBeDecayed(*p)) continue;

    GHepParticle&  mother = *p; // change the name now we know it will decay
    std::vector<DarkSectorDecayer::DecayChannel> dcs;
    int pdg_code = mother.Pdg();
    if(pdg_code == kPdgDNuMediator){
      dcs = DarkMediatorDecayChannels(mother, event);
    }
    else if(pdg_code == kPdgDarkNeutrino ||
            pdg_code == kPdgAntiDarkNeutrino){
      dcs = DarkNeutrinoDecayChannels(mother, event);
    }
    double total_amplitude = std::accumulate(dcs.begin(), dcs.end(), 0.,
                                             [](double total,
                                                const DarkSectorDecayer::DecayChannel& dc)
                                               {return total + dc.second;});

    int dcid = SelectDecayChannel(mother, event, dcs, total_amplitude);
    std::vector<GHepParticle> daughters = Decay(mother, dcs[dcid].first);
    SetSpaceTime(daughters, mother, total_amplitude);

    for(auto & daughter: daughters){
      daughter.SetFirstMother(ipos);
      event->AddParticle(daughter);
    }
  }
  
    // int pdg_code = p->Pdg();
    // GHepStatus_t status_code = p->Status();

    // //    std::cout << "Decaing particle " << ipos << " with PDG " << pdg_code << std::endl ; 

    // // if(!this->IsHandled  (pdg_code)) continue;


    // LOG("DarkSectorDecayer", pNOTICE)
    //       << "Decaying unstable particle: " << p->Name();


    // // TODO: here check pdg of particle
    // // if it's dark mediator call DecayDarkMediator
    // // if it's dark neutrino call DecayDarkNeutrino
    
    // //TODO: this block below
    // bool decayed = this->Decay(ipos, event);
    // if (!decayed) {
    //   LOG("DarkSectorDecayer", pWARN) << "Dark stuff not decayed!" ;
    //   LOG("DarkSectorDecayer", pWARN) << "Quitting the current event generation thread" ;

    //   event -> EventFlags() -> SetBitNumber(kDecayErr, true);

    //   genie::exceptions::EVGThreadException exception;
    //   exception.SetReason("Dark stuff not decayed"); // TODO
    //   exception.SwitchOnFastForward();
    //   throw exception;

    //   return ;
    // }


  LOG("DarkSectorDecayer", pNOTICE)
     << "Done finding & decaying dark sector particles";
}
//____________________________________________________________________________
std::vector<GHepParticle> DarkSectorDecayer::Decay(
  const GHepParticle & mother,
  const std::vector<int> & pdg_daughters) const
{
  TLorentzVector mother_p4 = *(mother.P4());
  LOG("DarkSectorDecayer", pINFO)
    << "Decaying a " << mother.GetName()
    << " with P4 = " << utils::print::P4AsString(&mother_p4);

  unsigned int nd = pdg_daughters.size();
  double mass[nd] = {0.};

  for(unsigned int iparticle = 0; iparticle < nd; iparticle++) {
    TParticlePDG * daughter = PDGLibrary::Instance()->Find(pdg_daughters[iparticle]);
    assert(daughter);

    mass[iparticle] = daughter->Mass();

    SLOG("DarkSectorDecayer", pINFO)
      << "+ daughter[" << iparticle << "]: "
      << daughter->GetName() << " (pdg-code = "
      << pdg_daughters[iparticle] << ", mass = " << mass[iparticle] << ")";
  }

  bool is_permitted = fPhaseSpaceGenerator.SetDecay(mother_p4, nd, mass);
  assert(is_permitted);

  // Find the maximum phase space decay weight
  double wmax = -1;
  for(int i=0; i<50; i++) {
    double w = fPhaseSpaceGenerator.Generate();
    wmax = TMath::Max(wmax,w);
  }
  assert(wmax>0);
  LOG("DarkSectorDecayer", pINFO)
    << "Max phase space gen. weight for current decay: " << wmax;

  // Generating un-weighted decays
  RandomGen * rnd = RandomGen::Instance();
  wmax *= 2;
  bool accept_decay=false;
  unsigned int itry=0;

  while(!accept_decay){
    itry++;
    assert(itry<kMaxUnweightDecayIterations);

    double w  = fPhaseSpaceGenerator.Generate();
    double gw = wmax * rnd->RndDec().Rndm();

    if(w>wmax) {
      LOG("DarkSectorDecayer", pWARN)
        << "Current decay weight = " << w << " > wmax = " << wmax;
    }
    LOG("DarkSectorDecayer", pINFO)
      << "Current decay weight = " << w << " / R = " << gw;

    accept_decay = (gw<=w);
  }

  // A decay was generated - Copy to the event record
  std::vector<GHepParticle> particles;
  // Loop over daughter list and add corresponding GHepParticles
  for(unsigned int id = 0; id < nd; id++) {
    TLorentzVector * daughter_p4 = fPhaseSpaceGenerator.GetDecay(id);
    LOG("DarkSectorDecayer", pDEBUG)
      << "Adding daughter particle with PDG code = " << pdg_daughters[id]
      << " and mass = " << mass[id] << " GeV";
    GHepStatus_t daughter_status_code = (pdg_daughters[id]==kPdgDNuMediator)
      ? kIStDecayedState : kIStStableFinalState;
    particles.push_back(GHepParticle(pdg_daughters[id], daughter_status_code,
                                     -1, -1, -1, -1,
                                     *daughter_p4, TLorentzVector()));
  }

  return particles;
}
//____________________________________________________________________________
int DarkSectorDecayer::SelectDecayChannel(
  const GHepParticle & mother,
  const GHepRecord * event,
  const std::vector<DecayChannel> & dcs,
  const double total_amplitude) const
{
  // Select a decay based on the amplitudes
  unsigned int ich = 0, sel_ich; // id of selected decay channel
  RandomGen * rnd = RandomGen::Instance();
  double x = total_amplitude * rnd->RndDec().Rndm();
  double partial_sum = 0. ;
  do {
    sel_ich = ich;
    partial_sum += dcs.at(ich++).second;
  } while (x > partial_sum );
  return sel_ich;
}
//____________________________________________________________________________
std::vector<DecayChannel> DarkSectorDecayer::DarkMediatorDecayChannels(
  const GHepParticle & mother,
  const GHepRecord * event) const
{
  // eq (4) and (5) and maybe some other higher order variations
  // TODO DNu: what alpha_D is?

  const double alpha_D = 0.25; // value on the paper
  std::array<int, 3> neutrinos = {kPdgNuE, kPdgNuMu, kPdgNuTau};
  std::array<int, 3> antineutrinos = {kPdgAntiNuE, kPdgAntiNuMu, kPdgAntiNuTau};
  std::vector<DarkSectorDecayer::DecayChannel> dcs;

  for(size_t i=0; i<neutrinos.size(); ++i){
    for(size_t j=0; j<antineutrinos.size(); ++j){// for antineutrinos
      const double decay_width = alpha_D/3. * fMixing2s[i]*fMixing2s[j] * fDMediatorMass;
      dcs.push_back(DecayChannel{{neutrinos[i], antineutrinos[j]}, decay_width});
    }
  }

  if(fDMediatorMass > 2.*PDGLibrary::Instance()->Find(kPdgElectron)->Mass()){
    const double decay_width = kAem*fEps2/3. * fDMediatorMass;
    dcs.push_back(DecayChannel{{kPdgElectron, kPdgPositron}, decay_width});
  }
  // In the future for the decay to muons
  // if(fDMediatorMass > 2.*PDGLibrary::Instance()->Find(kPdgMuon)->Mass()){
  //   const double decay_width = kAem*epsilon2/3. * fDMediatorMass;
  //   dcs.push_back(DecayChannel{{kPdgMuon, kPdgAntiMuon}, decay_width});
  // }
  return dcs;
}
//____________________________________________________________________________
  const GHepParticle & mother,
  const GHepRecord * event) const
std::vector<DarkSectorDecayer::DecayChannel> DarkSectorDecayer::DarkNeutrinoDecayChannels(
{
  // eq (3) and higher order variations 

  const double alpha_D = 0.25; // value on the paper
  std::array<int, 3> neutrinos = {kPdgNuE, kPdgNuMu, kPdgNuTau};
  std::array<int, 3> antineutrinos = {kPdgAntiNuE, kPdgAntiNuMu, kPdgAntiNuTau};
  std::vector<DarkSectorDecayer::DecayChannel> dcs;

  if(fDNuMass > fDMediatorMass){
    for(size_t i=0; i<neutrinos.size(); ++i){
      const double mass2ratio = fDMediatorMass2/fDNuMass2;
      const double p0 = 0.5*alpha_D * fMixing2s[3] * fMixing2s[i];
      const double p1 = fDNuMass*fDNuMass2/fDMediatorMass2;
      const double p2 = 1 - mass2ratio;
      const double p3 = 1 + mass2ratio - 2*mass2ratio*mass2ratio;
      const double decay_width = p0 * p1 * p2 * p3;
      dcs.push_back(DecayChannel{{neutrinos[i], kPdgDNuMediator}, decay_width});
      // TODO DNu: how to do antineutrinos decay?
      // dcs.push_back(DecayChannel{{antineutrinos[i], kPdgDNuMediator}, decay_width});
    }
  }
  return dcs;
}
//____________________________________________________________________________
void DarkSectorDecayer::SetSpaceTime(
  std::vector<GHepParticle> & pp,
  const GHepParticle & mother,
  double total_amplitude) const
{
  // TODO DNU: pay attention to units being used in GENIE!
  // convert decay amplitude into time
  const double lifetime =  units::second*1e24/total_amplitude; // units of 10ˆ-24 s

  RandomGen * rnd = RandomGen::Instance();
  double t = rnd->RndDec().Exp(lifetime);

  // get beta of decaying particle
  const TLorentzVector mother_X4 = *(mother.X4());
  TVector3 mother_boost = mother.P4()->BoostVector();

  // transport decay_particle with respect to their mother
  double speed_of_light = units::second/units::meter;
  TVector3 daughter_position = mother_X4.Vect() + mother_boost * (speed_of_light * t * 1e-9);// in fm

  for(auto & p : pp){
    p.SetPosition(TLorentzVector(daughter_position, mother_X4.T() + t));
  }
}
//____________________________________________________________________________


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
//          LOG("DarkSectorDecayer", pNOTICE)
//             << "Disabling decay channel containing resonance 1114";;
//          md = 999999999;
//      }
//      mass += md;
//   }
//   return mass;
// }
//____________________________________________________________________________
bool DarkSectorDecayer::ToBeDecayed(const GHepParticle & p) const
{
  GHepStatus_t status_code = p.Status();
  if(status_code != kIStDecayedState) return false;

  int pdg_code = p.Pdg();
  bool is_handled = false;
  if(pdg_code == kPdgDNuMediator ||
     pdg_code == kPdgDarkNeutrino ||
     pdg_code == kPdgAntiDarkNeutrino){
    is_handled = true;
  }

  LOG("DarkSectorDecayer", pDEBUG)
      << "Can decay particle with PDG code = " << pdg_code
      << "? " << ((is_handled)? "Yes" : "No");

  // // Find the particle in the PDG library & quit if it does not exist
  // TParticlePDG * mother = PDGLibrary::Instance()->Find(pdg_code);
  // if(!mother && is_handled) {
  //   LOG("DarkSectorDecayer", pERROR)
  //     << "\n *** The particle with PDG code = " << pdg_code
  //     << " was not found in PDGLibrary";
  //   //exit;
  // }

  return is_handled;
}
//____________________________________________________________________________
// void DarkSectorDecayer::InhibitDecay(int pdgc, TDecayChannel * dc) const
// {
//   if(! this->IsHandled(pdgc)) return;
//   if(!dc) return;

//   //
//   // Not implemented
//   //
// }
// //____________________________________________________________________________
// void DarkSectorDecayer::UnInhibitDecay(int pdgc, TDecayChannel * dc) const
// {
//   if(! this->IsHandled(pdgc)) return;
//   if(!dc) return;

//   //
//   // Not implemented
//   //
// }
// //____________________________________________________________________________
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
//      LOG("DarkSectorDecayer", pERROR)
//         << "\n *** The particle with PDG code = " << decay_particle_pdg_code
//          << " was not found in PDGLibrary";
//      return 0;
//   }
//   LOG("DarkSectorDecayer", pINFO)
//     << "Decaying a " << mother->GetName()
//     << " with P4 = " << utils::print::P4AsString(&decay_particle_p4);

//   // Get the resonance mass W (generally different from the mass associated
//   // with the input PDG code, since the it is produced off the mass shell)
//   double W = decay_particle_p4.M();
//   LOG("DarkSectorDecayer", pINFO) << "Available mass W = " << W;

//   // Get all decay channels
//   TObjArray * original_decay_list = mother->DecayList();

//   unsigned int nch = original_decay_list -> GetEntries();
//   LOG("DarkSectorDecayer", pINFO)
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

//       SLOG("DarkSectorDecayer", pDEBUG)
//                 << "Using channel: " << ich
//                 << " with final state mass = " << fsmass << " GeV";

//       tot_BR += ch->BranchingRatio();

//     } else {
//       SLOG("DarkSectorDecayer", pINFO)
//                 << "Suppresing channel: " << ich
//                 << " with final state mass = " << fsmass << " GeV";
//     } // final state mass

//     BR[ich] = tot_BR;
//   }//channel loop

//   if( tot_BR <= 0. ) {
//     SLOG("DarkSectorDecayer", pWARN)
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

//   LOG("DarkSectorDecayer", pINFO)
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

//      SLOG("DarkSectorDecayer", pINFO)
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
//   LOG("DarkSectorDecayer", pINFO)
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
//          LOG("DarkSectorDecayer", pWARN)
//             << "Current decay weight = " << w << " > wmax = " << wmax;
//       }
//       LOG("DarkSectorDecayer", pINFO)
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
//      LOG("DarkSectorDecayer", pDEBUG)
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
void DarkSectorDecayer::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DarkSectorDecayer::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DarkSectorDecayer::LoadConfig(void)
{
  double DKineticMixing = 0.;    // \varepsilon
  this->GetParam("Dark-KineticMixing", DKineticMixing);
  fEps2 = DKineticMixing * DKineticMixing;

  double DTheta = 0.;            // \theta
  this->GetParam("Dark-Theta", DTheta);
  fTheta2 = DTheta * DTheta;

  double DGaugeCoupling = 0.;   // g_D
  this->GetParam("Dark-GaugeCoupling", DGaugeCoupling);
  fgD2 = DGaugeCoupling * DGaugeCoupling;

  fDNuMass = 0.;
  this->GetParam("Dark-NeutrinoMass", fDNuMass);
  fDNuMass2 = fDNuMass * fDNuMass;

  fDMediatorMass = 0.;
  this->GetParam("Dark-MediatorMass", fDMediatorMass);
  fDMediatorMass2 = fDMediatorMass * fDMediatorMass;

  // Decayer::LoadConfig() ;

  // this -> GetParam( "DarkMediatorMass", fDarkMediatorMass ) ;

  // this -> GetParamDef( "Delta-ThetaOnly", fDeltaThetaOnly, true ) ;

  // this -> GetParamDef( "DeltaDecayMaximumTolerance", fMaxTolerance, 0.0005 ) ;

  // bool invalid_configuration = false ;

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
