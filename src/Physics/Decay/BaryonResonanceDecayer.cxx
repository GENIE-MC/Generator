//____________________________________________________________________________
/*
 Copyright (c) 2003-2018, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab
*/
//____________________________________________________________________________

#include <TClonesArray.h>
#include <TDecayChannel.h>
#include <TMath.h>

#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Utils/PrintUtils.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/Decay/BaryonResonanceDecayer.h"

using namespace genie;
using namespace genie::controls;
using namespace genie::constants;
//____________________________________________________________________________
BaryonResonanceDecayer::BaryonResonanceDecayer() :
Decayer("genie::BaryonResonanceDecayer")
{
  this->Initialize();
}
//____________________________________________________________________________
BaryonResonanceDecayer::BaryonResonanceDecayer(string config) :
Decayer("genie::BaryonResonanceDecayer", config)
{
  this->Initialize();
}
//____________________________________________________________________________
BaryonResonanceDecayer::~BaryonResonanceDecayer()
{

}
//____________________________________________________________________________
void BaryonResonanceDecayer::ProcessEventRecord(GHepRecord * event) const
{
  // Loop over particles, find unstable ones and decay them
  TObjArrayIter piter(event);
  GHepParticle * p = 0;
  int ipos = -1;

  while( (p = (GHepParticle *) piter.Next()) ) {

    ipos++;
    LOG("ResonanceDecay", pINFO) << "Checking: " << p->Name();

    int pdg_code = p->Pdg();
    GHepStatus_t status_code = p->Status();

    if(!this->IsHandled  (pdg_code)) continue;
    if(!this->ToBeDecayed(pdg_code, status_code)) continue;

    LOG("ResonanceDecay", pNOTICE)
          << "Decaying unstable particle: " << p->Name();

    bool decayed = this->Decay(ipos, event);
    assert(decayed); // handle this more graciously and throw an exception
  } // loop over particles

  LOG("ResonanceDecay", pNOTICE)
     << "Done finding & decaying baryon resonances";
}
//____________________________________________________________________________
bool BaryonResonanceDecayer::Decay(
  int decay_particle_id, GHepRecord * event) const
{
  // Reset previous decay weight
  fWeight = 1.;

  // Get particle to be decayed
  GHepParticle * decay_particle = event->Particle(decay_particle_id);
  if(!decay_particle) return false;

  // Select a decay channel
  TDecayChannel * selected_decay_channel =
     this->SelectDecayChannel(decay_particle_id, event);
  if(!selected_decay_channel) return false;

  // Decay the exclusive state and copy daughters in the event record
  this->DecayExclusive(decay_particle_id, event, selected_decay_channel);

  // Update the event weight for each weighted particle decay
  double weight = event->Weight() * fWeight;
  event->SetWeight(weight);

  // Mark input particle as a 'decayed state' & add its daughter links
  decay_particle->SetStatus(kIStDecayedState);

  return true;
}
//____________________________________________________________________________
TDecayChannel * BaryonResonanceDecayer::SelectDecayChannel(
   int decay_particle_id, GHepRecord * event) const
{
  // Get particle to be decayed
  GHepParticle * decay_particle = event->Particle(decay_particle_id);
  if(!decay_particle) return 0;

  // Get the particle 4-momentum and PDG code
  TLorentzVector decay_particle_p4 = *(decay_particle->P4());
  int decay_particle_pdg_code = decay_particle->Pdg();

  // Find the particle in the PDG library & quit if it does not exist
  TParticlePDG * mother =
     PDGLibrary::Instance()->Find(decay_particle_pdg_code);
  if(!mother) {
     LOG("ResonanceDecay", pERROR)
        << "\n *** The particle with PDG code = " << decay_particle_pdg_code
         << " was not found in PDGLibrary";
     return 0;
  }
  LOG("ResonanceDecay", pINFO)
    << "Decaying a " << mother->GetName()
    << " with P4 = " << utils::print::P4AsString(&decay_particle_p4);

  // Get the resonance mass W (generally different from the mass associated
  // with the input PDG code, since the it is produced off the mass shell)
  double W = decay_particle_p4.M();
  LOG("ResonanceDecay", pINFO) << "Available mass W = " << W;

  // Get all decay channels
  TObjArray * decay_list = mother->DecayList();
  unsigned int nch = decay_list->GetEntries();
  LOG("ResonanceDecay", pINFO)
    << mother->GetName() << " has: " << nch << " decay channels";

  // Loop over the decay channels (dc) and write down the branching
  // ratios to be used for selecting a decay channel.
  // Since a baryon resonance can be created at W < Mres, explicitly
  // check and inhibit decay channels for which W > final-state-mass

  double BR[nch], tot_BR = 0;

  bool isdelta =
    ( decay_particle_pdg_code ==  kPdgP33m1232_Delta0 ||
      decay_particle_pdg_code == -kPdgP33m1232_Delta0 ||
      decay_particle_pdg_code ==  kPdgP33m1232_DeltaP ||
      decay_particle_pdg_code == -kPdgP33m1232_DeltaP );

	for(unsigned int ich = 0; ich < nch; ich++) {

      TDecayChannel * ch = (TDecayChannel *) decay_list->At(ich);
      double fsmass = this->FinalStateMass(ch);

      if(fsmass < W) {
         SLOG("ResonanceDecay", pDEBUG)
            << "Using channel: " << ich
            << " with final state mass = " << fsmass << " GeV";
         double ch_BR = 0;
         if(isdelta) {
           ch_BR = this->DealsDeltaNGamma(decay_particle_pdg_code, ich, W);
         }
         else {
           ch_BR = ch->BranchingRatio();
         }
         tot_BR += ch_BR;
      } else {
         SLOG("ResonanceDecay", pINFO)
            << "Suppresing channel: " << ich
            << " with final state mass = " << fsmass << " GeV";
      } // final state mass
      BR[ich] = tot_BR;
  }//channel loop

  if(tot_BR==0) {
    SLOG("ResonanceDecay", pWARN)
      << "None of the " << nch << " decay channels is available @ W = " << W;
    return 0;
  }

  // Select a resonance based on the branching ratios
  unsigned int ich = 0, sel_ich; // id of selected decay channel
  RandomGen * rnd = RandomGen::Instance();
  double x = tot_BR * rnd->RndDec().Rndm();
  do {
    sel_ich = ich;
  } while (x > BR[ich++]);

  TDecayChannel * sel_ch = (TDecayChannel *) decay_list->At(sel_ich);

  LOG("ResonanceDecay", pINFO)
    << "Selected " << sel_ch->NDaughters() << "-particle decay channel ("
    << sel_ich << ") has BR = " << sel_ch->BranchingRatio();

  return sel_ch;
}
//____________________________________________________________________________
void BaryonResonanceDecayer::DecayExclusive(
  int decay_particle_id, GHepRecord * event, TDecayChannel * ch) const
     //int pdg_code, TLorentzVector & p, TDecayChannel * ch) const
{
  // Find the particle to be decayed in the event record
  GHepParticle * decay_particle = event->Particle(decay_particle_id);
  if(!decay_particle) return;

  // Get the decayed particle 4-momentum, 4-position and PDG code
  TLorentzVector decay_particle_p4 = *(decay_particle->P4());
  TLorentzVector decay_particle_x4 = *(decay_particle->X4());
  int decay_particle_pdg_code = decay_particle->Pdg();

  // Get the final state mass spectrum and the particle codes
  // for the selected decay channel
  unsigned int nd = ch->NDaughters();

  int    pdgc[nd];
  double mass[nd];

  for(unsigned int iparticle = 0; iparticle < nd; iparticle++) {
     int daughter_code = ch->DaughterPdgCode(iparticle);
     TParticlePDG * daughter = PDGLibrary::Instance()->Find(daughter_code);
     assert(daughter);

     pdgc[iparticle] = daughter_code;
     mass[iparticle] = daughter->Mass();

     SLOG("ResonanceDecay", pINFO)
         << "+ daughter[" << iparticle << "]: "
         << daughter->GetName() << " (pdg-code = "
         << pdgc[iparticle] << ", mass = " << mass[iparticle] << ")";
  }

  // Check whether the expected channel is Delta->pion+nucleon
  bool twobody =
    this->IsPiNDecayChannel(ch) &&
    this->IsDelta(decay_particle_pdg_code);

  // Decay the resonance using an N-body phase space generator
  // The particle will be decayed in its rest frame and then the daughters
  // will be boosted back to the original frame.

  bool is_permitted = fPhaseSpaceGenerator.SetDecay(decay_particle_p4, nd, mass);
  assert(is_permitted);

  // Find the maximum phase space decay weight
  // double wmax = fPhaseSpaceGenerator.GetWtMax();
  double wmax = -1;
  for(int i=0; i<50; i++) {
     double w = fPhaseSpaceGenerator.Generate();
     wmax = TMath::Max(wmax,w);
  }
  assert(wmax>0);
  LOG("ResonanceDecay", pINFO)
    << "Max phase space gen. weight for current decay: " << wmax;

  // Define variables for the Wtheta selection
  double aidrnd          = 0;
  double wthetacheck     = 0;
  double p32check        = 0.75;
  double p12check        = 1-p32check;
  double p2costhetacheck = 0;
  double costhetacheck   = 0;

  TLorentzVector vpioncheck;
  TLorentzVector vcheckdelta;

  if(fGenerateWeighted)
  {
    // Generating weighted decays
    // Do a single draw of momentum 4-vectors and then stop,
    // taking into account the weight for this particular draw
    double w = fPhaseSpaceGenerator.Generate();
    fWeight *= TMath::Max(w/wmax, 1.);
  }
  else
  {
    // Generating un-weighted decays
    RandomGen * rnd = RandomGen::Instance();
    wmax *= 2;
    bool accept_decay=false;
    unsigned int itry=0;

    while(!accept_decay)
    {
      itry++;
      assert(itry<kMaxUnweightDecayIterations);

      double w  = fPhaseSpaceGenerator.Generate();
      double gw = wmax * rnd->RndDec().Rndm();

      if(w>wmax) {
         LOG("ResonanceDecay", pWARN)
            << "Current decay weight = " << w << " > wmax = " << wmax;
      }
      LOG("ResonanceDecay", pINFO)
        << "Current decay weight = " << w << " / R = " << gw;

      // Extra logic that applies only for Delta -> N gamma
      if(twobody)
      {
      }

      accept_decay = (gw<=w);
    }//accept_decay
  }//fGenerateWeighted

  // A decay was generated - Copy to the event record

  // Check whether the interaction is off a nuclear target or free nucleon
  // Depending on whether this module is run before or after the hadron
  // transport module it would affect the daughters status code
  GHepParticle * target_nucleus = event->TargetNucleus();
  bool in_nucleus = (target_nucleus!=0);

  // Loop over daughter list and add corresponding GHepParticles
  for(unsigned int id = 0; id < nd; id++) {

     int daughter_pdg_code = pdgc[id];

     TLorentzVector * daughter_p4 = fPhaseSpaceGenerator.GetDecay(id);
     LOG("ResonanceDecay", pDEBUG)
        << "Adding daughter particle with PDG code = " << pdgc[id]
        << " and mass = " << mass[id] << " GeV";

     bool is_hadron = pdg::IsHadron(daughter_pdg_code);
     bool hadron_in_nuc = (in_nucleus && is_hadron && fRunBefHadroTransp);

     GHepStatus_t daughter_status_code = (hadron_in_nuc) ?
          kIStHadronInTheNucleus : kIStStableFinalState;

     event->AddParticle(
       daughter_pdg_code, daughter_status_code, decay_particle_id,-1,-1,-1,
       *daughter_p4, decay_particle_x4);
  }

  //
  // while(1){  //start a loop until break;
  //
  //   if(fGenerateWeighted)
  //   {
  //     // *** generating weighted decays ***
  //     double w = fPhaseSpaceGenerator.Generate();
  //     fWeight *= TMath::Max(w/wmax, 1.);
  //   }
  //   else
  //   {
  //     // *** generating un-weighted decays ***
  //     RandomGen * rnd = RandomGen::Instance();
  //     wmax *= 2;
  //     bool accept_decay=false;
  //     unsigned int itry=0;
  //
  //     while(!accept_decay)
  //     {
  //       itry++;
  //       assert(itry<kMaxUnweightDecayIterations);
  //
  //       double w  = fPhaseSpaceGenerator.Generate();
  //       double gw = wmax * rnd->RndDec().Rndm();
  //
  //       if(w>wmax) {
  //          LOG("Decay", pWARN)
  //             << "Current decay weight = " << w << " > wmax = " << wmax;
  //       }
  //       LOG("Decay", pINFO)
  //         << "Current decay weight = " << w << " / R = " << gw;
  //
  //       accept_decay = (gw<=w);
  //     }//accept_decay
  //   }//fGenerateWeighted
  //
  //   // // Add the mother particle to the event record (KS=11 as in PYTHIA)
  //   // TParticlePDG * mother = PDGLibrary::Instance()->Find(pdg_code);
  //   //
  //   // double px   = p.Px();
  //   // double py   = p.Py();
  //   // double pz   = p.Pz();
  //   // double E    = p.Energy();
  //   // double M    = mother->Mass();
  //
  //   if(twobody)
  //   {
  //     // // restore mother particle's 4-momentum.
  //     // vcheckdelta.SetPxPyPzE(px,py,pz,E);}
  //     //
  //     // // restore mother particle to the temp_particle list.
  //     // new ( (*temp_particle_list)[0] )
  //     //   TMCParticle(11,pdg_code,0,0,0,px,py,pz,E,M,0,0,0,0,0);
  //
  //     // Add the daughter particles to the event record
  //     for(unsigned int id = 0; id < nd; id++) {
  //
  //        TLorentzVector * p4 = fPhaseSpaceGenerator.GetDecay(id);
  //        LOG("Decay", pDEBUG)
  //           << "Adding final state particle PDGC = " << pdgc[id]
  //           << " with mass = " << mass[id] << " GeV";
  //
  //        px   = p4->Px();
  //        py   = p4->Py();
  //        pz   = p4->Pz();
  //        E    = p4->Energy();
  //        M    = mass[id];
  //
  //        // Set aidrnd and wthetacheck to control the angular distribution.
	//        if(twobody){
  //           if(pdgc[id]==kPdgPiP || pdgc[id]==kPdgPi0 ||pdgc[id]==kPdgPiM)
  //           {
	// 		       vpioncheck.SetPxPyPzE(px,py,pz,E);
	// 			     // Boost pion 4-vec from lab frame into CM frame
  //            vpioncheck.Boost(-vcheckdelta.BoostVector());
	// 	         costhetacheck=vpioncheck.Pz()/sqrt(
  //               vpioncheck.Px()*vpioncheck.Px()+
  //               vpioncheck.Py()*vpioncheck.Py()+
  //               vpioncheck.Pz()*vpioncheck.Pz());
  //            p2costhetacheck=0.5*(3*costhetacheck*costhetacheck-1);
  //            wthetacheck=1-p32check*(p2costhetacheck)+p12check*(p2costhetacheck);
  //            aidrnd=1.25*gRandom->Rndm();
	// 	       }  //end pion-selection
	//        }  //end twobody
  //
  //        new ( (*temp_particle_list)[1+id] )
  //           TMCParticle(1,pdgc[id],0,0,0,px,py,pz,E,M,0,0,0,0,0);
  //     }//end daughter particle loop
  //
  //     if(!twobody) break;
  //     if(twobody && wthetacheck>=aidrnd) break;
  //
  // }//end while(1)
  //
  // // particle_list=temp_particle_list;
  // // particle_list->SetOwner(true);
  // // return particle_list;
  //

}
//____________________________________________________________________________
double BaryonResonanceDecayer::DealsDeltaNGamma(
  int id_mother, int ichannel, double W) const
{
  //-- auxiliary parameters
  int DeltaFlag = 0;
	  if (id_mother == 2114 || id_mother==-2114) {
		  DeltaFlag = 1; // Delta0 or Delta0_bar
	  }
	  else if (id_mother == 2214 || id_mother==-2214) {
	      DeltaFlag = 2; // Delta+ or anti_Delta+
	  }
	  else  {
     // cout<<"Mother particle is not Delta+ or Delta0!!!"<<endl;
	  return 0;
	  }
  double mN  =   genie::constants::kNucleonMass;
  double mPi = genie::constants::kPi0Mass;
//  double mN    = kNucleonMass;
//  double mPi   = kPi0Mass;

  if (W<=mN+mPi) {
	  if (ichannel == 0) {return 0;} // ichannel =0,1,2 has to match
	                                // the channel order in genie_pdg_table.dat
	  if (ichannel == 1) {return 0;}
	  if (ichannel == 2) {return 1;}
  } else {

  double m  = 1.232;

  double width0= 0.12;

  double m_2   = TMath::Power(m, 2);
  double mN_2  = TMath::Power(mN,   2);
  double W_2   = TMath::Power(W,    2);
  double m_aux1= TMath::Power(mN+mPi, 2);
  double m_aux2= TMath::Power(mN-mPi, 2);

  double BRPi0    = 0.994;
  double BRPi01   = 0.667002;
  double BRPi02   = 0.332998;
  double BRgamma0 = 0.006;
  double widPi0   = width0*BRPi0;
  double widgamma0= width0*BRgamma0;

  double pPiW   = TMath::Sqrt((W_2-m_aux1)*(W_2-m_aux2))/(2*W);
  double pPim   = TMath::Sqrt((m_2-m_aux1)*(m_2-m_aux2))/(2*m);
  double EgammaW= (W_2-mN_2)/(2*W);
  double Egammam= (m_2-mN_2)/(2*m);
  double TPiW=TMath::Power(pPiW, 3);
  double TPim=TMath::Power(pPim, 3);
  double fgammaW= 1/(TMath::Power(1+EgammaW*EgammaW/0.706, 2));
  double fgammam= 1/(TMath::Power(1+Egammam*Egammam/0.706, 2));


  double Rinverse = widPi0*TMath::Power(Egammam, 3)*TMath::Power(fgammam, 2)*TPiW
	     /(widgamma0*TMath::Power(EgammaW, 3)*TMath::Power(fgammaW, 2)*TPim);
  double BRPi = Rinverse/(1+Rinverse);
  double BRgamma = 1/(1+Rinverse);

  if (DeltaFlag==1) {
  	  if (ichannel == 0) {return BRPi*BRPi02;}
	    if (ichannel == 1) {return BRPi*BRPi01;}
	    if (ichannel == 2) {return BRgamma;}
  }
  if (DeltaFlag==2) {
  	  if (ichannel == 0) {return BRPi*BRPi01;}
	    if (ichannel == 1) {return BRPi*BRPi02;}
	    if (ichannel == 2) {return BRgamma;}
  }
  }
 return 0;
}
//____________________________________________________________________________
double BaryonResonanceDecayer::Weight(void) const
{
  return fWeight;
}
//____________________________________________________________________________
double BaryonResonanceDecayer::FinalStateMass(TDecayChannel * ch) const
{
// Computes the total mass of the final state system

  double mass = 0;
  unsigned int nd = ch->NDaughters();

  for(unsigned int iparticle = 0; iparticle < nd; iparticle++) {

     int daughter_code = ch->DaughterPdgCode(iparticle);
     TParticlePDG * daughter = PDGLibrary::Instance()->Find(daughter_code);
     assert(daughter);

     double md = daughter->Mass();

     // hack to switch off channels giving rare  occurences of |1114| that has
     // no decay channels in the pdg table (08/2007)
     if(TMath::Abs(daughter_code) == 1114) {
         LOG("ResonanceDecay", pNOTICE)
            << "Disabling decay channel containing resonance 1114";;
         md = 999999999;
     }
     mass += md;
  }
  return mass;
}
//____________________________________________________________________________
bool BaryonResonanceDecayer::IsDelta(int pdg_code) const
{
  return (pdg_code == kPdgP33m1232_DeltaPP ||
          pdg_code == kPdgP33m1232_DeltaP  ||
          pdg_code == kPdgP33m1232_Delta0);
}
//____________________________________________________________________________
bool BaryonResonanceDecayer::IsPiNDecayChannel(TDecayChannel * ch) const
{
  if(!ch) return false;

  unsigned int nd = ch->NDaughters();
  if(nd != 2) return false;

  int  npi    = 0;
  int  nnucl  = 0;

  for(unsigned int iparticle = 0; iparticle < nd; iparticle++) {

     int daughter_code = ch->DaughterPdgCode(iparticle);

	   if(daughter_code == kPdgPiP ||
        daughter_code == kPdgPi0 ||
        daughter_code == kPdgPiM)
     {
         npi++;
     }
     else
     if(daughter_code == kPdgNeutron ||
        daughter_code == kPdgProton)
     {
         nnucl++;
     }
	 }//iparticle
   if(npi == 1 && nnucl == 1) return true;

   return false;
}
//____________________________________________________________________________
void BaryonResonanceDecayer::Initialize(void) const
{

}
//____________________________________________________________________________
bool BaryonResonanceDecayer::IsHandled(int pdg_code) const
{
  if( utils::res::IsBaryonResonance(pdg_code) ) return true;

  LOG("ResonanceDecay", pINFO)
      << "This algorithm can not decay particles with PDG code = "
      << pdg_code;

  return false;
}
//____________________________________________________________________________
void BaryonResonanceDecayer::InhibitDecay(int pdgc, TDecayChannel * dc) const
{
  if(! this->IsHandled(pdgc)) return;
  if(!dc) return;

  //
  // Not implemented
  //
}
//____________________________________________________________________________
void BaryonResonanceDecayer::UnInhibitDecay(int pdgc, TDecayChannel * dc) const
{
  if(! this->IsHandled(pdgc)) return;
  if(!dc) return;

  //
  // Not implemented
  //
}
//____________________________________________________________________________
