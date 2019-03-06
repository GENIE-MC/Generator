//____________________________________________________________________________
/*
 Copyright (c) 2003-2017, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - November 27, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Oct 02, 2009 - CA
   Add dummy `UnInhibitDecay(int,TDecayChannel*) const' and `InhibitDecay(int,
   TDecayChannel*) const' methods to conform to the DecayModelI interface.
   To implement soon.
*/
//____________________________________________________________________________

#include <TClonesArray.h>
#include <TDecayChannel.h>
#include <RVersion.h>
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,15,6)
#include <TMCParticle.h>
#else
#include <TMCParticle6.h>
#endif
#include <TMath.h>
#include <TF1.h>

#include "Framework/Conventions/Controls.h"
#include "Physics/Decay/RadiativeDecayer.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Utils/PrintUtils.h"
#include "Framework/Conventions/Constants.h"

using namespace genie;
using namespace genie::controls;
using namespace genie::constants;
//____________________________________________________________________________
RadiativeDecayer::RadiativeDecayer() :
DecayModelI("genie::RadiativeDecayer")
{

}
//____________________________________________________________________________
RadiativeDecayer::RadiativeDecayer(string config) :
DecayModelI("genie::RadiativeDecayer", config)
{

}
//____________________________________________________________________________
RadiativeDecayer::~RadiativeDecayer()
{

}
//____________________________________________________________________________
bool RadiativeDecayer::IsHandled(int code) const
{
// handles only requests to decay baryon resonances
  if( code == 11 || code == -11 ) return true;

  LOG("RadiativeDecay", pINFO) 
      << "This algorithm can not decay particles with PDG code = " << code;

  return false;
}
//____________________________________________________________________________
void RadiativeDecayer::Initialize(void) const
{

}
//____________________________________________________________________________
double RadiativeDecayer::findEnergyLoss(const DecayerInputs_t & inp) const 
{
	double e = inp.P4->E();
	double Q = 1.; // typical momentum scale,
	double a = (kAem/kPi)*(2*TMath::Log(Q/kElectronMass) - 1.);
	TF1 *f = new TF1("f","([0]/x)*TMath::Power(x/[1],[0])",0,e);
	f->SetParameter(0,a);
	f->SetParameter(1,e);
	double energyloss = f->GetRandom();
	LOG("RadiativeDecay", pNOTICE) << "energy loss " << energyloss;
	return energyloss;	
}
//____________________________________________________________________________
TClonesArray* RadiativeDecayer::Decay(const DecayerInputs_t & inp) const
{  
  if ( ! this->IsHandled(inp.PdgCode) ) return 0;

  //-- Find the particle in the PDG library & quit if it does not exist
  TParticlePDG * decayed_mother = PDGLibrary::Instance()->Find(inp.PdgCode);
  decayed_mother->Print();

  if(!decayed_mother) {
     LOG("RadiativeDecay", pERROR)
          << "\n *** The particle with PDG-Code = " << inp.PdgCode
                                         << " was not found in PDGLibrary";
     return 0;                               
  }  
  LOG("RadiativeDecay", pINFO)
       << "Decaying a " << decayed_mother->GetName()
                        << " with P4 = " << utils::print::P4AsString(inp.P4);
  
  //-- Reset previous weight
  fWeight = 1.;

//  double W = inp.P4->M();
//  LOG("Decay", pINFO) << "Available mass W = " << W;
//
//  //-- Get all decay channels
//  TObjArray * decay_list = decayed_mother->DecayList();
//  decayed_mother->Print();
//  unsigned int nch = decay_list->GetEntries();
//  LOG("Decay", pINFO)
//               << decayed_mother->GetName() << " has: " << nch << " decay channels";
//
//  //-- Loop over the decay channels (dc) and write down the branching
//  //  //   ratios to be used for selecting a decay channel.
//  //    //   Since a baryon resonance can be created at W < Mres, explicitly
//  //      //   check and inhibit decay channels for which W > final-state-mass
//  double BR[nch], tot_BR = 0;
//
//  if(inp.PdgCode== 2114 || inp.PdgCode==-2114 || inp.PdgCode==2214||inp.PdgCode==-2214){
//          for(unsigned int ich = 0; ich < nch; ich++) {
//
//     TDecayChannel * ch = (TDecayChannel *) decay_list->At(ich);
//     double fsmass = this->FinalStateMass(ch);
//
//     if(fsmass < W) {
//       SLOG("Decay", pDEBUG)
//               << "Using channel: " << ich
//                        << " with final state mass = " << fsmass << " GeV";
//     } else {
//       SLOG("Decay", pINFO)
//               << "Suppresing channel: " << ich
//                        << " with final state mass = " << fsmass << " GeV";
//     }
//     BR[ich] = tot_BR;
//  }
//  }  else {
//  for(unsigned int ich = 0; ich < nch; ich++) {
//
//     TDecayChannel * ch = (TDecayChannel *) decay_list->At(ich);
//     double fsmass = this->FinalStateMass(ch);
//
//     if(fsmass < W) {
//       SLOG("Decay", pDEBUG)
//               << "Using channel: " << ich
//                        << " with final state mass = " << fsmass << " GeV";
//       tot_BR += ch->BranchingRatio();
//     } else {
//       SLOG("Decay", pINFO)
//               << "Suppresing channel: " << ich
//                        << " with final state mass = " << fsmass << " GeV";
//     }
//     BR[ich] = tot_BR;
//  }
//  }
//  if(tot_BR==0) {
//    SLOG("Decay", pWARN)
//      << "None of the " << nch << " decay chans is available @ W = " << W;
//    return 0;
//  }
//
//  //-- Select a resonance based on the branching ratios
//  unsigned int ich = 0, sel_ich; // id of selected decay channel
//  RandomGen * rnd = RandomGen::Instance();
//  double x = tot_BR * rnd->RndDec().Rndm();
//  do {
//    sel_ich = ich;
//
//  } while (x > BR[ich++]);
//
//  TDecayChannel * ch = (TDecayChannel *) decay_list->At(sel_ich);
//
//  LOG("Decay", pINFO)
//    << "Selected " << ch->NDaughters() << "-particle decay chan ("
//    << sel_ich << ") has BR = " << ch->BranchingRatio();
//
  //-- Decay the exclusive state and return the particle list
  TLorentzVector p4(*inp.P4);
  int pdg_code = inp.PdgCode;
  //-- Get the final state mass spectrum and the particle codes
  unsigned int nd = 2; 

  int    pdgc[nd];
  double mass[nd];

// Customized part
  bool twobody=false;      // flag of expected channel Delta->pion+nucleon
///////////////////////////////////////

  for(unsigned int iparticle = 0; iparticle < nd; iparticle++) {

     int daughter_code = 0;
     if (iparticle == 0) daughter_code = 11; 
     else if (iparticle == 1) daughter_code = 22;
     TParticlePDG * daughter = PDGLibrary::Instance()->Find(daughter_code);
     assert(daughter);

     pdgc[iparticle] = daughter_code;
     mass[iparticle] = daughter->Mass();

     SLOG("RadiativeDecay", pINFO)
       << "+ daughter[" << iparticle << "]: "
        << daughter->GetName() << " (pdg-code = "
          << pdgc[iparticle] << ", mass = " << mass[iparticle] << ")";
  }

  //-- Decay the resonance using an N-body phase space generator
  //   The particle will be decayed in its rest frame and then the daughters
  //   will be boosted back to the original frame.

  //bool is_permitted = fPhaseSpaceGenerator.SetDecay(p4, nd, mass);
  //assert(is_permitted);
  //SLOG("RadiativeDecay", pINFO)
  //     << "is permitted "<< is_permitted;

// Customized part--define variables for the Wtheta selection---------------
  double aidrnd=0;
  double wthetacheck=0;
  double p32check=0.75; 
  double p12check=1-p32check;
  double p2costhetacheck=0;
  double costhetacheck=0;

  TLorentzVector vpioncheck;
  TLorentzVector vcheckdelta;


  //-- Create the event record
  TClonesArray * particle_list = new TClonesArray("TMCParticle", 1+nd);
  TClonesArray * temp_particle_list = new TClonesArray("TMCParticle", 1+nd);//A temprary record.

  ////double wmax = fPhaseSpaceGenerator.GetWtMax();
  //double wmax = -1;
  //for(int i=0; i<50; i++) {
  //   double w = fPhaseSpaceGenerator.Generate();
  //   wmax = TMath::Max(wmax,w);
  //}
  //assert(wmax>0);
  //LOG("RadiativeDecay", pINFO)
  //   << "Max phase space gen. weight for current decay: " << wmax;

  //if(fGenerateWeighted)
  //{
  //   // *** generating weighted decays ***
  //   double w = fPhaseSpaceGenerator.Generate();
  //   fWeight *= TMath::Max(w/wmax, 1.);
  //}
  //else
  //{
  //   // *** generating un-weighted decays ***
  //   RandomGen * rnd1 = RandomGen::Instance();
  //   double wmax = 0.000000001;
  //   wmax *= 2;
  //   bool accept_decay=false;
  //   register unsigned int itry=0;

  //   while(!accept_decay)
  //   {
  //     itry++;
  //     assert(itry<kMaxUnweightDecayIterations);

  //     double w  = fPhaseSpaceGenerator.Generate();
  //     double gw = wmax * rnd1->RndDec().Rndm();

  //     if(w>wmax) {
  //        LOG("RadiativeDecay", pWARN) 
  //           << "Current decay weight = " << w << " > wmax = " << wmax;
  //     }
  //     LOG("RadiativeDecay", pINFO) 
  //        << "Current decay weight = " << w << " / R = " << gw;

  //     accept_decay = (gw<=w);
  //   }
  //}

  SLOG("RadiativeDecay", pINFO) <<"after Create the event record";
  //-- Add the mother particle to the event record (KS=11 as in PYTHIA)
  TParticlePDG * mother = PDGLibrary::Instance()->Find(pdg_code);

  double px   = p4.Px();
  double py   = p4.Py();
  double pz   = p4.Pz();
  double E    = p4.Energy();
  double M    = mother->Mass();

  //if(twobody){vcheckdelta.SetPxPyPzE(px,py,pz,E);}  // restore mother particle's 4-momentum.

  new ( (*temp_particle_list)[0] )
                    TMCParticle(0,pdg_code,0,0,0,px,py,pz,E,M,0,0,0,0,0); //restore mother particle to the temp_particle list.

  SLOG("RadiativeDecay", pINFO) <<"after Add the mother particle to the event record";
  //-- Add the daughter particles to the event record
  for(unsigned int id = 0; id < nd; id++) {

       //TLorentzVector * p4d = fPhaseSpaceGenerator.GetDecay(id);
       //LOG("RadiativeDecay", pINFO)
       //        << "Adding final state particle PDGC = " << pdgc[id]
       //                            << " with mass = " << mass[id] << " GeV";
       //px   = p4d->Px();
       //py   = p4d->Py();
       //pz   = p4d->Pz();
       //E    = p4d->Energy();
       //M    = mass[id];

       LOG("RadiativeDecay", pINFO)
               << "Adding final state particle PDGC = " << pdgc[id]
                                   << " with mass = " << mass[id] << " GeV";

  double x = 0.000005;
  int ks = 0; // particle state 
  if (pdgc[id] == 11 ) {
       ks = 2; // the electron carries on to be an intermidiate state
       px   = p4.Px();
       py   = p4.Py();
       pz   = p4.Pz()*(1-x);
       M    = mass[id];
       E    = sqrt(pow(pz,2)+pow(M,2));
  }
  else if (pdgc[id] == 22 ) {
       ks = 1; // the photon is a final state stable particle;
       px   = 0.;
       py   = 0.;
       pz   = p4.Pz()*x;
       E    = pz;
       M    = mass[id];
  }
  new ( (*temp_particle_list)[1+id] )
                         TMCParticle(ks,pdgc[id],0,0,0,px,py,pz,E,M,0,0,0,0,0);
  LOG("RadiativeDecay", pINFO) <<"after Add the daughter particles to the event record";
  }//end daughter particle loop

  LOG("RadiativeDecay", pINFO) <<"after daughter particle loop"; 

  particle_list=temp_particle_list;
  //-- Set owner and return
  particle_list->SetOwner(true);
  return particle_list;
}
//____________________________________________________________________________
double RadiativeDecayer::Weight(void) const
{
  return fWeight;
}
//____________________________________________________________________________
void RadiativeDecayer::InhibitDecay(int pdgc, TDecayChannel * dc) const
{
  if(! this->IsHandled(pdgc)) return;

  if(!dc) return;

}
//____________________________________________________________________________
void RadiativeDecayer::UnInhibitDecay(int pdgc, TDecayChannel * dc) const
{
  if(! this->IsHandled(pdgc)) return;

  if(!dc) return;

}
//____________________________________________________________________________
double RadiativeDecayer::FinalStateMass(TDecayChannel * ch) const
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
         LOG("RadiativeDecay", pNOTICE)
                  << "Disabling decay channel containing resonance 1114";;
         md = 999999999;
     }
     mass += md;
  }  
  return mass;
}
//____________________________________________________________________________
void RadiativeDecayer::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void RadiativeDecayer::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void RadiativeDecayer::LoadConfig(void)
{
// Read configuration options or set defaults

  //-- Generated weighted or un-weighted hadronic systems
  //fGenerateWeighted = fConfig->GetBoolDef("generate-weighted", false);
  //GetParam( "generate-weighted",fGenerateWeighted);
}
//____________________________________________________________________________
