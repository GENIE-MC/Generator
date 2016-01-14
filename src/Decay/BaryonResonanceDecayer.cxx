//____________________________________________________________________________
/*
 Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
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

#include "BaryonResonance/BaryonResUtils.h"
#include "Conventions/Controls.h"
#include "Decay/BaryonResonanceDecayer.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGLibrary.h"
#include "Numerical/RandomGen.h"
#include "Utils/PrintUtils.h"

using namespace genie;
using namespace genie::controls;

//____________________________________________________________________________
BaryonResonanceDecayer::BaryonResonanceDecayer() :
DecayModelI("genie::BaryonResonanceDecayer")
{

}
//____________________________________________________________________________
BaryonResonanceDecayer::BaryonResonanceDecayer(string config) :
DecayModelI("genie::BaryonResonanceDecayer", config)
{

}
//____________________________________________________________________________
BaryonResonanceDecayer::~BaryonResonanceDecayer()
{

}
//____________________________________________________________________________
bool BaryonResonanceDecayer::IsHandled(int code) const
{
// handles only requests to decay baryon resonances

  if( utils::res::IsBaryonResonance(code) ) return true;

  LOG("Decay", pINFO) 
      << "This algorithm can not decay particles with PDG code = " << code;

  return false;
}
//____________________________________________________________________________
TClonesArray* BaryonResonanceDecayer::Decay(const DecayerInputs_t & inp) const
{  
  if ( ! this->IsHandled(inp.PdgCode) ) return 0;

  //-- Find the particle in the PDG library & quit if it does not exist
  TParticlePDG * mother = PDGLibrary::Instance()->Find(inp.PdgCode);

  if(!mother) {
     LOG("Decay", pERROR)
          << "\n *** The particle with PDG-Code = " << inp.PdgCode
                                         << " was not found in PDGLibrary";
     return 0;                               
  }  
  LOG("Decay", pINFO)
       << "Decaying a " << mother->GetName()
                        << " with P4 = " << utils::print::P4AsString(inp.P4);
  
  //-- Reset previous weight
  fWeight = 1.;

  //-- Get the resonance mass W (generally different from the mass associated
  //   with the input pdg_code, since the it is produced off the mass shell)
  double W = inp.P4->M();
  LOG("Decay", pINFO) << "Available mass W = " << W;
  
  //-- Get all decay channels
  TObjArray * decay_list = mother->DecayList();
  unsigned int nch = decay_list->GetEntries();
  LOG("Decay", pINFO)
               << mother->GetName() << " has: " << nch << " decay channels";

  //-- Loop over the decay channels (dc) and write down the branching
  //   ratios to be used for selecting a decay channel.
  //   Since a baryon resonance can be created at W < Mres, explicitly
  //   check and inhibit decay channels for which W > final-state-mass
  
  double BR[nch], tot_BR = 0;    
  
  for(unsigned int ich = 0; ich < nch; ich++) {

     TDecayChannel * ch = (TDecayChannel *) decay_list->At(ich);
     double fsmass = this->FinalStateMass(ch);

     if(fsmass < W) {
       SLOG("Decay", pDEBUG)
               << "Using channel: " << ich 
                        << " with final state mass = " << fsmass << " GeV";         
       tot_BR += ch->BranchingRatio();
     } else {       
       SLOG("Decay", pINFO)
               << "Suppresing channel: " << ich 
                        << " with final state mass = " << fsmass << " GeV";         
     }
     BR[ich] = tot_BR;
  }

  if(tot_BR==0) {
    SLOG("Decay", pWARN) 
      << "None of the " << nch << " decay chans is available @ W = " << W;
    return 0;    
  }

  //-- Select a resonance based on the branching ratios
  unsigned int ich = 0, sel_ich; // id of selected decay channel
  RandomGen * rnd = RandomGen::Instance();
  double x = tot_BR * rnd->RndDec().Rndm();
  do { 
    sel_ich = ich;
    
  } while (x > BR[ich++]);

  TDecayChannel * ch = (TDecayChannel *) decay_list->At(sel_ich);

  LOG("Decay", pINFO) 
    << "Selected " << ch->NDaughters() << "-particle decay chan (" 
    << sel_ich << ") has BR = " << ch->BranchingRatio();

  //-- Decay the exclusive state and return the particle list
  TLorentzVector p4(*inp.P4);
  return ( this->DecayExclusive(inp.PdgCode, p4, ch) );
}
//____________________________________________________________________________
void BaryonResonanceDecayer::Initialize(void) const
{

}
//____________________________________________________________________________
TClonesArray * BaryonResonanceDecayer::DecayExclusive(
     int pdg_code, TLorentzVector & p, TDecayChannel * ch) const
{
  //-- Get the final state mass spectrum and the particle codes
  unsigned int nd = ch->NDaughters();

  int    pdgc[nd];
  double mass[nd];

  for(unsigned int iparticle = 0; iparticle < nd; iparticle++) {

     int daughter_code = ch->DaughterPdgCode(iparticle);
     TParticlePDG * daughter = PDGLibrary::Instance()->Find(daughter_code);
     assert(daughter);

     pdgc[iparticle] = daughter_code;
     mass[iparticle] = daughter->Mass();

     SLOG("Decay", pINFO)
       << "+ daughter[" << iparticle << "]: "
        << daughter->GetName() << " (pdg-code = "
          << pdgc[iparticle] << ", mass = " << mass[iparticle] << ")";
  }

  //-- Decay the resonance using an N-body phase space generator
  //   The particle will be decayed in its rest frame and then the daughters
  //   will be boosted back to the original frame.

  bool is_permitted = fPhaseSpaceGenerator.SetDecay(p, nd, mass);
  assert(is_permitted);

  //double wmax = fPhaseSpaceGenerator.GetWtMax();
  double wmax = -1;
  for(int i=0; i<50; i++) {
     double w = fPhaseSpaceGenerator.Generate();
     wmax = TMath::Max(wmax,w);
  }
  assert(wmax>0);
  LOG("Decay", pINFO)
     << "Max phase space gen. weight for current decay: " << wmax;

  if(fGenerateWeighted)
  {
     // *** generating weighted decays ***
     double w = fPhaseSpaceGenerator.Generate();
     fWeight *= TMath::Max(w/wmax, 1.);
  }
  else
  {
     // *** generating un-weighted decays ***
     RandomGen * rnd = RandomGen::Instance();
     wmax *= 2;
     bool accept_decay=false;
     register unsigned int itry=0;

     while(!accept_decay)
     {
       itry++;
       assert(itry<kMaxUnweightDecayIterations);

       double w  = fPhaseSpaceGenerator.Generate();
       double gw = wmax * rnd->RndDec().Rndm();

       if(w>wmax) {
          LOG("Decay", pWARN) 
             << "Current decay weight = " << w << " > wmax = " << wmax;
       }
       LOG("Decay", pINFO) 
          << "Current decay weight = " << w << " / R = " << gw;

       accept_decay = (gw<=w);
     }
  }

  //-- Create the event record
  TClonesArray * particle_list = new TClonesArray("TMCParticle", 1+nd);

  //-- Add the mother particle to the event record (KS=11 as in PYTHIA)
  TParticlePDG * mother = PDGLibrary::Instance()->Find(pdg_code);

  double px   = p.Px();
  double py   = p.Py();
  double pz   = p.Pz();
  double E    = p.Energy();
  double M    = mother->Mass();

  new ( (*particle_list)[0] )
                    TMCParticle(11,pdg_code,0,0,0,px,py,pz,E,M,0,0,0,0,0);

  //-- Add the daughter particles to the event record
  for(unsigned int id = 0; id < nd; id++) {

       TLorentzVector * p4 = fPhaseSpaceGenerator.GetDecay(id);
       LOG("Decay", pDEBUG)
               << "Adding final state particle PDGC = " << pdgc[id]
                                   << " with mass = " << mass[id] << " GeV";
       px   = p4->Px();
       py   = p4->Py();
       pz   = p4->Pz();
       E    = p4->Energy();
       M    = mass[id];

       new ( (*particle_list)[1+id] )
                         TMCParticle(1,pdgc[id],0,0,0,px,py,pz,E,M,0,0,0,0,0);
  }

  //-- Set owner and return
  particle_list->SetOwner(true);
  return particle_list;
}
//____________________________________________________________________________
double BaryonResonanceDecayer::Weight(void) const
{
  return fWeight;
}
//____________________________________________________________________________
void BaryonResonanceDecayer::InhibitDecay(int pdgc, TDecayChannel * dc) const
{
  if(! this->IsHandled(pdgc)) return;

  if(!dc) return;

}
//____________________________________________________________________________
void BaryonResonanceDecayer::UnInhibitDecay(int pdgc, TDecayChannel * dc) const
{
  if(! this->IsHandled(pdgc)) return;

  if(!dc) return;

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
         LOG("Decay", pNOTICE)
                  << "Disabling decay channel containing resonance 1114";;
         md = 999999999;
     }
     mass += md;
  }  
  return mass;
}
//____________________________________________________________________________
void BaryonResonanceDecayer::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void BaryonResonanceDecayer::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void BaryonResonanceDecayer::LoadConfig(void)
{
// Read configuration options or set defaults

  //-- Generated weighted or un-weighted hadronic systems
  fGenerateWeighted = fConfig->GetBoolDef("generate-weighted", false);
}
//____________________________________________________________________________
