//____________________________________________________________________________
/*!

\class    genie::BaryonResonanceDecayer

\brief    Baryon resonance decayer.

          A simple decayer based on resonance's branching fractions (BRs) and
          an N-body phase space generator. Since the resonance can be produced
          off-shell, decay channels with total-mass > W are suppressed. \n

          Unlike PythiaDecayer, which relies on PYTHIA, this algorithm is not
          based on any rich underlying model. Should be used for decaying
          baryon resonances only. \n

          The baryon resonance PDG codes follow the MINOS extensions to
          PDG tables. \n

          The BaryonResonanceDecayer is a concrete implementation of the
          DecayModelI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  November 27, 2004

*/
//____________________________________________________________________________

#include <TGenPhaseSpace.h>
#include <TMCParticle6.h>

#include "BaryonResonance/BaryonResUtils.h"
#include "Decay/BaryonResonanceDecayer.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGLibrary.h"
#include "Numerical/RandomGen.h"
#include "Utils/PrintUtils.h"

using namespace genie;

//____________________________________________________________________________
BaryonResonanceDecayer::BaryonResonanceDecayer() :
DecayModelI()
{
  fName     = "genie::BaryonResonanceDecayer";
}
//____________________________________________________________________________
BaryonResonanceDecayer::BaryonResonanceDecayer(const char * param_set) :
DecayModelI(param_set)
{
  fName = "genie::BaryonResonanceDecayer";

  FindConfig();
}
//____________________________________________________________________________
BaryonResonanceDecayer::~BaryonResonanceDecayer()
{

}
//____________________________________________________________________________
bool BaryonResonanceDecayer::IsHandled(int pdg_code) const
{
// handles only requests to decay baryon resonances

  if( res_utils::IsBaryonResonance(pdg_code) ) {

    LOG("Decay", pINFO) 
         << "\n *** The particle with PDG-Code = "
                           << pdg_code << " is not decayed by this algorithm";
    return true;
    
  } else return false;
}
//____________________________________________________________________________
TClonesArray* BaryonResonanceDecayer::Decay(const DecayerInputs_t & inp) const
{  
  if ( ! this->IsHandled(inp.PdgCode) ) return 0;

  //-- Find the particle in the PDG library

  TParticlePDG * mother = PDGLibrary::Instance()->Find(inp.PdgCode);

  if(!mother) {

     LOG("Decay", pERROR)
          << "\n *** The particle with PDG-Code = " << inp.PdgCode
                                         << " was not found in PDGLibrary";
     return 0;                               
  }  

  LOG("Decay", pINFO)
           << "Decaying resonance = " << mother->GetName()
                         << " with P4 = " << print_utils::P4AsString(inp.P4);
  
  //-- Get the resonance mass W (generally different from the mass associated
  //   with the input pdg_code)

  double W = inp.P4->M();
  
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
     
     if(fsmass < W) tot_BR += ch->BranchingRatio();
     else {       
       tot_BR += 0;

       SLOG("Decay", pINFO)
               << "Suppresing channel: " << ich 
                        << " with final state mass = " << fsmass << " GeV";         
     }
     BR[ich] = tot_BR;
  }

  //-- Select a resonance based on the branching ratios
  
  unsigned int ich = 0, sel_ich; // id of selected decay channel

  RandomGen * rnd = RandomGen::Instance();

  double x = tot_BR * rnd->Random2().Rndm();

  do {
    
    sel_ich = ich;
    
  } while (x > BR[ich++]);

  LOG("Decay", pDEBUG)
       << "\n Selected decay channel: " << sel_ich
               << "\n rnd = " << x << ", SumBR(RES : 0 -> "
                                       << sel_ich << ") = " << BR[sel_ich];

  TDecayChannel * ch = (TDecayChannel *) decay_list->At(sel_ich);

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
           << "Adding daughter[" << iparticle << "]: PDGC = "
                     << pdgc[iparticle] << ", mass = " << mass[iparticle];
  }

  //-- Decay the resonance using an N-body phase space generator
  //   The particle will be decayed in its rest frame and then the daughters
  //   will be boosted back to the original frame.

  TGenPhaseSpace phase_space_generator;

  bool is_permitted = phase_space_generator.SetDecay(p, nd, mass);

  assert(is_permitted);

  phase_space_generator.Generate();

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

       TLorentzVector * p4 = phase_space_generator.GetDecay(id);

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
double BaryonResonanceDecayer::FinalStateMass(TDecayChannel * channel) const
{
// Computes the total mass of the final state system

  double mass = 0;
  
  unsigned int nd = channel->NDaughters();

  for(unsigned int iparticle = 0; iparticle < nd; iparticle++) {

     int daughter_code = channel->DaughterPdgCode(iparticle);

     TParticlePDG * daughter = PDGLibrary::Instance()->Find(daughter_code);

     assert(daughter);

     mass += ( daughter->Mass() );
  }  
  return mass;
}
//____________________________________________________________________________
