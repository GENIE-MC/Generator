//____________________________________________________________________________
/*!

\class   genie::RSPPHadronicSystemGenerator

\brief   Generates the 'final state' hadronic system in v SPP interactions.
         
         It adds the remnant nucleus (if any) and the baryon resonance decay 
         products at the GHEP record. The resonance decay products are pre-
         determined since in this thread we generate exclusive SPP reactions.
         The module uses a simple phase space decay.
         Is a concrete implementation of the EventRecordVisitorI interface.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created November 23, 2004

*/
//____________________________________________________________________________

#include "BaryonResonance/BaryonResonance.h"
#include "BaryonResonance/BaryonResUtils.h"
#include "Conventions/Constants.h"
#include "EVGModules/RSPPHadronicSystemGenerator.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepParticle.h"
#include "Interaction/Interaction.h"
#include "Interaction/SppChannel.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "Utils/PrintUtils.h"

using namespace genie;
using namespace genie::constants;

//___________________________________________________________________________
RSPPHadronicSystemGenerator::RSPPHadronicSystemGenerator() :
HadronicSystemGenerator("genie::RSPPHadronicSystemGenerator")
{

}
//___________________________________________________________________________
RSPPHadronicSystemGenerator::RSPPHadronicSystemGenerator(string config):
HadronicSystemGenerator("genie::RSPPHadronicSystemGenerator", config)
{

}
//___________________________________________________________________________
RSPPHadronicSystemGenerator::~RSPPHadronicSystemGenerator()
{

}
//___________________________________________________________________________
void RSPPHadronicSystemGenerator::ProcessEventRecord(GHepRecord* evrec) const
{
// This method generates the final state hadronic system

  //-- If the struck nucleon was within a nucleus, then add the final state
  //   nucleus at the EventRecord
  this->AddTargetNucleusRemnant(evrec);

  //-- Add the baryon resonance decay products at the event record
  this->AddResonanceDecayProducts(evrec);
}
//___________________________________________________________________________
void RSPPHadronicSystemGenerator::AddResonanceDecayProducts(
                                                   GHepRecord * evrec) const
{
// generate momenta for the baryon resonance decay products and add them at
// the event record

  //-- find out which SPP channel we are generating
  Interaction * interaction = evrec->GetInteraction();
  SppChannel_t spp_channel = SppChannel::FromInteraction(interaction);

  //-- get the final state nucleon and pion
  int nuc_pdgc = SppChannel::FinStateNucleon (spp_channel);
  int pi_pdgc  = SppChannel::FinStatePion    (spp_channel);

  //-- get the total 4-p for the two-hadron system (= parent resonance 4-p)

  const InitialState & init_state = interaction->GetInitialState();
  bool is_nucleus = init_state.GetTarget().IsNucleus();

  //-- access the resonance entry at the GHEP record
  int res_pos = 0;
  if(is_nucleus) res_pos = 4;
  else           res_pos = 3;

  GHepParticle * res  = evrec->Particle(res_pos);
  
  //-- mark the resonance as decayed
  res->SetStatus(kIstDecayedState);

  //-- generate 4-p for the two-hadron system
  double mnuc = PDGLibrary::Instance() -> Find(nuc_pdgc) -> Mass();
  double mpi  = PDGLibrary::Instance() -> Find(pi_pdgc)  -> Mass();

  double mass[2] = { mnuc, mpi };

  TLorentzVector * p4 = res->GetP4();

  LOG("RESHadronicVtx", pINFO)
                 << "\n RES 4-P = " << utils::print::P4AsString(p4);

  bool is_permitted = fPhaseSpaceGenerator.SetDecay(*p4, 2, mass);
  assert(is_permitted);

  fPhaseSpaceGenerator.Generate();

  //-- add the two hadrons at the event record
  TLorentzVector & p4_nuc = *fPhaseSpaceGenerator.GetDecay(0);
  TLorentzVector & p4_pi  = *fPhaseSpaceGenerator.GetDecay(1);
  TLorentzVector vdummy(0,0,0,0); // dummy 'vertex'

  int mom = res_pos;
  evrec->AddParticle(
               nuc_pdgc,kIStStableFinalState, mom,-1,-1,-1, p4_nuc, vdummy);
  evrec->AddParticle(
               pi_pdgc, kIStStableFinalState, mom,-1,-1,-1, p4_pi,  vdummy);
  delete p4;
}
//___________________________________________________________________________

