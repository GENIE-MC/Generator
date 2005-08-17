//____________________________________________________________________________
/*!

\class   genie::RESHadronicSystemGenerator

\brief   Generates the 'final state' hadronic system in v RES interactions.

         It creates the GHepParticle entries for the target nucleus (if any)
         and the res. decay products and they are added to the GHEP record. \n

         The resonance decay products should be known at the time this visitor
         acts on th event record (this visitor should run on event generation
         threads initiated by selecting an exlusive RES channel).

         It does not handle the propagation of generated hadrons out of the
         nuclear medium and it does not handle decays of unstable particles
         (these would be handled by other event record visitors later on the
         event generation threads). So the 'final state' might not be final
         after all.

         Is a concrete implementation of the EventRecordVisitorI interface.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created November 23, 2004

*/
//____________________________________________________________________________

#include <TGenPhaseSpace.h>

#include "BaryonResonance/BaryonResonance.h"
#include "BaryonResonance/BaryonResUtils.h"
#include "Conventions/Constants.h"
#include "EVGModules/RESHadronicSystemGenerator.h"
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
RESHadronicSystemGenerator::RESHadronicSystemGenerator() :
HadronicSystemGenerator()
{
  fName = "genie::RESHadronicSystemGenerator";
}
//___________________________________________________________________________
RESHadronicSystemGenerator::RESHadronicSystemGenerator(
                                                    const char * param_set) :
HadronicSystemGenerator(param_set)
{
  fName = "genie::RESHadronicSystemGenerator";

  this->FindConfig();
}
//___________________________________________________________________________
RESHadronicSystemGenerator::~RESHadronicSystemGenerator()
{

}
//___________________________________________________________________________
void RESHadronicSystemGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
// This method generates the final state hadronic system

  //-- If the struck nucleon was within a nucleus, then add the final state
  //   nucleus at the EventRecord
  this->AddTargetNucleusRemnant(evrec);

  //-- Add the baryon resonance decay products at the event record
  this->AddResonanceDecayProducts(evrec);
}
//___________________________________________________________________________
void RESHadronicSystemGenerator::AddResonanceDecayProducts(
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

  int res_pos = 0;
  if(is_nucleus) res_pos = 4;
  else           res_pos = 3;

  GHepParticle * res  = evrec->GetParticle(res_pos);
  TLorentzVector * p4 = res->GetP4();

  LOG("RESHadronicVtx", pINFO)
                 << "\n RES 4-P = " << print_utils::P4AsString(p4);

  //-- generate 4-p for the two-hadron system
  double mnuc = PDGLibrary::Instance() -> Find(nuc_pdgc) -> Mass();
  double mpi  = PDGLibrary::Instance() -> Find(pi_pdgc)  -> Mass();

  double mass[2] = { mnuc, mpi };

  TGenPhaseSpace phase_space_generator;

  bool is_permitted = phase_space_generator.SetDecay(*p4, 2, mass);
  assert(is_permitted);

  phase_space_generator.Generate();

  //-- add the two hadrons at the event record
  TLorentzVector & p4_nuc = *phase_space_generator.GetDecay(0);
  TLorentzVector & p4_pi  = *phase_space_generator.GetDecay(1);
  TLorentzVector vdummy(0,0,0,0); // dummy 'vertex'

  int mom = res_pos;
  evrec->AddParticle(
               nuc_pdgc,kIStStableFinalState, mom,-1,-1,-1, p4_nuc, vdummy);
  evrec->AddParticle(
               pi_pdgc, kIStStableFinalState, mom,-1,-1,-1, p4_pi,  vdummy);

  delete p4;
}
//___________________________________________________________________________

