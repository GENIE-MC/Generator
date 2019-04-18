//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - November 23, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Feb 07, 2009 - CA
   Removed call to AddTargetNucleusRemnant(). This simulation step is now
   performed further upstream in the processing chain.
 @ Mar 03, 2009 - CA
   Moved into the new RES package from its previous location (EVGModules).

*/
//____________________________________________________________________________

#include "Framework/Conventions/Constants.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Interaction/SppChannel.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/BaryonResonance.h"
#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Physics/Resonance/EventGen/RSPPHadronicSystemGenerator.h"

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
  Interaction * interaction = evrec->Summary();
  SppChannel_t spp_channel = SppChannel::FromInteraction(interaction);

  //-- get the final state nucleon and pion
  int nuc_pdgc = SppChannel::FinStateNucleon (spp_channel);
  int pi_pdgc  = SppChannel::FinStatePion    (spp_channel);

  //-- get the total 4-p for the two-hadron system (= parent resonance 4-p)

  const InitialState & init_state = interaction->InitState();
  bool is_nucleus = init_state.Tgt().IsNucleus();

  //-- access the resonance entry at the GHEP record
  int res_pos = 0;
  if(is_nucleus) res_pos = 4;
  else           res_pos = 3;

  GHepParticle * res  = evrec->Particle(res_pos);
  const TLorentzVector & x4 = *(res->X4());  

  //-- mark the resonance as decayed
  res->SetStatus(kIStDecayedState);

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

  // decide the particle status
  GHepStatus_t ist = (is_nucleus) ? 
                              kIStHadronInTheNucleus : kIStStableFinalState;
  int mom = res_pos;
  evrec->AddParticle(nuc_pdgc, ist, mom,-1,-1,-1, p4_nuc, x4);
  evrec->AddParticle(pi_pdgc,  ist, mom,-1,-1,-1, p4_pi,  x4);
  delete p4;
}
//___________________________________________________________________________

