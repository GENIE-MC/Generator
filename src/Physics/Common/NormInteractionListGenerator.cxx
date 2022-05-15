//____________________________________________________________________________
/*
 Copyright (c) 2003-2022, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Igor Kakorin <kakorin@jinr.ru>
 Joint Institute for Nuclear Research
*/
//____________________________________________________________________________

#include <TMath.h>

#include "Physics/Common/NormInteractionListGenerator.h"
#include "Framework/EventGen/InteractionList.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/ParticleData/PDGCodes.h"

using namespace genie;

//___________________________________________________________________________
NormInteractionListGenerator::NormInteractionListGenerator() :
InteractionListGeneratorI("genie::NormInteractionListGenerator")
{

}
//___________________________________________________________________________
NormInteractionListGenerator::NormInteractionListGenerator(string config) :
InteractionListGeneratorI("genie::NormInteractionListGenerator", config)
{

}
//___________________________________________________________________________
NormInteractionListGenerator::~NormInteractionListGenerator()
{

}
//___________________________________________________________________________
InteractionList * NormInteractionListGenerator::CreateInteractionList(
                                      const InitialState & init_state) const
{
  InteractionList * intlist = new InteractionList;
  bool isNC = init_state.ProbePdg()==kPdgNuE      || init_state.ProbePdg()==kPdgAntiNuE  || init_state.ProbePdg()==kPdgNuMu || init_state.ProbePdg()==kPdgAntiNuMu || init_state.ProbePdg()==kPdgNuTau || init_state.ProbePdg()==kPdgAntiNuTau;
  bool isEM = init_state.ProbePdg()==kPdgElectron || init_state.ProbePdg()==kPdgPositron || init_state.ProbePdg()==kPdgMuon || init_state.ProbePdg()==kPdgAntiMuon || init_state.ProbePdg()==kPdgTau   || init_state.ProbePdg()==kPdgAntiTau;
  if (!isNC && !isEM)
  {
      LOG("IntLst", pWARN)
        << "Unknown InteractionType! Returning NULL InteractionList "
        << "for init-state: " << init_state.AsString();
      return 0;
  }
  
  ProcessInfo proc_info(kScNorm,  isNC?kIntWeakNC:kIntEM);
  Interaction * interaction = new Interaction(init_state, proc_info);

 
  intlist->push_back(interaction);

  return intlist;
}
//___________________________________________________________________________
void NormInteractionListGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void NormInteractionListGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void NormInteractionListGenerator::LoadConfigData(void)
{

}
//____________________________________________________________________________
