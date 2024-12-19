//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool 
*/
//____________________________________________________________________________

#include "Framework/EventGen/InteractionList.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/AnomalyMediatedNuGamma/EventGen/AMNuGammaInteractionListGenerator.h"
#include "Framework/ParticleData/PDGCodes.h"

using namespace genie;

//___________________________________________________________________________
AMNuGammaInteractionListGenerator::AMNuGammaInteractionListGenerator() :
InteractionListGeneratorI(
	"genie::AMNuGammaInteractionListGenerator")
{

}
//___________________________________________________________________________
AMNuGammaInteractionListGenerator::AMNuGammaInteractionListGenerator(
 string config):
InteractionListGeneratorI(
	"genie::AMNuGammaInteractionListGenerator", config)
{

}
//___________________________________________________________________________
AMNuGammaInteractionListGenerator::~AMNuGammaInteractionListGenerator()
{

}
//___________________________________________________________________________
InteractionList * AMNuGammaInteractionListGenerator::CreateInteractionList(
                                       const InitialState & init_state) const
{
  LOG("IntLst", pINFO)
     << "InitialState = " << init_state.AsString();

  int nupdg   = init_state.ProbePdg();
  int tgtpdg  = init_state.Tgt().Pdg();

  InteractionList * intlist = new InteractionList;

  const Target & target = init_state.Tgt();

  if(target.Z()>0) {
    Interaction * interaction = Interaction::AMNuGamma(tgtpdg,kPdgProton,nupdg,0);
    intlist->push_back(interaction);
  }
  if(target.N()>0) {
    Interaction * interaction = Interaction::AMNuGamma(tgtpdg,kPdgNeutron,nupdg,0);
    intlist->push_back(interaction);
  }

  return intlist;
}
//___________________________________________________________________________
