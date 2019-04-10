//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - May 13, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Feb 25, 2008 - CA
   This interction list generator was first added in version 2.3.1 as part of
   the new event generation thread handling amonaly-mediated single gamma
   interactions.
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
