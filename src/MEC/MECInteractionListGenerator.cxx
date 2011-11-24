//____________________________________________________________________________
/*
 Copyright (c) 2003-2011, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 22, 2008 - CA
   This interction list generator was first added in version 2.5.1 as part of
   the new event generation thread handling MEC interactions.
*/
//____________________________________________________________________________

#include "EVGCore/InteractionList.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "MEC/MECInteractionListGenerator.h"

using namespace genie;

//___________________________________________________________________________
MECInteractionListGenerator::MECInteractionListGenerator() :
InteractionListGeneratorI("genie::MECInteractionListGenerator")
{

}
//___________________________________________________________________________
MECInteractionListGenerator::MECInteractionListGenerator(string config):
InteractionListGeneratorI("genie::MECInteractionListGenerator", config)
{

}
//___________________________________________________________________________
MECInteractionListGenerator::~MECInteractionListGenerator()
{

}
//___________________________________________________________________________
InteractionList * 
  MECInteractionListGenerator::CreateInteractionList(
      const InitialState & init_state) const
{
  LOG("IntLst", pINFO) << "InitialState = " << init_state.AsString();

  int nupdg  = init_state.ProbePdg();
  int tgtpdg = init_state.Tgt().Pdg();

  const Target & target = init_state.Tgt();
  if(target.A() < 4) return 0;

  InteractionList * intlist = new InteractionList;

  Interaction * interaction = Interaction::MECCC(tgtpdg,nupdg,0);
  intlist->push_back(interaction);

  return intlist;
}
//___________________________________________________________________________
