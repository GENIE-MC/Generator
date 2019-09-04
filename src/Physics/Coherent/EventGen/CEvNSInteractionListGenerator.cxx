//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Author: Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
         University of Liverpool & STFC Rutherford Appleton Laboratory

 For the class documentation see the corresponding header file.
*/
//____________________________________________________________________________

#include "Framework/EventGen/InteractionList.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/Coherent/EventGen/CEvNSInteractionListGenerator.h"

using namespace genie;

//___________________________________________________________________________
CEvNSInteractionListGenerator::CEvNSInteractionListGenerator() :
InteractionListGeneratorI("genie::CEvNSInteractionListGenerator")
{

}
//___________________________________________________________________________
CEvNSInteractionListGenerator::CEvNSInteractionListGenerator(string config) :
InteractionListGeneratorI("genie::CEvNSInteractionListGenerator", config)
{

}
//___________________________________________________________________________
CEvNSInteractionListGenerator::~CEvNSInteractionListGenerator()
{

}
//___________________________________________________________________________
InteractionList * CEvNSInteractionListGenerator::CreateInteractionList(
  const InitialState & init_state) const
{
  LOG("IntLst", pINFO)
      << "InitialState = " << init_state.AsString();

  const Target & target = init_state.Tgt();
  if(!target.IsNucleus()) {
    LOG("IntLst", pINFO)
      << "Not a nuclear target! Returning NULL CEvNS interaction list "
      << "for init-state: " << init_state.AsString();
    return 0;
  }

  InteractionList * intlist = new InteractionList;

  ProcessInfo proc_info(kScCoherentElastic, kIntWeakNC);
  Interaction * interaction = new Interaction(init_state, proc_info);

  intlist->push_back(interaction);

  return intlist;
}
//___________________________________________________________________________
