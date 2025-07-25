//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Jeremy Hewes, Georgia Karagiorgi
 University of Manchester
*/
//____________________________________________________________________________

#include "Framework/EventGen/InteractionList.h"
#include "Framework/Interaction/Interaction.h"
#include "Physics/NNBarOscillation/NNBarOscDummyInteractionListGenerator.h"

using namespace genie;

//___________________________________________________________________________
NNBarOscDummyInteractionListGenerator::NNBarOscDummyInteractionListGenerator() :
InteractionListGeneratorI("genie::NNBarOscDummyInteractionListGenerator")
{

}
//___________________________________________________________________________
NNBarOscDummyInteractionListGenerator::NNBarOscDummyInteractionListGenerator(string config):
InteractionListGeneratorI("genie::NNBarOscDummyInteractionListGenerator", config)
{

}
//___________________________________________________________________________
NNBarOscDummyInteractionListGenerator::~NNBarOscDummyInteractionListGenerator()
{

}
//___________________________________________________________________________
InteractionList * NNBarOscDummyInteractionListGenerator::CreateInteractionList(
    const InitialState & /*init_state*/) const
{
  return 0;
}
//___________________________________________________________________________
