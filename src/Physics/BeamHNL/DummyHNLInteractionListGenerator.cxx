//____________________________________________________________________________
/*
  Implementation file for DummyHNLInteractionGenerator.h
*/
//____________________________________________________________________________

#include "Framework/EventGen/InteractionList.h"
#include "Framework/Interaction/Interaction.h"
#include "DummyHNLInteractionListGenerator.h"

using namespace genie;

//___________________________________________________________________________
DummyHNLInteractionListGenerator::DummyHNLInteractionListGenerator() :
InteractionListGeneratorI("genie::DummyHNLInteractionListGenerator")
{

}
//___________________________________________________________________________
DummyHNLInteractionListGenerator::DummyHNLInteractionListGenerator(string config):
InteractionListGeneratorI("genie::DummyHNLInteractionListGenerator", config)
{

}
//___________________________________________________________________________
DummyHNLInteractionListGenerator::~DummyHNLInteractionListGenerator()
{

}
//___________________________________________________________________________
InteractionList * DummyHNLInteractionListGenerator::CreateInteractionList(
    const InitialState & /*init_state*/) const
{
  return 0;
}
//___________________________________________________________________________
