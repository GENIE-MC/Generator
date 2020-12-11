//____________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory
*/
//____________________________________________________________________________

#include "Framework/EventGen/InteractionList.h"
#include "Framework/Interaction/Interaction.h"
#include "Physics/NucleonDecay/DummyInteractionListGenerator.h"

using namespace genie;

//___________________________________________________________________________
DummyInteractionListGenerator::DummyInteractionListGenerator() :
InteractionListGeneratorI("genie::DummyInteractionListGenerator")
{

}
//___________________________________________________________________________
DummyInteractionListGenerator::DummyInteractionListGenerator(string config):
InteractionListGeneratorI("genie::DummyInteractionListGenerator", config)
{

}
//___________________________________________________________________________
DummyInteractionListGenerator::~DummyInteractionListGenerator()
{

}
//___________________________________________________________________________
InteractionList * DummyInteractionListGenerator::CreateInteractionList(
    const InitialState & /*init_state*/) const
{
  return 0;
}
//___________________________________________________________________________
