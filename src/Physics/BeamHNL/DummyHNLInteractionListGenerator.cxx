//____________________________________________________________________________
/*
  Copyright (c) 2003-2023, The GENIE Collaboration
  For the full text of the license visit http://copyright.genie-mc.org

  Author: John Plows <komninos-john.plows \at physics.ox.ac.uk>
          University of Oxford

          Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
	  University of Liverpool & STFC Rutherford Appleton Laboratory
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
