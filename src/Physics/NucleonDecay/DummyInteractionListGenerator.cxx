//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - Sep 03, 2008

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Nov 10, 2011 - CA
   First added in v2.7.1.

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
