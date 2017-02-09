//____________________________________________________________________________
/*
 Copyright (c) 2003-2017, GENIE Neutrino MC Generator Collaboration
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

#include "EVGCore/InteractionList.h"
#include "Interaction/Interaction.h"
#include "NeutronOsc/NOscDummyInteractionListGenerator.h"

using namespace genie;

//___________________________________________________________________________
NOscDummyInteractionListGenerator::NOscDummyInteractionListGenerator() :
InteractionListGeneratorI("genie::NOscDummyInteractionListGenerator")
{

}
//___________________________________________________________________________
NOscDummyInteractionListGenerator::NOscDummyInteractionListGenerator(string config):
InteractionListGeneratorI("genie::NOscDummyInteractionListGenerator", config)
{

}
//___________________________________________________________________________
NOscDummyInteractionListGenerator::~NOscDummyInteractionListGenerator()
{

}
//___________________________________________________________________________
InteractionList * NOscDummyInteractionListGenerator::CreateInteractionList(
    const InitialState & /*init_state*/) const
{
  return 0;
}
//___________________________________________________________________________
