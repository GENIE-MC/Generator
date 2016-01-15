//____________________________________________________________________________
/*
 Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Dec 14, 2009 - CA
   Was first added in v2.5.1

*/
//____________________________________________________________________________

#include "EVGCore/InteractionList.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "VHE/GLRESInteractionListGenerator.h"

using namespace genie;

//___________________________________________________________________________
GLRESInteractionListGenerator::GLRESInteractionListGenerator() :
InteractionListGeneratorI("genie::GLRESInteractionListGenerator")
{

}
//___________________________________________________________________________
GLRESInteractionListGenerator::GLRESInteractionListGenerator(string config) :
InteractionListGeneratorI("genie::GLRESInteractionListGenerator", config)
{

}
//___________________________________________________________________________
GLRESInteractionListGenerator::~GLRESInteractionListGenerator()
{

}
//___________________________________________________________________________
InteractionList * 
   GLRESInteractionListGenerator::CreateInteractionList(
                                       const InitialState & init_state) const
{
// channels:
// nuebar + e- -> W-

  if(init_state.ProbePdg() != kPdgAntiNuE) {
     LOG("IntLst", pDEBUG) 
          << "Return *null* interaction list";
     return 0;
  }

  int target = init_state.Tgt().Pdg();

  InteractionList * intlist = new InteractionList;

  Interaction * interaction = Interaction::GLR(target);
  intlist->push_back(interaction);

  return intlist;
}
//___________________________________________________________________________
