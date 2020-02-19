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
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Physics/GlashowResonance/EventGen/GLRESInteractionListGenerator.h"

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
