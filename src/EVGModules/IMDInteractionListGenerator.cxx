//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - July 13, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "EVGModules/IMDInteractionListGenerator.h"
#include "EVGCore/InteractionList.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"

using namespace genie;

//___________________________________________________________________________
IMDInteractionListGenerator::IMDInteractionListGenerator() :
InteractionListGeneratorI("genie::IMDInteractionListGenerator")
{

}
//___________________________________________________________________________
IMDInteractionListGenerator::IMDInteractionListGenerator(string config) :
InteractionListGeneratorI("genie::IMDInteractionListGenerator", config)
{

}
//___________________________________________________________________________
IMDInteractionListGenerator::~IMDInteractionListGenerator()
{

}
//___________________________________________________________________________
InteractionList * IMDInteractionListGenerator::CreateInteractionList(
                                       const InitialState & init_state) const
{
  LOG("InteractionList", pINFO) << "InitialState = " << init_state.AsString();

  if(init_state.ProbePdg() != kPdgNuMu) {
     LOG("InteractionList", pDEBUG) 
          << "Return *null* interaction list (non nu_mu probe in IMD thread)";
     return 0;
  }

  InteractionList * intlist = new InteractionList;

  // clone init state and de-activate the struck nucleon info
  InitialState init(init_state);
  init_state.TgtPtr()->SetHitNucPdg(0);

  ProcessInfo   proc_info(kScInverseMuDecay, kIntWeakCC);
  Interaction * interaction = new Interaction(init, proc_info);

  intlist->push_back(interaction);

  return intlist;
}
//___________________________________________________________________________
