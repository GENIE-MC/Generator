//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - December 19, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "EVGModules/COHInteractionListGenerator.h"
#include "EVGCore/InteractionList.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"

using namespace genie;

//___________________________________________________________________________
COHInteractionListGenerator::COHInteractionListGenerator() :
InteractionListGeneratorI("genie::COHInteractionListGenerator")
{

}
//___________________________________________________________________________
COHInteractionListGenerator::COHInteractionListGenerator(string config) :
InteractionListGeneratorI("genie::COHInteractionListGenerator", config)
{

}
//___________________________________________________________________________
COHInteractionListGenerator::~COHInteractionListGenerator()
{

}
//___________________________________________________________________________
InteractionList * COHInteractionListGenerator::CreateInteractionList(
                                      const InitialState & init_state) const
{
  LOG("InteractionList", pINFO)
                           << "InitialState = " << init_state.AsString();

  InteractionType_t inttype;
  if      (fIsCC) inttype = kIntWeakCC;
  else if (fIsNC) inttype = kIntWeakNC;
  else {
     LOG("InteractionList", pWARN)
       << "Unknown InteractionType! Returning NULL InteractionList "
                         << "for init-state: " << init_state.AsString();
     return 0;
  }

  int nupdg = init_state.ProbePdg();
  if( !pdg::IsNeutrino(nupdg) && !pdg::IsAntiNeutrino(nupdg) ) {
     LOG("InteractionList", pWARN)
       << "Can not handle probe! Returning NULL InteractionList "
                         << "for init-state: " << init_state.AsString();
     return 0;
  }
  const Target & target = init_state.Tgt();
  if(!target.IsNucleus()) {
     LOG("InteractionList", pWARN)
       << "Not a nuclear target! Returning NULL InteractionList "
                         << "for init-state: " << init_state.AsString();
     return 0;
  }

  InteractionList * intlist = new InteractionList;

  ProcessInfo proc_info(kScCoherent, inttype);
  Interaction * interaction = new Interaction(init_state, proc_info);

  if(fIsCC) {
    if(pdg::IsNeutrino(nupdg)) {
        //v A -> l- A pi+
        interaction->ExclTagPtr()->SetNPions(1,0,0);  
    } else {
        //vbar A -> l+ A pi-
        interaction->ExclTagPtr()->SetNPions(0,0,1); 
    }
  }
  else {interaction->ExclTagPtr()->SetNPions(0,1,0);} //v A -> v A pi0

  intlist->push_back(interaction);

  return intlist;
}
//___________________________________________________________________________
void COHInteractionListGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void COHInteractionListGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void COHInteractionListGenerator::LoadConfigData(void)
{
  fIsCC = fConfig->GetBoolDef("is-CC", false);
  fIsNC = fConfig->GetBoolDef("is-NC", false);
}
//____________________________________________________________________________

