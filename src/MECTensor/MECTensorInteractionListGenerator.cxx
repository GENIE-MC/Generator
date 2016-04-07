//____________________________________________________________________________
/*
 Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory 
         Jackie Schwehr <jackie.schwehr \at colostate.edu>
         Colorado State University

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 22, 2008 - CA
   This interction list generator was first added in version 2.5.1 as part of
   the new event generation thread handling MEC interactions.
 @ Sep 14, 2014 - JS
   Adapted for the "Valencia" MEC model.
*/
//____________________________________________________________________________

#include "EVGCore/InteractionList.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "MECTensor/MECTensorInteractionListGenerator.h"

using namespace genie;

//___________________________________________________________________________
MECTensorInteractionListGenerator::MECTensorInteractionListGenerator() :
InteractionListGeneratorI("genie::MECTensorInteractionListGenerator")
{

}
//___________________________________________________________________________
MECTensorInteractionListGenerator::MECTensorInteractionListGenerator(string config):
InteractionListGeneratorI("genie::MECTensorInteractionListGenerator", config)
{

}
//___________________________________________________________________________
MECTensorInteractionListGenerator::~MECTensorInteractionListGenerator()
{

}
//___________________________________________________________________________
InteractionList * 
  MECTensorInteractionListGenerator::CreateInteractionList(
      const InitialState & init_state) const
{
  LOG("IntLst", pINFO) << "InitialState = " << init_state.AsString();

  int nupdg  = init_state.ProbePdg();
  int tgtpdg = init_state.Tgt().Pdg();

  const Target & target = init_state.Tgt();

  if(target.A() < 4) return NULL;
  // All other nuclei pass this spot.
  // If the nucleus is known, then we can generate events and splines.
  // If the nucleus is not known, cross sections and splines return zero.
  
  InteractionList * intlist = new InteractionList;

  Interaction * interaction = Interaction::MECCC(tgtpdg, nupdg, 0.0);
  intlist->push_back(interaction);

  return intlist;
}
//___________________________________________________________________________
void MECTensorInteractionListGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void MECTensorInteractionListGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void MECTensorInteractionListGenerator::LoadConfigData(void)
{
  fIsCC = fConfig->GetBoolDef("is-CC",  false);
  fIsNC = fConfig->GetBoolDef("is-NC",  false);
  fIsEM = fConfig->GetBoolDef("is-EM",  false);
}
//____________________________________________________________________________

