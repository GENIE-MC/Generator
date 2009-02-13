//____________________________________________________________________________
/*
 Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - May 13, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
*/
//____________________________________________________________________________

#include "EVGCore/InteractionList.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "Diffractive/DfrcInteractionListGenerator.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"

using namespace genie;

//___________________________________________________________________________
DfrcInteractionListGenerator::DfrcInteractionListGenerator() :
InteractionListGeneratorI("genie::DfrcInteractionListGenerator")
{

}
//___________________________________________________________________________
DfrcInteractionListGenerator::DfrcInteractionListGenerator(string config):
InteractionListGeneratorI("genie::DfrcInteractionListGenerator", config)
{

}
//___________________________________________________________________________
DfrcInteractionListGenerator::~DfrcInteractionListGenerator()
{

}
//___________________________________________________________________________
InteractionList * DfrcInteractionListGenerator::CreateInteractionList(
                                       const InitialState & init_state) const
{
  LOG("IntLst", pINFO)
     << "InitialState = " << init_state.AsString();

  int  nupdg   = init_state.ProbePdg();
  int  tgtpdg  = init_state.Tgt().Pdg();
  bool hasP    = (init_state.Tgt().Z() > 0);
  bool hasN    = (init_state.Tgt().N() > 0);
  bool isnu    = pdg::IsNeutrino     (nupdg);
  bool isnubar = pdg::IsAntiNeutrino (nupdg);

  InteractionList * intlist = new InteractionList;

  if (isnu && hasN) {
     Interaction * interaction = Interaction::DFRCC(tgtpdg,kPdgNeutron,nupdg);
     intlist->push_back(interaction);
  } else if (isnubar && hasP) {
  
  }

  return intlist;
}
//___________________________________________________________________________
