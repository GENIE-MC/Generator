//____________________________________________________________________________
/*
 Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Corey Reed <cjreed \at nikhef.nl>
         Nikhef

 For the class documentation see the corresponding header file.
*/
//____________________________________________________________________________

#include "EVGCore/InteractionList.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "VLE/IBDInteractionListGenerator.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"

using namespace genie;

//___________________________________________________________________________
IBDInteractionListGenerator::IBDInteractionListGenerator() :
   InteractionListGeneratorI(
      "genie::IBDInteractionListGenerator")
{

}
//___________________________________________________________________________
IBDInteractionListGenerator::IBDInteractionListGenerator(
   string config):
   InteractionListGeneratorI(
      "genie::IBDInteractionListGenerator", config)
{

}
//___________________________________________________________________________
IBDInteractionListGenerator::~IBDInteractionListGenerator()
{

}
//___________________________________________________________________________
InteractionList * IBDInteractionListGenerator::CreateInteractionList(
   const InitialState & init_state) const
{
   LOG("IBD", pINFO) << "InitialState = " << init_state.AsString();

   const int  nupdg   = init_state.ProbePdg();
   const bool isnu    = pdg::IsNeutrino     (nupdg);
   const bool isnubar = pdg::IsAntiNeutrino (nupdg);
   
   const Target& target = init_state.Tgt();
   const int     tgtpdg = target.Pdg();
   const bool    hasP   = (target.Z() > 0);
   const bool    hasN   = (target.N() > 0);
   
   InteractionList * intlist = 0;

   if (isnu && hasN) {
      intlist = new InteractionList;
      intlist->push_back( Interaction::IBD(tgtpdg,kPdgNeutron, nupdg,0) );
   } else if (isnubar && hasP) {
      intlist = new InteractionList;
      intlist->push_back( Interaction::IBD(tgtpdg,kPdgProton,nupdg,0) );
   } else {
      LOG("IBD", pWARN)
	 << "Returning NULL InteractionList for init-state: "
	 << init_state.AsString();
   }

   return intlist;
}
//____________________________________________________________________________
void IBDInteractionListGenerator::Configure(const Registry & config)
{
   Algorithm::Configure(config);
   this->LoadConfigData();
}
//____________________________________________________________________________
void IBDInteractionListGenerator::Configure(string config)
{
   Algorithm::Configure(config);
   this->LoadConfigData();
}
//____________________________________________________________________________
void IBDInteractionListGenerator::LoadConfigData(void)
{
   // (no data to be loaded)
}
//____________________________________________________________________________
