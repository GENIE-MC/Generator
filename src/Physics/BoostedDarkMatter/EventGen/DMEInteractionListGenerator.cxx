//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Feb 09, 2009 - CA
   Moved into the DME package from its previous location (EVGModules package)
 @ Feb 12, 2013 - CA (code from Rosen Matev)
   Handle the IMD annihilation channel.

*/
//____________________________________________________________________________

#include "Physics/BoostedDarkMatter/EventGen/DMEInteractionListGenerator.h"
#include "Framework/EventGen/InteractionList.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"

using namespace genie;

//___________________________________________________________________________
DMEInteractionListGenerator::DMEInteractionListGenerator() :
InteractionListGeneratorI("genie::DMEInteractionListGenerator")
{

}
//___________________________________________________________________________
DMEInteractionListGenerator::DMEInteractionListGenerator(string config) :
InteractionListGeneratorI("genie::DMEInteractionListGenerator", config)
{

}
//___________________________________________________________________________
DMEInteractionListGenerator::~DMEInteractionListGenerator()
{

}
//___________________________________________________________________________
InteractionList * DMEInteractionListGenerator::CreateInteractionList(
                                       const InitialState & init_state) const
{
  LOG("IntLst", pINFO) << "InitialState = " << init_state.AsString();

  return this -> DMEELInteractionList (init_state);
}
//___________________________________________________________________________
InteractionList * DMEInteractionListGenerator::DMEELInteractionList(
                                       const InitialState & init_state) const
{
// channels:
// DM       + e- -> DM       + e-   
// DMbar    + e- -> DMbar   + e-   [CC + NC + interference]

  int nupdg = init_state.ProbePdg();
  InteractionList * intlist = new InteractionList;

  // clone init state and de-activate the struck nucleon info
  InitialState init(init_state);
  init_state.TgtPtr()->SetHitNucPdg(0);

  if(nupdg == kPdgDarkMatter  || nupdg == kPdgAntiDarkMatter) {
     ProcessInfo   proc_info(kScDarkMatterElectron,  kIntDarkMatter);
     Interaction * interaction = new Interaction(init, proc_info);
     intlist->push_back(interaction);
  }

  return intlist;
}
//___________________________________________________________________________
void DMEInteractionListGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DMEInteractionListGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DMEInteractionListGenerator::LoadConfig(void)
{
  
}
//____________________________________________________________________________

