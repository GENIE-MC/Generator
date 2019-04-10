//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Joshua Berger <jberger \at physics.wisc.edu>
         University of Wisconsin-Madison

         Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab
*/
//____________________________________________________________________________

#include "Framework/EventGen/InteractionList.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Physics/BoostedDarkMatter/EventGen/DMELInteractionListGenerator.h"

using namespace genie;

//___________________________________________________________________________
DMELInteractionListGenerator::DMELInteractionListGenerator() :
InteractionListGeneratorI("genie::DMELInteractionListGenerator")
{

}
//___________________________________________________________________________
DMELInteractionListGenerator::DMELInteractionListGenerator(string config):
InteractionListGeneratorI("genie::DMELInteractionListGenerator",  config)
{

}
//___________________________________________________________________________
DMELInteractionListGenerator::~DMELInteractionListGenerator()
{

}
//___________________________________________________________________________
InteractionList * DMELInteractionListGenerator::CreateInteractionList(
   const InitialState & init_state) const
{
  LOG("IntLst", pINFO)
     << "InitialState = " << init_state.AsString();

  if (fIsDM) 
     return this->CreateInteractionListDM(init_state);
  else {
     LOG("IntLst", pWARN)
       << "Unknown InteractionType! Returning NULL InteractionList "
       << "for init-state: " << init_state.AsString();
     return 0;
  }
  return 0;
}
//___________________________________________________________________________
InteractionList * DMELInteractionListGenerator::CreateInteractionListDM(
   const InitialState & init_state) const
{
  InteractionList * intlist = new InteractionList;

  int nuclpdg[2] = { kPdgProton, kPdgNeutron };

  int      dmpdg   = init_state.ProbePdg();
  bool     isdm    = pdg::IsDarkMatter     (dmpdg);

  if(!isdm) {
    LOG("IntLst", pWARN)
      << "Can not handle probe! Returning NULL InteractionList "
      << "for init-state: " << init_state.AsString();
    delete intlist;
    return 0;
  }

  for(int i=0; i<2; i++) {

     ProcessInfo   proc_info(kScDarkMatterElastic, kIntDarkMatter);
     Interaction * interaction = new Interaction(init_state, proc_info);

     Target * target  = interaction->InitStatePtr()->TgtPtr();
     bool     hasP    = (target->Z() > 0);
     bool     hasN    = (target->N() > 0);

     if(nuclpdg[i] == kPdgProton  && !hasP) {
       delete interaction;
       continue;
     }
     if(nuclpdg[i] == kPdgNeutron  && !hasN) {
       delete interaction;
       continue;
     }
     target->SetHitNucPdg(nuclpdg[i]);
     intlist->push_back(interaction);
  }

  if(intlist->size() == 0) {
     LOG("IntLst", pINFO)
       << "Returning NULL InteractionList for init-state: "
       << init_state.AsString();
     delete intlist;
     return 0;
  }
  return intlist;
}
//____________________________________________________________________________
void DMELInteractionListGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void DMELInteractionListGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void DMELInteractionListGenerator::LoadConfigData(void)
{
    GetParamDef( "is-DM", fIsDM, false  );
}
//____________________________________________________________________________

