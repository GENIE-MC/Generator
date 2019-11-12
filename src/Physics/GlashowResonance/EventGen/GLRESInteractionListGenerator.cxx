//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Dec 14, 2009 - CA
   Was first added in v2.5.1
 @ November 08, 2019 - Alfonso Garcia
   Modified to generate the kinematics of outgoing lepton properly.
   Phys. Rev. D 22, 2122 – Published 1 November 1980

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
// nuebar + e- -> W- -> nuebar + e-
// nuebar + e- -> W- -> nuebar + mu-
// nuebar + e- -> W- -> nuebar + tau-
// nuebar + e- -> W- -> hadrons

  if(init_state.ProbePdg() != kPdgAntiNuE) {
     LOG("IntLst", pDEBUG) 
          << "Return *null* interaction list";
     return 0;
  }

  InitialState init(init_state);
  init_state.TgtPtr()->SetHitNucPdg(0);  

  ProcessInfo   proc_info(kScGlashowResonance, kIntWeakCC);

  InteractionList * intlist = new InteractionList;

  if (fIsMu) {
    Interaction * interaction = new Interaction(init_state, proc_info);
    XclsTag exclusive_tag;
    exclusive_tag.SetFinalLepton(kPdgMuon);
    interaction->SetExclTag(exclusive_tag);
    intlist->push_back(interaction);  
  }
  else if (fIsTau) {
    Interaction * interaction = new Interaction(init_state, proc_info);
    XclsTag exclusive_tag;
    exclusive_tag.SetFinalLepton(kPdgTau);
    interaction->SetExclTag(exclusive_tag);
    intlist->push_back(interaction);  
  }
  else if (fIsEle) {
    Interaction * interaction = new Interaction(init_state, proc_info);
    XclsTag exclusive_tag;
    exclusive_tag.SetFinalLepton(kPdgElectron);
    interaction->SetExclTag(exclusive_tag);
    intlist->push_back(interaction);  
  }
  else if (fIsHad) {
    Interaction * interaction = new Interaction(init_state, proc_info);
    XclsTag exclusive_tag;
    exclusive_tag.SetFinalLepton(kPdgPiP);
    interaction->SetExclTag(exclusive_tag);
    intlist->push_back(interaction);  
  }

  if(intlist->size() == 0) {
     LOG("IntLst", pERROR)
         << "Returning NULL InteractionList for init-state: "
                                                  << init_state.AsString();
     delete intlist;
     return 0;
  }

  return intlist;
}
//___________________________________________________________________________
void GLRESInteractionListGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void GLRESInteractionListGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void GLRESInteractionListGenerator::LoadConfigData(void)
{

  GetParamDef("is-Mu",  fIsMu,  false ) ;
  GetParamDef("is-Tau", fIsTau, false ) ;
  GetParamDef("is-Ele", fIsEle, false ) ;
  GetParamDef("is-Had", fIsHad, false ) ;

}