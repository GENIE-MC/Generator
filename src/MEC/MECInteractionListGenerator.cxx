//____________________________________________________________________________
/*
 Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 22, 2008 - CA
   This interction list generator was first added in version 2.5.1 as part of
   the new event generation thread handling MEC interactions.
*/
//____________________________________________________________________________

#include "EVGCore/InteractionList.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "MEC/MECInteractionListGenerator.h"

using namespace genie;

//___________________________________________________________________________
MECInteractionListGenerator::MECInteractionListGenerator() :
InteractionListGeneratorI("genie::MECInteractionListGenerator")
{

}
//___________________________________________________________________________
MECInteractionListGenerator::MECInteractionListGenerator(string config):
InteractionListGeneratorI("genie::MECInteractionListGenerator", config)
{

}
//___________________________________________________________________________
MECInteractionListGenerator::~MECInteractionListGenerator()
{

}
//___________________________________________________________________________
InteractionList * 
  MECInteractionListGenerator::CreateInteractionList(
      const InitialState & init_state) const
{
  LOG("IntLst", pINFO) << "InitialState = " << init_state.AsString();

  int nupdg  = init_state.ProbePdg();
  int tgtpdg = init_state.Tgt().Pdg();

  const Target & target = init_state.Tgt();
  if(target.A() < 4) return 0;

  InteractionList * intlist = new InteractionList;

  const int nc = 3;
  const int nucleon_cluster[nc] = { 
    kPdgClusterNN, kPdgClusterNP, kPdgClusterPP };

  for(int ic = 0; ic < nc; ic++) {
     int ncpdg = nucleon_cluster[ic];
     if(fIsCC) {
       bool allowed = false;
       if(pdg::IsNeutrino(nupdg)) {
         // neutrino CC => final state primary lepton is -1
         // therefore the nucleon-cluster charge needs to be incremented by +1.
         if(ncpdg == kPdgClusterNN || ncpdg == kPdgClusterNP) {
            allowed = true;
         }
       }
       else 
       if(pdg::IsAntiNeutrino(nupdg)) {
         // anti-neutrino CC => final state primary lepton is +1
         // therefore the nucleon-cluster charge needs to be incremented by -1.
         if(ncpdg == kPdgClusterNP || ncpdg == kPdgClusterPP) {
            allowed = true;
         }
       }
       if(allowed) {
         Interaction * interaction = 
             Interaction::MECCC(tgtpdg,ncpdg,nupdg,0);
         intlist->push_back(interaction);
       }
     }//CC?

     else
     if(fIsNC) {
       Interaction * interaction = 
           Interaction::MECNC(tgtpdg,ncpdg,nupdg,0);
       intlist->push_back(interaction);
     }//NC?

     else
     if(fIsEM) {
       Interaction * interaction = 
           Interaction::MECEM(tgtpdg,ncpdg,nupdg,0);
       intlist->push_back(interaction);
     }//EM?

  }

  return intlist;
}
//___________________________________________________________________________
void MECInteractionListGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void MECInteractionListGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void MECInteractionListGenerator::LoadConfigData(void)
{
  fIsCC = fConfig->GetBoolDef("is-CC",  false);
  fIsNC = fConfig->GetBoolDef("is-NC",  false);
  fIsEM = fConfig->GetBoolDef("is-EM",  false);
}
//____________________________________________________________________________

