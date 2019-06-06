//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Authors: Chris Marshall <marshall \at pas.rochester.edu>
          University of Rochester

          Martti Nirkko
          University of Berne

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include "Framework/EventGen/InteractionList.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Physics/Strange/EventGen/SKInteractionListGenerator.h"

using namespace genie;

//___________________________________________________________________________
SKInteractionListGenerator::SKInteractionListGenerator() :
InteractionListGeneratorI("genie::SKInteractionListGenerator")
{

}
//___________________________________________________________________________
SKInteractionListGenerator::SKInteractionListGenerator(string config) :
InteractionListGeneratorI("genie::SKInteractionListGenerator", config)
{

}
//___________________________________________________________________________
SKInteractionListGenerator::~SKInteractionListGenerator()
{

}
//___________________________________________________________________________
InteractionList * SKInteractionListGenerator::CreateInteractionList(
                                      const InitialState & init_state) const
{
  LOG("IntLst", pINFO)
      << "InitialState = " << init_state.AsString();

  if (fIsNC) {
    // deltaS = deltaQ and deltaS = 1 for this process -- no NC
    LOG("IntLst", pWARN)
      << "Interaction type is NC for deltaS = 1 process! Returning NULL InteractionList "
       << "for init-state: " << init_state.AsString();
     return 0;
  }
  else if (!fIsCC) {
     // shouldn't happen... warn
     LOG("IntLst", pWARN)
       << "Unknown InteractionType! Returning NULL InteractionList "
       << "for init-state: " << init_state.AsString();
     return 0;
  }

  int probe_pdg = init_state.ProbePdg();
  bool isnu = pdg::IsNeutrino(probe_pdg);
  if( !isnu ) {
     // shouldn't happen... warn
     LOG("IntLst", pWARN)
       << "Can not handle probe! Returning NULL InteractionList "
       << "for init-state: " << init_state.AsString();
     return 0;
  }

  InteractionList * intlist = new InteractionList;

  const int nch = 3;
  int inuclpdg[nch] = {0}; // hit nucleon pdg
  int fnuclpdg[nch] = {0}; // FS nucleon pdg
  int kaonpdg[nch] = {0}; // FS kaon pdg
  if( pdg::IsNeutrino(probe_pdg) ) {
    inuclpdg[0] = kPdgProton; inuclpdg[1] = kPdgNeutron; inuclpdg[2] = kPdgNeutron;
    fnuclpdg[0] = kPdgProton; fnuclpdg[1] = kPdgNeutron; fnuclpdg[2] = kPdgProton;
    kaonpdg[0]  = kPdgKP;     kaonpdg[1]  = kPdgKP;      kaonpdg[2]  = kPdgK0;
  } else {
    inuclpdg[0] = kPdgProton; inuclpdg[1] = kPdgNeutron; inuclpdg[2] = kPdgProton;
    fnuclpdg[0] = kPdgProton; fnuclpdg[1] = kPdgNeutron; fnuclpdg[2] = kPdgNeutron;
    kaonpdg[0]  = kPdgKM;     kaonpdg[1]  = kPdgKM;      kaonpdg[2]  = kPdgK0;
  }

  for(int i=0; i<nch; i++) {

    ProcessInfo   proc_info(kScSingleKaon, kIntWeakCC);
    Interaction * interaction = new Interaction(init_state, proc_info);

    Target * target  = interaction->InitStatePtr()->TgtPtr();
    bool hasP = (target->Z() > 0);
    bool hasN = (target->N() > 0);

    XclsTag * xcls = interaction->ExclTagPtr();

    if(inuclpdg[i] == kPdgProton  && !hasP) {
      delete interaction;
      continue;
    }
    if(inuclpdg[i] == kPdgNeutron && !hasN) {
      delete interaction;
      continue;
    }
    target->SetHitNucPdg(inuclpdg[i]);
    xcls->SetStrange(kaonpdg[i]);
    if( fnuclpdg[i] == kPdgProton ) xcls->SetNProtons(1);
    else                            xcls->SetNNeutrons(1);

    intlist->push_back(interaction);
  }

  return intlist;
}
//___________________________________________________________________________
void SKInteractionListGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void SKInteractionListGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void SKInteractionListGenerator::LoadConfigData(void)
{
  this->GetParamDef("is-CC", fIsCC, false);
  this->GetParamDef("is-NC", fIsNC, false);
}
//____________________________________________________________________________
