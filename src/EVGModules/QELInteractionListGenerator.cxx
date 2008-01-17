//____________________________________________________________________________
/*
 Copyright (c) 2003-2008, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory - May 13, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "EVGModules/QELInteractionListGenerator.h"
#include "EVGCore/InteractionList.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"

using namespace genie;

//___________________________________________________________________________
QELInteractionListGenerator::QELInteractionListGenerator() :
InteractionListGeneratorI("genie::QELInteractionListGenerator")
{

}
//___________________________________________________________________________
QELInteractionListGenerator::QELInteractionListGenerator(string config):
InteractionListGeneratorI("genie::QELInteractionListGenerator",  config)
{

}
//___________________________________________________________________________
QELInteractionListGenerator::~QELInteractionListGenerator()
{

}
//___________________________________________________________________________
InteractionList * QELInteractionListGenerator::CreateInteractionList(
                                       const InitialState & init_state) const
{
  LOG("InteractionList", pINFO)
                          << "InitialState = " << init_state.AsString();

  if      (fIsCC && !fIsCharm) 
                return this->CreateInteractionListCC(init_state);
  else if (fIsNC && !fIsCharm) 
                return this->CreateInteractionListNC(init_state);
  else if (fIsCC &&  fIsCharm) 
                return this->CreateInteractionListCharmCC(init_state);
  else {
     LOG("InteractionList", pWARN)
       << "Unknown InteractionType! Returning NULL InteractionList "
                         << "for init-state: " << init_state.AsString();
     return 0;
  }
  return 0;
}
//___________________________________________________________________________
InteractionList * QELInteractionListGenerator::CreateInteractionListCC(
                                       const InitialState & init_state) const
{
  InteractionList * intlist = new InteractionList;

  ProcessInfo   proc_info(kScQuasiElastic, kIntWeakCC);
  Interaction * interaction = new Interaction(init_state, proc_info);

  int      nupdg   = init_state.ProbePdg();
  bool     isnu    = pdg::IsNeutrino     (nupdg);
  bool     isnubar = pdg::IsAntiNeutrino (nupdg);

  Target * target  = interaction->InitStatePtr()->TgtPtr();
  bool     hasP    = (target->Z() > 0);
  bool     hasN    = (target->N() > 0);

  if (isnu && hasN) {
     target->SetHitNucPdg(kPdgNeutron);
     intlist->push_back(interaction);

  } else if (isnubar && hasP) {
     target->SetHitNucPdg(kPdgProton);
     intlist->push_back(interaction);

  } else {
     LOG("InteractionList", pWARN)
       << "Returning NULL InteractionList for init-state: "
                                                  << init_state.AsString();
     delete interaction;
     delete intlist;
     return 0;
  }
  return intlist;
}
//___________________________________________________________________________
InteractionList * QELInteractionListGenerator::CreateInteractionListNC(
                                       const InitialState & init_state) const
{
  InteractionList * intlist = new InteractionList;

  int nuclpdg[2] = { kPdgProton, kPdgNeutron };

  int      nupdg   = init_state.ProbePdg();
  bool     isnu    = pdg::IsNeutrino     (nupdg);
  bool     isnubar = pdg::IsAntiNeutrino (nupdg);

  if(!isnu && !isnubar) {
     LOG("InteractionList", pWARN)
       << "Can not handle probe! Returning NULL InteractionList "
                         << "for init-state: " << init_state.AsString();
     delete intlist;
     return 0;
  }

  for(int i=0; i<2; i++) {

     ProcessInfo   proc_info(kScQuasiElastic, kIntWeakNC);
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
     LOG("InteractionList", pERROR)
       << "Returning NULL InteractionList for init-state: "
                                                  << init_state.AsString();
     delete intlist;
     return 0;
  }
  return intlist;
}
//___________________________________________________________________________
InteractionList * 
  QELInteractionListGenerator::CreateInteractionListCharmCC(
                                      const InitialState & init_state) const
{
  //   vl + n --> l- + Lambda_{c}^{+} (2285)
  //   vl + n --> l- + Sigma_{c}^{+}  (2455)
  //   vl + p --> l- + Sigma_{c}^{++} (2455)

  int  nupdg = init_state.ProbePdg();
  bool isnu = pdg::IsNeutrino(nupdg);
  if(!isnu) {
     LOG("InteractionList", pERROR)
       << "Returning NULL InteractionList for init-state: "
                                                  << init_state.AsString();
     return 0;
  }

  const int nch = 3;
  int nuclpdg [nch] = { kPdgNeutron,  kPdgNeutron, kPdgProton   };
  int charmpdg[nch] = { kPdgLambdaPc, kPdgSigmaPc, kPdgSigmaPPc };

  InteractionList * intlist = new InteractionList;

  for(int i=0; i<nch; i++) {

     ProcessInfo   proc_info(kScQuasiElastic, kIntWeakCC);
     Interaction * interaction = new Interaction(init_state, proc_info);

     Target * target  = interaction->InitStatePtr()->TgtPtr();
     bool hasP = (target->Z() > 0);
     bool hasN = (target->N() > 0);

     XclsTag * xcls = interaction->ExclTagPtr();

     if(nuclpdg[i] == kPdgProton  && !hasP) {
       delete interaction;
       continue;
     }
     if(nuclpdg[i] == kPdgNeutron  && !hasN) {
       delete interaction;
       continue;
     }
     target->SetHitNucPdg(nuclpdg[i]);
     xcls->SetCharm(charmpdg[i]);

     intlist->push_back(interaction);
  }
  return intlist;
}
//____________________________________________________________________________
void QELInteractionListGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void QELInteractionListGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void QELInteractionListGenerator::LoadConfigData(void)
{
  fIsCC    = fConfig->GetBoolDef("is-CC",    false);
  fIsNC    = fConfig->GetBoolDef("is-NC",    false);
  fIsCharm = fConfig->GetBoolDef("is-Charm", false);
}
//____________________________________________________________________________

