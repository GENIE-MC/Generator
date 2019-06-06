//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Mar 03, 2009 - CA
   Moved into the new QEL package from its previous location (EVGModules)
 @ Sep 15, 2009 - CA
   Generate interaction lists for charged lepton scattering too.

*/
//____________________________________________________________________________

#include "Framework/EventGen/InteractionList.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Physics/QuasiElastic/EventGen/QELInteractionListGenerator.h"

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
  LOG("IntLst", pINFO)
     << "InitialState = " << init_state.AsString();

  if (fIsCC && !fIsCharm && !fIsStrange) 
     return this->CreateInteractionListCC(init_state);
  else 
  if (fIsNC && !fIsCharm && !fIsStrange) 
     return this->CreateInteractionListNC(init_state);
  else 
  if (fIsEM) 
     return this->CreateInteractionListEM(init_state);
  else 
  if (fIsCC &&  fIsCharm) 
     return this->CreateInteractionListCharmCC(init_state);
  else 
  if (fIsCC &&  fIsStrange) 
     return this->CreateInteractionListStrangeCC(init_state);
  else {
     LOG("IntLst", pWARN)
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

  if(!isnu && !isnubar) {
     LOG("IntLst", pWARN)
       << "Can not handle probe! Returning NULL InteractionList "
       << "for init-state: " << init_state.AsString();
     delete intlist;
     return 0;
  }

  if (isnu && hasN) {
     target->SetHitNucPdg(kPdgNeutron);
     intlist->push_back(interaction);

  } else if (isnubar && hasP) {
     target->SetHitNucPdg(kPdgProton);
     intlist->push_back(interaction);

  } else {
     LOG("IntLst", pINFO)
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
     LOG("IntLst", pWARN)
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
     LOG("IntLst", pINFO)
       << "Returning NULL InteractionList for init-state: "
       << init_state.AsString();
     delete intlist;
     return 0;
  }
  return intlist;
}
//___________________________________________________________________________
InteractionList * QELInteractionListGenerator::CreateInteractionListEM(
   const InitialState & init_state) const
{
  InteractionList * intlist = new InteractionList;

  int tgtpdg = init_state.Tgt().Pdg();
  int ppdg   = init_state.ProbePdg();

  bool ischgl = pdg::IsChargedLepton(ppdg);
  if(!ischgl) {
     LOG("IntLst", pWARN)
       << "Can not handle probe! Returning NULL InteractionList "
       << "for init-state: " << init_state.AsString();
     delete intlist;
     return 0;
  }

  bool hasP = (init_state.Tgt().Z() > 0);
  if(hasP) {
    Interaction * interaction = Interaction::QELEM(tgtpdg,kPdgProton,ppdg);
    intlist->push_back(interaction);
  }
  bool hasN = (init_state.Tgt().N() > 0);
  if(hasN) {
    Interaction * interaction = Interaction::QELEM(tgtpdg,kPdgNeutron,ppdg);
    intlist->push_back(interaction);
  }

  if(intlist->size() == 0) {
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
     LOG("IntLst", pINFO)
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
//___________________________________________________________________________
InteractionList * 
  QELInteractionListGenerator::CreateInteractionListStrangeCC(
    const InitialState & init_state) const
{
  // v + n --> mu+ + Sigma^{-} 
  // v + p --> mu+ + Lambda^{0}
  // v + p --> mu+ + Sigma^{0}

  int  nupdg = init_state.ProbePdg();
  bool isnu = pdg::IsNeutrino(nupdg);
  if(isnu) {
     LOG("IntLst", pERROR)
       << "Returning NULL InteractionList for init-state: "
       << init_state.AsString();
     return 0;
  }

  const int nch = 3;
  int nuclpdg [nch]   = { kPdgNeutron,  kPdgProton, kPdgProton   };
  int strangepdg[nch] = { kPdgSigmaM ,  kPdgLambda, kPdgSigma0   };

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
     xcls->SetStrange(strangepdg[i]);

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
  GetParamDef( "is-CC", fIsCC, false ) ;
  GetParamDef( "is-NC", fIsNC, false ) ;
  GetParamDef( "is-EM", fIsEM, false ) ;
  GetParamDef( "is-Charm", fIsCharm, false  ) ;
  GetParamDef( "is-Strange", fIsStrange, false  ) ;
}
//____________________________________________________________________________

