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
   Moved into the NuE package from its previous location (EVGModules package)
 @ Feb 12, 2013 - CA (code from Rosen Matev)
   Handle the IMD annihilation channel.

*/
//____________________________________________________________________________

#include "Physics/NuElectron/EventGen/NuEInteractionListGenerator.h"
#include "Framework/EventGen/InteractionList.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"

using namespace genie;

//___________________________________________________________________________
NuEInteractionListGenerator::NuEInteractionListGenerator() :
InteractionListGeneratorI("genie::NuEInteractionListGenerator")
{

}
//___________________________________________________________________________
NuEInteractionListGenerator::NuEInteractionListGenerator(string config) :
InteractionListGeneratorI("genie::NuEInteractionListGenerator", config)
{

}
//___________________________________________________________________________
NuEInteractionListGenerator::~NuEInteractionListGenerator()
{

}
//___________________________________________________________________________
InteractionList * NuEInteractionListGenerator::CreateInteractionList(
                                       const InitialState & init_state) const
{
  LOG("IntLst", pINFO) << "InitialState = " << init_state.AsString();

  if(fIsIMD)  return this -> IMDInteractionList   (init_state);
  else if(fIsIMDAnh)  return this -> IMDAnnihilationInteractionList (init_state);
  else        return this -> NuEELInteractionList (init_state);
}
//___________________________________________________________________________
InteractionList * NuEInteractionListGenerator::IMDInteractionList(
                                       const InitialState & init_state) const
{
// channels:
// numu + e- -> mu- + nu_e [CC] -- 'inverse muon decay'

  if(init_state.ProbePdg() != kPdgNuMu) {
     LOG("IntLst", pDEBUG) 
          << "Return *null* interaction list (non nu_mu probe in IMD thread)";
     return 0;
  }

  InteractionList * intlist = new InteractionList;

  // clone init state and de-activate the struck nucleon info
  InitialState init(init_state);
  init_state.TgtPtr()->SetHitNucPdg(0);

  ProcessInfo   proc_info(kScInverseMuDecay, kIntWeakCC);
  Interaction * interaction = new Interaction(init, proc_info);

  intlist->push_back(interaction);

  return intlist;
}
//___________________________________________________________________________
InteractionList * NuEInteractionListGenerator::IMDAnnihilationInteractionList(
                                       const InitialState & init_state) const
{
// channels:
// nuebar + e- -> mu- + nu_e [CC] -- 'inverse muon decay annihilation channel'

  if(init_state.ProbePdg() != kPdgAntiNuE) {
     LOG("IntLst", pDEBUG) 
          << "Return *null* interaction list (non anti_nu_e probe in IMDAnnihilation thread)";
     return 0;
  }

  InteractionList * intlist = new InteractionList;

  // clone init state and de-activate the struck nucleon info
  InitialState init(init_state);
  init_state.TgtPtr()->SetHitNucPdg(0);

  ProcessInfo   proc_info(kScIMDAnnihilation, kIntWeakCC);
  Interaction * interaction = new Interaction(init, proc_info);

  intlist->push_back(interaction);

  return intlist;
}
//___________________________________________________________________________
InteractionList * NuEInteractionListGenerator::NuEELInteractionList(
                                       const InitialState & init_state) const
{
// channels:
// nue      + e- -> nue      + e-   [CC + NC + interference]
// nuebar   + e- -> nuebar   + e-   [CC + NC + interference]
// numu     + e- -> numu     + e-   [NC]
// numu     + e- -> mu-      + nu_e [CC] -- handled by the IMD thread
// nutau    + e- -> nutau    + e-   [NC]
// nutau    + e- -> tau    - + nu_e [CC] -- neglected
// numubar  + e- -> numubar  + e-   [NC]
// nutaubar + e- -> nutaubar + e-   [NC]

  int nupdg = init_state.ProbePdg();
  InteractionList * intlist = new InteractionList;

  // clone init state and de-activate the struck nucleon info
  InitialState init(init_state);
  init_state.TgtPtr()->SetHitNucPdg(0);

  // NC
  if(nupdg == kPdgNuMu  || nupdg == kPdgAntiNuMu || 
     nupdg == kPdgNuTau || nupdg == kPdgAntiNuTau) {
     ProcessInfo   proc_info(kScNuElectronElastic,  kIntWeakNC);
     Interaction * interaction = new Interaction(init, proc_info);
     intlist->push_back(interaction);
  }

  // CC+NC+interference
  if(nupdg == kPdgNuE  || nupdg == kPdgAntiNuE) { 
     ProcessInfo   proc_info(kScNuElectronElastic,  kIntWeakMix); 
     Interaction * interaction = new Interaction(init, proc_info);
     intlist->push_back(interaction);
  }

  return intlist;
}
//___________________________________________________________________________
void NuEInteractionListGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void NuEInteractionListGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void NuEInteractionListGenerator::LoadConfig(void)
{
	GetParamDef( "is-IMD", fIsIMD, false ) ;
	GetParamDef( "is-IMD-ANH", fIsIMDAnh, false ) ;
}
//____________________________________________________________________________

