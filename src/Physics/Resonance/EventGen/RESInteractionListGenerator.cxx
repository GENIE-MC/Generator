//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - May 13, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Mar 03, 2009 - CA
   Moved into the new RES package from its previous location (EVGModules)
 @ Sep 21, 2009 - CA
   Generate interaction lists for charge lepton scattering

*/
//____________________________________________________________________________

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/EventGen/InteractionList.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Physics/Resonance/EventGen/RESInteractionListGenerator.h"

using namespace genie;

//___________________________________________________________________________
RESInteractionListGenerator::RESInteractionListGenerator() :
InteractionListGeneratorI("genie::RESInteractionListGenerator")
{

}
//___________________________________________________________________________
RESInteractionListGenerator::RESInteractionListGenerator(string config) :
InteractionListGeneratorI("genie::RESInteractionListGenerator", config)
{

}
//___________________________________________________________________________
RESInteractionListGenerator::~RESInteractionListGenerator()
{

}
//___________________________________________________________________________
InteractionList * RESInteractionListGenerator::CreateInteractionList(
                                       const InitialState & init_state) const
{
  LOG("IntLst", pINFO) << "InitialState = " << init_state.AsString();

  // In the thread generating interactions from the list produced here (RES),
  // we simulate (for free and nuclear targets) semi-inclusive resonance
  // interactions: v + N -> v(l) + R -> v(l) + X
  // Specifically, the RES thread generates:
  //
  //  CC:
  //    nu       + p (A) -> l-       R (A), for all resonances with Q=+2
  //    nu       + n (A) -> l-       R (A), for all resonances with Q=+1
  //    \bar{nu} + p (A) -> l+       R (A), for all resonances with Q= 0
  //    \bar{nu} + n (A) -> l+       R (A), for all resonances with Q=-1
  //  NC:
  //    nu       + p (A) -> nu       R (A), for all resonances with Q=+1
  //    nu       + n (A) -> nu       R (A), for all resonances with Q= 0
  //    \bar{nu} + p (A) -> \bar{nu} R (A), for all resonances with Q=+1
  //    \bar{nu} + n (A) -> \bar{nu} R (A), for all resonances with Q= 0
  //
  // and then the resonance R should be allowed to decay to get the full
  // hadronic final state X. All decay channels kinematically accessible
  // to the (off the mass-shell produced) resonance can be allowed.

  // specify the requested interaction type
  InteractionType_t inttype;
  if      (fIsCC) inttype = kIntWeakCC;
  else if (fIsNC) inttype = kIntWeakNC;
  else if (fIsEM) inttype = kIntEM;
  else {
     LOG("IntLst", pWARN)
       << "Unknown InteractionType! Returning NULL InteractionList "
       << "for init-state: " << init_state.AsString();
     return 0;
  }

  // create a process information object
  ProcessInfo proc_info(kScResonant, inttype);

  // learn whether the input nuclear or free target has avail. p and n
  const Target & inp_target = init_state.Tgt();
  bool hasP = (inp_target.Z() > 0);
  bool hasN = (inp_target.N() > 0);

  // possible hit nucleons
  const int hit_nucleon[2] = {kPdgProton, kPdgNeutron};

  // create an interaction list
  InteractionList * intlist = new InteractionList;

  // loop over all baryon resonances considered in current MC job
  unsigned int nres = fResList.NResonances();
  for(unsigned int ires = 0; ires < nres; ires++) {

     //get current resonance
     Resonance_t res = fResList.ResonanceId(ires);

     // loop over hit nucleons
     for(int i=0; i<2; i++) {

       // proceed only if the hit nucleon exists in the current init state
       if(hit_nucleon[i]==kPdgProton  && !hasP) continue;
       if(hit_nucleon[i]==kPdgNeutron && !hasN) continue;

       // proceed only if the current resonance conserves charge
       // (the only problematic case is when the RES charge has to be +2
       //  because then only Delta resonances are possible)
       bool skip_res =  proc_info.IsWeakCC() &&
                        pdg::IsNeutrino(init_state.ProbePdg()) &&
                        (hit_nucleon[i]==kPdgProton) &&
                        (!utils::res::IsDelta(res));
       if(skip_res) continue;

       // create an interaction
       Interaction * interaction = new Interaction(init_state, proc_info);

       // add the struck nucleon
       Target * target = interaction->InitStatePtr()->TgtPtr();
       target->SetHitNucPdg(hit_nucleon[i]);

       // add the baryon resonance in the exclusive tag
       XclsTag * xcls = interaction->ExclTagPtr();
       xcls->SetResonance(res);

       // add the interaction at the interaction list
       intlist->push_back(interaction);

     }//hit nucleons
  } //resonances

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
void RESInteractionListGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void RESInteractionListGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void RESInteractionListGenerator::LoadConfigData(void)
{
  string resonances = "";
  this->GetParam("ResonanceNameList", resonances);
  SLOG("IntLst", pDEBUG) << "Resonance list: " << resonances;

  fResList.Clear();
  fResList.DecodeFromNameList(resonances);
  LOG("IntLst", pDEBUG) << fResList;

  this->GetParamDef("is-CC", fIsCC, false);
  this->GetParamDef("is-NC", fIsNC, false);
  this->GetParamDef("is-EM", fIsEM, false);
}
//____________________________________________________________________________
