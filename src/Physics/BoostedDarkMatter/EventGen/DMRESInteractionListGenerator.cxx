//____________________________________________________________________________
/*
 Copyright (c) 2003-2023, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Zachary W. Orr <zachorrt@colostate.edu>
 Colorado State University

 based on code by
 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory
*/
//____________________________________________________________________________

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/EventGen/InteractionList.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Physics/BoostedDarkMatter/EventGen/DMRESInteractionListGenerator.h"

using namespace genie;

//___________________________________________________________________________
DMRESInteractionListGenerator::DMRESInteractionListGenerator() :
InteractionListGeneratorI("genie::DMRESInteractionListGenerator")
{

}
//___________________________________________________________________________
DMRESInteractionListGenerator::DMRESInteractionListGenerator(string config) :
InteractionListGeneratorI("genie::DMRESInteractionListGenerator", config)
{

}
//___________________________________________________________________________
DMRESInteractionListGenerator::~DMRESInteractionListGenerator()
{

}
//___________________________________________________________________________
InteractionList * DMRESInteractionListGenerator::CreateInteractionList(
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
  // Only accept DM scattering here
  if      (fIsDM) inttype = kIntDarkMatter;
  else {
     LOG("IntLst", pWARN)
       << "Unknown InteractionType! Returning NULL InteractionList "
       << "for init-state: " << init_state.AsString();
     return 0;
  }

  // create a process information object
  ProcessInfo proc_info(kScDarkMatterResonant, inttype);

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
//       bool skip_res =  proc_info.IsWeakCC() &&
//                        pdg::IsNeutrino(init_state.ProbePdg()) &&
//                        (hit_nucleon[i]==kPdgProton) &&
//                        (!utils::res::IsDelta(res));
//       if(skip_res) continue;

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
void DMRESInteractionListGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void DMRESInteractionListGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void DMRESInteractionListGenerator::LoadConfigData(void)
{
  string resonances = "";
  this->GetParam("ResonanceNameList", resonances);
  SLOG("IntLst", pDEBUG) << "Resonance list: " << resonances;

  fResList.Clear();
  fResList.DecodeFromNameList(resonances);
  LOG("IntLst", pDEBUG) << fResList;

  this->GetParamDef("is-DM", fIsDM, false);
}
//____________________________________________________________________________
