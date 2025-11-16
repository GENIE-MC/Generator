//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Alfonso Garcia <aagarciasoto \at km3net.de>
 IFIC & Harvard University
*/
//____________________________________________________________________________

#include "Framework/EventGen/InteractionList.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Physics/HELepton/EventGen/HELeptonInteractionListGenerator.h"

using namespace genie;

//___________________________________________________________________________
HELeptonInteractionListGenerator::HELeptonInteractionListGenerator() :
InteractionListGeneratorI("genie::HELeptonInteractionListGenerator")
{

}
//___________________________________________________________________________
HELeptonInteractionListGenerator::HELeptonInteractionListGenerator(string config) :
InteractionListGeneratorI("genie::HELeptonInteractionListGenerator", config)
{

}
//___________________________________________________________________________
HELeptonInteractionListGenerator::~HELeptonInteractionListGenerator()
{

}
//___________________________________________________________________________
InteractionList *
   HELeptonInteractionListGenerator::GLRESInteraction(
                                       const InitialState & init_state) const
{


  InteractionList * intlist = new InteractionList;

  ProcessInfo proc_info(kScGlashowResonance, kIntWeakCC);

  int probepdg = init_state.ProbePdg();

  if(probepdg == kPdgAntiNuE) {
    InitialState init(init_state);
    init_state.TgtPtr()->SetHitNucPdg(0);  
    Interaction * interaction = new Interaction(init_state, proc_info);
    XclsTag exclusive_tag;
    if      (fIsGLRESMu)  exclusive_tag.SetFinalLepton(kPdgMuon);
    else if (fIsGLRESTau) exclusive_tag.SetFinalLepton(kPdgTau);
    else if (fIsGLRESEle) exclusive_tag.SetFinalLepton(kPdgElectron);
    else if (fIsGLRESHad) exclusive_tag.SetFinalLepton(kPdgPiP);
    interaction->SetExclTag(exclusive_tag);
    intlist->push_back(interaction);
  }

  return intlist;

}
//___________________________________________________________________________
InteractionList *
   HELeptonInteractionListGenerator::HENuElectronInteraction(
                                       const InitialState & init_state) const
{


  InteractionList * intlist = new InteractionList;

  int probepdg = init_state.ProbePdg();

  if (fIsHENuElCC) {
    ProcessInfo proc_info(kScGlashowResonance, kIntWeakCC);
    InitialState init(init_state);
    init_state.TgtPtr()->SetHitNucPdg(0);  
    Interaction * interaction = new Interaction(init_state, proc_info);
    XclsTag exclusive_tag; //charged lepton
    if      ( pdg::IsNuMu(probepdg)      ) exclusive_tag.SetFinalLepton(kPdgMuon);
    else if ( pdg::IsNuTau(probepdg)     ) exclusive_tag.SetFinalLepton(kPdgTau);
    else if ( pdg::IsNuE(probepdg)       ) exclusive_tag.SetFinalLepton(kPdgElectron);
    else return intlist;
    interaction->SetExclTag(exclusive_tag);
    intlist->push_back(interaction);
  }
  else if (fIsHENuElNC) {
    ProcessInfo proc_info(kScGlashowResonance, kIntWeakNC);
    InitialState init(init_state);
    init_state.TgtPtr()->SetHitNucPdg(0);  
    Interaction * interaction = new Interaction(init_state, proc_info);
    XclsTag exclusive_tag; //charged lepton
    if      ( pdg::IsNuMu(probepdg)      ) exclusive_tag.SetFinalLepton(kPdgElectron);
    else if ( pdg::IsNuTau(probepdg)     ) exclusive_tag.SetFinalLepton(kPdgElectron);
    else if ( pdg::IsAntiNuMu(probepdg)  ) exclusive_tag.SetFinalLepton(kPdgElectron);
    else if ( pdg::IsAntiNuTau(probepdg) ) exclusive_tag.SetFinalLepton(kPdgElectron);
    else return intlist;
    interaction->SetExclTag(exclusive_tag);
    intlist->push_back(interaction);
  }

  return intlist;

}
//___________________________________________________________________________
InteractionList *
   HELeptonInteractionListGenerator::PhotonRESInteraction(
                                       const InitialState & init_state) const
{

  InteractionList * intlist = new InteractionList;

  ProcessInfo proc_info(kScPhotonResonance, kIntWeakCC);

  int probepdg = init_state.ProbePdg();
  bool hasP = (init_state.Tgt().Z() > 0);
  bool hasN = (init_state.Tgt().N() > 0);

  int nuclpdg[2] = { kPdgProton, kPdgNeutron };
  for(int inucl=0; inucl<2; inucl++) {
    int struck_nucleon = nuclpdg[inucl];
    if( (struck_nucleon == kPdgProton  && hasP) || (struck_nucleon == kPdgNeutron && hasN) ) {
      Interaction * interaction = new Interaction(init_state, proc_info);
      Target * target = interaction->InitStatePtr()->TgtPtr();
      target->SetHitNucPdg(struck_nucleon);
      XclsTag exclusive_tag;
      if      (fIsPhotonRESMu)  exclusive_tag.SetFinalLepton( (probepdg>0) ? kPdgAntiMuon : kPdgMuon     );
      else if (fIsPhotonRESTau) exclusive_tag.SetFinalLepton( (probepdg>0) ? kPdgAntiTau  : kPdgTau      );
      else if (fIsPhotonRESEle) exclusive_tag.SetFinalLepton( (probepdg>0) ? kPdgPositron : kPdgElectron );
      else if (fIsPhotonRESHad) exclusive_tag.SetFinalLepton( (probepdg>0) ? kPdgPiP      : kPdgPiM      );
      interaction->SetExclTag(exclusive_tag);
      intlist->push_back(interaction);
    }
  }

  return intlist;

}
//___________________________________________________________________________
InteractionList *
   HELeptonInteractionListGenerator::PhotonCOHInteraction(
                                       const InitialState & init_state) const
{

  InteractionList * intlist = new InteractionList;
  ProcessInfo   proc_info(kScPhotonCoherent, kIntWeakCC);
  Interaction * interaction = new Interaction(init_state, proc_info);
  intlist->push_back(interaction);
  return intlist;

}
//___________________________________________________________________________
InteractionList *
   HELeptonInteractionListGenerator::CreateInteractionList(
                                       const InitialState & init_state) const
{
// channels:
// nuebar   + e-   -> W- -> nuebar + e-   [CC+NC]
// nuebar   + e-   -> W- -> nuebar + mu-  [CC]
// nuebar   + e-   -> W- -> nuebar + tau- [CC]
// nuebar   + e-   -> W- -> hadrons       [CC]
// nue      + e-   -> e     + nue         [CC+NC]
// numu     + e-   -> mu    + nue         [CC]
// nutau    + e-   -> tau   + nue         [CC]
// numu     + e-   -> numu  + e           [NC]
// nutau    + e-   -> nutau + e           [NC]
// numubar  + e-   -> numubar  + e        [NC]
// nutaubar + e-   -> nutaubar + e        [NC]
// nu     + gamma* -> l- + W+ (coherent & resonant)
// nubar  + gamma* -> l+ + W- (coherent & resonant)

  int ppdg = init_state.ProbePdg();
  if( !pdg::IsNeutralLepton(ppdg) ) {
     LOG("IntLst", pWARN)
       << "Can not handle probe! Returning NULL InteractionList "
       << "for init-state: " << init_state.AsString();
     return 0;
  }

  if      (fIsGLRESMu)      return GLRESInteraction(init_state);
  else if (fIsGLRESTau)     return GLRESInteraction(init_state);
  else if (fIsGLRESEle)     return GLRESInteraction(init_state);
  else if (fIsGLRESHad)     return GLRESInteraction(init_state);
  else if (fIsHENuElCC)     return HENuElectronInteraction(init_state);
  else if (fIsHENuElNC)     return HENuElectronInteraction(init_state);
  else if (fIsPhotonRESMu)  return PhotonRESInteraction(init_state);
  else if (fIsPhotonRESTau) return PhotonRESInteraction(init_state);
  else if (fIsPhotonRESEle) return PhotonRESInteraction(init_state);
  else if (fIsPhotonRESHad) return PhotonRESInteraction(init_state);
  else if (fIsPhotonCOH)    return PhotonCOHInteraction(init_state);
  else {
     LOG("IntLst", pERROR)
         << "Returning NULL InteractionList for init-state: " << init_state.AsString();
     return 0;
  }

}
//___________________________________________________________________________
void HELeptonInteractionListGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void HELeptonInteractionListGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void HELeptonInteractionListGenerator::LoadConfigData(void)
{

  GetParamDef("is-GLRES-Mu",      fIsGLRESMu,      false ) ;
  GetParamDef("is-GLRES-Tau",     fIsGLRESTau,     false ) ;
  GetParamDef("is-GLRES-Ele",     fIsGLRESEle,     false ) ;
  GetParamDef("is-GLRES-Had",     fIsGLRESHad,     false ) ;
  GetParamDef("is-HENuEl-CC",     fIsHENuElCC,     false ) ;
  GetParamDef("is-HENuEl-NC",     fIsHENuElNC,     false ) ;
  GetParamDef("is-PhotonRES-Mu",  fIsPhotonRESMu,  false ) ;
  GetParamDef("is-PhotonRES-Tau", fIsPhotonRESTau, false ) ;
  GetParamDef("is-PhotonRES-Ele", fIsPhotonRESEle, false ) ;
  GetParamDef("is-PhotonRES-Had", fIsPhotonRESHad, false ) ;
  GetParamDef("is-PhotonCOH",     fIsPhotonCOH,    false ) ;

}