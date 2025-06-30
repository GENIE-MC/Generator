//____________________________________________________________________________
/*!

\class    genie::HELeptonInteractionListGenerator

\brief    Interaction list generator in HELepton

\author   Alfonso Garcia <aagarciasoto \at km3net.de>
          IFIC & Harvard University

\created  Dec 8, 2021

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _HE_LEPTON_INTERACTION_GENERATOR_H_
#define _HE_LEPTON_INTERACTION_GENERATOR_H_

#include "Framework/EventGen/InteractionListGeneratorI.h"

namespace genie {

class HELeptonInteractionListGenerator : public InteractionListGeneratorI {

public :
  HELeptonInteractionListGenerator();
  HELeptonInteractionListGenerator(string config);
 ~HELeptonInteractionListGenerator();

  // implement the InteractionListGeneratorI interface
  InteractionList * CreateInteractionList(const InitialState & init) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  InteractionList * GLRESInteraction        (const InitialState & init_state) const;
  InteractionList * HENuElectronInteraction (const InitialState & init_state) const;
  InteractionList * PhotonRESInteraction    (const InitialState & init_state) const;
  InteractionList * PhotonCOHInteraction    (const InitialState & init_state) const;

  void LoadConfigData(void);
  
  bool fIsGLRESMu;
  bool fIsGLRESTau;
  bool fIsGLRESEle;
  bool fIsGLRESHad;
  bool fIsHENuElCC;
  bool fIsHENuElNC;
  bool fIsPhotonRESMu;
  bool fIsPhotonRESEle;
  bool fIsPhotonRESTau;
  bool fIsPhotonRESHad;
  bool fIsPhotonCOH;

};

}      // genie namespace

#endif // _HE_LEPTON_INTERACTION_GENERATOR_H_
