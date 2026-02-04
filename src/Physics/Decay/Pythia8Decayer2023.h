//____________________________________________________________________________
/*!

\class    genie::Pythia8Decayer2023

\brief    Interface to PYTHIA particle decayer. \n
          The Pythia8Decayer2023 is a concrete implementation of the Decayer
          interface.

\author   Robert Hatcher <rhatcher \at fnal.gov>
          Fermilab

\created  December 21, 2023

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _PYTHIA8DECAYER2023_H_
#define _PYTHIA8DECAYER2023_H_

#ifdef __GENIE_PYTHIA8_ENABLED__
#include "Framework/Utils/Pythia8Singleton.h"
#endif

#include "Physics/Decay/Decayer.h"

namespace genie {

class GHepParticle;
class Pythia8Decayer2023: protected Decayer {

public:
  Pythia8Decayer2023();
  Pythia8Decayer2023(string config);
  virtual ~Pythia8Decayer2023();

  // Implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event) const;

private:

  void   Initialize             (void)                           const;
  bool   IsHandled              (int pdgc)                       const;
  void   InhibitDecay           (int pdgc, TDecayChannel * ch=0) const;
  void   UnInhibitDecay         (int pdgc, TDecayChannel * ch=0) const;
  bool   Decay                  (int ip, GHepRecord * event)     const;
  double SumOfBranchingRatios   (int kc)                         const;
  int    FindPythiaDecayChannel (int kc, TDecayChannel * ch)     const;
  bool   MatchDecayChannels     (int ic, TDecayChannel * ch)     const;

  mutable double     fWeight;
};

}         // genie namespace
#endif    // _PYTHIA8DECAYER2023_H_
