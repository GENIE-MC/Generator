//____________________________________________________________________________
/*!

\class    genie::Pythia6Decayer2023

\brief    Interface to PYTHIA particle decayer. \n
          The Pythia6Decayer2023 is a concrete implementation of the Decayer
          interface.

\author   Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  June 20, 2004

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _PYTHIA6_DECAYER_I_H_
#define _PYTHIA6_DECAYER_I_H_

#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_PYTHIA6_ENABLED__
#include <TPythia6.h>
#endif

#include "Physics/Decay/Decayer.h"

namespace genie {

class GHepParticle;
class Pythia6Decayer2023: protected Decayer {

public:
  Pythia6Decayer2023();
  Pythia6Decayer2023(string config);
  virtual ~Pythia6Decayer2023();

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

#ifdef __GENIE_PYTHIA6_ENABLED__
  mutable TPythia6 * fPythia;  ///< PYTHIA6 wrapper class
#endif
  mutable double     fWeight;
};

}         // genie namespace
#endif    // _PYTHIA6_DECAYER_H_
