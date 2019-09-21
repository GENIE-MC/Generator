//____________________________________________________________________________
/*!

\class    genie::Pythia6Decayer

\brief    Interface to PYTHIA6 particle decayer. \n
          The Pythia6Decayer is a concrete implementation of the Decayer
          interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  June 20, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _PYTHIA6_DECAYER_I_H_
#define _PYTHIA6_DECAYER_I_H_

#ifdef __GENIE_PYTHIA6_ENABLED__
#include <TPythia6.h>
#endif

#include "Physics/Decay/Decayer.h"

namespace genie {

class GHepParticle;
class Pythia6Decayer: protected Decayer {

public:
  Pythia6Decayer();
  Pythia6Decayer(string config);
  virtual ~Pythia6Decayer();

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
