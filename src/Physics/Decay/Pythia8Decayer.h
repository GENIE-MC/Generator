//____________________________________________________________________________
/*!

\class    genie::Pythia8Decayer

\brief    Interface to PYTHIA8 particle decayer. \n
          The Pythia8Decayer is a concrete implementation of the Decayer
          interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

          Shivesh Mandalia <s.p.mandalia@qmul.ac.uk>
          Queen Mary University of London

\created  June 20, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _PYTHIA8_DECAYER_I_H_
#define _PYTHIA8_DECAYER_I_H_

#include "Physics/Decay/Decayer.h"
#include "Physics/Hadronization/Pythia8Singleton.h"

namespace genie {

class GHepParticle;
class Pythia8Decayer: protected Decayer {

public:
  Pythia8Decayer();
  Pythia8Decayer(string config);
  virtual ~Pythia8Decayer();

  // Implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event) const;

private:

  void   Initialize             (void)                                 const;
  bool   IsHandled              (int pdgc)                             const;
  void   InhibitDecay           (int pdgc, TDecayChannel * ch=0)       const;
  void   UnInhibitDecay         (int pdgc, TDecayChannel * ch=0)       const;
  bool   Decay                  (int ip, GHepRecord * event)           const;
  double SumOfBranchingRatios   (int pdgc)                             const;
  int    FindPythiaDecayChannel (int pdgc, TDecayChannel * ch)         const;
  bool   MatchDecayChannels     (int pdgc, int ic, TDecayChannel * ch) const;

  mutable Pythia8Singleton * fPythia;  ///< PYTHIA8 wrapper class
  mutable double     fWeight;
};

}         // genie namespace
#endif    // _PYTHIA8_DECAYER_H_
