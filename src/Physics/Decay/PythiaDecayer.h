//____________________________________________________________________________
/*!

\class    genie::PythiaDecayer

\brief    Interface to PYTHIA particle decayer. \n
          The PythiaDecayer is a concrete implementation of the Decayer
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

#include <TPythia6.h>

#include "Physics/Decay/Decayer.h"

namespace genie {

class GHepParticle;
class PythiaDecayer: protected Decayer {

public:
  PythiaDecayer();
  PythiaDecayer(string config);
  virtual ~PythiaDecayer();

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

  mutable TPythia6 * fPythia;  ///< PYTHIA6 wrapper class
  mutable double     fWeight;
};

}         // genie namespace
#endif    // _PYTHIA6_DECAYER_H_
