//____________________________________________________________________________
/*!

\class    genie::CEvNSEventGenerator

\brief    Generates complete CEvNS events.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
          University of Liverpool & STFC Rutherford Appleton Laboratory

\created  July 16, 2019

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _CEVNS_EVENT_GENERATOR_H_
#define _CEVNS_EVENT_GENERATOR_H_

#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Framework/Utils/Range1.h"

namespace genie {

class XSecAlgorithmI;

class CEvNSEventGenerator : public EventRecordVisitorI {

public :
  CEvNSEventGenerator();
  CEvNSEventGenerator(string config);
 ~CEvNSEventGenerator();

  // Implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event) const;

  // Overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

public:

  // Methods to load sub-algorithms and config data from the Registry
  void LoadConfig (void);

  // Event generation methods
  void GenerateKinematics    (GHepRecord * event) const;
  void AddFinalStateNeutrino (GHepRecord * event) const;
  void AddRecoilNucleus      (GHepRecord * event) const;

  mutable const XSecAlgorithmI * fXSecModel; ///<

  bool   fGenerateUniformly;    ///<
  double fSafetyFactor;         ///<
  double fMaxXSecDiffTolerance; ///<
};

}      // genie namespace

#endif // _CEVNS_EVENT_GENERATOR_H_
