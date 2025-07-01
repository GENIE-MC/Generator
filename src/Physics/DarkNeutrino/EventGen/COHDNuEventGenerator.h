//____________________________________________________________________________
/*!
\class    genie::COHDNuEventGenerator

\brief    Generates complete COHDNu events.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Iker de Icaza <i.de-icaza-astiz \at sussex.ac.uk>
          University of Sussex

          Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  June 30, 2020

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _COHDNu_EVENT_GENERATOR_H_
#define _COHDNu_EVENT_GENERATOR_H_

#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Framework/Utils/Range1.h"

namespace genie {

class XSecAlgorithmI;

class COHDNuEventGenerator : public EventRecordVisitorI {

public :
  COHDNuEventGenerator();
  COHDNuEventGenerator(string config);
 ~COHDNuEventGenerator();

  // Implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event) const;

  // Overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  // Methods to load sub-algorithms and config data from the Registry
  void LoadConfig (void);

  // Event generation methods
  void GenerateKinematics        (GHepRecord * event) const;
  void AddFinalStateDarkNeutrino (GHepRecord * event) const;
  void AddRecoilNucleus          (GHepRecord * event) const;

  mutable const XSecAlgorithmI * fXSecModel; ///<

  bool   fGenerateUniformly;    ///<
  double fSafetyFactor;         ///<
  double fMaxXSecDiffTolerance; ///<

  double fDNuMass, fDNuMass2;

};

} // genie namespace

#endif // _COHDNu_EVENT_GENERATOR_H_
