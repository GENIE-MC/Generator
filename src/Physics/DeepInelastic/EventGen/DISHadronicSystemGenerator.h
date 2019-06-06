//____________________________________________________________________________
/*!

\class    genie::DISHadronicSystemGenerator

\brief    Generates the final state hadronic system in v DIS interactions.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  October 03, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _DIS_HADRONIC_SYSTEM_GENERATOR_H_
#define _DIS_HADRONIC_SYSTEM_GENERATOR_H_

#include "Physics/Common/HadronicSystemGenerator.h"

namespace genie {

class HadronizationModelI;

class DISHadronicSystemGenerator : public HadronicSystemGenerator {

public :
  DISHadronicSystemGenerator();
  DISHadronicSystemGenerator(string config);
 ~DISHadronicSystemGenerator();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void AddFragmentationProducts (GHepRecord * event_rec) const;
  void SimulateFormationZone    (GHepRecord * event_rec) const;

  void LoadConfig (void);

  const HadronizationModelI * fHadronizationModel;

  bool   fFilterPreFragmEntries;
  double fR0;          ///< param controling nuclear size
  double fNR;          ///< how far beyond the nuclear boundary does the particle tracker goes?
  double fct0pion;     ///< formation zone (c * formation time) - for pions
  double fct0nucleon;  ///< formation zone (c * formation time) - for nucleons
  double fK;           ///< param multiplying pT^2 in formation zone calculation
};

}      // genie namespace

#endif // _DIS_HADRONIC_SYSTEM_GENERATOR_H_
