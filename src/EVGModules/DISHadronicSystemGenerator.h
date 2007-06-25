//____________________________________________________________________________
/*!

\class    genie::DISHadronicSystemGenerator

\brief    Generates the final state hadronic system in v DIS interactions.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  October 03, 2004

\cpright  Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _DIS_HADRONIC_SYSTEM_GENERATOR_H_
#define _DIS_HADRONIC_SYSTEM_GENERATOR_H_

#include "EVGModules/HadronicSystemGenerator.h"

namespace genie {

class HadronizationModelI;

class DISHadronicSystemGenerator : public HadronicSystemGenerator {

public :

  DISHadronicSystemGenerator();
  DISHadronicSystemGenerator(string config);
  ~DISHadronicSystemGenerator();

  //-- implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void AddFragmentationProducts (GHepRecord * event_rec) const;
  void SimulateFormationZone    (GHepRecord * event_rec) const;

  void LoadConfig (void);

  const HadronizationModelI * fHadronizationModel;

  bool   fFilterPreFragmEntries;
  double fct0;        ///< formation zone (c * formation time)
  double fK;          ///< param multiplying pT^2 in formation zone calculation
};

}      // genie namespace

#endif // _DIS_HADRONIC_SYSTEM_GENERATOR_H_
