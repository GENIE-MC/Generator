//____________________________________________________________________________
/*!

\class    genie::HEDISGenerator

\brief    Generates the final state leptonic and hadronic system in v HEDIS 
          interactions.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Alfonso Garcia <alfonsog \at nikhef.nl>
          NIKHEF

\created  August 28, 2019

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _HEDIS_GENERATOR_H_
#define _HEDIS_GENERATOR_H_

#include "Physics/Common/HadronicSystemGenerator.h"


namespace genie {

class HEDISGenerator : public HadronicSystemGenerator {

public :
  HEDISGenerator();
  HEDISGenerator(string config);
 ~HEDISGenerator();

  // implement the EventRecordVisitorI interface
  void Initialize        (void)               const;
  void ProcessEventRecord(GHepRecord * evrec) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void AddPrimaryLepton         (GHepRecord * evrec) const;

  void LoadConfig (void);

  const EventRecordVisitorI * fHadronizationModel;

};

}      // genie namespace

#endif // _HEDIS_HADRONIC_SYSTEM_GENERATOR_H_
