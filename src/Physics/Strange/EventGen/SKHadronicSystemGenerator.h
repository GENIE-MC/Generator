//____________________________________________________________________________
/*!

\class    genie::SKHadronicSystemGenerator

\brief    Generates the f/s hadronic system in single-Kaon production interactions.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  March 20, 2014

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _SK_HADRONIC_SYSTEM_GENERATOR_H_
#define _SK_HADRONIC_SYSTEM_GENERATOR_H_

#include "Physics/Common/HadronicSystemGenerator.h"

namespace genie {

class SKHadronicSystemGenerator : public HadronicSystemGenerator {

public :
  SKHadronicSystemGenerator();
  SKHadronicSystemGenerator(string config);
 ~SKHadronicSystemGenerator();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;
  void CalculateHadronicSystem_AtharSingleKaon(GHepRecord * event_rec) const;
};

}      // genie namespace
#endif // _SK_HADRONIC_SYSTEM_GENERATOR_H_

