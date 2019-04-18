//____________________________________________________________________________
/*!

\class    genie::COHPiHadronicSystemGenerator

\brief    Generates the f/s hadronic system in v COH elastic interactions.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  November 29, 2007

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _COHEL_HADRONIC_SYSTEM_GENERATOR_H_
#define _COHEL_HADRONIC_SYSTEM_GENERATOR_H_

#include "Physics/Common/HadronicSystemGenerator.h"

namespace genie {

class COHElHadronicSystemGenerator : public HadronicSystemGenerator {

public :
  COHElHadronicSystemGenerator();
  COHElHadronicSystemGenerator(string config);
 ~COHElHadronicSystemGenerator();

  //-- implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;
};

}      // genie namespace
#endif // _COHEL_HADRONIC_SYSTEM_GENERATOR_H_

