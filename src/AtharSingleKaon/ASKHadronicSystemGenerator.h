//____________________________________________________________________________
/*!

\class    genie::ASKHadronicSystemGenerator

\brief    Generates the f/s hadronic system in v ASK pi production interactions.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  October 03, 2004

\cpright  Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _ASK_HADRONIC_SYSTEM_GENERATOR_H_
#define _ASK_HADRONIC_SYSTEM_GENERATOR_H_

#include "EVGModules/HadronicSystemGenerator.h"

namespace genie {

class ASKHadronicSystemGenerator : public HadronicSystemGenerator {

public :
  ASKHadronicSystemGenerator();
  ASKHadronicSystemGenerator(string config);
 ~ASKHadronicSystemGenerator();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;
  void CalculateHadronicSystem_AtharSingleKaon(GHepRecord * event_rec) const;
};

}      // genie namespace
#endif // _ASK_HADRONIC_SYSTEM_GENERATOR_H_

