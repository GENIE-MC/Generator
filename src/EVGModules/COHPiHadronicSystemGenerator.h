//____________________________________________________________________________
/*!

\class    genie::COHPiHadronicSystemGenerator

\brief    Generates the f/s hadronic system in v COH pi production interactions.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  October 03, 2004

\cpright  Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _COHPI_HADRONIC_SYSTEM_GENERATOR_H_
#define _COHPI_HADRONIC_SYSTEM_GENERATOR_H_

#include "EVGModules/HadronicSystemGenerator.h"

namespace genie {

class COHPiHadronicSystemGenerator : public HadronicSystemGenerator {

public :
  COHPiHadronicSystemGenerator();
  COHPiHadronicSystemGenerator(string config);
  ~COHPiHadronicSystemGenerator();

  //-- implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;
};

}      // genie namespace
#endif // _COHPI_HADRONIC_SYSTEM_GENERATOR_H_

