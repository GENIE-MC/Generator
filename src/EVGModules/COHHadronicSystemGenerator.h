//____________________________________________________________________________
/*!

\class   genie::COHHadronicSystemGenerator

\brief   Generates the final state hadronic system in v COH interactions.

         Is a concrete implementation of the EventRecordVisitorI interface.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 03, 2004

*/
//____________________________________________________________________________

#ifndef _COH_HADRONIC_SYSTEM_GENERATOR_H_
#define _COH_HADRONIC_SYSTEM_GENERATOR_H_

#include "EVGModules/HadronicSystemGenerator.h"

namespace genie {

class COHHadronicSystemGenerator : public HadronicSystemGenerator {

public :
  COHHadronicSystemGenerator();
  COHHadronicSystemGenerator(string config);
  ~COHHadronicSystemGenerator();

  //-- implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;
};

}      // genie namespace
#endif // _COH_HADRONIC_SYSTEM_GENERATOR_H_

