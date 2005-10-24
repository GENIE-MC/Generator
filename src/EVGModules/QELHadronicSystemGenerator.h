//____________________________________________________________________________
/*!

\class   genie::QELHadronicSystemGenerator

\brief   Generates the final state hadronic system in v QEL interactions.

         Is a concrete implementation of the EventRecordVisitorI interface.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 03, 2004

*/
//____________________________________________________________________________

#ifndef _QEL_HADRONIC_SYSTEM_GENERATOR_H_
#define _QEL_HADRONIC_SYSTEM_GENERATOR_H_

#include "EVGModules/HadronicSystemGenerator.h"

namespace genie {

class QELHadronicSystemGenerator : public HadronicSystemGenerator {

public :

  QELHadronicSystemGenerator();
  QELHadronicSystemGenerator(string config);
  ~QELHadronicSystemGenerator();

  //-- implement the EventRecordVisitorI interface

  void ProcessEventRecord (GHepRecord * event_rec) const;

private:

  void AddRecoilNucleon (GHepRecord * event_rec) const;

};

}      // genie namespace

#endif // _QEL_HADRONIC_SYSTEM_GENERATOR_H_
