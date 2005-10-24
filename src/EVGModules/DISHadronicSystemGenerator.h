//____________________________________________________________________________
/*!

\class   genie::DISHadronicSystemGenerator

\brief   Generates the final state hadronic system in v DIS interactions.

         Is a concrete implementation of the EventRecordVisitorI interface.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 03, 2004

*/
//____________________________________________________________________________

#ifndef _DIS_HADRONIC_SYSTEM_GENERATOR_H_
#define _DIS_HADRONIC_SYSTEM_GENERATOR_H_

#include <TVector3.h>
#include "EVGModules/HadronicSystemGenerator.h"

namespace genie {

class DISHadronicSystemGenerator : public HadronicSystemGenerator {

public :

  DISHadronicSystemGenerator();
  DISHadronicSystemGenerator(string config);
  ~DISHadronicSystemGenerator();

  //-- implement the EventRecordVisitorI interface

  void ProcessEventRecord(GHepRecord * event_rec) const;

private:

  void AddFragmentationProducts (GHepRecord * event_rec) const;

  TVector3 HCM2LAB (GHepRecord * event_rec) const;
};

}      // genie namespace

#endif // _DIS_HADRONIC_SYSTEM_GENERATOR_H_
