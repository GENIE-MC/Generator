//____________________________________________________________________________
/*!

\class   genie::DISPrimaryLeptonGenerator

\brief   Generates the final state primary lepton in v DIS interactions.

         Is a concrete implementation of the EventRecordVisitorI interface.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 03, 2004

*/
//____________________________________________________________________________

#ifndef _DIS_PRIMARY_LEPTON_GENERATOR_H_
#define _DIS_PRIMARY_LEPTON_GENERATOR_H_

#include "EVGModules/PrimaryLeptonGenerator.h"

namespace genie {

class DISPrimaryLeptonGenerator : public PrimaryLeptonGenerator {

public :

  DISPrimaryLeptonGenerator();
  DISPrimaryLeptonGenerator(const char * param_set);
  ~DISPrimaryLeptonGenerator();

  //-- implement the EventRecordVisitorI interface

  void ProcessEventRecord(GHepRecord * event_rec) const;
};

}      // genie namespace

#endif // _DIS_PRIMARY_LEPTON_GENERATOR_H_
