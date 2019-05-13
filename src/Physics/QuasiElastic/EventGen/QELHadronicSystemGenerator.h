//____________________________________________________________________________
/*!

\class    genie::QELHadronicSystemGenerator

\brief    Generates the final state hadronic system in v QEL interactions.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  October 03, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _QEL_HADRONIC_SYSTEM_GENERATOR_H_
#define _QEL_HADRONIC_SYSTEM_GENERATOR_H_

#include "Physics/Common/HadronicSystemGenerator.h"

namespace genie {

class QELHadronicSystemGenerator : public HadronicSystemGenerator {

public :
  QELHadronicSystemGenerator();
  QELHadronicSystemGenerator(string config);
 ~QELHadronicSystemGenerator();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord (GHepRecord * event_rec) const;

private:
  void AddRecoilBaryon    (GHepRecord * event_rec) const;
};

}      // genie namespace
#endif // _QEL_HADRONIC_SYSTEM_GENERATOR_H_
