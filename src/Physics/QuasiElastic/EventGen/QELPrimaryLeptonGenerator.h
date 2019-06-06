//____________________________________________________________________________
/*!

\class    genie::QELPrimaryLeptonGenerator

\brief    Generates the final state primary lepton in v QEL interactions.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  October 03, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _QEL_PRIMARY_LEPTON_GENERATOR_H_
#define _QEL_PRIMARY_LEPTON_GENERATOR_H_

#include "Physics/Common/PrimaryLeptonGenerator.h"

namespace genie {

class QELPrimaryLeptonGenerator : public PrimaryLeptonGenerator {

public :
  QELPrimaryLeptonGenerator();
  QELPrimaryLeptonGenerator(string config);
 ~QELPrimaryLeptonGenerator();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;
};

}      // genie namespace
#endif // _QEL_PRIMARY_LEPTON_GENERATOR_H_
