//____________________________________________________________________________
/*!

\class    genie::QELPrimaryLeptonGenerator

\brief    Generates the final state primary lepton in v QEL interactions.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  October 03, 2004

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _QEL_PRIMARY_LEPTON_GENERATOR_H_
#define _QEL_PRIMARY_LEPTON_GENERATOR_H_

#include "EVGModules/PrimaryLeptonGenerator.h"

namespace genie {

class QELPrimaryLeptonGenerator : public PrimaryLeptonGenerator {

public :

  QELPrimaryLeptonGenerator();
  QELPrimaryLeptonGenerator(string config);
  ~QELPrimaryLeptonGenerator();

  //-- implement the EventRecordVisitorI interface

  void ProcessEventRecord(GHepRecord * event_rec) const;
};

}      // genie namespace

#endif // _QEL_PRIMARY_LEPTON_GENERATOR_H_
