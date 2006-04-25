//____________________________________________________________________________
/*!

\class    genie::RESPrimaryLeptonGenerator

\brief    Generates the final state primary lepton in v RES interactions.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  October 03, 2004

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _RES_PRIMARY_LEPTON_GENERATOR_H_
#define _RES_PRIMARY_LEPTON_GENERATOR_H_

#include "EVGModules/PrimaryLeptonGenerator.h"

namespace genie {

class RESPrimaryLeptonGenerator : public PrimaryLeptonGenerator {

public :

  RESPrimaryLeptonGenerator();
  RESPrimaryLeptonGenerator(string config);
  ~RESPrimaryLeptonGenerator();

  //-- implement the EventRecordVisitorI interface

  void ProcessEventRecord(GHepRecord * event_rec) const;
};

}      // genie namespace

#endif // _RES_PRIMARY_LEPTON_GENERATOR_H_
