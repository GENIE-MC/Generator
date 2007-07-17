//____________________________________________________________________________
/*!

\class    genie::DISPrimaryLeptonGenerator

\brief    Generates the final state primary lepton in v DIS interactions.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  October 03, 2004

\cpright  Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _DIS_PRIMARY_LEPTON_GENERATOR_H_
#define _DIS_PRIMARY_LEPTON_GENERATOR_H_

#include "EVGModules/PrimaryLeptonGenerator.h"

namespace genie {

class DISPrimaryLeptonGenerator : public PrimaryLeptonGenerator {

public :

  DISPrimaryLeptonGenerator();
  DISPrimaryLeptonGenerator(string config);
  ~DISPrimaryLeptonGenerator();

  //-- implement the EventRecordVisitorI interface

  void ProcessEventRecord(GHepRecord * event_rec) const;
};

}      // genie namespace

#endif // _DIS_PRIMARY_LEPTON_GENERATOR_H_
