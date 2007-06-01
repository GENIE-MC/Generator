//____________________________________________________________________________
/*!

\class    genie::NuEPrimaryLeptonGenerator

\brief    Generates the final state primary lepton in neutrino-electron events.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  July 13, 2005

\cpright  Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _NUE_PRIMARY_LEPTON_GENERATOR_H_
#define _NUE_PRIMARY_LEPTON_GENERATOR_H_

#include "EVGModules/PrimaryLeptonGenerator.h"

namespace genie {

class NuEPrimaryLeptonGenerator : public PrimaryLeptonGenerator {

public :
  NuEPrimaryLeptonGenerator();
  NuEPrimaryLeptonGenerator(string config);
 ~NuEPrimaryLeptonGenerator();

  //-- implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;
};

}      // genie namespace
#endif // _NUE_PRIMARY_LEPTON_GENERATOR_H_
