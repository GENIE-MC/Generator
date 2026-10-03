//____________________________________________________________________________
/*!

\class    genie::DMRESPrimaryLeptonGenerator

\brief    Generates the final state primary lepton in DM RES interactions.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Zach Orr Colorado State

\created  Jan 27 2024

\cpright  Copyright (c) 2003-2023, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _DMRES_OUTGOING_DARK_GENERATOR_H_
#define _DMRES_OUTGOING_DARK_GENERATOR_H_

#include "Physics/Common/OutgoingDarkGenerator.h"

namespace genie {

  class DMRESOutgoingDarkGenerator : public OutgoingDarkGenerator {

public :
  DMRESOutgoingDarkGenerator();
  DMRESOutgoingDarkGenerator(string config);
 ~DMRESOutgoingDarkGenerator();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;
};

}      // genie namespace
#endif // _DMRES_OUTGOING_DARK_GENERATOR_H_
