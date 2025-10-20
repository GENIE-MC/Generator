//____________________________________________________________________________
/*!

\class    genie::DMEOutgoingDarkGenerator

\brief    Generates the final state primary lepton in neutrino-electron events.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  July 13, 2005

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          
*/
//____________________________________________________________________________

#ifndef _DME_OUTGOING_DARK_GENERATOR_H_
#define _DME_OUTGOING_DARK_GENERATOR_H_

#include "Physics/Common/OutgoingDarkGenerator.h"

namespace genie {

class DMEOutgoingDarkGenerator : public OutgoingDarkGenerator {

public :
  DMEOutgoingDarkGenerator();
  DMEOutgoingDarkGenerator(string config);
 ~DMEOutgoingDarkGenerator();

  //-- implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;
};

}      // genie namespace
#endif // _DME_OUTGOING_DARK_GENERATOR_H_
