//____________________________________________________________________________
/*!

\class    genie::DMDISOutgoingDarkGenerator

\brief    Generates the final state dark matter in DM DIS interactions.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Joshua Berger <jberger \at physics.wisc.edu
          University of Wisconsin-Madison
          Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  September 4, 2017

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          
*/
//____________________________________________________________________________

#ifndef _DMDIS_OUTGOING_DARK_GENERATOR_H_
#define _DMDIS_OUTGOING_DARK_GENERATOR_H_

#include "Physics/Common/OutgoingDarkGenerator.h"

namespace genie {

class DMDISOutgoingDarkGenerator : public OutgoingDarkGenerator {

public :
  DMDISOutgoingDarkGenerator();
  DMDISOutgoingDarkGenerator(string config);
  ~DMDISOutgoingDarkGenerator();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;
};

}      // genie namespace
#endif // _DMDIS_OUTGOING_DARK_GENERATOR_H_
