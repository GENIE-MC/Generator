//____________________________________________________________________________
/*!

\class    genie::DMELOutgoingDarkGenerator

\brief    Generates the final state primary lepton in v DMEL interactions.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Joshua Berger <jberger \at physics.wisc.edu>
          University of Wisconsin-Madison
          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  September 4, 2017

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _DMEL_PRIMARY_LEPTON_GENERATOR_H_
#define _DMEL_PRIMARY_LEPTON_GENERATOR_H_

#include "Physics/Common/OutgoingDarkGenerator.h"

namespace genie {

class DMELOutgoingDarkGenerator : public OutgoingDarkGenerator {

public :
  DMELOutgoingDarkGenerator();
  DMELOutgoingDarkGenerator(string config);
 ~DMELOutgoingDarkGenerator();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;
};

}      // genie namespace
#endif // _DMEL_PRIMARY_LEPTON_GENERATOR_H_
