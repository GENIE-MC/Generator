//____________________________________________________________________________
/*!

\class    genie::DFRPrimaryLeptonGenerator

\brief    Generates the final state primary lepton in diffractive reactions.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
          University of Liverpool & STFC Rutherford Appleton Laboratory

\created  Feb 15th, 2009

\cpright  Copyright (c) 2003-2020, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _DIFFRACTIVE_PRIMARY_LEPTON_GENERATOR_H_
#define _DIFFRACTIVE_PRIMARY_LEPTON_GENERATOR_H_

#include "Physics/Common/PrimaryLeptonGenerator.h"

namespace genie {

class DFRPrimaryLeptonGenerator : public PrimaryLeptonGenerator {

public :
  DFRPrimaryLeptonGenerator();
  DFRPrimaryLeptonGenerator(string config);
 ~DFRPrimaryLeptonGenerator();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;
};

}      // genie namespace
#endif // _DIFFRACTIVE_PRIMARY_LEPTON_GENERATOR_H_
