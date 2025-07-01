//____________________________________________________________________________
/*!

\class    genie::SKPrimaryLeptonGenerator

\brief    Generates the final state primary lepton in single-Kaon interactions.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  March 20, 2014

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org        
*/
//____________________________________________________________________________

#ifndef _SK_PRIMARY_LEPTON_GENERATOR_H_
#define _SK_PRIMARY_LEPTON_GENERATOR_H_

#include "Physics/Common/PrimaryLeptonGenerator.h"

namespace genie {

class SKPrimaryLeptonGenerator : public PrimaryLeptonGenerator {

public :

  SKPrimaryLeptonGenerator();
  SKPrimaryLeptonGenerator(string config);
 ~SKPrimaryLeptonGenerator();

  // Implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

//void CalculatePrimaryLepton(GHepRecord * event_rec) const;
};

}      // genie namespace

#endif // _SK_PRIMARY_LEPTON_GENERATOR_H_
