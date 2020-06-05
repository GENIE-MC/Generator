//____________________________________________________________________________
/*!

\class    genie::NuEPrimaryLeptonGenerator

\brief    Generates the final state primary lepton in neutrino-electron events.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
          University of Liverpool & STFC Rutherford Appleton Laboratory

\created  July 13, 2005

\cpright  Copyright (c) 2003-2020, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _NUE_PRIMARY_LEPTON_GENERATOR_H_
#define _NUE_PRIMARY_LEPTON_GENERATOR_H_

#include "Physics/Common/PrimaryLeptonGenerator.h"

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
