//____________________________________________________________________________
/*!

\class    genie::ASKPrimaryLeptonGenerator

\brief    Generates the final state primary lepton in v ASK NC interactions.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  September 26, 2005

\cpright  Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _ASK_PRIMARY_LEPTON_GENERATOR_H_
#define _ASK_PRIMARY_LEPTON_GENERATOR_H_

#include "EVGModules/PrimaryLeptonGenerator.h"

namespace genie {

class ASKPrimaryLeptonGenerator : public PrimaryLeptonGenerator {

public :

  ASKPrimaryLeptonGenerator();
  ASKPrimaryLeptonGenerator(string config);
  ~ASKPrimaryLeptonGenerator();

  //-- implement the EventRecordVisitorI interface

  void ProcessEventRecord(GHepRecord * event_rec) const;
  void CalculatePrimaryLepton(GHepRecord * event_rec) const;
};

}      // genie namespace

#endif // _ASK_PRIMARY_LEPTON_GENERATOR_H_
