//____________________________________________________________________________
/*!

\class    genie::COHPrimaryLeptonGenerator

\brief    Generates the final state primary lepton in v COH NC interactions.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  September 26, 2005

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _COH_PRIMARY_LEPTON_GENERATOR_H_
#define _COH_PRIMARY_LEPTON_GENERATOR_H_

#include "Physics/Common/PrimaryLeptonGenerator.h"

namespace genie {

  class COHPrimaryLeptonGenerator : public PrimaryLeptonGenerator {

  public :

    COHPrimaryLeptonGenerator();
    COHPrimaryLeptonGenerator(string config);
    ~COHPrimaryLeptonGenerator();

    //-- implement the EventRecordVisitorI interface

    void ProcessEventRecord(GHepRecord * event_rec) const;

  private :
    
    void CalculatePrimaryLepton_AlvarezRuso(GHepRecord * event_rec) const;

  };

}      // genie namespace

#endif // _COH_PRIMARY_LEPTON_GENERATOR_H_
