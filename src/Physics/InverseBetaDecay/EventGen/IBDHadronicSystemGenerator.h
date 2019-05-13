//____________________________________________________________________________
/*!

\class    genie::IBDHadronicSystemGenerator

\brief    Generates the final state hadronic system in v IBD interactions.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Corey Reed <cjreed \at nikhef.nl> - October 29, 2009
          using code from the QELKinematicGenerator written by
          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  October 29, 2009

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _IBD_HADRONIC_SYSTEM_GENERATOR_H_
#define _IBD_HADRONIC_SYSTEM_GENERATOR_H_

#include "Physics/Common/HadronicSystemGenerator.h"

namespace genie {

class IBDHadronicSystemGenerator : public HadronicSystemGenerator {

public :
  IBDHadronicSystemGenerator();
  IBDHadronicSystemGenerator(string config);
  virtual ~IBDHadronicSystemGenerator();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord (GHepRecord * event_rec) const;

private:
  void AddRecoilBaryon    (GHepRecord * event_rec) const;
};

}      // genie namespace
#endif // _IBD_HADRONIC_SYSTEM_GENERATOR_H_
