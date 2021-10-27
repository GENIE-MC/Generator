//____________________________________________________________________________
/*!

\class    genie::HENuElGenerator

\brief    Glashow resonance event generator

\author   Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
          University of Liverpool & STFC Rutherford Appleton Laboratory

\created  Feb 15, 2008

\cpright  Copyright (c) 2003-2020, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _HE_NUEL_GENERATOR_H_
#define _HE_NUEL_GENERATOR_H_

#define __GENIE_PYTHIA6_ENABLED__

#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Physics/HELepton/XSection/Born.h"

namespace genie {

class HENuElGenerator : public EventRecordVisitorI {

public :
  HENuElGenerator();
  HENuElGenerator(string config);
 ~HENuElGenerator();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord (GHepRecord * event) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void LoadConfig(void);

  Born * born;

};

}      // genie namespace
#endif // _GLASHOW_RESONANCE_GENERATOR_H_
