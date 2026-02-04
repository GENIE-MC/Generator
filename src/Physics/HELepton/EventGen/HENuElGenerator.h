//____________________________________________________________________________
/*!

\class    genie::HENuElGenerator

\brief    Generator for high energy neutrino-electron scattering.

\author   Alfonso Garcia <aagarciasoto \at km3net.de>
          IFIC & Harvard University

\created  Dec 8, 2021

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _HE_NUEL_GENERATOR_H_
#define _HE_NUEL_GENERATOR_H_

#include "Framework/Conventions/GBuild.h"
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
