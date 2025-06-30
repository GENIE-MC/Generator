//____________________________________________________________________________
/*!

\class    genie::HELeptonKinematicsGenerator

\brief    Kinematics generator for HELepton.

\author   Alfonso Garcia <aagarciasoto \at km3net.de>
          IFIC & Harvard University

\created  Dec 8, 2021

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _HE_LEPTON_KINEMATICS_GENERATOR_H_
#define _HE_LEPTON_KINEMATICS_GENERATOR_H_

#include "Physics/Common/KineGeneratorWithCache.h"

namespace genie {

class HELeptonKinematicsGenerator : public KineGeneratorWithCache {

public :
  HELeptonKinematicsGenerator();
  HELeptonKinematicsGenerator(string config);
 ~HELeptonKinematicsGenerator();

  //-- implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  //-- methods to load sub-algorithms and config data from the Registry
  void   LoadConfig (void);

  //-- overload KineGeneratorWithCache methods
  double ComputeMaxXSec (const Interaction * in) const;
  double Energy         (const Interaction * in) const;

  double fDeltaN1NuE;
  double fDeltaN1NuMu;
  double fDeltaN1NuTau;

};

}      // genie namespace
#endif // _HE_LEPTON_KINEMATICS_GENERATOR_H_
