//____________________________________________________________________________
/*!

\class    genie::SKKinematicsGenerator

\brief    Generates values for the kinematic variables describing neutrino-nucleus 
          single kaon production events.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Chris Marshall <marshall \at pas.rochester.edu>
          University of Rochester

\created  October 03, 2014

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _SK_KINEMATICS_GENERATOR_H_
#define _SK_KINEMATICS_GENERATOR_H_

#include "Framework/Utils/Range1.h"
#include "Physics/Common/KineGeneratorWithCache.h"

namespace genie {

class SKKinematicsGenerator : public KineGeneratorWithCache {

public :
  SKKinematicsGenerator();
  SKKinematicsGenerator(string config);
 ~SKKinematicsGenerator();

  // Implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  // Overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

public:
  // Methods to load sub-algorithms and config data from the Registry
  void LoadConfig (void);

  // Different kinematics calculators for different models
  void CalculateKin_AtharSingleKaon(GHepRecord * event_rec) const;

  double ComputeMaxXSec (const Interaction * in) const;

  // Overload KineGeneratorWithCache method to get energy
  double Energy (const Interaction * in) const;

private:

  // In computeMaxXSec method, scan log(1-cos(theta)) from this value up to log(2)
  double fMinLog1MinusCosTheta;

};

}      // genie namespace
#endif // _SK_KINEMATICS_GENERATOR_H_
