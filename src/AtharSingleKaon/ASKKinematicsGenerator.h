//____________________________________________________________________________
/*!

\class    genie::ASKKinematicsGenerator

\brief    Generates values for the kinematic variables describing neutrino-nucleus 
          single kaon production events.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Chris Marshall <marshall \at pas.rochester.edu>
          University of Rochester

\created  October 03, 2014

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _ASK_KINEMATICS_GENERATOR_H_
#define _ASK_KINEMATICS_GENERATOR_H_

#include "EVGModules/KineGeneratorWithCache.h"
#include "Utils/Range1.h"

class TF3;

namespace genie {

class ASKKinematicsGenerator : public KineGeneratorWithCache {

public :
  ASKKinematicsGenerator();
  ASKKinematicsGenerator(string config);
 ~ASKKinematicsGenerator();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

public:
  // methods to load sub-algorithms and config data from the Registry
  void LoadConfig (void);

  // different kinematics calculators for different models
  void   CalculateKin_AtharSingleKaon(GHepRecord * event_rec) const;

  double ComputeMaxXSec (const Interaction * in) const;

  // overload KineGeneratorWithCache method to get energy
  double Energy         (const Interaction * in) const;

private:

  // In computeMaxXSec method, scan log(1-cos(theta)) from this value up to log(2)
  double fMinLog1MinusCosTheta;

};

}      // genie namespace
#endif // _ASK_KINEMATICS_GENERATOR_H_
