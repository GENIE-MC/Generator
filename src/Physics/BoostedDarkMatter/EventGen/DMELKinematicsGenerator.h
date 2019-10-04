//____________________________________________________________________________
/*!

\class    genie::DMELKinematicsGenerator

\brief    Generates values for the kinematic variables describing DM elastic
          interaction events.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Joshua Berger <jberger \at physics.wisc.edu>
          University of Wisconsin-Madison
          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  September 4, 2017

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _DM_KINEMATICS_GENERATOR_H_
#define _DM_KINEMATICS_GENERATOR_H_

#include "Physics/Common/KineGeneratorWithCache.h"
#include "Framework/Utils/Range1.h"

namespace genie {

class DMELKinematicsGenerator : public KineGeneratorWithCache {

public :
  DMELKinematicsGenerator();
  DMELKinematicsGenerator(string config);
 ~DMELKinematicsGenerator();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void SpectralFuncExperimentalCode(GHepRecord * event_rec) const;

  void   LoadConfig     (void);
  double ComputeMaxXSec (const Interaction * in) const;
};

}      // genie namespace
#endif // _DM_KINEMATICS_GENERATOR_H_
