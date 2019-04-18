//____________________________________________________________________________
/*!

\class    genie::DMDISKinematicsGenerator

\brief    Generates values for the kinematic variables describing DIS DM
          interaction events.
          Is a concrete implementation of the EventRecordVisitorI interface.

          Part of its implementation, related with the caching and retrieval of
          previously computed values, is inherited from the KineGeneratorWithCache
          abstract class.

\author   Joshua Berger <jberger \at physics.wisc.edu>
          University of Wisconsin-Madison
          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  September 1, 2017

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _DMDIS_KINEMATICS_GENERATOR_H_
#define _DMDIS_KINEMATICS_GENERATOR_H_

#include "Physics/Common/KineGeneratorWithCache.h"
#include "Framework/Utils/Range1.h"

namespace genie {

class DMDISKinematicsGenerator : public KineGeneratorWithCache {

public :
  DMDISKinematicsGenerator();
  DMDISKinematicsGenerator(string config);
  ~DMDISKinematicsGenerator();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void   LoadConfig      (void);
  double ComputeMaxXSec  (const Interaction * interaction) const;
};

}      // genie namespace

#endif // _DMDIS_KINEMATICS_GENERATOR_H_
