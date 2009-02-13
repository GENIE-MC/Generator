//____________________________________________________________________________
/*!

\class    genie::DfrcKinematicsGenerator

\brief    Generates values for the kinematic variables describing DIS v
          interaction events.
          Is a concrete implementation of the EventRecordVisitorI interface.

          Part of its implementation, related with the caching and retrieval of
          previously computed values, is inherited from the KineGeneratorWithCache
          abstract class.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  October 03, 2004

\cpright  Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _DIFFRACTIVE_KINEMATICS_GENERATOR_H_
#define _DIFFRACTIVE_KINEMATICS_GENERATOR_H_

#include "EVGModules/KineGeneratorWithCache.h"
#include "Utils/Range1.h"

namespace genie {

class DfrcKinematicsGenerator : public KineGeneratorWithCache {

public :
  DfrcKinematicsGenerator();
  DfrcKinematicsGenerator(string config);
  ~DfrcKinematicsGenerator();

  //! implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  //! overload the Algorithm::Configure() methods to load private data
  //! members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void   LoadConfig      (void);
  double ComputeMaxXSec  (const Interaction * interaction) const;
};

}      // genie namespace

#endif // _DIFFRACTIVE_KINEMATICS_GENERATOR_H_
