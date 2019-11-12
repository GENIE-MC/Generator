//____________________________________________________________________________
/*!
\class    genie::GLRESKinematicsGenerator

\brief    Generates values for the kinematic variables describing Glashow resonance.
          Is a concrete implementation of the EventRecordVisitorI interface.
          Part of its implementation, related with the caching and retrieval of
          previously computed values, is inherited from the KineGeneratorWithCache
          abstract class.

\author   Alfonso Garcia <alfonsog \at nikhef.nl>
          NIKHEF (Amsterdam)

\created  November 8, 2019

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _GLASHOW_RESONANCE_KINEMATICS_GENERATOR_H_
#define _GLASHOW_RESONANCE_KINEMATICS_GENERATOR_H_

#include "Physics/Common/KineGeneratorWithCache.h"

namespace genie {

class GLRESKinematicsGenerator : public KineGeneratorWithCache {

public :
  GLRESKinematicsGenerator();
  GLRESKinematicsGenerator(string config);
 ~GLRESKinematicsGenerator();

  //-- implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

public:

  //-- methods to load sub-algorithms and config data from the Registry
  void   LoadConfig (void);

  //-- overload KineGeneratorWithCache methods
  double ComputeMaxXSec (const Interaction * in) const;
  double Energy         (const Interaction * in) const;
};

}      // genie namespace
#endif // _GLASHOW_RESONANCE_KINEMATICS_GENERATOR_H_
