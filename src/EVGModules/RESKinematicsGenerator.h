//____________________________________________________________________________
/*!

\class    genie::RESKinematicsGenerator

\brief    Generates resonance event (v+N->l+Resonance) kinematics.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  November 18, 2004

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _RES_KINEMATICS_GENERATOR_H_
#define _RES_KINEMATICS_GENERATOR_H_

#include "EVGModules/KineGeneratorWithCache.h"
#include "Utils/Range1.h"

class TF2;

namespace genie {

class RESKinematicsGenerator : public KineGeneratorWithCache {

public :
  RESKinematicsGenerator();
  RESKinematicsGenerator(string config);
  ~RESKinematicsGenerator();

  //-- implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void   LoadConfig      (void);
  double ComputeMaxXSec  (const Interaction * interaction) const;

  mutable TF2 * fEnvelope; ///< 2-D envelope used for importance sampling
  double fWcut;            ///< Wcut parameter in DIS/RES join scheme
};

}      // genie namespace
#endif // _RES_KINEMATICS_GENERATOR_H_
