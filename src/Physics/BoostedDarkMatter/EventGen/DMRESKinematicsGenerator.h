//____________________________________________________________________________
/*!

\class    genie::DMRESKinematicsGenerator

\brief    Generates resonance event (DM+N->l+Resonance) kinematics.
          Is a concrete implementation of the EventRecordVisitorI interface.
          Follows same procedure as RESKinematicsGenerator.h

\author   Zach Orr
          Colorado State University

\created  November 18, 2004

\cpright  Copyright (c) 2003-2023, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _DMRES_KINEMATICS_GENERATOR_H_
#define _DMRES_KINEMATICS_GENERATOR_H_

#include "Framework/Utils/Range1.h"
#include "Physics/Common/KineGeneratorWithCache.h"

class TF2;

namespace genie {

class DMRESKinematicsGenerator : public KineGeneratorWithCache {

public :
  DMRESKinematicsGenerator();
  DMRESKinematicsGenerator(string config);
 ~DMRESKinematicsGenerator();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void   LoadConfig      (void);
  double ComputeMaxXSec  (const Interaction * interaction) const;

  mutable TF2 * fEnvelope; ///< 2-D envelope used for importance sampling
  double fWcut;            ///< Wcut parameter in DIS/RES join scheme
};

}      // genie namespace
#endif // _DMRES_KINEMATICS_GENERATOR_H_
