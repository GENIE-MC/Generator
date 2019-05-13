//____________________________________________________________________________
/*!

\class    genie::DFRKinematicsGenerator

\brief    Generates kinematics for diffractive reactions.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  Feb 15, 2009

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _DIFFRACTIVE_KINEMATICS_GENERATOR_H_
#define _DIFFRACTIVE_KINEMATICS_GENERATOR_H_

#include "Physics/Common/KineGeneratorWithCache.h"
#include "Framework/Utils/Range1.h"

namespace genie {

class DFRKinematicsGenerator : public KineGeneratorWithCache {

public :
  DFRKinematicsGenerator();
  DFRKinematicsGenerator(string config);
 ~DFRKinematicsGenerator();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void   LoadConfig      (void);
  double ComputeMaxXSec  (const Interaction * interaction) const;

  double fBeta;
};

}      // genie namespace
#endif // _DIFFRACTIVE_KINEMATICS_GENERATOR_H_
