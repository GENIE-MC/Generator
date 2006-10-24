//____________________________________________________________________________
/*!

\class    genie::IMDKinematicsGenerator

\brief    Generates values for the kinematic variables describing inverse muon
          decay (IMD) neutrino interaction events.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  July 13, 2005

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _IMD_KINEMATICS_GENERATOR_H_
#define _IMD_KINEMATICS_GENERATOR_H_

#include "EVGModules/KineGeneratorWithCache.h"
#include "Utils/Range1.h"

namespace genie {

class IMDKinematicsGenerator : public KineGeneratorWithCache {

public :

  IMDKinematicsGenerator();
  IMDKinematicsGenerator(string config);
  ~IMDKinematicsGenerator();

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

#endif // _IMD_KINEMATICS_GENERATOR_H_
