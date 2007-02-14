//____________________________________________________________________________
/*!

\class    genie::QELKinematicsGenerator

\brief    Generates values for the kinematic variables describing QEL neutrino
          interaction events.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  October 03, 2004

\cpright  Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _QEL_KINEMATICS_GENERATOR_H_
#define _QEL_KINEMATICS_GENERATOR_H_

#include "EVGModules/KineGeneratorWithCache.h"
#include "Utils/Range1.h"

namespace genie {

class QELKinematicsGenerator : public KineGeneratorWithCache {

public :
  QELKinematicsGenerator();
  QELKinematicsGenerator(string config);
  ~QELKinematicsGenerator();

  //-- implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void   LoadConfig     (void);
  double ComputeMaxXSec (const Interaction * in) const;
};

}      // genie namespace

#endif // _QEL_KINEMATICS_GENERATOR_H_
