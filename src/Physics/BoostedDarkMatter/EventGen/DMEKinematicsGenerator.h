//____________________________________________________________________________
/*!

\class    genie::DMEKinematicsGenerator

\brief    Generates kinematics for neutrino-electron events.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  July 13, 2005

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _DME_KINEMATICS_GENERATOR_H_
#define _DME_KINEMATICS_GENERATOR_H_

#include "Physics/Common/KineGeneratorWithCache.h"
#include "Framework/Utils/Range1.h"

namespace genie {

class DMEKinematicsGenerator : public KineGeneratorWithCache {

public :
  DMEKinematicsGenerator();
  DMEKinematicsGenerator(string config);
 ~DMEKinematicsGenerator();

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
#endif // _DME_KINEMATICS_GENERATOR_H_
