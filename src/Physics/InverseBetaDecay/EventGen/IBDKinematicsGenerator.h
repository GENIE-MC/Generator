//____________________________________________________________________________
/*!

\class    genie::IBDKinematicsGenerator

\brief    Generates values for the kinematic variables describing IBD neutrino
          interaction events.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Corey Reed <cjreed \at nikhef.nl> - October 29, 2009
          using code from the QELKinematicGenerator written by
          Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  October 29, 2009

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          
*/
//____________________________________________________________________________

#ifndef _IBD_KINEMATICS_GENERATOR_H_
#define _IBD_KINEMATICS_GENERATOR_H_

#include "Physics/Common/KineGeneratorWithCache.h"
#include "Framework/Utils/Range1.h"

namespace genie {

class IBDKinematicsGenerator : public KineGeneratorWithCache {

public :
  IBDKinematicsGenerator();
  IBDKinematicsGenerator(string config);
  virtual ~IBDKinematicsGenerator();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void   LoadConfig     (void);
  double ComputeMaxXSec (const Interaction * in) const;
};

}      // genie namespace
#endif // _IBD_KINEMATICS_GENERATOR_H_
