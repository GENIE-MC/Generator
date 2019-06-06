//____________________________________________________________________________
/*!

\class    genie::COHElKinematicsGenerator

\brief    Generates values for the kinematic variables describing coherent 
          neutrino-nucleus elastic scattering events.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  November 22, 2007

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _COHEL_KINEMATICS_GENERATOR_H_
#define _COHEL_KINEMATICS_GENERATOR_H_

#include "Physics/Common/KineGeneratorWithCache.h"
#include "Framework/Utils/Range1.h"

namespace genie {

class COHElKinematicsGenerator : public KineGeneratorWithCache {

public :
  COHElKinematicsGenerator();
  COHElKinematicsGenerator(string config);
 ~COHElKinematicsGenerator();

  //-- implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  //-- overload the Algorithm::Configure() methods to load private data
  //-- members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

public:
  //-- methods to load sub-algorithms and config data from the Registry
  void LoadConfig (void);

  //-- overload KineGeneratorWithCache method to compute max xsec
  double ComputeMaxXSec (const Interaction * in) const;

  //-- overload KineGeneratorWithCache method to get energy
  double Energy         (const Interaction * in) const;
};

}      // genie namespace
#endif // _COHEL_KINEMATICS_GENERATOR_H_
