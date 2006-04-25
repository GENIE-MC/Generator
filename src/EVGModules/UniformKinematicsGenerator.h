//____________________________________________________________________________
/*!

\class    genie::UniformKinematicsGenerator

\brief    Generates values for the kinematic variables describing QEL,RES,DIS
          or COH neutrino interaction events (depending on the interaction
          already found at the event record it visits) using a flat probability
          distribution.
          Only use this one if you want to generate samples with many events
          from improbable phase space regions without having to generate very
          large statistics.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  September 29, 2005

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _UNIFORM_KINEMATICS_GENERATOR_H_
#define _UNIFORM_KINEMATICS_GENERATOR_H_

#include "EVGCore/EventRecordVisitorI.h"

namespace genie {

class UniformKinematicsGenerator : public EventRecordVisitorI {

public :

  UniformKinematicsGenerator();
  UniformKinematicsGenerator(string config);
  ~UniformKinematicsGenerator();

  //-- implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

private:
  void GenerateUnifQELKinematics (GHepRecord * evrec) const;
  void GenerateUnifRESKinematics (GHepRecord * evrec) const;
  void GenerateUnifDISKinematics (GHepRecord * evrec) const;
  void GenerateUnifCOHKinematics (GHepRecord * evrec) const;
};

}      // genie namespace

#endif // _UNIFORM_KINEMATICS_GENERATOR_H_

