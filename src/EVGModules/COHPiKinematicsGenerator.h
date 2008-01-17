//____________________________________________________________________________
/*!

\class    genie::COHPiKinematicsGenerator

\brief    Generates values for the kinematic variables describing coherent 
          neutrino-nucleus pion production events.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  October 03, 2004

\cpright  Copyright (c) 2003-2008, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _COHPI_KINEMATICS_GENERATOR_H_
#define _COHPI_KINEMATICS_GENERATOR_H_

#include "EVGModules/KineGeneratorWithCache.h"
#include "Utils/Range1.h"

class TF2;

namespace genie {

class COHPiKinematicsGenerator : public KineGeneratorWithCache {

public :
  COHPiKinematicsGenerator();
  COHPiKinematicsGenerator(string config);
 ~COHPiKinematicsGenerator();

  //! implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  //! overload the Algorithm::Configure() methods to load private data
  //! members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

public:
  //! methods to load sub-algorithms and config data from the Registry
  void LoadConfig (void);

  //! overload KineGeneratorWithCache method to compute max xsec
  double ComputeMaxXSec (const Interaction * in) const;

  //! overload KineGeneratorWithCache method to get energy
  double Energy         (const Interaction * in) const;

  mutable TF2 * fEnvelope; ///< 2-D envelope used for importance sampling
  double fRo;              ///< nuclear scale parameter
};

}      // genie namespace
#endif // _COHPI_KINEMATICS_GENERATOR_H_
