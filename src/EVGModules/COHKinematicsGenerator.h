//____________________________________________________________________________
/*!

\class   genie::COHKinematicsGenerator

\brief   Generates values for the kinematic variables describing coherent NC
         neutrino interaction events.

         Is a concrete implementation of the EventRecordVisitorI interface.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 03, 2004

*/
//____________________________________________________________________________

#ifndef _COH_KINEMATICS_GENERATOR_H_
#define _COH_KINEMATICS_GENERATOR_H_

#include "EVGModules/KineGeneratorWithCache.h"
#include "Utils/Range1.h"

namespace genie {

class XSecAlgorithmI;

class COHKinematicsGenerator : public KineGeneratorWithCache {

public :

  COHKinematicsGenerator();
  COHKinematicsGenerator(string config);
  ~COHKinematicsGenerator();

  //-- implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

public:

  //-- methods to load sub-algorithms and config data from the Registry
  void LoadSubAlg     (void);
  void LoadConfigData (void);

  //-- compute kinematical limits
  Range1D_t yRange (const Interaction * in) const;

  //-- overload KineGeneratorWithCache methods
  double ComputeMaxXSec (const Interaction * in) const;
  double Energy         (const Interaction * in) const;

  //-- private data members
  double                 fSafetyFactor;
  const XSecAlgorithmI * fXSecModel;
};

}      // genie namespace

#endif // _COH_KINEMATICS_GENERATOR_H_
