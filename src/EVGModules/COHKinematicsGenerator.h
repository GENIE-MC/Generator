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

class COHKinematicsGenerator : public KineGeneratorWithCache {

public :

  COHKinematicsGenerator();
  COHKinematicsGenerator(const char * param_set);
  ~COHKinematicsGenerator();

  //-- implement the EventRecordVisitorI interface

  void ProcessEventRecord(GHepRecord * event_rec) const;

public:

  Range1D_t yRange         (const Interaction * in) const;
  double    ComputeMaxXSec (const Interaction * in) const;
};

}      // genie namespace

#endif // _COH_KINEMATICS_GENERATOR_H_
