//____________________________________________________________________________
/*!

\class   genie::IMDKinematicsGenerator

\brief   Generates values for the kinematic variables describing inverse muon
         decay (IMD) neutrino interaction events.

         Is a concrete implementation of the EventRecordVisitorI interface.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created July 13, 2005

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
  IMDKinematicsGenerator(const char * param_set);
  ~IMDKinematicsGenerator();

  //-- implement the EventRecordVisitorI interface

  void ProcessEventRecord(GHepRecord * event_rec) const;

public:

  double  ComputeMaxXSec (const Interaction * in) const;
};

}      // genie namespace

#endif // _IMD_KINEMATICS_GENERATOR_H_
