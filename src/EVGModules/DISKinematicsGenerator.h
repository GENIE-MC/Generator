//____________________________________________________________________________
/*!

\class   genie::DISKinematicsGenerator

\brief   Generates values for the kinematic variables describing DIS v
         interaction events.

         Is a concrete implementation of the EventRecordVisitorI interface.

         Part of its implementation, related with the caching and retrieval of
         previously computed values, is inherited from the KineGeneratorWithCache
         abstract class.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 03, 2004

*/
//____________________________________________________________________________

#ifndef _DIS_KINEMATICS_GENERATOR_H_
#define _DIS_KINEMATICS_GENERATOR_H_

#include "EVGModules/KineGeneratorWithCache.h"
#include "Utils/Range1.h"

namespace genie {

class DISKinematicsGenerator : public KineGeneratorWithCache {

public :

  DISKinematicsGenerator();
  DISKinematicsGenerator(const char * param_set);
  ~DISKinematicsGenerator();

  //-- implement the EventRecordVisitorI interface

  void ProcessEventRecord(GHepRecord * event_rec) const;

private:

  Range1D_t XRange          (void) const;
  Range1D_t YRange          (void) const;
  Range1D_t WRange          (const Interaction * interaction) const;
  Range1D_t Q2Range         (const Interaction * interaction) const;
  bool      ValidKinematics (const Interaction * interaction) const;
  double    ComputeMaxXSec  (const Interaction * interaction) const;
};

}      // genie namespace

#endif // _DIS_KINEMATICS_GENERATOR_H_
