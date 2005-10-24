//____________________________________________________________________________
/*!

\class   genie::RESKinematicsGenerator

\brief   Generates resonance event (v+N->l+Resonance) kinematics.

         Is a concrete implementation of the EventRecordVisitorI interface.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created November 18, 2004

*/
//____________________________________________________________________________

#ifndef _RES_KINEMATICS_GENERATOR_H_
#define _RES_KINEMATICS_GENERATOR_H_

#include "EVGModules/KineGeneratorWithCache.h"
#include "Utils/Range1.h"

namespace genie {

class RESKinematicsGenerator : public KineGeneratorWithCache {

public :

  RESKinematicsGenerator();
  RESKinematicsGenerator(string config);
  ~RESKinematicsGenerator();

  //-- implement the EventRecordVisitorI interface

  void ProcessEventRecord(GHepRecord * event_rec) const;

private:

  Range1D_t WRange  (const Interaction * interaction) const;
  Range1D_t Q2Range (const Interaction * interaction) const;

  double ComputeMaxXSec  (const Interaction * interaction) const;
};

}      // genie namespace

#endif // _RES_KINEMATICS_GENERATOR_H_
