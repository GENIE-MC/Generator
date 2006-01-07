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

class XSecAlgorithmI;

class DISKinematicsGenerator : public KineGeneratorWithCache {

public :

  DISKinematicsGenerator();
  DISKinematicsGenerator(string config);
  ~DISKinematicsGenerator();

  //-- implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void      LoadSubAlg      (void);
  void      LoadConfigData  (void);
  Range1D_t XRange          (void) const;
  Range1D_t YRange          (void) const;
  Range1D_t WRange          (const Interaction * interaction) const;
  Range1D_t Q2Range         (const Interaction * interaction) const;
  void      SetKineXY       (const Interaction * interaction) const;
  double    ComputeMaxXSec  (const Interaction * interaction) const;

  double fWmin;
  double fWmax;
  double fQ2min;
  double fQ2max;

  const XSecAlgorithmI * fXSecModel;
};

}      // genie namespace

#endif // _DIS_KINEMATICS_GENERATOR_H_
