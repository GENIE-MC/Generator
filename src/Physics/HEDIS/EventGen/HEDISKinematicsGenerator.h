//____________________________________________________________________________
/*!

\class    genie::HEDISKinematicsGenerator

\brief    Generates values for the kinematic variables describing HEDIS v
          interaction events.
          Is a concrete implementation of the EventRecordVisitorI interface.

          Max Xsec are precomputed as the total xsec and they are stored in 
          ascii files. This saves a lot of computational time.

\author   Alfonso Garcia <alfonsog \at nikhef.nl>
          NIKHEF

\created  August 28, 2019

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _HEDIS_KINEMATICS_GENERATOR_H_
#define _HEDIS_KINEMATICS_GENERATOR_H_

#include "Physics/Common/KineGeneratorWithCache.h"

namespace genie {

class HEDISKinematicsGenerator : public KineGeneratorWithCache {

public :
  HEDISKinematicsGenerator();
  HEDISKinematicsGenerator(string config);
  ~HEDISKinematicsGenerator();

  double ComputeMaxXSec       (const Interaction * interaction) const;

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  double Scan(Interaction * interaction, Range1D_t xrange,Range1D_t Q2range, int NKnotsQ2, int NKnotsX, double ME2, double & x_scan, double & Q2_scan) const;

  void   LoadConfig           (void);

  int    fWideNKnotsX;
  int    fWideNKnotsQ2;
  double fWideRange;
  int    fFineNKnotsX;
  int    fFineNKnotsQ2;
  double fFineRange;

  double fSFXmin;   ///< minimum value of x for which SF tables are computed
  double fSFQ2min;  ///< minimum value of Q2 for which SF tables are computed 
  double fSFQ2max;  ///< maximum value of Q2 for which SF tables are computed 

};

}      // genie namespace

#endif // _HEDIS_KINEMATICS_GENERATOR_H_
