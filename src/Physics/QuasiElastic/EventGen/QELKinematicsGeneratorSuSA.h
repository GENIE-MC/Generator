//____________________________________________________________________________
/*!

\class    genie::QELKinematicsGeneratorSuSA

\brief    Event generator for SuSAv2 1p1h interactions

\author   Stephen Dolan <stephen.joseph.dolan \at cern.ch>
          European Organization for Nuclear Research (CERN)

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _QEL_KINEMATICS_GENERATOR_SUSA_H_
#define _QEL_KINEMATICS_GENERATOR_SUSA_H_

#include <TGenPhaseSpace.h>

#include "Physics/NuclearState/NuclearModelI.h"
#include "Framework/ParticleData/PDGCodeList.h"
#include "Physics/Common/KineGeneratorWithCache.h"
#include "Framework/Utils/Range1.h"

namespace genie {

class QELKinematicsGeneratorSuSA : public KineGeneratorWithCache {

public :
  QELKinematicsGeneratorSuSA();
  QELKinematicsGeneratorSuSA(string config);
 ~QELKinematicsGeneratorSuSA();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord (GHepRecord * event) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void    LoadConfig              (void);
  void    AddTargetNucleusRemnant (GHepRecord * event) const;
  void    SelectLeptonKinematics  (GHepRecord * event) const;
  void    GenerateNucleon         (GHepRecord * event) const;
  double  ComputeMaxXSec(const Interaction * in) const;
  double  ComputeMaxXSec_elas(const Interaction * in) const;


  const NuclearModelI *          fNuclModel;

  double fQ3Max;
  bool fForceBound;
  bool fForceEbFromModel;
  bool fForceFixEb;
  double fEbOR;

  // Carbon removal energy - used for scaling
  double fEbC;

};

}      // genie namespace
#endif // _QEL_KINEMATICS_GENERATOR_SUSA_H_
