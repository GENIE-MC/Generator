//____________________________________________________________________________
/*!

\class    genie::QELKinematicsGeneratorSuSA

\brief    Event generator for SuSAv2 1p1h interactions

\author   Stephen Dolan

\created  Sep. 22, 2008

\cpright  Copyright (c) 2003-2018, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _QEL_KINEMATICS_GENERATOR_SUSA_H_
#define _QEL_KINEMATICS_GENERATOR_SUSA_H_

#include <TGenPhaseSpace.h>

#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Framework/ParticleData/PDGCodeList.h"

namespace genie {

class XSecAlgorithmI;
class NuclearModelI;

class QELKinematicsGeneratorSuSA : public EventRecordVisitorI {

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
  
  mutable const XSecAlgorithmI * fXSecModel;
  const NuclearModelI *          fNuclModel;

  double fQ3Max;
};

}      // genie namespace
#endif // _QEL_KINEMATICS_GENERATOR_SUSA_H_
