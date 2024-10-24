//____________________________________________________________________________
/*!

\class    genie::QELEventGeneratorSuSA

\brief    Event generator for SuSAv2 1p1h interactions

\author   Stephen Dolan <stephen.joseph.dolan \at cern.ch>
          European Organization for Nuclear Research (CERN)

          Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory

\cpright  Copyright (c) 2003-2023, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _QEL_KINEMATICS_GENERATOR_MARTINI_H_
#define _QEL_KINEMATICS_GENERATOR_MARTINI_H_

#include <TGenPhaseSpace.h>

#include "Physics/NuclearState/NuclearModelI.h"
#include "Framework/ParticleData/PDGCodeList.h"
#include "Physics/Common/KineGeneratorWithCache.h"
#include "Framework/Utils/Range1.h"

namespace genie {

class QELEventGeneratorMartini : public KineGeneratorWithCache {

public :
  QELEventGeneratorMartini();
  QELEventGeneratorMartini(string config);
 ~QELEventGeneratorMartini();

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

  const NuclearModelI *          fNuclModel;

  double fQ3Max;
  bool fForceBound;
  bool fForceEbFromModel;
  bool fForceFixEb;
  double fEbOR;

  // Carbon removal energy - used for scaling
  double fEbC;

  /// Delegate event generation for free nucleon targets (which are
  /// not handled by the SuSAv2 calculation) to a different event
  /// generator
  const EventRecordVisitorI* fFreeNucleonEventGenerator;
};

}      // genie namespace
#endif // _QEL_KINEMATICS_GENERATOR_MARTINI_H_
