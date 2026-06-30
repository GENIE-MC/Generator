//____________________________________________________________________________
/*!

\class    genie::QELEventGeneratorINCL

\brief    Generates values for the kinematic variables describing QEL neutrino
          interaction events.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Andrew Furmanski

\created  August 04, 2014

\cpright  Copyright (c) 2003-2024, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________
#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_INCL_ENABLED__

#ifndef _QEL_EVENT_GENERATOR_INCL_H_
#define _QEL_EVENT_GENERATOR_INCL_H_

#include "Physics/NuclearState/NucleusGenI.h"
#include "Physics/Common/KineGeneratorWithCache.h"
#include "Physics/QuasiElastic/XSection/QELUtils.h"
#include "Framework/Utils/Range1.h"
#include "Framework/Conventions/Controls.h"

namespace genie {

class QELEventGeneratorINCL: public KineGeneratorWithCache {

public :
  QELEventGeneratorINCL();
  QELEventGeneratorINCL(string config);
 ~QELEventGeneratorINCL();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  mutable double fEb; // Binding energy

  void   LoadConfig     (void);
  double ComputeMaxXSec(const Interaction* in) const;

  void AddTargetNucleusRemnant (GHepRecord * evrec) const; ///< add a recoiled nucleus remnant

  const NucleusGenI   *  fNucleusGen;  ///< nucleus generator

  mutable double fMinAngleEM;

  /// Enum that indicates which approach should be used to handle the binding
  /// energy of the struck nucleon
  QELEvGen_BindingMode_t fHitNucleonBindingMode;

  /// The number of nucleons to sample from the nuclear model when choosing a maximum
  /// momentum to use in ComputeMaxXSec()
  int fMaxXSecNucleonThrows;

}; // class definition

} // genie namespace

#endif // _QEL_EVENT_GENERATOR_INCL_H_

#endif // __GENIE_INCL_ENABLED__
