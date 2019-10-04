//____________________________________________________________________________
/*!

\class    genie::DMELEventGenerator

\brief    Generates values for the kinematic variables describing DMEL neutrino
          interaction events.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Andrew Furmanski

\created  August 04, 2014

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________


#ifndef _DMEL_EVENT_GENERATOR_H_
#define _DMEL_EVENT_GENERATOR_H_


#include "Physics/NuclearState/NuclearModelI.h"
#include "Physics/Common/KineGeneratorWithCache.h"
#include "Physics/BoostedDarkMatter/XSection/DMELUtils.h"
#include "Framework/Utils/Range1.h"
#include "Framework/Conventions/Controls.h"


namespace genie {

class DMELEventGenerator: public KineGeneratorWithCache {

public :
  DMELEventGenerator();
  DMELEventGenerator(string config);
 ~DMELEventGenerator();

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

  const NuclearModelI *  fNuclModel;   ///< nuclear model

  mutable double fMinAngleEM;

  /// Enum that indicates which approach should be used to handle the binding
  /// energy of the struck nucleon
  DMELEvGen_BindingMode_t fHitNucleonBindingMode;

  /// The number of nucleons to sample from the nuclear model when choosing a maximum
  /// momentum to use in ComputeMaxXSec()
  int fMaxXSecNucleonThrows;

}; // class definition

} // genie namespace

#endif // _DMEL_EVENT_GENERATOR_H_
