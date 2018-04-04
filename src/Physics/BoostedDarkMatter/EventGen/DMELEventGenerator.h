//____________________________________________________________________________
/*!

\class    genie::DMELKinematicsGenerator

\brief    Generates values for the kinematic variables describing EL DM
          interaction events.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Joshua Berger <jberger \at physics.wisc.edu>
          University of Wisconsin-Madison
          Andrew Furmanski

\created  September 4, 2017

\cpright  Copyright (c) 2003-2018, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _DMEL_EVENT_GENERATOR_H_
#define _DMEL_EVENT_GENERATOR_H_

#include "Physics/NuclearState/NuclearModelI.h"
#include "Physics/Common/KineGeneratorWithCache.h"
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

  double ComputeXSec (Interaction * in, double costheta, double phi) const;
  TVector3 COMframe2Lab(InitialState initialState) const;

  double COMJacobian(TLorentzVector lepton, TLorentzVector leptonCOM, TLorentzVector outNucleon, TVector3 beta) const;
  
  // unused // double fQ2min;
  mutable double fEb; // Binding energy

  void   LoadConfig     (void);
  double  ComputeMaxXSec(const Interaction * in) const;

  void AddTargetNucleusRemnant (GHepRecord * evrec) const; ///< add a recoiled nucleus remnant


  const NuclearModelI *  fNuclModel;   ///< nuclear model

  //mutable double fQ2min;
  //mutable double fQ2max;
  //
  mutable double fMinAngleEM;


}; // class definition

} // genie namespace

#endif // _DMEL_EVENT_GENERATOR_H_
