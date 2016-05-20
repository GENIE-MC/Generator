//____________________________________________________________________________
/*!

\class    genie::QELKinematicsGenerator

\brief    Generates values for the kinematic variables describing QEL neutrino
          interaction events.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Andrew Furmanski

\created  August 04, 2014

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________


#ifndef _QEL_Event_GENERATOR_H_
#define _QEL_Event_GENERATOR_H_


#include "EVGModules/KineGeneratorWithCache.h"
#include "Nuclear/NuclearModelI.h"
#include "Utils/Range1.h"
#include "Conventions/Controls.h"


namespace genie {

class QELEventGenerator: public KineGeneratorWithCache {

public :
  QELEventGenerator();
  QELEventGenerator(string config);
 ~QELEventGenerator();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

  double ComputeXSec (Interaction * in, double costheta, double phi) const;
  TVector3 COMframe2Lab(InitialState initialState) const;

  double COMJacobian(TLorentzVector lepton, TLorentzVector leptonCOM, TLorentzVector outNucleon, TVector3 beta) const;

private:
  
  double fQ2min;
  mutable double fEb; // Binding energy

  void   LoadConfig     (void);
  double  ComputeMaxXSec(const Interaction * in) const;

  void AddTargetNucleusRemnant (GHepRecord * evrec) const; ///< add a recoiled nucleus remnant


  const NuclearModelI *  fNuclModel;   ///< nuclear model

  mutable double fXSecMax;
  //mutable double fQ2min;
  //mutable double fQ2max;
  //
  mutable double fMinAngleEM;


}; // class definition

} // genie namespace

#endif // _QEL_Event_GENERATOR_H_
