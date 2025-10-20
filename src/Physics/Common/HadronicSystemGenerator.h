//____________________________________________________________________________
/*!

\class    genie::HadronicSystemGenerator

\brief    Abstract class. Is used to pass some commonly recurring methods to
          all concrete implementations of the EventRecordVisitorI interface
          generating the hadronic system for a specific processes (QEL,DIS,
          RES,...)

\author   Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  July 16, 2005

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _HADRONIC_SYSTEM_GENERATOR_H_
#define _HADRONIC_SYSTEM_GENERATOR_H_

#include <TLorentzVector.h>
#include <TVector3.h>

#include "Framework/EventGen/EventRecordVisitorI.h"

namespace genie {

class HadronicSystemGenerator : public EventRecordVisitorI {

public :

  // Do not implement the EventRecordVisitorI interface.
  // Leave it for the concrete subclasses


  // Common methods for all concrete subclasses

  void AddTargetNucleusRemnant  (GHepRecord * event_rec) const;
  void AddFinalHadronicSyst     (GHepRecord * event_rec) const;
  void PreHadronTransportDecays (GHepRecord * event_rec) const;

  TLorentzVector Hadronic4pLAB       (GHepRecord * event_rec) const;
  TLorentzVector MomentumTransferLAB (GHepRecord * event_rec) const;
  TVector3       HCM2LAB             (GHepRecord * event_rec) const;
  int            HadronShowerCharge  (GHepRecord * event_rec) const;
  int            ResonanceCharge     (GHepRecord * event_rec) const;

protected:

  // Abstract base class
  HadronicSystemGenerator();
  HadronicSystemGenerator(string name);
  HadronicSystemGenerator(string name, string config);
 ~HadronicSystemGenerator();

  const EventRecordVisitorI * fPreINukeDecayer;
};

}      // genie namespace

#endif // _HADRONIC_SYSTEM_GENERATOR_H_
