//____________________________________________________________________________
/*!

\class    genie::OutgoingDarkGenerator

\brief    Abstract class. Is used to pass common implementation to concrete
          implementations of the EventRecordVisitorI interface generating the
          primary lepton for a specific processes (QEL,DIS,RES,IMD,...)

\author   Joshua Berger <jberger \at physics.wisc.edu>
          University of Wisconsin-Madison
          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  October 03, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _OUTGOING_DARK_GENERATOR_H_
#define _OUTGOING_DARK_GENERATOR_H_

class TVector3;
class TLorentzVector;

#include "Framework/EventGen/EventRecordVisitorI.h"

namespace genie {

class OutgoingDarkGenerator : public EventRecordVisitorI {

public :

  //-- Standard EventRecordVisitorI interface implementation
  virtual void ProcessEventRecord(GHepRecord * evrec) const;

  //-- Common methods for all concrete OutgoingDarkGenerator-type
  //   EventRecordVisitors
  virtual void     SetPolarization  (GHepRecord * ev) const;
  virtual TVector3 NucRestFrame2Lab (GHepRecord * ev) const;
  virtual void     AddToEventRecord (
              GHepRecord * ev, int pdgc, const TLorentzVector & p4) const;

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string config);

protected:

  //-- Abstract class - Can only be instantiated by its children.
  OutgoingDarkGenerator();
  OutgoingDarkGenerator(string name);
  OutgoingDarkGenerator(string name, string config);
  virtual ~OutgoingDarkGenerator();

  void LoadConfig(void);

  bool fApplyCoulombCorrection;
};

}      // genie namespace

#endif // _OUTGOING_DARK_GENERATOR_H_
