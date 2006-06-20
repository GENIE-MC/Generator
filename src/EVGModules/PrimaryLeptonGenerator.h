//____________________________________________________________________________
/*!

\class    genie::PrimaryLeptonGenerator

\brief    Abstract class. Is used to pass common implementation to concrete
          implementations of the EventRecordVisitorI interface generating the
          primary lepton for a specific processes (QEL,DIS,RES,IMD,...)

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  October 03, 2004

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _PRIMARY_LEPTON_GENERATOR_H_
#define _PRIMARY_LEPTON_GENERATOR_H_

class TVector3;
class TLorentzVector;

#include "EVGCore/EventRecordVisitorI.h"

namespace genie {

class PrimaryLeptonGenerator : public EventRecordVisitorI {

public :

  //! Standard EventRecordVisitorI interface implementation
  virtual void ProcessEventRecord(GHepRecord * evrec) const;

  //! Common methods for all concrete PrimaryLeptonGenerator-type
  //! EventRecordVisitors
  virtual void     SetPolarization  (GHepRecord * ev) const;
  virtual TVector3 NucRestFrame2Lab (GHepRecord * ev) const;
  virtual void     AddToEventRecord (
              GHepRecord * ev, int pdgc, const TLorentzVector & p4) const;
protected:

  //! Abstract class - Can only be instantiated by its children.
  PrimaryLeptonGenerator();
  PrimaryLeptonGenerator(string name);
  PrimaryLeptonGenerator(string name, string config);
  virtual ~PrimaryLeptonGenerator();
};

}      // genie namespace

#endif // _PRIMARY_LEPTON_GENERATOR_H_
