//____________________________________________________________________________
/*!

\class   genie::PrimaryLeptonGenerator

\brief   Abstract class. Is used to pass common implementation to concrete
         implementations of the EventRecordVisitorI interface generating the
         primary lepton for a specific processes (QEL,DIS,RES,IMD,...)

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 03, 2004

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

  //-- Do not implement the EventRecordVisitorI interface.
  //   Leave it for its children that are EventRecordVisitors too

  //-- Common methods for all concrete PrimaryLeptonGenerator-type
  //   EventRecordVisitors

  TVector3 *
      NucRestFrame2Lab (GHepRecord * ev) const;
  TLorentzVector *
      P4InNucRestFrame (GHepRecord * ev, double costh, double El) const;
  TLorentzVector *
      Rotate4P (TLorentzVector * p4nu, int pdgc, double cThSc, double El) const;
  void
     AddToEventRecord (GHepRecord * ev, int pdgc, const TLorentzVector * p4) const;

protected:

  //-- Abstract class - Can only be instantiated by its children.

  PrimaryLeptonGenerator();
  PrimaryLeptonGenerator(string name);
  PrimaryLeptonGenerator(string name, string config);
  ~PrimaryLeptonGenerator();
};

}      // genie namespace

#endif // _PRIMARY_LEPTON_GENERATOR_H_
