//____________________________________________________________________________
/*!

\class   genie::Intranuke

\brief   The INTRANUKE cascading MC for intranuclear rescattering.

         add description here
         this is a EventRecordVisitorI template

         Is a concrete implementation of the EventRecordVisitorI interface.

\author

\created Month xx, yyyy

*/
//____________________________________________________________________________

#ifndef _INTRANUKE_H_
#define _INTRANUKE_H_

#include "EVGCore/EventRecordVisitorI.h"
#include "EVGModules/InukeInt.h"

class TLorentzVector;

namespace genie {

class Target;

class Intranuke : public EventRecordVisitorI {

public :
  Intranuke();
  Intranuke(string config);
  ~Intranuke();

  //-- implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

private:

  TLorentzVector StepParticle (const TLorentzVector & x4,
                               const TLorentzVector & p4, double step) const;
  double         MeanFreePath (const Target & target,
                               const TLorentzVector & p4, int pdgc) const;
  InukeInt_t     PionFate     (void) const;
};

}      // genie namespace

#endif // _INTRANUKE_H_
