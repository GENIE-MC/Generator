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
class TVector3;

namespace genie {

// Cummulative interaction probabilities for pi+Fe in 50 MeV bins
// Data from NeuGEN's Intranuke
static const int kPNDataPoints = 12;
static const double kPElastic[kPNDataPoints] = {
       .97,.94,.93,.92,.91,.90,.90,.90,.90,.90,.90,.90
};
static const double kPInelastic[kPNDataPoints] = {
       .51,.53,.55,.56,.56,.57,.57,.57,.57,.57,.57,.57
};
static const double kPAbsorption[kPNDataPoints] = {
       .13,.18,.23,.23,.18,.13,.10,.08,.07,.06,.05,.05
};

class Target;

class Intranuke : public EventRecordVisitorI {

public :
  Intranuke();
  Intranuke(string config);
  ~Intranuke();

  //-- implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

private:

  TVector3       Hadronic3P   (GHepRecord * event_rec) const;
  TLorentzVector StepParticle (const TLorentzVector & x4,
                               const TLorentzVector & p4, double step) const;
  double         MeanFreePath (const Target & target,
                               const TLorentzVector & p4, int pdgc) const;
  InukeInt_t     PionFate     (const Target & target,
                               const TLorentzVector & p4, int pdgc) const;
};

}      // genie namespace

#endif // _INTRANUKE_H_
