#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_GEANT4_INTERFACE_ENABLED__
//____________________________________________________________________________
/*!

\class    genie::HG4BertCascIntranuke

\brief    Interface to the Geant4 Bertini intranuclear cascade
          A concrete implementation of the EventRecordVisitorI interface

\ref      D.H. Wright and M.H. Kelsey, "The Geant4 Bertini Cascade",
          Nucl. Inst. & Meth. A804 (2015) 175.

\author   Dennis Wright <dwright@slac.stanford.edu>

\created  31 January 2017

*/
//____________________________________________________________________________

#ifndef _HG4BERTCASCINTERNUKE_H_
#define _HG4BERTCASCINTERNUKE_H_

#include "Framework/EventGen/EventRecordVisitorI.h"
//rwh//#include "Physics/HadronTransport/INukeHadroFates.h"
#include "Physics/NuclearState/NuclearModelI.h"
#include "Framework/Conventions/GMode.h"

class TLorentzVector;

class G4ParticleDefinition;
class G4KineticTrackVector;

namespace genie {

class AlgFactory;
class GHepParticle;
class INukeHadroData;


class HG4BertCascIntranuke : public EventRecordVisitorI {

public :
  HG4BertCascIntranuke();
  HG4BertCascIntranuke(string config);
  int G4BertCascade(GHepRecord * event_rec) const;
 ~HG4BertCascIntranuke();

  void ProcessEventRecord(GHepRecord* event_rec) const;

  void Configure(const Registry & config);
  void Configure(string param_set);

private:

  void LoadConfig(void);

  void InitG4Particles() const;
  void TransportHadrons(GHepRecord* ev) const;
  G4ParticleDefinition* PDGtoG4Particle(int pdg) const;
  G4KineticTrackVector* ConvertGenieSecondariesToG4(GHepRecord* evrec) const;

  bool Conserve4Momentum(GHepRecord* ev) const;

  void   GenerateVertex     (GHepRecord * ev) const;
  bool   IsInNucleus        (const GHepParticle* p) const;
  void SetTrackingRadius(const GHepParticle* p) const;

  // utility objects & params
  mutable double fTrackingRadius;  // tracking radius for nucleus current event
  //rwh//INukeHadroData* fHadroData;      // a collection of h+N,h+A data & calculations
  //rwh//AlgFactory* fAlgf;               // algorithm factory instance
  const NuclearModelI* fNuclmodel; // nuclear model used to generate fermi momentum
  mutable int fRemnA;              // remnant nucleus A
  mutable int fRemnZ;              // remnant nucleus Z
  mutable GEvGenMode_t   fGMode;
  // configuration parameters
  double fR0;                      // effective nuclear size param
  double fNR;                      // param multiplying the nuclear radius,
                                   // determining how far to track hadrons
                                   //  beyond the "nuclear boundary"
};

}      // genie namespace

#endif // _HG4BERTCASCINTERNUKE_H_
#endif // __GENIE_GEANT4_INTERFACE_ENABLED__
