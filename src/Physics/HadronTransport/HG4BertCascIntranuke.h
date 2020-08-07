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
#include "Framework/Conventions/GBuild.h"

#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Framework/Conventions/GMode.h"
#include "Physics/NuclearState/NuclearModelI.h"


#include <TLorentzVector.h>
class TVector3;

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
  virtual string GetINukeMode() const {return "hG4BertCasc";};


  void Configure(const Registry & config);
  void Configure(string param_set);

private:

  void LoadConfig (void);

  void InitG4Particles() const;
  void TransportHadrons(GHepRecord* ev) const;
  const G4ParticleDefinition* PDGtoG4Particle(int pdg) const;
  G4KineticTrackVector* ConvertGenieSecondariesToG4(GHepRecord* evrec) const;
  G4KineticTrackVector* ConvertGenieSecondariesToG4(std::vector<GHepParticle> partList) const;

  bool   Conserve4Momentum  (GHepRecord* ev) const;
  bool   CanRescatter       (const GHepParticle* p) const;
  bool   IsBaryon           (const GHepParticle* p) const;
  void   GenerateVertex     (GHepRecord * ev) const;
  bool   IsInNucleus        (const GHepParticle* p) const;
  void   SetTrackingRadius  (const GHepParticle* p) const;
  double GenerateStep       (GHepRecord* ev, GHepParticle* p) const;
  bool   NeedsRescattering  (const GHepParticle * p) const;

  // utility objects & params
  mutable double fTrackingRadius;  // tracking radius for nucleus current event

  const NuclearModelI* fNuclmodel; // nuclear model used to generate fermi momentum
  mutable int fRemnA;              // remnant nucleus A
  mutable int fRemnZ;              // remnant nucleus Z
  mutable GEvGenMode_t   fGMode;
  // configuration parameters
  double fR0;                      // effective nuclear size param
  double fNR;                      // param multiplying the nuclear radius,
                                   // determining how far to track hadrons
                                   //  beyond the "nuclear boundary"
  double       fNucRmvE;      ///< binding energy to subtract from cascade nucleons
  double       fDelRPion;     ///< factor by which Pion Compton wavelength gets multiplied to become nuclear size enhancement
  double       fDelRNucleon;  ///< factor by which Nucleon Compton wavelength gets multiplied to become nuclear size enhancement
  double       fHadStep;      ///< step size for intranuclear hadron transport
  double       fNucAbsFac;    ///< absorption xsec correction factor (hN Mode)
  double       fNucCEXFac;    ///< charge exchange xsec correction factor (hN Mode)
  double       fEPreEq;       ///< threshold for pre-equilibrium reaction
  double       fFermiFac;     ///< testing parameter to modify fermi momentum
  double       fFermiMomentum;     ///< whether or not particle collision is pauli blocked
  bool         fUseOset;      ///< Oset model for low energy pion in hN
  bool         fAltOset;      ///< NuWro's table-based implementation (not recommended)
  bool         fXsecNNCorr;   ///< use nuclear medium correction for NN cross section
  bool         fDoFermi;
  double       fPionMFPScale;
  double       fNucleonMFPScale;
};

}      // genie namespace

#endif // _HG4BERTCASCINTERNUKE_H_
#endif // __GENIE_GEANT4_INTERFACE_ENABLED__
