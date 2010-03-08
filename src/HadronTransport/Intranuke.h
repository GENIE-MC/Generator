//____________________________________________________________________________
/*!

\class    genie::Intranuke

\brief    The INTRANUKE intranuclear hadron transport MC.
          Is a concrete implementation of the EventRecordVisitorI interface.

\ref      R.Merenyi et al., Phys.Rev.D45 (1992)
          R.D.Ransome, Nucl.Phys.B 139 (2005)

          Current INTRANUKE development is led by S.Dytman and H.Gallagher.
          The original INTRANUKE cascade MC was developed (in fortran) for the
          NeuGEN MC by R.Edgecock, G.F.Pearce, W.A.Mann, R.Merenyi and others.

\author   Steve Dytman <dytman+@pitt.edu>, Pittsburgh University
          Hugh Gallagher <gallag@minos.phy.tufts.edu>, Tufts University
          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk> STFC, Rutherford Lab

\created  September 20, 2005

\cpright  Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _INTRANUKE_H_
#define _INTRANUKE_H_

#include <TGenPhaseSpace.h>

#include "EVGCore/EventRecordVisitorI.h"
#include "HadronTransport/INukeMode.h"
#include "HadronTransport/INukeHadroFates.h"

class TLorentzVector;
class TVector3;

namespace genie {

class GHepParticle;
class INukeHadroData;
class PDGCodeList;

class Intranuke : public EventRecordVisitorI {

friend class IntranukeTester;

public :
  Intranuke();
  Intranuke(string config);
 ~Intranuke();

  //-- implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

private:

  //-- private methods:

  // methods for loading configuration
  void LoadConfig (void);

  // general methods for the cascade mc structure
  void   TransportHadrons   (GHepRecord * ev) const;
  void   GenerateVertex     (GHepRecord * ev) const;
  bool   NeedsRescattering  (const GHepParticle* p) const;
  bool   CanRescatter       (const GHepParticle* p) const;
  bool   IsInNucleus        (const GHepParticle* p) const;
  void   SetNuclearRadius   (const GHepParticle* p) const;
  void   StepParticle       (GHepParticle * p, double dr) const;
  void   SimHadroProc       (GHepRecord* ev, GHepParticle* p) const;
  double GenerateStep       (GHepRecord* ev, GHepParticle* p) const;
  double MeanFreePath       (GHepRecord* ev, GHepParticle* p) const;

  // methods specific to intranuke HA-mode
  INukeFateHA_t HadronFateHA    (const GHepParticle* p) const;
  double        FateWeight      (int pdgc, INukeFateHA_t fate) const;
  void          SimHadroProcHA  (GHepRecord* ev, GHepParticle* p) const;
  void          Inelastic       (GHepRecord* ev, GHepParticle* p, INukeFateHA_t fate) const;
  void          PiSlam          (GHepRecord* ev, GHepParticle* p, INukeFateHA_t fate) const;
  void          PnSlam          (GHepRecord* ev, GHepParticle* p, INukeFateHA_t fate) const;
  double        PiBounce        (void) const;
  double        PnBounce        (void) const;
  double        Degrade         (double ke) const;

  // methods specific to intranuke HN-mode
  void          SimHadroProcHN  (GHepRecord* ev, GHepParticle* p) const;

  // general phase space decay method
  bool PhaseSpaceDecay (GHepRecord* ev, GHepParticle* p, const PDGCodeList & pdgv) const;

  //-- utility objects & params
  mutable double         fNuclRadius;
  mutable TGenPhaseSpace fGenPhaseSpace; ///< a phase space generator
  INukeHadroData *       fHadroData;     ///< a collection of h+N,h+A data & calculations
  mutable int            fRemnA;         ///< remnant nucleus A
  mutable int            fRemnZ;         ///< remnant nucleus Z
  mutable TLorentzVector fRemnP4;        ///< P4 of remnant system

  //-- configuration parameters
  INukeMode_t  fMode;         ///< intranuke mode (h+A, h+N)
  bool         fInTestMode;   ///<
  double       fR0;           ///< effective nuclear size param
  double       fNR;           ///< param multiplying the nuclear radius, determining how far to track hadrons beyond the "nuclear boundary"
  double       fNucRmvE;      ///< binding energy to subtract from cascade nucleons
  double       fDelRPion;     ///< factor by which Pion Compton wavelength gets multiplied to become nuclear size enhancement 
  double       fDelRNucleon;  ///< factor by which Nucleon Compton wavelength gets multiplied to become nuclear size enhancement 
  double       fHadStep;      ///< step size for intranuclear hadron transport
};

}      // genie namespace

#endif // _INTRANUKE_H_
