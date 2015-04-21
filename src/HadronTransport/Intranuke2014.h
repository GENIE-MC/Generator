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
          Aaron Meyer <asm58@pitt.edu>, Pittsburgh University
	  Alex Bell, Pittsburgh University
          Hugh Gallagher <gallag@minos.phy.tufts.edu>, Tufts University
          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk> STFC, Rutherford Lab

\created  September 20, 2005

\cpright  Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _INTRANUKE_2014_H_
#define _INTRANUKE_2014_H_

#include <TGenPhaseSpace.h>


#include "Algorithm/AlgFactory.h"
#include "EVGCore/EventRecordVisitorI.h"
#include "Conventions/GMode.h"
#include "HadronTransport/INukeMode.h"
#include "HadronTransport/INukeHadroFates.h"
#include "Nuclear/NuclearModelI.h"

class TLorentzVector;
class TVector3;

namespace genie {

class GHepParticle;
class INukeHadroData2014;
class PDGCodeList;
class HNIntranuke2014;
class HAIntranuke2014;

class Intranuke2014 : public EventRecordVisitorI {

friend class IntranukeTester;

public :
  Intranuke2014();
  Intranuke2014(string name);
  Intranuke2014(string name, string config);
 ~Intranuke2014();

  // implement the EventRecordVisitorI interface 
  virtual void ProcessEventRecord(GHepRecord * event_rec) const;

  // override the Algorithm::Configure methods to load configuration
  // data to protected data members
  void Configure (const Registry & config);
  void Configure (string param_set);

protected:

  // methods for loading configuration
  virtual void LoadConfig (void)=0;

  // general methods for the cascade mc structure
  void   TransportHadrons   (GHepRecord * ev) const;
  void   GenerateVertex     (GHepRecord * ev) const;
  bool   NeedsRescattering  (const GHepParticle* p) const;
  bool   CanRescatter       (const GHepParticle* p) const;
  bool   IsInNucleus        (const GHepParticle* p) const;
  void   SetTrackingRadius  (const GHepParticle* p) const;
  double GenerateStep       (GHepRecord* ev, GHepParticle* p) const;

  // virtual functions for individual modes
  virtual void SimulateHadronicFinalState(GHepRecord* ev, GHepParticle* p) const = 0;
  virtual int HandleCompoundNucleus(GHepRecord* ev, GHepParticle* p, int mom) const = 0;

  // utility objects & params
  mutable double         fTrackingRadius;///< tracking radius for the nucleus in the current event
  mutable TGenPhaseSpace fGenPhaseSpace; ///< a phase space generator
  INukeHadroData2014 *       fHadroData2014;     ///< a collection of h+N,h+A data & calculations
  AlgFactory *           fAlgf;          ///< algorithm factory instance
  const NuclearModelI *  fNuclmodel;     ///< nuclear model used to generate fermi momentum
  mutable int            fRemnA;         ///< remnant nucleus A
  mutable int            fRemnZ;         ///< remnant nucleus Z
  mutable TLorentzVector fRemnP4;        ///< P4 of remnant system
  mutable GEvGenMode_t   fGMode;         ///< event generation mode (lepton+A, hadron+A, ...)

  // configuration parameters
  double       fR0;           ///< effective nuclear size param
  double       fNR;           ///< param multiplying the nuclear radius, determining how far to track hadrons beyond the "nuclear boundary"
  double       fNucRmvE;      ///< binding energy to subtract from cascade nucleons
  double       fDelRPion;     ///< factor by which Pion Compton wavelength gets multiplied to become nuclear size enhancement 
  double       fDelRNucleon;  ///< factor by which Nucleon Compton wavelength gets multiplied to become nuclear size enhancement 
  double       fHadStep;      ///< step size for intranuclear hadron transport
  double       fNucAbsFac;    ///< absorption xsec correction factor (hN Mode)
  double       fNucCEXFac;    ///< charge exchange xsec correction factor (hN Mode)
  double       fEPreEq;       ///< threshold for pre-equilibrium reaction
  double       fFermiFac;     ///< testing parameter to modify fermi momentum
  double       fFreeStep;     ///< produced particle free stem, in fm
  double       fFermiMomentum;     ///< whether or not particle collision is pauli blocked
  bool         fDoFermi;      ///< whether or not to do fermi mom. 
  bool         fDoMassDiff;   ///< whether or not to do mass diff. mode
  bool         fDoCompoundNucleus; ///< whether or not to do compound nucleus considerations
};

}      // genie namespace

#endif // _INTRANUKE_2014_H_
