//____________________________________________________________________________
/*!

\class    genie::Intranuke

\brief    The INTRANUKE cascading MC for intranuclear rescattering.
          Is a concrete implementation of the EventRecordVisitorI interface.

\ref      R.Merenyi et al., Phys.Rev.D45 (1992)
          R.D.Ransome, Nucl.Phys.B 139 (2005)

          The original INTRANUKE cascade MC was developed (in fortran) for the
          NeuGEN MC by G.F.Pearce, R.Edgecock, W.A.Mann and H.Gallagher.

          Current INTRANUKE development is led by S.Dytman and H.Gallagher.

\author   Steve Dytman <dytman+@pitt.edu>, Pittsburgh University
          Hugh Gallagher <gallag@minos.phy.tufts.edu>, Tufts University
          Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk> CCLRC, Rutherford Lab

\created  September 20, 2005

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _INTRANUKE_H_
#define _INTRANUKE_H_

#include <TGenPhaseSpace.h>

#include "EVGCore/EventRecordVisitorI.h"
#include "HadronTransport/INukeProc.h"
#include "HadronTransport/INukeMode.h"

class TLorentzVector;
class TVector3;

namespace genie {

class GHepParticle;
class INukeHadroData;

class Intranuke : public EventRecordVisitorI {

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

  //-- private methods
  void        LoadConfig             (void);
  void        TransportHadrons       (GHepRecord * ev) const;
  void        GenerateVertex         (GHepRecord * ev) const;
  bool        NeedsRescattering      (const GHepParticle* p) const;
  bool        CanRescatter           (const GHepParticle* p) const;
  bool        IsInNucleus            (const GHepParticle* p) const;
  void        SetNuclearRadius       (const GHepParticle* p) const;
  INukeProc_t HadronFate             (const GHepParticle* p) const;
  void        StepParticle           (GHepParticle * p, double dr) const;
  bool        IsFreshHadron          (GHepRecord* ev, GHepParticle* p) const;
  double      FormationZone          (GHepRecord* ev, GHepParticle* p) const;
  void        AdvanceFreshHadron     (GHepRecord* ev, GHepParticle* p) const;
  double      GenerateStep           (GHepRecord* ev, GHepParticle* p) const;
  double      MeanFreePath           (GHepRecord* ev, GHepParticle* p) const;
  void        SimHadronicInteraction (GHepRecord* ev, GHepParticle* p) const;
  void        SimAbsorption          (GHepRecord* ev, GHepParticle* p) const;
  void        SimChargeExchange      (GHepRecord* ev, GHepParticle* p) const;
  void        SimInelasticScattering (GHepRecord* ev, GHepParticle* p) const;
  void        SimElasticScattering   (GHepRecord* ev, GHepParticle* p) const;

  mutable double fNuclRadius;

  //-- utility objects
  TGenPhaseSpace   fGenPhaseSp; ///< a phase space generator
  INukeHadroData * fHadroData;  ///< a collection of h+N,h+A data & calculations

  //-- configuration parameters
  INukeMode_t  fMode;  ///< h+A, h+N
  double       fct0;   ///< formation zone (c * formation time)
  double       fK;
  double       fR0;
};

}      // genie namespace

#endif // _INTRANUKE_H_
