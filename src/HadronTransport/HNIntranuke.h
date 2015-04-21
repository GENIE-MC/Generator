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

#ifndef _HN_INTRANUKE_H_
#define _HN_INTRANUKE_H_

#include <TGenPhaseSpace.h>

#include "Algorithm/AlgFactory.h"
#include "EVGCore/EventRecordVisitorI.h"
#include "HadronTransport/INukeMode.h"
#include "HadronTransport/INukeHadroFates.h"
#include "HadronTransport/Intranuke.h"
#include "Nuclear/NuclearModelI.h"

class TLorentzVector;
class TVector3;

namespace genie {

class GHepParticle;
class INukeHadroData;
class PDGCodeList;

class HNIntranuke : public Intranuke { 

friend class IntranukeTester;

public :
  HNIntranuke();
  HNIntranuke(string config);
 ~HNIntranuke();

  void ProcessEventRecord(GHepRecord * event_rec) const;

private:

  void LoadConfig (void);

  // methods specific to intranuke HN-mode
  void          SimulateHadronicFinalState (GHepRecord* ev, GHepParticle* p) const;
  INukeFateHN_t HadronFateHN      (const GHepParticle* p) const;
  double        FateWeight        (int pdgc, INukeFateHN_t fate) const;
  void          ElasHN	          (GHepRecord* ev, GHepParticle* p, INukeFateHN_t fate) const;
  void          AbsorbHN	  (GHepRecord* ev, GHepParticle* p, INukeFateHN_t fate) const;
  void          InelasticHN	  (GHepRecord* ev, GHepParticle* p) const;
  void          GammaInelasticHN  (GHepRecord* ev, GHepParticle* p, INukeFateHN_t fate) const; 
  bool          HandleCompoundNucleus(GHepRecord* ev, GHepParticle* p, int mom) const;           

  // data members specific to intranuke HN-mode
  double fNucQEFac;

};

}      // genie namespace

#endif // _HN_INTRANUKE_H_

