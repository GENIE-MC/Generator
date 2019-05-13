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

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _HN_INTRANUKE_2018_H_
#define _HN_INTRANUKE_2018_H_

#include <TGenPhaseSpace.h>

#include "Physics/NuclearState/NuclearModelI.h"
#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Physics/HadronTransport/INukeMode.h"
#include "Physics/HadronTransport/INukeHadroFates.h"
#include "Physics/HadronTransport/Intranuke2018.h"

class TLorentzVector;
class TVector3;

namespace genie {

class GHepParticle;
class INukeHadroData;
class PDGCodeList;

class HNIntranuke2018 : public Intranuke2018 { 

friend class IntranukeTester;

public :
  HNIntranuke2018();
  HNIntranuke2018(string config);
 ~HNIntranuke2018();

  void ProcessEventRecord(GHepRecord * event_rec) const;

  virtual string GetINukeMode() const {return "hN2018";};
  virtual string GetGenINukeMode() const {return "hN";};

private:

  void LoadConfig (void);

  // methods specific to intranuke HN-mode
  void          SimulateHadronicFinalState (GHepRecord* ev, GHepParticle* p) const;
  INukeFateHN_t HadronFateHN      (const GHepParticle* p) const;
  INukeFateHN_t HadronFateOset () const;
  double        FateWeight        (int pdgc, INukeFateHN_t fate) const;
  void          ElasHN	          (GHepRecord* ev, GHepParticle* p, INukeFateHN_t fate) const;
  void          AbsorbHN	  (GHepRecord* ev, GHepParticle* p, INukeFateHN_t fate) const;
  void          InelasticHN	  (GHepRecord* ev, GHepParticle* p) const;
  void          GammaInelasticHN  (GHepRecord* ev, GHepParticle* p, INukeFateHN_t fate) const; 
  bool          HandleCompoundNucleusHN (GHepRecord* ev, GHepParticle* p) const;
  int           HandleCompoundNucleus(GHepRecord* ev, GHepParticle* p, int mom) const;           

   mutable int nuclA;     ///< value of A for the target nucleus in hA mode

  // data members specific to intranuke HN-mode
  double fNucQEFac;

};

}      // genie namespace

#endif // _HN_INTRANUKE_ALT_H_

