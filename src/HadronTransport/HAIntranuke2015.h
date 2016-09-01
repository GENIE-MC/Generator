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

#ifndef _HA_INTRANUKE_2015_H_
#define _HA_INTRANUKE_2015_H_

#include <TGenPhaseSpace.h>

#include "Interfaces/NuclearModelI.h"
#include "Algorithm/AlgFactory.h"
#include "EVGCore/EventRecordVisitorI.h"
#include "HadronTransport/INukeMode.h"
#include "HadronTransport/INukeHadroFates.h"
#include "HadronTransport/Intranuke2015.h"

class TLorentzVector;
class TVector3;

namespace genie {

class GHepParticle;
class INukeHadroData2015;
class PDGCodeList;

class HAIntranuke2015 : public Intranuke2015 {

friend class IntranukeTester;

public :
  HAIntranuke2015();
  HAIntranuke2015(string config);
 ~HAIntranuke2015();

  void ProcessEventRecord(GHepRecord * event_rec) const;

private:

  void LoadConfig (void);

  void  SimulateHadronicFinalState           (GHepRecord* ev, GHepParticle* p) const;
  void  SimulateHadronicFinalStateKinematics (GHepRecord* ev, GHepParticle* p) const;

  INukeFateHA_t HadronFateHA     (const GHepParticle* p) const;
  INukeFateHA_t HadronFateOset   (void) const;
  void          Inelastic        (GHepRecord* ev, GHepParticle* p, INukeFateHA_t fate) const;
  void          ElasHA           (GHepRecord* ev, GHepParticle* p, INukeFateHA_t fate) const;
  void          InelasticHA      (GHepRecord* ev, GHepParticle* p, INukeFateHA_t fate) const;
  double        PiBounce         (void) const;
  double        PnBounce         (void) const;
  int           HandleCompoundNucleus(GHepRecord* ev, GHepParticle* p, int mom) const;           

  mutable int nuclA;     ///< value of A for the target nucleus in hA mode
  mutable unsigned int fNumIterations;
};

}      // genie namespace

#endif // _HA_INTRANUKE_2015_H_
