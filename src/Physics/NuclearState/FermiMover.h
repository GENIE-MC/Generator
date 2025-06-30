//____________________________________________________________________________
/*!

\class    genie::FermiMover

\brief    It visits the event record & computes a Fermi motion momentum for
          initial state nucleons bound in nuclei.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  October 08, 2004

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          
*/
//____________________________________________________________________________

#ifndef _FERMI_MOVER_H_
#define _FERMI_MOVER_H_

#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Physics/NuclearState/FermiMomentumTable.h"
#include "Framework/Interaction/Target.h"
#include "Physics/NuclearState/SRCNuclearRecoil.h"
#include "Physics/NuclearState/SecondNucleonEmissionI.h"

namespace genie {

class NuclearModelI;

class FermiMover : public EventRecordVisitorI {

public :
  FermiMover();
  FermiMover(string config);
 ~FermiMover();

  //-- implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void KickHitNucleon          (GHepRecord * evrec) const; ///< give hit nucleon a momentum

  void AddTargetNucleusRemnant (GHepRecord * evrec) const; ///< add a recoiled nucleus remnant

  void LoadConfig (void);

  bool  fKeepNuclOnMassShell;          ///< keep hit bound nucleon on the mass shell?
  bool  fMomDepErmv;                   ///< use momentum dependent calculation of Ermv
  const NuclearModelI *  fNuclModel;   ///< nuclear model

  const SecondNucleonEmissionI *  fSecondEmitter ; 

};

}      // genie namespace
#endif // _FERMI_MOVER_H_
