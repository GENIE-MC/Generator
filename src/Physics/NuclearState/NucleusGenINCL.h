//____________________________________________________________________________
/*!

\class    genie::NucleusGenINCL

\brief    It visits the event record & computes a Fermi motion momentum for
          initial state nucleons bound in nuclei.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Liang Liu <liangliu \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  October 17, 2024

\cpright  Copyright (c) 2003-2024, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          
*/
//____________________________________________________________________________

#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_INCL_ENABLED__

#ifndef _NUCLEUS_GEN_INCL_H_
#define _NUCLEUS_GEN_INCL_H_

#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Physics/NuclearState/FermiMomentumTable.h"
#include "Framework/Interaction/Target.h"
#include "Physics/NuclearState/SRCNuclearRecoil.h"
#include "Physics/NuclearState/SecondNucleonEmissionI.h"

namespace genie {


class NucleusGenINCL : public EventRecordVisitorI {

public :
  NucleusGenINCL();
  NucleusGenINCL(string config);
 ~NucleusGenINCL();

  //-- implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void setInitialStateVertex   (GHepRecord * evrec) const; ///< give hit nucleon a position

  void setInitialStateMomentum (GHepRecord * evrec) const; ///< give hit nucleon a momentum

  void setTargetNucleusRemnant (GHepRecord * evrec) const; ///< add a recoiled nucleus remnant

  // TODO; check the energy-momentum conservation

  void LoadConfig (void);

  // bool  fKeepNuclOnMassShell;          ///< keep hit bound nucleon on the mass shell?
  // bool  fMomDepErmv;                   ///< use momentum dependent calculation of Ermv
  // const NuclearModelI *  fNuclModel;   ///< nuclear model

  // const SecondNucleonEmissionI *  fSecondEmitter ; 

};

}      // genie namespace



#endif // _NUCLEUS_GEN_INCL_H_
#endif // end  __GENIE_INCL_ENABLED__
