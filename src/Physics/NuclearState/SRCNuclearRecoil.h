//____________________________________________________________________________
/*!

\class    genie::SRCNuclearRecoil

<<<<<<< 25269eb370b20cd2d9b07b10212aa503c67430cd
\brief       Created this new module that controls the addition of the recoil nucleon in the event record 
             and extracts its kinematics

\author   Afroditi Papadopoulou <apapadop \at mit.edu>
          Massachusetts Institute of Technology - October 04, 2019

\created  October 04, 2019
=======
\brief    It visits the event record & computes a Fermi motion momentum for
          initial state nucleons bound in nuclei.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  October 08, 2004
>>>>>>> refactorization of FermiMover to accomodate the possibility of an SRC recoil nucleon

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _SRC_NUCLEAR_RECOIL_H_
#define _SRC_NUCLEAR_RECOIL_H_

#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Physics/NuclearState/FermiMomentumTable.h"
#include "Framework/Interaction/Target.h"

namespace genie {

class NuclearModelI;

class SRCNuclearRecoil : public EventRecordVisitorI {

public :
  SRCNuclearRecoil();
  SRCNuclearRecoil(string config);
 ~SRCNuclearRecoil();

  //-- implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  int SRCRecoilPDG(GHepParticle * nucleon, GHepParticle * nucleus, Target* tgt, double pF2) const; // determine the PDG code of the SRC pair

<<<<<<< 25269eb370b20cd2d9b07b10212aa503c67430cd
  void LoadConfig (void);

=======
//  void Emit2ndNucleonFromSRC   (GHepRecord * evrec,
//                                const int eject_nucleon_pdg) const;
//                                ///^ emit a 2nd nucleon due to short range corellations

//  void KickHitNucleon          (GHepRecord * evrec) const; ///< give hit nucleon a momentum

//  void AddTargetNucleusRemnant (GHepRecord * evrec) const; ///< add a recoiled nucleus remnant

  void LoadConfig (void);

//  bool  fKeepNuclOnMassShell;          ///< keep hit bound nucleon on the mass shell?
//  bool  fSRCRecoilNucleon;             ///< simulate recoil nucleon due to short range corellation?

>>>>>>> refactorization of FermiMover to accomodate the possibility of an SRC recoil nucleon
  const NuclearModelI *  fNuclModel;   ///< nuclear model

  double fPPPairPercentage;
  double fPNPairPercentage;
  const FermiMomentumTable * fKFTable;
  string fKFTableName;
};

}      // genie namespace
#endif // _SRC_NUCLEAR_RECOIL_H_
