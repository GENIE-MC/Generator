//____________________________________________________________________________
/*!

\class    genie::FermiMover

\brief    It visits the event record & computes a Fermi motion momentum for
          initial state nucleons bound in nuclei.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  October 08, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _FERMI_MOVER_H_
#define _FERMI_MOVER_H_

#include "Framework/EventGen/EventRecordVisitorI.h"

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
  void Emit2ndNucleonFromSRC   (GHepRecord * evrec,
                                const int eject_nucleon_pdg) const;
                                ///^ emit a 2nd nucleon due to short range corellations
  void AddTargetNucleusRemnant (GHepRecord * evrec) const; ///< add a recoiled nucleus remnant

  void LoadConfig (void);

  bool  fKeepNuclOnMassShell;          ///< keep hit bound nucleon on the mass shell?
  bool  fSRCRecoilNucleon;             ///< simulate recoil nucleon due to short range corellation?
  const NuclearModelI *  fNuclModel;   ///< nuclear model
};

}      // genie namespace
#endif // _FERMI_MOVER_H_
