//____________________________________________________________________________
/*!

\class    genie::SRCNuclearRecoil

\brief       Created this new module that controls the addition of the recoil nucleon in the event record 
             and extracts its kinematics

\author   Afroditi Papadopoulou <apapadop \at mit.edu>
          Massachusetts Institute of Technology - October 04, 2019

\created  October 04, 2019

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
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
#include "Physics/NuclearState/SecondNucleonEmissionI.h"

namespace genie {

class SRCNuclearRecoil : public SecondNucleonEmissionI {

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


 protected:
  void LoadConfig (void);  
  
  int SRCRecoilPDG( const GHepParticle & nucleon, const Target & tgt) const; // determine the PDG code of the SRC pair

private:


  double fPPPairPercentage;
  double fPNPairPercentage;

};

}      // genie namespace
#endif // _SRC_NUCLEAR_RECOIL_H_
