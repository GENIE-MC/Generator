//____________________________________________________________________________
/*!

\class    genie::NucDeExcitationSim

\brief    Generates nuclear de-excitation gamma rays

\author   Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\ref      16O:
           H.Ejiri,Phys.Rev.C48,1442(1993);
           K.Kobayashi et al., Nucl.Phys.B (proc Suppl) 139 (2005)

\created  March 05, 2008

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _NUCLEAR_DEEXCITATION_H_
#define _NUCLEAR_DEEXCITATION_H_

#include <TLorentzVector.h>

#include "Framework/EventGen/EventRecordVisitorI.h"

namespace genie {

class NucDeExcitationSim : public EventRecordVisitorI {

public :
  NucDeExcitationSim();
  NucDeExcitationSim(string config);
 ~NucDeExcitationSim();

  void Configure( const Registry& config ) override;
  void Configure( std::string config ) override;

  //-- implement the EventRecordVisitorI interface
  void ProcessEventRecord (GHepRecord * evrec) const override;

private:
  void           OxygenTargetSim      (GHepRecord * evrec) const;
  void           CarbonTargetSim      (GHepRecord * evrec) const;
  void           ArgonTargetSim       (GHepRecord * evrec) const;
  void           AddPhoton            (GHepRecord * evrec, double E0, double t) const;
  double         PhotonEnergySmearing (double E0, double t) const;
  TLorentzVector Photon4P             (double E) const;

  void           LoadConfig();

  // Configuration flags that enable/disable de-excitations for specific nuclei
  bool fDoCarbon = false;
  bool fDoArgon = false;
};

}      // genie namespace
#endif // _NUCLEAR_DEEXCITATION_H_
