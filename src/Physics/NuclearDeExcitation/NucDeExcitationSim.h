//____________________________________________________________________________
/*!

\class    genie::NucDeExcitationSim

\brief    Generates nuclear de-excitation gamma rays

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\ref      16O: 
           H.Ejiri,Phys.Rev.C48,1442(1993); 
           K.Kobayashi et al., Nucl.Phys.B (proc Suppl) 139 (2005)

\created  March 05, 2008

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
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

  //-- implement the EventRecordVisitorI interface
  void ProcessEventRecord (GHepRecord * evrec) const;

private:
  void           OxygenTargetSim      (GHepRecord * evrec) const;
  void           AddPhoton            (GHepRecord * evrec, double E0, double t) const;
  double         PhotonEnergySmearing (double E0, double t) const;
  TLorentzVector Photon4P             (double E) const;
};

}      // genie namespace
#endif // _NUCLEAR_DEEXCITATION_H_
