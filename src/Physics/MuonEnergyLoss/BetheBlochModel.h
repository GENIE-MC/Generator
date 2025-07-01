//____________________________________________________________________________
/*!

\class    genie::mueloss::BetheBlochModel

\brief    Bethe-Bloch model for muon energy loss due to Ionization
          Concrete implementation of the MuELossI interface.

\ref      W.Lohmann, R.Kopp and R.Voss,
          Energy Loss of Muons in the Energy Range 1-10000 GeV, CERN 85-03

\author   Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  December 10, 2003

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org        
*/
//____________________________________________________________________________

#ifndef _BETHE_BLOCH_MODEL_H_
#define _BETHE_BLOCH_MODEL_H_

#include "Physics/MuonEnergyLoss/MuELossI.h"

namespace genie   {
namespace mueloss {

class BetheBlochModel : public MuELossI
{
public:
  BetheBlochModel();
  BetheBlochModel(string config);
  virtual ~BetheBlochModel();

  //! implement the MuELossI interface
  double        dE_dx   (double E, MuELMaterial_t material) const;
  MuELProcess_t Process (void) const { return eMupIonization; }
};

}      // mueloss namespace
}      // genie   namespace
#endif // _BETHE_BLOCH_MODEL_H_
