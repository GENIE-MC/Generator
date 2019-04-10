//____________________________________________________________________________
/*!

\class    genie::MuELossI

\brief    Cross Section Calculation Interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  December 10, 2003

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _MUELOSS_I_H_
#define _MUELOSS_I_H_

#include "Framework/Algorithm/Algorithm.h"
#include "Physics/MuonEnergyLoss/MuELMaterial.h"
#include "Physics/MuonEnergyLoss/MuELProcess.h"

namespace genie   {
namespace mueloss {

const double kMaxMuE = 10000; // 10 TeV

class MuELossI : public Algorithm {

public:
  virtual ~MuELossI();

  virtual double        dE_dx   (double E, MuELMaterial_t m) const = 0;
  virtual MuELProcess_t Process (void) const = 0;

protected:
  MuELossI();
  MuELossI(string name);
  MuELossI(string name, string config);
};

}       // mueloss namespace
}       // genie   namespace

#endif  // _MUELOSS_I_H_
