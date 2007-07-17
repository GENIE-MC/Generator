//____________________________________________________________________________
/*!

\class    genie::MuELossI

\brief    Cross Section Calculation Interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  December 10, 2003

\cpright  Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _MUELOSS_I_H_
#define _MUELOSS_I_H_

#include "Algorithm/Algorithm.h"
#include "MuELoss/MuELMaterial.h"
#include "MuELoss/MuELProcess.h"

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
