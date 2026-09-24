//____________________________________________________________________________
/*!

\class    genie::DISNuclearModelI

\brief    Pure abstract base class.
          Defines the EMCModelI interface to be implemented by any physics
          model describing the ratio of l-A DIS cross-sections with respect to
	  the free nucleon calculation

\author   J. Tena Vidal <julia.tena-vidal@ific.uv.es>
          Universitat de Valencia

\created  February, 2026

\cpright  Copyright (c) 2003-2027, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _DIS_NUCLEAR_MODEL_I_H_
#define _DIS_NUCLEAR_MODEL_I_H_

#include <string>
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Interaction/Target.h"
#include "Framework/Interaction/Kinematics.h"

using std::ostringstream;
using namespace genie;
using namespace genie::constants;
using namespace genie::controls;

namespace genie {

class DISNuclearModelI : public Algorithm {

public:
  virtual ~DISNuclearModelI() {};
  virtual double DISACorrection (const Interaction * interaction) const = 0;

protected:

  DISNuclearModelI();
  DISNuclearModelI(string name);
  DISNuclearModelI(string name, string config);

};

}         // genie namespace
#endif    // _DIS_NUCLEAR_MODEL_I_H_
