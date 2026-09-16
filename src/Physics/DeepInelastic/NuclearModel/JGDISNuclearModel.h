//____________________________________________________________________________
/*!

\class    genie::BYDISNuclearModel

\brief    Uses DIS Nuclear Model Correction from J.Gomez et. al.
          Reference: journals.aps.org/prd/pdf/10.1103/PhysRevD.49.4348

\author   J. Tena Vidal <julia.tena-vidal@ific.uv.es>
          Universitat de Valencia

\created  February, 2026

\cpright  Copyright (c) 2003-2027, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _JG_DIS_NUCLEAR_MODEL_H_
#define _JG_DIS_NUCLEAR_MODEL_H_

#include <string>

#include "Physics/DeepInelastic/NuclearModel/DISNuclearModelI.h"

namespace genie {

class JGDISNuclearModel : public DISNuclearModelI {

public:
  JGDISNuclearModel();
  JGDISNuclearModel(string config);
  virtual ~JGDISNuclearModel() {};

  double DISACorrection (const Interaction * interaction) const ;
  void Configure (const Registry & config);
  void Configure (string config);

private:
  void LoadConfig (void);

  double fAc0; //> JG-NuclModel-c0
  double fAc1; //> JG-NuclModel-c1
  double fAc2; //> JG-NuclModel-c2

  double fAa0; //> JG-NuclModel-a0
  double fAa1; //> JG-NuclModel-a1
  double fAa2; //> JG-NuclModel-a2
  double fAa3; //> JG-NuclModel-a3
  double fAa4; //> JG-NuclModel-a4
  double fAa5; //> JG-NuclModel-a5
  double fAa6; //> JG-NuclModel-a6
  double fAa7; //> JG-NuclModel-a7
  double fAa8; //> JG-NuclModel-a8
};

}         // genie namespace
#endif    // _JG_DIS_NUCLEAR_MODEL_H_
