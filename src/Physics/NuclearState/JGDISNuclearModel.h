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

#include "Physics/NuclearState/DISNuclearModelI.h"

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
};

}         // genie namespace
#endif    // _JG_DIS_NUCLEAR_MODEL_H_
