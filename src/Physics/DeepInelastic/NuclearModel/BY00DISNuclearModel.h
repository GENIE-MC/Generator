//____________________________________________________________________________
/*!

\class    genie::BY00DISNuclearModel

\brief    Uses DIS Nuclear Model Correction from Bodek-Yang

\author   J. Tena Vidal <julia.tena-vidal@ific.uv.es>
          Universitat de Valencia

\created  February, 2026

\cpright  Copyright (c) 2003-2027, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _BY00_DIS_NUCLEAR_MODEL_H_
#define _BY00_DIS_NUCLEAR_MODEL_H_

#include <string>

#include "Physics/DeepInelastic/NuclearModel/DISNuclearModelI.h"

namespace genie {

class BY00DISNuclearModel : public DISNuclearModelI {

public:
  BY00DISNuclearModel();
  BY00DISNuclearModel(string config);
  virtual ~BY00DISNuclearModel() {};

  double DISACorrection (const Interaction * interaction) const ;
  void Configure (const Registry & config);
  void Configure (string config);

private:
  void LoadConfig (void);

  double f2HScale ;  //> BY00-NuclModel-2H-Scale
  double f2Hf0 ;     //> BY00-NuclModel-2H-f0
  double f2Hf1 ;     //> BY00-NuclModel-2H-f1
  double f2Hf2 ;     //> BY00-NuclModel-2H-f2
  double f2Hf3 ;     //> BY00-NuclModel-2H-f3
  double f2Hf4 ;     //> BY00-NuclModel-2H-f4
  double f2Hf5 ;     //> BY00-NuclModel-2H-f5

  double fFef0 ;     //> BY00-NuclModel-Fe-f0
  double fFef1 ;     //> BY00-NuclModel-Fe-f1
  double fFef2 ;     //> BY00-NuclModel-Fe-f2
  double fFef3 ;     //> BY00-NuclModel-Fe-f3
  double fFef4 ;     //> BY00-NuclModel-Fe-f4

};

}         // genie namespace
#endif    // _BY_DIS_NUCLEAR_MODEL_H_
