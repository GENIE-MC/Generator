//____________________________________________________________________________
/*!

\class    genie::BY21DISNuclearModel

\brief    Uses DIS Nuclear Model Correction from the Updated Bodek-Yang Model
          Reference: arxiv.org/pdf/2108.09240 (Section 9)

\author   J. Tena Vidal <julia.tena-vidal@ific.uv.es>
          Universitat de Valencia

\created  February, 2026

\cpright  Copyright (c) 2003-2027, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _BY_DIS_NUCLEAR_MODEL_21_H_
#define _BY_DIS_NUCLEAR_MODEL_21_H_

#include <string>

#include "Physics/DeepInelastic/NuclearModel/DISNuclearModelI.h"

namespace genie {

class BY21DISNuclearModel : public DISNuclearModelI {

public:
  BY21DISNuclearModel();
  BY21DISNuclearModel(string config);
  virtual ~BY21DISNuclearModel() {};

  double DISACorrection (const Interaction * interaction) const ;
  void Configure (const Registry & config);
  void Configure (string config);

private:
  void LoadConfig (void);

  double f2HScale ;  //> BY21-NuclModel-2H-Scale
  double f2Hf0 ;     //> BY21-NuclModel-2H-f0
  double f2Hf1 ;     //> BY21-NuclModel-2H-f1
  double f2Hf2 ;     //> BY21-NuclModel-2H-f2
  double f2Hf3 ;     //> BY21-NuclModel-2H-f3
  double f2Hf4 ;     //> BY21-NuclModel-2H-f4
  double f2Hf5 ;     //> BY21-NuclModel-2H-f5

  double fFef0 ;     //> BY21-NuclModel-Fe-f0
  double fFef1 ;     //> BY21-NuclModel-Fe-f1
  double fFef2 ;     //> BY21-NuclModel-Fe-f2
  double fFef3 ;     //> BY21-NuclModel-Fe-f3
  double fFef4 ;     //> BY21-NuclModel-Fe-f4

  double fCf0 ;     //> BY21-NuclModel-C-f0
  double fCf1 ;     //> BY21-NuclModel-C-f1
  double fCf2 ;     //> BY21-NuclModel-C-f2
  double fCf3 ;     //> BY21-NuclModel-C-f3
  double fCf4 ;     //> BY21-NuclModel-C-f4
  double fCf5 ;     //> BY21-NuclModel-C-f5

  double fAuf0 ;     //> BY21-NuclModel-Au-f0
  double fAuf1 ;     //> BY21-NuclModel-Au-f1
  double fAuf2 ;     //> BY21-NuclModel-Au-f2
  double fAuf3 ;     //> BY21-NuclModel-Au-f3
  double fAuf4 ;     //> BY21-NuclModel-Au-f4
  double fAuf5 ;     //> BY21-NuclModel-Au-f5
  double fAuf6 ;     //> BY21-NuclModel-Au-f6

};

}         // genie namespace
#endif    // _BY_DIS_NUCLEAR_MODEL_21_H_
