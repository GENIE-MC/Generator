//____________________________________________________________________________
/*!

\class    genie::BYDISNuclearModel

\brief    Uses DIS Nuclear Model Correction from Bodek-Yang 
          Reference: arxiv.org/pdf/2108.09240 (Section 9)

\author   J. Tena Vidal <julia.tena-vidal@ific.uv.es>
          Universitat de Valencia

\created  February, 2026

\cpright  Copyright (c) 2003-2027, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _BY_DIS_NUCLEAR_MODEL_I_H_
#define _BY_DIS_NUCLEAR_MODEL_I_H_

#include <string>

#include "Physics/NuclearState/DISNuclearModelI.h"

namespace genie {

class BYDISNuclearModel : public DISNuclearModelI {

public:
  BYDISNuclearModel();
  BYDISNuclearModel(string config);
  virtual ~BYDISNuclearModel() {};

  double DISACorrection (const Interaction & interaction) const ;
  void Configure (const Registry & config);
  void Configure (string config);

private:
  void LoadConfig (void);
};

}         // genie namespace
#endif    // _BY_DIS_NUCLEAR_MODEL_I_H_
