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

#include "Physics/NuclearState/DISNuclearModelI.h"

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
};

}         // genie namespace
#endif    // _BY_DIS_NUCLEAR_MODEL_H_
