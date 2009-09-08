//____________________________________________________________________________
/*!

\class    genie::DISStructureFuncModelCC

\brief    Structure functions for neutrino DIS CC.
          Is a concrete implementation of the DISStructureFuncModelI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  May 03, 2004

\cpright  Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _DIS_STRUCTURE_FUNC_MODEL_CC_H_
#define _DIS_STRUCTURE_FUNC_MODEL_CC_H_

#include "PartonModel/DISStructureFuncModel.h"

namespace genie {

class DISStructureFuncModelCC : public DISStructureFuncModel {

public:
  DISStructureFuncModelCC();
  DISStructureFuncModelCC(string config);
  virtual ~DISStructureFuncModelCC();
};

}         // genie namespace
#endif    // _DIS_STRUCTURE_FUNC_MODEL_CC_H_

