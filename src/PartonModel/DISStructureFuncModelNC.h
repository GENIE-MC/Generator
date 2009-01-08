//____________________________________________________________________________
/*!

\class    genie::DISStructureFuncModelNC

\brief    Form Factors for neutrino - free nucleon DIS NC interactions.
          Is a concrete implementation of the DISStructureFuncModelI interface.

\ref      E.A.Paschos and J.Y.Yu, Phys.Rev.D 65.033002

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  May 03, 2004

\cpright  Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _DIS_STRUCTURE_FUNC_MODEL_NC_H_
#define _DIS_STRUCTURE_FUNC_MODEL_NC_H_

#include "PartonModel/DISStructureFuncModel.h"

namespace genie {

class DISStructureFuncModelNC : public DISStructureFuncModel {

public:
  DISStructureFuncModelNC();
  DISStructureFuncModelNC(string config);
  virtual ~DISStructureFuncModelNC();
};

}         // genie namespace
#endif    // _DIS_STRUCTURE_FUNC_MODEL_NC_H_

