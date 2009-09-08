//____________________________________________________________________________
/*!

\class    genie::DISStructureFuncModelEM

\brief    Structure functions for charged lepton scattering.
          Is a concrete implementation of the DISStructureFuncModelI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  Sep 09, 2009

\cpright  Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _DIS_STRUCTURE_FUNC_MODEL_EM_H_
#define _DIS_STRUCTURE_FUNC_MODEL_EM_H_

#include "PartonModel/DISStructureFuncModel.h"

namespace genie {

class DISStructureFuncModelEM : public DISStructureFuncModel {

public:
  DISStructureFuncModelEM();
  DISStructureFuncModelEM(string config);
  virtual ~DISStructureFuncModelEM();
};

}         // genie namespace
#endif    // _DIS_STRUCTURE_FUNC_MODEL_EM_H_

