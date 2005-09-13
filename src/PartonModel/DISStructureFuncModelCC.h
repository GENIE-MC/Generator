//____________________________________________________________________________
/*!

\class    genie::DISStructureFuncModelCC

\brief    Form Factors for neutrino - free nucleon DIS CC interactions.
          Is a concrete implementation of the DISStructureFuncModelI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004

*/
//____________________________________________________________________________

#ifndef _DIS_STRUCTURE_FUNC_MODEL_CC_H_
#define _DIS_STRUCTURE_FUNC_MODEL_CC_H_

#include "PartonModel/DISStructureFuncModel.h"

namespace genie {

class DISStructureFuncModelCC : public DISStructureFuncModel {

public:

  DISStructureFuncModelCC();
  DISStructureFuncModelCC(const char * param_set);
  virtual ~DISStructureFuncModelCC();

  //-- DISStructureFuncModelI interface implementation

  // override just this interface method and take any other implementation
  // from DISStructureFuncModel
  void Calculate(const Interaction * interaction) const;
};

}         // genie namespace

#endif    // _DIS_STRUCTURE_FUNC_MODEL_CC_H_

