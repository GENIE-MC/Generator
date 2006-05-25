//____________________________________________________________________________
/*!

\class    genie::BYStructureFuncModelCC

\brief    Computes CC vN DIS Form Factors according to the Bodek-Yang model.
          Inherits part of its implemenation from the DISBodekYangFormFactors
          abstract class.

          Check out DISBodekYangFormFactors for comments and references.

          Is a concrete implementation of the DISFormFactorsModelI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  September 28, 2004

*/
//____________________________________________________________________________

#ifndef _BODEK_YANG_STRUCTURE_FUNCTIONS_MODEL_CC_H_
#define _BODEK_YANG_STRUCTURE_FUNCTIONS_MODEL_CC_H_

#include "BodekYang/BYStructureFuncModel.h"

namespace genie {

class BYStructureFuncModelCC : public BYStructureFuncModel {

public:
  BYStructureFuncModelCC();
  BYStructureFuncModelCC(string config);
  ~BYStructureFuncModelCC();
};

}         // genie namespace

#endif    // _BODEK_YANG_STRUCTURE_FUNCTIONS_MODEL_CC_H_
