//____________________________________________________________________________
/*!

\class    genie::BYStructureFuncModelNC

\brief    Computes NC vN DIS Form Factors according to the Bodek-Yang model.
          Inherits part of its implemenation from the DISBodekYangFormFactors
          abstract class.

          Check out DISBodekYangFormFactors for comments and references.

          Is a concrete implementation of the DISFormFactorsModelI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  September 28, 2004

*/
//____________________________________________________________________________

#ifndef _BODEK_YANG_STRUCTURE_FUNCTIONS_MODEL_NC_H_
#define _BODEK_YANG_STRUCTURE_FUNCTIONS_MODEL_NC_H_

#include "BodekYang/BYStructureFuncModel.h"

namespace genie {

class BYStructureFuncModelNC : public BYStructureFuncModel {

public:
  BYStructureFuncModelNC();
  BYStructureFuncModelNC(string config);
  ~BYStructureFuncModelNC();
};

}         // genie namespace
#endif    // _BODEK_YANG_STRUCTURE_FUNCTIONS_MODEL_NC_H_
