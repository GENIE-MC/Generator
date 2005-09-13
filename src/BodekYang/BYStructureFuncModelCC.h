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

#ifndef _BODEK_YANG_FORM_FACTORS_CC_H_
#define _BODEK_YANG_FORM_FACTORS_CC_H_

#include "BodekYang/BYStructureFuncModel.h"

namespace genie {

class BYStructureFuncModelCC : public BYStructureFuncModel {

public:

  BYStructureFuncModelCC();
  BYStructureFuncModelCC(const char * param_set);
  ~BYStructureFuncModelCC();

  //-- DISFormFactorsModelI interface implementation

  double xF1 (const Interaction * interaction) const;
  double F2  (const Interaction * interaction) const;
  double xF3 (const Interaction * interaction) const;
  double F4  (const Interaction * interaction) const;
  double xF5 (const Interaction * interaction) const;
  double F6  (const Interaction * interaction) const;
};

}         // genie namespace

#endif    // _BODEK_YANG_FORM_FACTORS_CC_H_
