//____________________________________________________________________________
/*!

\class    genie::DISStructureFuncModelNC

\brief    Form Factors for neutrino - free nucleon DIS NC interactions.
          Is a concrete implementation of the DISStructureFuncModelI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004

*/
//____________________________________________________________________________

#ifndef _DIS_STRUCTURE_FUNC_MODEL_NC_H_
#define _DIS_STRUCTURE_FUNC_MODEL_NC_H_

#include "PartonModel/DISStructureFuncModel.h"

namespace genie {

class DISStructureFuncModelNC : public DISStructureFuncModel {

public:

  DISStructureFuncModelNC();
  DISStructureFuncModelNC(const char * param_set);
  virtual ~DISStructureFuncModelNC();

  //-- DISStructureFuncModelI interface implementation

  double xF1 (const Interaction * interaction) const;
  double F2  (const Interaction * interaction) const;
  double xF3 (const Interaction * interaction) const;
  double F4  (const Interaction * interaction) const;
  double xF5 (const Interaction * interaction) const;
  double F6  (const Interaction * interaction) const;
};

}         // genie namespace

#endif    // _DIS_STRUCTURE_FUNC_MODEL_NC_H_

