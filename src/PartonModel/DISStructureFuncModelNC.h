//____________________________________________________________________________
/*!

\class    genie::DISStructureFuncModelNC

\brief    Form Factors for neutrino - free nucleon DIS NC interactions.
          Is a concrete implementation of the DISStructureFuncModelI interface.

\ref      E.A.Paschos and J.Y.Yu, Phys.Rev.D 65.033002

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
  DISStructureFuncModelNC(string config);
  virtual ~DISStructureFuncModelNC();

  //-- DISStructureFuncModelI interface implementation

  // override just this interface method and take any other implementation
  // from DISStructureFuncModel
  void Calculate(const Interaction * interaction) const;

private:

  void LoadConfig();

  double fGL;
  double fGR;
  double fGLprime;
  double fGRprime;
  double fGL2;
  double fGR2;
  double fGLprime2;
  double fGRprime2;
};

}         // genie namespace

#endif    // _DIS_STRUCTURE_FUNC_MODEL_NC_H_

