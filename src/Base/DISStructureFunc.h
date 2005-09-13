//____________________________________________________________________________
/*!

\class    genie::DISStructureFunc

\brief    A class holding Deep Inelastic Scattering (DIS) Form Factors
          (invariant structure funstions)

          This class is using the \b Strategy Pattern. \n
          It can accept requests to calculate itself, for a given interaction,
          that it then delegates to the algorithmic object, implementing the
          DISStructureFuncModelI interface, that it finds attached to itself.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 05, 2004

*/
//____________________________________________________________________________

#ifndef _DIS_STRUCTURE_FUNCTIONS_H_
#define _DIS_STRUCTURE_FUNCTIONS_H_

#include <iostream>

#include "Base/DISStructureFuncModelI.h"
#include "Interaction/Interaction.h"

using std::ostream;

namespace genie {

class DISStructureFunc {

public:

  DISStructureFunc();
  DISStructureFunc(const DISStructureFunc & form_factors);
  virtual ~DISStructureFunc() { }

  void   SetModel  (const DISStructureFuncModelI * model);
  void   Calculate (const Interaction * interaction);

  double xF1 (void) const { return fxF1; }
  double F2  (void) const { return fF2;  }
  double xF3 (void) const { return fxF3; }
  double F4  (void) const { return fF4;  }
  double xF5 (void) const { return fxF5; }
  double F6  (void) const { return fF6;  }

  void Print(ostream & stream) const;

  friend ostream & operator << (ostream & stream, const DISStructureFunc & sf);

private:

  void   InitFormFactors(void);

  double fxF1;
  double fF2;
  double fxF3;
  double fF4;
  double fxF5;
  double fF6;

  const DISStructureFuncModelI * fModel;
};

}       // genie namespace

#endif  // _DIS_STRUCTURE_FUNCTIONS_H_
