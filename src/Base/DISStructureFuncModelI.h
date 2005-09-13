//____________________________________________________________________________
/*!

\class    genie::DISStructureFuncModelI

\brief    Pure Abstract Base Class. Defines the DISStructureFuncModelI
          interface to be implemented by any algorithmic class computing DIS
          structure functions.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004

*/
//____________________________________________________________________________

#ifndef _DIS_STRUCTURE_FUNCTIONS_MODEL_I_H_
#define _DIS_STRUCTURE_FUNCTIONS_MODEL_I_H_

#include "Algorithm/Algorithm.h"
#include "Interaction/Interaction.h"

namespace genie {

class DISStructureFuncModelI : public Algorithm {

public:

  virtual ~DISStructureFuncModelI();

  //-- define DISStructureFuncModelI interface

  virtual void Calculate (const Interaction * interaction) const = 0;

  virtual double F1 (void) const = 0;
  virtual double F2 (void) const = 0;
  virtual double F3 (void) const = 0;
  virtual double F4 (void) const = 0;
  virtual double F5 (void) const = 0;
  virtual double F6 (void) const = 0;

protected:

  DISStructureFuncModelI();
  DISStructureFuncModelI(const char * param_set);
};

}         // genie namespace

#endif    // _DIS_STRUCTURE_FUNCTIONS_MODEL_I_H_
