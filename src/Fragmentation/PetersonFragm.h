//____________________________________________________________________________
/*!

\class    genie::PetersonFragm

\brief    The Peterson fragmentation function.

          Is a concrete implementation of the FragmentationFunctionI interface.

\ref      C.Peterson et al., Phys.Rev.D23, 56 (1981)

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  June 15, 2004

*/
//____________________________________________________________________________

#ifndef _PETERSON_FRAGM_H_
#define _PETERSON_FRAGM_H_

#include <TF1.h>

#include "Fragmentation/FragmentationFunctionI.h"

namespace genie {

class PetersonFragm : public FragmentationFunctionI {

public:

  PetersonFragm();
  PetersonFragm(string config);
  ~PetersonFragm();

  //-- implements the FragmentationFunctionI interface

  double Value     (double z) const;
  double GenerateZ (void)     const;

private:

  void BuildFunction (void);

  TF1 * fFunc;
};

}      // genie namespace

#endif // _PETERSON_FRAGM_H_
