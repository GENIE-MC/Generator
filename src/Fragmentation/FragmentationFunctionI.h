//____________________________________________________________________________
/*!

\class    genie::FragmentationFunctionI

\brief    Pure abstract base class.
          Defines the FragmentationFunctionI interface to be implemented by
          any algorithmic class implementing a fragmentation function.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  June 15, 2004

*/
//____________________________________________________________________________

#ifndef _FRAGMENTATION_FUNCTION_I_H_
#define _FRAGMENTATION_FUNCTION_I_H_

#include "Algorithm/Algorithm.h"

namespace genie {

class FragmentationFunctionI : public Algorithm {

public:

  virtual ~FragmentationFunctionI();

  //-- define FragmentationFunctionI interface

  virtual double Value     (double z) const = 0;
  virtual double GenerateZ (void)     const = 0;

protected:

  FragmentationFunctionI();
  FragmentationFunctionI(const char * param_set);
};

}      // genie namespace

#endif // _FRAGMENTATION_FUNCTION_I_H_
