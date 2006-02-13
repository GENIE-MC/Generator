//____________________________________________________________________________
/*!

\class    genie::Simpson1D

\brief    The 1-dimensional extended Simpson rule (an open integration formula)
          of order 1/N^4.

\ref      Numerical Recipes in C, Cambridge Univ. Press, 2002, page 134

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004

*/
//____________________________________________________________________________

#ifndef _SIMPSON_1D_H_
#define _SIMPSON_1D_H_

#include "Numerical/IntegratorI.h"

namespace genie {

class Simpson1D: public IntegratorI
{
public:

  Simpson1D();
  virtual ~Simpson1D();

  double Integrate(FunctionMap & func_map) const;
  double EvalError(FunctionMap & func_map) const { return 0; }
};

}        // namespace

#endif   // _SIMPSON_1D_H_
