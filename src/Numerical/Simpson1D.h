//____________________________________________________________________________
/*!

\class    genie::Simpson1D

\brief

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
};

}        // namespace

#endif   // _SIMPSON_1D_H_
