//____________________________________________________________________________
/*!

\class    genie::Simpson2D

\brief

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004

*/
//____________________________________________________________________________

#ifndef _SIMPSON_2D_H_
#define _SIMPSON_2D_H_

#include "Numerical/IntegratorI.h"

namespace genie {

class Simpson2D: public IntegratorI
{
public:

  Simpson2D();
  virtual ~Simpson2D();

  double Integrate(FunctionMap & func_map) const;
};

}        // namespace

#endif   // _SIMPSON_2D_H_

