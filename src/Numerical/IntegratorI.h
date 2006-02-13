//____________________________________________________________________________
/*!

\class    genie::IntegratorI

\brief

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  July 01, 2004

*/
//____________________________________________________________________________

#ifndef _INTEGRATOR_I_H_
#define _INTEGRATOR_I_H_

#include "Algorithm/Algorithm.h"
#include "Numerical/FunctionMap.h"

namespace genie {

class IntegratorI : public Algorithm
{
public:

  virtual double Integrate(FunctionMap & func_map) const = 0;
  virtual double EvalError(FunctionMap & func_map) const = 0;

protected:

  IntegratorI();
  IntegratorI(string name);
  IntegratorI(string name, string config);
  virtual ~IntegratorI();
};

}        // namespace
#endif   // _INTEGRATOR_I_H_
