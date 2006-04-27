//____________________________________________________________________________
/*!

\class    genie::IntegratorI

\brief    Numerical integration algorithm ABC

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  July 01, 2004

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _INTEGRATOR_I_H_
#define _INTEGRATOR_I_H_

#include "Algorithm/Algorithm.h"

namespace genie {

class GSFunc;

class IntegratorI : public Algorithm
{
public:
  virtual ~IntegratorI();

  virtual double Integrate(GSFunc & gsfunc) const = 0;

protected:
  IntegratorI();
  IntegratorI(string name);
  IntegratorI(string name, string config);
};

}        // namespace
#endif   // _INTEGRATOR_I_H_
