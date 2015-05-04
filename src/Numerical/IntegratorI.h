//____________________________________________________________________________
/*!

\class    genie::IntegratorI

\brief    Numerical integration algorithm ABC

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  July 01, 2004

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
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
