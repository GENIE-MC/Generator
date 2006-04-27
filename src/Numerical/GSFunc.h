//____________________________________________________________________________
/*!

\class    genie::GSFunc

\brief    GENIE scalar function ABC

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  July 01, 2004

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _GENIE_SCALAR_FUNCTION_H_
#define _GENIE_SCALAR_FUNCTION_H_

#include <vector>

#include "Numerical/GFunc.h"

using std::vector;

namespace genie {

class GSFunc : public GFunc
{
public:
  virtual ~GSFunc();
  virtual double operator () (const vector<double> & x) = 0;

protected:
  GSFunc(unsigned int nds);
};

}        // genie namespace
#endif   // _GENIE_SCALAR_FUNCTION_H_
