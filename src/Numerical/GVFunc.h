//____________________________________________________________________________
/*!

\class    genie::GVFunc

\brief    GENIE vector function ABC

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  July 01, 2004

\cpright  Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _GENIE_VECTOR_FUNCTION_H_
#define _GENIE_VECTOR_FUNCTION_H_

#include <vector>

#include "Numerical/GFunc.h"

using std::vector;

namespace genie {

class GVFunc : public GFunc
{
public:
  virtual ~GVFunc();
  virtual const vector<double> & operator () (const vector<double> & x) = 0;

protected:
  GVFunc(unsigned int nds);

  vector<double> fOutV;
};

}        // genie namespace
#endif   // _GENIE_VECTOR_FUNCTION_H_

