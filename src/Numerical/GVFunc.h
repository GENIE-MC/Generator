//____________________________________________________________________________
/*!

\class    genie::GVFunc

\brief    GENIE vector function ABC

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  July 01, 2004

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
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

