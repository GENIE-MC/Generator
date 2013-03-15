//____________________________________________________________________________
/*!

\class    genie::GBFunc

\brief    GENIE boolean function ABC

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  July 01, 2004

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _GENIE_BOOL_FUNCTION_H_
#define _GENIE_BOOL_FUNCTION_H_

#include <vector>

#include "Numerical/GFunc.h"

using std::vector;

namespace genie {

class GBFunc : public GFunc
{
public:
  virtual ~GBFunc();
  virtual bool operator () (const vector<double> & x) = 0;

protected:
  GBFunc(unsigned int nds);
};

}        // genie namespace
#endif   // _GENIE_BOOL_FUNCTION_H_
