//____________________________________________________________________________
/*!

\class    genie::GFunc

\brief    GENIE function ABC

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  July 01, 2004

\cpright  Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _GENIE_FUNCTION_H_
#define _GENIE_FUNCTION_H_

#include <vector>
#include <string>

using std::vector;
using std::string;

#include "Utils/Range1.h"

namespace genie {

class Range1D_t;

class GFunc
{
public:
  virtual ~GFunc();

  virtual void         SetParam    (unsigned int i, string name, double min, double max);
  virtual void         SetParam    (unsigned int i, string name, Range1D_t & l);
  virtual unsigned int NParams     (void) const { return fNDim; }
  virtual Range1D_t    ParamLimits (unsigned int i) const;
  virtual string       ParamName   (unsigned int i) const;

protected:
  GFunc();
  GFunc(unsigned int nds);
  GFunc(const GFunc & func);

  virtual void Init    (unsigned int nds);
  virtual void CleanUp (void);

  unsigned int     fNDim;
  vector<double> * fVMin;
  vector<double> * fVMax;
  vector<string> * fVName;
};

}        // namespace

#endif   // _GENIE_FUNCTION_H_
