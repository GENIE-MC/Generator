//____________________________________________________________________________
/*!

\class    genie::GFunc

\brief    GENIE function ABC

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  July 01, 2004

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
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
