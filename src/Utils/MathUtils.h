//____________________________________________________________________________
/*!

\namespace  genie::utils::math

\brief      Simple mathematical utilities not found in ROOT's TMath

\author     Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
            University of Liverpool & STFC Rutherford Appleton Lab

\created    May 06, 2004

\cpright    Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
            For the full text of the license visit http://copyright.genie-mc.org
            or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _MATH_UTILS_H_
#define _MATH_UTILS_H_

#include <vector>

#include "Utils/Range1.h"

using std::vector;

namespace genie {
namespace utils {

namespace math
{
  double KahanSummation (double x[], unsigned int n);
  double KahanSummation (const vector<double> & x);
  bool   AreEqual       (double x1, double x2);
  bool   AreEqual       (float  x1, float  x2);
  bool   IsWithinLimits (double x, Range1D_t range);
  bool   IsWithinLimits (float  x, Range1F_t range);
  bool   IsWithinLimits (int    i, Range1I_t range);
  double NonNegative    (double x);
  double NonNegative    (float  x);

} // math  namespace
} // utils namespace
} // genie namespace

#endif // _MATH_UTILS_H_
