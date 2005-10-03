//____________________________________________________________________________
/*!

\namespace  genie::utils::math

\brief      Simple mathematical utilities not found in ROOT's TMath

\author     Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
            CCLRC, Rutherford Appleton Laboratory

\created    May 06, 2004

*/
//____________________________________________________________________________

#ifndef _MATH_UTILS_H_
#define _MATH_UTILS_H_

#include "Utils/Range1.h"

namespace genie {
namespace utils {

namespace math
{
  bool   AreEqual (double x1, double x2);
  bool   AreEqual (float  x1, float  x2);

  bool   IsWithinLimits (double x, Range1D_t range);
  bool   IsWithinLimits (float  x, Range1F_t range);
  bool   IsWithinLimits (int    i, Range1I_t range);

  double NonNegative (double x);
  double NonNegative (float  x);

} // math  namespace
} // utils namespace
} // genie namespace

#endif // _MATH_UTILS_H_
