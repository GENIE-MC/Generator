//____________________________________________________________________________
/*!

\namespace  genie::math_utils

\brief      Simple mathematical utilities not found in ROOT's TMath
          
\author     Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
            CCLRC, Rutherford Appleton Laboratory

\created    May 06, 2004

*/ 
//____________________________________________________________________________

#include <float.h>

#include <TMath.h>

#include "Utils/MathUtils.h"

//____________________________________________________________________________
bool genie::math_utils::AreEqual(double x1, double x2)
{
  double err = 10. * DBL_EPSILON;

  double dx  = TMath::Abs(x1-x2);

  return (dx < err);
}
//____________________________________________________________________________
bool genie::math_utils::AreEqual(float x1, float x2)
{
  float err = 10. * FLT_EPSILON;

  float dx  = TMath::Abs(x1-x2);

  return (dx < err);
}
//____________________________________________________________________________
bool genie::math_utils::IsWithinLimits(double x, Range1D_t range)
{
  return ( x >= range.min && x <= range.max );
}
//____________________________________________________________________________
bool genie::math_utils::IsWithinLimits(float x, Range1F_t range)
{
  return ( x >= range.min && x <= range.max );
}
//____________________________________________________________________________
bool genie::math_utils::IsWithinLimits(int i, Range1I_t range)
{
  return ( i >= range.min && i <= range.max );
}
//____________________________________________________________________________
double genie::math_utils::NonNegative(double x)
{
// this is used to handle very small numbers in sqrts 

  return TMath::Max(0., x);
}
//____________________________________________________________________________
double genie::math_utils::NonNegative(float x)
{
// this is used to handle very small numbers in sqrts 

  return TMath::Max( (float)0., x);
}
//____________________________________________________________________________

