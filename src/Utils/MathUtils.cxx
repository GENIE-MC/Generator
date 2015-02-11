//____________________________________________________________________________
/*
 Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - May 06, 2004

 For documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <float.h>

#include <TMath.h>

#include "Messenger/Messenger.h"
#include "Utils/MathUtils.h"

//____________________________________________________________________________
double genie::utils::math::KahanSummation(double x[], unsigned int n)
{
// the Kahan summation algorithm - minimizes the error when adding a sequence
// of finite precision floating point numbers (compensated summation)

  double sum = x[0];
  double c   = 0.0;
  for(unsigned int i=1; i<n; i++) {
    double y = x[i]-c;
    double t = sum+y;
    c   = (t-sum) - y;
    sum = t;
  }
  return sum;
}
//____________________________________________________________________________
double genie::utils::math::KahanSummation(const vector<double> & x)
{
// the Kahan summation algorithm - minimizes the error when adding a sequence
// of finite precision floating point numbers (compensated summation)

  double sum = x[0];
  double c   = 0.0;
  for(unsigned int i=1; i<x.size(); i++) {
    double y = x[i]-c;
    double t = sum+y;
    c   = (t-sum) - y;
    sum = t;
  }
  return sum;
}
//____________________________________________________________________________
bool genie::utils::math::AreEqual(double x1, double x2)
{
  double err = 0.001*DBL_EPSILON;
  double dx  = TMath::Abs(x1-x2);
  if(dx<err) {
    LOG("Math", pINFO) << x1 << " := " << x2;
    return true;
  }
  return false;;
}
//____________________________________________________________________________
bool genie::utils::math::AreEqual(float x1, float x2)
{
  float err = FLT_EPSILON;
  float dx  = TMath::Abs(x1-x2);
  if(dx<err) {
    LOG("Math", pINFO) << x1 << " := " << x2;
    return true;
  }
  return false;;
}
//____________________________________________________________________________
bool genie::utils::math::IsWithinLimits(double x, Range1D_t range)
{
  return ( x >= range.min && x <= range.max );
}
//____________________________________________________________________________
bool genie::utils::math::IsWithinLimits(float x, Range1F_t range)
{
  return ( x >= range.min && x <= range.max );
}
//____________________________________________________________________________
bool genie::utils::math::IsWithinLimits(int i, Range1I_t range)
{
  return ( i >= range.min && i <= range.max );
}
//____________________________________________________________________________
double genie::utils::math::NonNegative(double x)
{
// this is used to handle very small numbers in sqrts

  return TMath::Max(0., x);
}
//____________________________________________________________________________
double genie::utils::math::NonNegative(float x)
{
// this is used to handle very small numbers in sqrts

  return TMath::Max( (float)0., x);
}
//____________________________________________________________________________

