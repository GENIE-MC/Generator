//____________________________________________________________________________
/*!

\class    genie::Interpolator2D

\brief    A 2D interpolator using the GSL spline type
          If GSL version is not sufficient, does an inefficient version using TGraph2D.

\author   Steve Dennis <s.r.dennis \at liverpool.ac.uk>
          University of Liverpool

\created  November, 2017

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <cstdlib>

#ifndef GENIE_INTERPOLATOR2D_H_
#define GENIE_INTERPOLATOR2D_H_

namespace genie {

class Interpolator2D
{
  public:
    Interpolator2D(const size_t & size_x, const double * grid_x,
                   const size_t & size_y, const double * grid_y,
                   const double * knots);
    ~Interpolator2D();
    
    double Eval    (const double & x, const double & y) const;
    double DerivX  (const double & x, const double & y) const;
    double DerivY  (const double & x, const double & y) const;
    double DerivXX (const double & x, const double & y) const;
    double DerivXY (const double & x, const double & y) const;
    double DerivYY (const double & x, const double & y) const;
    
  private:
    // Done using PIMPL to avoid GSL vs ROOT mess in libraries
    // Struct type declarations will be done in object code
    struct spline2d_container    ; // stores type gsl_spline2d
    struct interp_accel_container; // stores type gsl_interp_accel
    // And these are our actual members
    spline2d_container             * fSpline;
    mutable interp_accel_container * fAcc_x;
    mutable interp_accel_container * fAcc_y;
};

} // namespace genie

#endif //GENIE_INTERPOLATOR2D_H_
