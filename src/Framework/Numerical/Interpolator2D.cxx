//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Steve Dennis <s.r.dennis \at liverpool.ac.uk>
*/
//____________________________________________________________________________

#include <cassert>
#include <cstdlib>

#include "Framework/Numerical/Interpolator2D.h"
#include "gsl/gsl_version.h"

using namespace genie;

// Implemented primarily with the GSL v2 spline2d class.
// However, not all sites have this, so there's a backup implementation using TGraph2D
#if GSL_MAJOR_VERSION >= 2
#include "gsl/gsl_spline2d.h"

// define our Pimpl structs
struct Interpolator2D::spline2d_container
{
  spline2d_container() : spl(NULL) {};
  ~spline2d_container() { if (spl) delete spl; };
  gsl_spline2d * spl;
};
//____________________________________________________________________________
struct Interpolator2D::interp_accel_container
{
  interp_accel_container() : acc(NULL) {};
  ~interp_accel_container() { if (acc) delete acc; };
  gsl_interp_accel * acc;
};
//____________________________________________________________________________

// Now define our actual code (GSL)
//____________________________________________________________________________
Interpolator2D::Interpolator2D(
  const size_t & size_x, const double * grid_x,
  const size_t & size_y, const double * grid_y,
  const double * knots) :
  fSpline (new Interpolator2D::spline2d_container()     ),
  fAcc_x  (new Interpolator2D::interp_accel_container() ),
  fAcc_y  (new Interpolator2D::interp_accel_container() )
{
  fSpline->spl = gsl_spline2d_alloc(gsl_interp2d_bilinear,size_x,size_y);
  gsl_spline2d_init(fSpline->spl,grid_x,grid_y,knots,size_x,size_y);
  fAcc_x->acc = gsl_interp_accel_alloc();
  fAcc_y->acc = gsl_interp_accel_alloc();
}
//____________________________________________________________________________
Interpolator2D::~Interpolator2D()
{
  if (fSpline) delete fSpline;
  if (fAcc_x ) delete fAcc_x;
  if (fAcc_y ) delete fAcc_y;
}
//____________________________________________________________________________
double Interpolator2D::Eval(const double & x, const double & y) const
{
  return gsl_spline2d_eval(
    fSpline->spl, x, y,
    fAcc_x->acc,
    fAcc_y->acc);
}
//____________________________________________________________________________
double Interpolator2D::DerivX(const double & x, const double & y) const
{
  return gsl_spline2d_eval_deriv_x(
    fSpline->spl, x, y,
    fAcc_x->acc,
    fAcc_y->acc);
}
//____________________________________________________________________________
double Interpolator2D::DerivY(const double & x, const double & y) const
{
  return gsl_spline2d_eval_deriv_y(
    fSpline->spl, x, y,
    fAcc_x->acc,
    fAcc_y->acc);
}
//____________________________________________________________________________
double Interpolator2D::DerivXX(const double & x, const double & y) const
{
  return gsl_spline2d_eval_deriv_xx(
    fSpline->spl, x, y,
    fAcc_x->acc,
    fAcc_y->acc);
}
//____________________________________________________________________________
double Interpolator2D::DerivXY(const double & x, const double & y) const
{
  return gsl_spline2d_eval_deriv_xy(
    fSpline->spl, x, y,
    fAcc_x->acc,
    fAcc_y->acc);
}
//____________________________________________________________________________
double Interpolator2D::DerivYY(const double & x, const double & y) const
{
  return gsl_spline2d_eval_deriv_yy(
    fSpline->spl, x, y,
    fAcc_x->acc,
    fAcc_y->acc);
}
//____________________________________________________________________________
//____________________________________________________________________________
//____________________________________________________________________________
// And now the TGraph2D version
#else
#include "TGraph2D.h"
// define our Pimpl structs for TGraph
struct Interpolator2D::spline2d_container
{
  spline2d_container() : spl(NULL) {};
  ~spline2d_container() { if (spl) delete spl; };
  TGraph2D * spl;
};
//____________________________________________________________________________
struct Interpolator2D::interp_accel_container
{
  interp_accel_container() {};
  ~interp_accel_container() { };
};
//____________________________________________________________________________
// Now define our actual code (TGraph2D)
//____________________________________________________________________________
Interpolator2D::Interpolator2D(
  const size_t & size_x, const double * grid_x,
  const size_t & size_y, const double * grid_y,
  const double * knots) :
  fSpline (new spline2d_container()),
  fAcc_x  (NULL),
  fAcc_y  (NULL)
{
  fSpline->spl = new TGraph2D(size_x*size_y);
  fSpline->spl->SetDirectory(0);
  size_t iz = 0;
  for (size_t iy = 0 ; iy < size_y ; iy++) {
    for (size_t ix = 0 ; ix < size_x ; ix++) {
      fSpline->spl->SetPoint(iz,grid_x[ix],grid_y[iy],knots[iz]);
      iz++;
    }
  }
}
//____________________________________________________________________________
Interpolator2D::~Interpolator2D()
{
  if (fSpline) delete fSpline;
  if (fAcc_x ) delete fAcc_x;
  if (fAcc_y ) delete fAcc_y;
}
//____________________________________________________________________________
double Interpolator2D::Eval(const double & x, const double & y) const
{
  return fSpline->spl->Interpolate(x,y);
}
//____________________________________________________________________________
double Interpolator2D::DerivX(const double & x, const double & y) const
{
  assert(!"Method requires GSL version 2 or higher.");
  return -999;
}
//____________________________________________________________________________
double Interpolator2D::DerivY(const double & x, const double & y) const
{
  assert(!"Method requires GSL version 2 or higher.");
  return -999;
}
//____________________________________________________________________________
double Interpolator2D::DerivXX(const double & x, const double & y) const
{
  assert(!"Method requires GSL version 2 or higher.");
  return -999;
}
//____________________________________________________________________________
double Interpolator2D::DerivXY(const double & x, const double & y) const
{
  assert(!"Method requires GSL version 2 or higher.");
  return -999;
}
//____________________________________________________________________________
double Interpolator2D::DerivYY(const double & x, const double & y) const
{
  assert(!"Method requires GSL version 2 or higher.");
  return -999;
}
//____________________________________________________________________________

#endif // GSL_MAJOR_VERSION
