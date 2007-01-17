//____________________________________________________________________________
/*!

\program testNumerical

\brief   

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created May 4, 2004

\cpright Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
         All rights reserved.
         For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#include "Algorithm/AlgFactory.h"
#include "Numerical/GSFunc.h"
#include "Numerical/Spline.h"
#include "Numerical/IntegratorI.h"
#include "Messenger/Messenger.h"
#include "Utils/Range1.h"

using namespace genie;

class GExampleFunc1D : public GSFunc
{
public:
  GExampleFunc1D() : GSFunc(1) {}
  ~GExampleFunc1D() {}

  double operator() (const vector<double> & p)
  {
    double x = p[0];
    double f = 3*x*x+2*x+7; // f(x)
    return f;
  }
};

class GExampleFunc2D : public GSFunc
{
public:
  GExampleFunc2D() : GSFunc(2) {}
  ~GExampleFunc2D() {}

  double operator() (const vector<double> & p)
  {
    double x2 = TMath::Power(p[0],2);
    double y2 = TMath::Power(p[1],2);
    double f  = TMath::Exp(-x2+9)*TMath::Exp(-y2+16); // f(x,y)
    return f;
  }
};

int main(int /*argc*/, char ** /*argv*/)
{
  //****** NUMERICAL INTEGRATION

  // get integrator
  AlgFactory * algf = AlgFactory::Instance();
  const IntegratorI * integrator1 = 
        dynamic_cast<const IntegratorI *> (
           algf->GetAlgorithm("genie::Simpson1D","Default-Linear"));
  const IntegratorI * integrator2 = 
        dynamic_cast<const IntegratorI *> (
           algf->GetAlgorithm("genie::Simpson2D","Default-Linear"));

  // create test functions & set limits

  // 1-D:
  GExampleFunc1D myfunc1d;

  Range1D_t lim(0,10);
  myfunc1d.SetParam(0,"x", lim);

  // 2-D:
  GExampleFunc2D myfunc2d;

  myfunc2d.SetParam(0,"x",-12.,12.);
  myfunc2d.SetParam(1,"y",-12.,12.);

  // integrate functions
  double I1 = integrator1->Integrate(myfunc1d);
  double I2 = integrator2->Integrate(myfunc2d);

  LOG("Main",pINFO) 
         << "1-D Integral <for given numerical accuracy> = " << I1;
  LOG("Main",pINFO) 
         << "2-D Integral <for given numerical accuracy> = " << I2;

  //****** INTERPOLATION

  //evaluate the example function at various points

  const int N  = 40;
  double x[N], y[N];
  double dx = (lim.max - lim.min)/(N-1);
  vector<double> p(1);

  for(int i=0; i<N; i++) {
    p[0] = lim.min + i*dx;
    x[i] = p[0];
    y[i] = myfunc1d(p);
  }

  // create an 1-D spline
  Spline * spline = new Spline(N,x,y);

  // evaluate the original function and the spline at various points

  int M = 100; // != N, otherwise we test only at the spline knots
  dx = (lim.max - lim.min)/(M-1);

  for(int i=0; i<M; i++) {

    p[0] = lim.min + i*dx;
    double yfunc = myfunc1d(p);
    double yspln = spline->Evaluate(p[0]);

    LOG("Main",pINFO) << "x = " << p[0] 
               << ", f(x) = " << yfunc << ", spline(x) = " << yspln;
  }
  delete spline;

  return 0;
}
