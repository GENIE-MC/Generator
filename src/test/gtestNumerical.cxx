//____________________________________________________________________________
/*!

\program gtestNumerical

\brief   Program testing GENIE's numerical utility classes

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created May 4, 2004

\cpright Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <TMath.h>

#include "Algorithm/AlgFactory.h"
#include "Numerical/GSFunc.h"
#include "Numerical/Spline.h"
#include "Numerical/IntegratorI.h"
#include "Messenger/Messenger.h"
#include "Utils/Range1.h"

using namespace genie;

//____________________________________________________________________________
// example 1-D GSFunc
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
//____________________________________________________________________________
// example 2-D GSFunc
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
//____________________________________________________________________________

void testNumericalIntegration (void);
void testSplineInterpolation  (void);
void testSplineOperators      (void);

int main(int /*argc*/, char ** /*argv*/)
{
  testNumericalIntegration();
  //testSplineInterpolation();
  //testSplineOperators();

  return 0;
}
//____________________________________________________________________________
void testNumericalIntegration(void)
{
// test numerial integration

  // get integrator
  AlgFactory * algf = AlgFactory::Instance();
  const IntegratorI * integrator1a = 
        dynamic_cast<const IntegratorI *> (
           algf->GetAlgorithm("genie::Simpson1D","Default-Linear"));

  const IntegratorI * integrator2a = 
        dynamic_cast<const IntegratorI *> (
           algf->GetAlgorithm("genie::Simpson2D","Default-Linear"));
  const IntegratorI * integrator2b = 
        dynamic_cast<const IntegratorI *> (
           algf->GetAlgorithm("genie::Simpson2D","Default-Logarithmic"));
  const IntegratorI * integrator2c = 
        dynamic_cast<const IntegratorI *> (
           algf->GetAlgorithm("genie::Simpson2D","Default-Linear-FixedNBins-201"));
  const IntegratorI * integrator2d = 
        dynamic_cast<const IntegratorI *> (
           algf->GetAlgorithm("genie::Simpson2D","Default-Logarithmic-FixedNBins-201"));
  const IntegratorI * integrator2e = 
        dynamic_cast<const IntegratorI *> (
           algf->GetAlgorithm("genie::Simpson2DFixN","LogLog-201"));

  // create test functions & set limits

  // 1-D:
  GExampleFunc1D myfunc1d;
  Range1D_t lim(1,10);
  myfunc1d.SetParam(0,"x", lim);

  // 2-D:
  GExampleFunc2D myfunc2d;
  myfunc2d.SetParam(0,"x",0.001,20.);
  myfunc2d.SetParam(1,"y",1.,10.);

  // integrate functions
  double I1a = integrator1a->Integrate(myfunc1d);
  double I2a = integrator2a->Integrate(myfunc2d);
  double I2b = integrator2b->Integrate(myfunc2d);
  double I2c = integrator2c->Integrate(myfunc2d);
  double I2d = integrator2d->Integrate(myfunc2d);
  double I2e = integrator2e->Integrate(myfunc2d);

  LOG("test",pINFO) 
         << "1-D Integral <for given numerical accuracy> = " << I1a;
  LOG("test",pINFO) 
         << "2-D Integral <for given numerical accuracy> a= " << I2a;
  LOG("test",pINFO) 
         << "2-D Integral <for given numerical accuracy> b= " << I2b;
  LOG("test",pINFO) 
         << "2-D Integral <for given numerical accuracy> c= " << I2c;
  LOG("test",pINFO) 
         << "2-D Integral <for given numerical accuracy> d= " << I2d;
  LOG("test",pINFO) 
         << "2-D Integral <for given numerical accuracy> e= " << I2e;
}
//____________________________________________________________________________
void testSplineInterpolation(void)
{
// test spline interpolation

  //test func
  GExampleFunc1D myfunc1d;
  Range1D_t lim(0,10);
  myfunc1d.SetParam(0,"x", lim);

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
    LOG("test",pINFO) << "x = " << p[0] 
               << ", f(x) = " << yfunc << ", spline(x) = " << yspln;
  }
  delete spline;
}
//____________________________________________________________________________
void testSplineOperators(void)
{
  const int N = 10;
  double x1[N], y1[N];
  double x2[N], y2[N];

  double xmin = 0.;
  double xmax = 1.;
  double dx   = (xmax-xmin)/(N-1);

  for(int i=0; i<N; i++) {
    x1[i] = xmin + i*dx;
    x2[i] = xmin + i*dx;
    y1[i] = TMath::Sin(x1[i]) + 0.2;
    y2[i] = TMath::Power(TMath::Cos(3*x2[i]), 2) + 0.5;
  }

  Spline spl1(N,x1,y1);
  spl1.SetName("test spline #1");

  Spline spl2(N,x1,y1);
  spl2.SetName("test spline #2");

  LOG("test",pINFO) << "Printing spline_1";
  LOG("test",pINFO) << spl1;

  LOG("test",pINFO) << "Printing spline_2";
  LOG("test",pINFO) << spl2;

  LOG("test",pINFO) << "...adding 3 to spline_1";
  spl1.Add(3.);
  LOG("test",pINFO) << spl1;

  LOG("test",pINFO) << "...multiplying spline_1 with 3";
  spl1.Multiply(3.);
  LOG("test",pINFO) << spl1;

  LOG("test",pINFO) << "...multiplying spline_1 with 2*spline_2";
  spl1.Multiply(spl2, 2.);
  LOG("test",pINFO) << spl1;
}
//____________________________________________________________________________


