//____________________________________________________________________________
/*!

\program testNumerical

\brief   

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created May 4, 2004
*/
//____________________________________________________________________________

#include "Algorithm/AlgFactory.h"
#include "Numerical/GSFunc.h"
#include "Numerical/Spline.h"
#include "Numerical/IntegratorI.h"
#include "Messenger/Messenger.h"
#include "Utils/Range1.h"

using namespace genie;

class GExampleFunc : public GSFunc
{
public:
  GExampleFunc() : GSFunc(1) {}
  ~GExampleFunc() {}

  double operator() (const vector<double> & p)
  {
    double x = p[0];
    double y = 3*x*x+2*x+7;
    return y;
  }
};

int main(int argc, char ** argv)
{
  //****** NUMERICAL INTEGRATION

  // get integrator
  AlgFactory * algf = AlgFactory::Instance();
  const IntegratorI * integrator = 
        dynamic_cast<const IntegratorI *> (
           algf->GetAlgorithm("genie::Simpson1D","Default-Linear"));

  // create function
  GExampleFunc myfunction;

  // set function limits
  Range1D_t xlimits(0.,10.);
  myfunction.SetParam(0,"x",xlimits);

  // integrate function
  double I = integrator->Integrate(myfunction);

  LOG("Main",pINFO) 
              << "Integral <for given numerical accuracy> = " << I;

  //****** INTERPOLATION

  //evaluate the example function at various points

  const int N  = 40;
  double x[N], y[N];
  double dx = (xlimits.max - xlimits.min)/(N-1);
  vector<double> p(1);

  for(int i=0; i<N; i++) {

    p[0] = xlimits.min + i*dx;
    x[i] = p[0];
    y[i] = myfunction(p);
  }

  // create an 1-D spline
  Spline * spline = new Spline(N,x,y);

  // evaluate the original function and the spline at various points


  int M = 100; // != N, otherwise we test only at the spline knots
  dx = (xlimits.max - xlimits.min)/(M-1);

  for(int i=0; i<M; i++) {

    p[0] = xlimits.min + i*dx;
    double yfunc = myfunction(p);
    double yspln = spline->Evaluate(p[0]);

    LOG("Main",pINFO) << "x = " << p[0] 
               << ", f(x) = " << yfunc << ", spline(x) = " << yspln;
  }
  delete spline;

  return 0;
}
