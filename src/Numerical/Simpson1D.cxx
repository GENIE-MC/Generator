//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - May 03, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Conventions/GBuild.h"
#include "Numerical/Simpson1D.h"
#include "Numerical/GSFunc.h"
#include "Numerical/FunctionMap.h"
#include "Numerical/UnifGrid.h"
#include "Numerical/UnifGridDimension.h"
#include "Messenger/Messenger.h"

using namespace genie;

//____________________________________________________________________________
Simpson1D::Simpson1D():
IntegratorI("genie::Simpson1D")
{

}
//____________________________________________________________________________
Simpson1D::Simpson1D(string config) :
IntegratorI("genie::Simpson1D", config)
{

}
//____________________________________________________________________________
Simpson1D::~Simpson1D()
{

}
//____________________________________________________________________________
double Simpson1D::Integrate(GSFunc & gsfunc) const
{
  unsigned int ndim = gsfunc.NParams();
  assert(ndim==1); // Simpson1D requires an 1-D function

  UnifGrid init_grid(gsfunc, fSpacing); // initial grid for the function
  FunctionMap fmap(init_grid);          // a function map for this grid

  vector<double> x(ndim);    // input param vector for scalar function
  double y        = 0;       // scalar function output
  double sum      = 0;       // computed integral - current step
  double sum_old  = 9999999; // computed integral - previous step
  double err      = 0;       // evaluated numerical error
  unsigned int n  = fNo;     // param controling num of integration steps
  unsigned int np = 0;       // number of integration steps

  // Increase the number of integration steps (2**N+1) until the computed
  // integral value converges to the real one within the required accuracy
  for(unsigned int iter=0; iter<fIMaxConv; iter++) {

    np = (unsigned int) TMath::Power(2,(int)n) + 1;
    n++;

    fmap.IncreaseGridDensity(np);
    const UnifGrid & curr_grid = fmap.GetGrid();

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    LOG("Simpson1D", pINFO)
      << "Integration: iter = " << iter << ", using grid: " << curr_grid;
#endif

    // populate the function map with the values of the input function
    // computed on the grid points
    for(unsigned int i=0; i<np; i++) {

       //-- fill input vector and evaluate the scalar function
       x[0] = curr_grid(0, i); // x

       if(!fmap.ValueIsSet(x)) {
          y = gsfunc(x); // f(x)

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
          LOG("Simpson1D", pDEBUG)
               << "grid point...." << i << "/" << np << " : "
                                << "func(x = " << x[0] << " ) = " << y;
#endif
          //-- note that if the grid points are distributed uniformly in
          //   ln(x) then the scalar function has to be multiplied by x:
          //   integral { f(x)dx } = integral { x*f(x) dln(x) }
          if(fSpacing==kGSpLoge) y *= x[0];

          //-- add the computed point at the function map
          fmap.SetValue(y, x);
       } else {
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
           LOG("Simpson1D", pDEBUG)
              << "grid point...." << i << "/" << np << " : "
              << "func at x = " << x[0] << " computed at previous step";
#endif
       }
    }

    // compute the integral using the Simpson rule and evaluate the error
    sum = this->SimpsonRule(fmap);
    if(sum+sum_old == 0) return 0;
    err = 200*TMath::Abs( (sum-sum_old)/(sum+sum_old) ); // in %

    LOG("Simpson1D", pINFO)
       << "Integral = " << sum << " (prev = " << sum_old
                                << ") / Estimated err = " << err << " %";

    if(err < fMaxPcntErr) {
       LOG("Simpson1D", pNOTICE)
          << "Integral = " << sum << " / Estimated err = " << err << " %";
       return sum;
    } else {
      sum_old = sum;
    }
  }
  LOG("Simpson1D", pERROR)
            << "Maximum numerical error allowed = " << fMaxPcntErr << " %";
  LOG("Simpson1D", pFATAL)
              << "Integral didn't converge to required numerical accuracy";
  LOG("Simpson1D", pFATAL)
            << "Estimated Error = " << err
                       << " % - Aborting @ " << np << " integration steps";
  abort();
}
//____________________________________________________________________________
double Simpson1D::SimpsonRule(FunctionMap & func_map) const
{
  const UnifGrid & grid = func_map.GetGrid();
  const UnifGridDimension & gd = *grid[0];

  unsigned int N    = gd.NPoints();
  double       step = gd.Step();

  LOG("Simpson1D", pDEBUG) << "N-Points = " << N << ", Step-Size = " << step;

  vector<unsigned int> pos(1);

  double sum = 0;

  pos[0] = 0;
  sum += (0.5*func_map.Value(pos));
  pos[0] = N-1;
  sum += (0.5*func_map.Value(pos));
  for(unsigned int i = 1; i< N-1; i++)  {
       pos[0] = i;
       sum    += (func_map.Value(pos) * (i%2 + 1));
  }
  sum *= (2.*step/3.);

  return sum;
}
//____________________________________________________________________________
void Simpson1D::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void Simpson1D::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfigData();
}
//____________________________________________________________________________
void Simpson1D::LoadConfigData(void)
{
  fIMaxConv   = (unsigned int) fConfig->GetInt("MaxNIter");
  fNo         = (unsigned int) fConfig->GetInt("InitNStep");
  fMaxPcntErr = fConfig->GetDouble("MaxErr");

  bool inloge = fConfig->GetBool("InLoge");
  if(inloge) fSpacing = kGSpLoge;
  else       fSpacing = kGSpLinear;
}
//____________________________________________________________________________

