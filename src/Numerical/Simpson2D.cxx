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
#include "Numerical/Simpson2D.h"
#include "Numerical/FunctionMap.h"
#include "Numerical/GSFunc.h"
#include "Numerical/UnifGrid.h"
#include "Numerical/UnifGridDimension.h"
#include "Messenger/Messenger.h"

using namespace genie;

//____________________________________________________________________________
Simpson2D::Simpson2D():
IntegratorI("genie::Simpson2D")
{

}
//____________________________________________________________________________
Simpson2D::Simpson2D(string config) :
IntegratorI("genie::Simpson2D", config)
{

}
//____________________________________________________________________________
Simpson2D::~Simpson2D()
{

}
//____________________________________________________________________________
double Simpson2D::Integrate(GSFunc & gsfunc) const
{
  unsigned int ndim = gsfunc.NParams();
  assert(ndim==2); // Simpson2D requires an 2-D function

  UnifGrid init_grid(gsfunc, fSpacing); // a grid for the input function
  FunctionMap fmap(init_grid);          // a function map for this grid

  vector<double> x(ndim);    // input param vector for scalar function
  double y        = 0;       // scalar function output
  double sum      = 0;       // computed integral - current step
  double sum_old  = 9999999; // computed integral - previous step
  double err      = 0;       // evaluated numerical error
  unsigned int n  = fNo;     // param controling num of integration steps
  unsigned int np = 0;       // number of integration steps

  // Check whether the use prefers to use a fixed number of integration  
  // steps. Note that if you do this I won;t give any guarantee for the
  // convergence of the numerical integration

  if(fForceFixedNBins) {  /* ok, enter in cheat mode */

    fmap.IncreaseGridDensity(fNBinsD0, 0);
    fmap.IncreaseGridDensity(fNBinsD1, 1);
    const UnifGrid & curr_grid = fmap.GetGrid();
    // populate the function map with the values of the input function
    // computed on the grid points
    for(unsigned int i=0; i < curr_grid[0]->NPoints(); i++) {
      x[0] = curr_grid(0, i);
      for(unsigned int j=0; j < curr_grid[1]->NPoints(); j++) {
         x[1] = curr_grid(1, j);

         y = gsfunc(x); // f(x)
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
         LOG("Simpson2D", pDEBUG)
               << "grid point...." << i << "," << j
                 << "/" << np << "," << np << " : "
                    << "func(x = " << x[0] << ", " << x[1] << ") = " << y;
#endif
         if(fSpacing==kGSpLoge) y *= (x[0]*x[1]);
         fmap.SetValue(y, x);
      } //x1
    } //x0

    // compute the sum using the Simpson method & return
    sum = this->SimpsonRule(fmap);
    LOG("Simpson2D", pINFO)
        << "Integral = " << sum 
                  << " / Estimated err = *** check disabled by user ***";
    return sum; 

  } // end-cheat-mode


  // Perform integration without fixing the number of integration steps:

  // Increase the number of integration steps (2**N+1) until the computed
  // integral value converges to the real one within the required accuracy
  for(unsigned int iter=0; iter<fIMaxConv; iter++) {

    int idim=-1;
    if(fFastDensityIncrease) {
      // increase the grid density fast - all dimensions simultaneously
      np = (unsigned int) TMath::Power(2,(int)n) + 1;
      n++;
      fmap.IncreaseGridDensity(np);
    } else {
      // increase the grid density slowly - 1 dimension at a time...
      if(iter%ndim==0) {
        np = (unsigned int) TMath::Power(2,(int)n) + 1;
        n++;
        idim = 0;
      }
      fmap.IncreaseGridDensity(np, idim++);
    }

    const UnifGrid & curr_grid = fmap.GetGrid();
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    LOG("Simpson2D", pINFO)
      << "Integration: iter = " << iter << ", using grid: " << curr_grid;
#endif
    // populate the function map with the values of the input function
    // computed on the grid points
    for(unsigned int i=0; i < curr_grid[0]->NPoints(); i++) {
      x[0] = curr_grid(0, i);
      for(unsigned int j=0; j < curr_grid[1]->NPoints(); j++) {
         x[1] = curr_grid(1, j);

         if(!fmap.ValueIsSet(x)) {
            y = gsfunc(x); // f(x)
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
            LOG("Simpson2D", pDEBUG)
               << "grid point...." << i << "," << j
                 << "/" << np << "," << np << " : "
                    << "func(x = " << x[0] << ", " << x[1] << ") = " << y;
#endif
            //-- note that if the grid points are distributed uniformly in
            //   ln(x) then the scalar function has to be multiplied by x:
            //   integral { f(x)dx } = integral { x*f(x) dln(x) }
            if(fSpacing==kGSpLoge) y *= (x[0]*x[1]);

            //-- add the computed point at the function map
            fmap.SetValue(y, x);
         } else {
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
           LOG("Simpson2D", pDEBUG)
               << "grid point...." << i << "," << j
                << "/" << np << "," << np << " : " << "func at (x = " 
                 << x[0] << ", " << x[1] << ") computed at previous step";
#endif
         }
      }
    }

    // compute the integral using the Simpson rule and evaluate the error
    sum = this->SimpsonRule(fmap);
    if(sum+sum_old == 0) return 0;
    err = 200*TMath::Abs( (sum-sum_old)/(sum+sum_old) ); // in %

    LOG("Simpson2D", pINFO)
       << "Integral = " << sum << " (prev = " << sum_old
                                  << ") / Estimated err = " << err << " %";
    if(err < fMaxPcntErr) {
       LOG("Simpson2D", pNOTICE)
           << "Integral = " << sum << " / Estimated err = " << err << " %";
        return sum;
    } else {
      sum_old = sum;
    }
  }
  LOG("Simpson2D", pERROR)
            << "Maximum numerical error allowed = " << fMaxPcntErr << " %";
  LOG("Simpson2D", pFATAL)
              << "Integral didn't converge to required numerical accuracy";
  LOG("Simpson2D", pFATAL)
            << "Estimated Error = " << err
                       << " % - Aborting @ " << np << " integration steps";
  abort();
}
//____________________________________________________________________________
double Simpson2D::SimpsonRule(FunctionMap & func_map) const
{
  const UnifGrid & grid = func_map.GetGrid();

  unsigned int N[2];
  double       step[2];

  for(unsigned int idim=0; idim<2; idim++) {
    N[idim]    = grid[idim]->NPoints();
    step[idim] = grid[idim]->Step();
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    LOG("Simpson2D", pDEBUG) << "DIM: " << idim
            << " -> N = " << N[idim] << ", dx = " << step[idim];
#endif
  }

  vector<unsigned int>  pos(2);
  vector<double> sum1d(N[0]);

  for(unsigned int i=0; i<N[0]; i++) {
    sum1d[i] = 0;

    pos[0] = i;
    pos[1] = 0;      sum1d[i] += (0.5*func_map.Value(pos));
    pos[1] = N[1]-1; sum1d[i] += (0.5*func_map.Value(pos));
    for(unsigned int j = 1; j < N[1]-1; j++)  {
       pos[1] = j; sum1d[i] += (func_map.Value(pos) * (j%2 + 1));
    }
    sum1d[i] *= (2.*step[1]/3.);
  }

  double sum2d = (sum1d[0]+sum1d[N[0]-1])/2.;
  for(unsigned int i=1; i<N[0]-1; i++) {
    sum2d += (sum1d[i] * (i%2 + 1));
  }
  sum2d *= (2.*step[0]/3.);

  return sum2d;
}
//____________________________________________________________________________
void Simpson2D::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void Simpson2D::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfigData();
}
//____________________________________________________________________________
void Simpson2D::LoadConfigData(void)
{
  fIMaxConv   = (unsigned int) fConfig->GetIntDef("MaxNIter", 20);
  fNo         = (unsigned int) fConfig->GetIntDef("InitNStep", 3);
  fMaxPcntErr = fConfig->GetDoubleDef("MaxErr", 0.1); //%

  bool inloge = fConfig->GetBoolDef("InLoge", false);
  if(inloge) fSpacing = kGSpLoge;
  else       fSpacing = kGSpLinear;

  // check the preferred grid density increase rate method
  fFastDensityIncrease = fConfig->GetBoolDef("FastDensityIncrease", false);

  // check whether the user wants to use a fixed number of bins
  // *** notice that if this is used there is no guarantee for the convergence
  // *** of the numerical integration

  fForceFixedNBins = fConfig->GetBoolDef("ForceFixedNBins", false);
  fNBinsD0 = 0;
  fNBinsD1 = 0;
  if(fForceFixedNBins) {
     fNBinsD0 = fConfig->GetIntDef("NBinsDim0", 401);
     fNBinsD1 = fConfig->GetIntDef("NBinsDim1", 401);
  }
}
//____________________________________________________________________________

