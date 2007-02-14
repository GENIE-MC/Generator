//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - May 03, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Numerical/SimMCIntegrator.h"
#include "Numerical/GSFunc.h"
#include "Numerical/RandomGen.h"
#include "Messenger/Messenger.h"

using namespace genie;

//____________________________________________________________________________
SimMCIntegrator::SimMCIntegrator():
IntegratorI("genie::SimMCIntegrator")
{

}
//____________________________________________________________________________
SimMCIntegrator::SimMCIntegrator(string config) :
IntegratorI("genie::SimMCIntegrator", config)
{

}
//____________________________________________________________________________
SimMCIntegrator::~SimMCIntegrator()
{

}
//____________________________________________________________________________
double SimMCIntegrator::Integrate(GSFunc & gsfunc) const
{
  unsigned int ndim = gsfunc.NParams();

  vector<double> x    (ndim); // input param vector for scalar function
  vector<double> dx   (ndim); // input param (max-min) ranges
  vector<double> xmin (ndim); // input paran min
  double y        = 0;        // estimated scalar function output
  double sum_y2   = 0;        // sum{y^2} over the sample
  double sum_y    = 0;        // sum{y} over the sample
  double sum_y_2  = 0;        // sum{y}^2 over the sample
  double av_y2    = 0;        // <y^2> over the sample
  double av_y     = 0;        // <y> over the sample
  double av_y_2   = 0;        // <y>^2 over the sample
  double sum      = 0;        // computed integral - current step
  double err      = 99999;    // evaluated numerical error

  // Get the random number generators
  RandomGen * rnd = RandomGen::Instance();

  // Init
  for(unsigned int id=0; id<ndim; id++) {
    Range1D_t l = gsfunc.ParamLimits(id);
    dx[id]      = l.max-l.min;
    xmin[id]    = l.min;
  }

  // Compute "hyper"volume
  double vol = 1.;
  for(unsigned int id=0; id<ndim; id++) {
    vol *= dx[id];
  }
  LOG("SimMCIntegrator", pNOTICE) << "Hypervolume = " << vol;
  if(vol==0) return 0;

  // Increase the number of function evalations until the computed
  // integral value converges to the real one within the required accuracy
  unsigned int ntot = 0;

  while(err>fMaxPcntErr) {
    // Generate an input parameter vector and evaluate the function
    for(unsigned int id=0; id<ndim; id++) {
          x[id] = xmin[id] + dx[id] * rnd->RndNum().Rndm();
    }
    y = gsfunc(x);

    ntot++;
    sum_y  += y;
    sum_y2 += (y*y);

    // do this every fNMin evaluations
    if(ntot>0 && ntot%fNMin==0) {

      sum_y_2  = sum_y*sum_y;
      av_y2    = sum_y2/ (double)ntot;
      av_y     = sum_y / (double)ntot;
      av_y_2   = av_y*av_y;

      // compute the integral
      sum = av_y*vol;
      if(sum==0) return 0;

      // evaluate the error
      err = 100*vol * TMath::Sqrt((av_y2-av_y_2)/(double)ntot) / sum;

      LOG("SimMCIntegrator", pINFO)
        << "Integral = " << sum << " / Estimated err = " << err << " %";

      if(ntot>fNMax) {
        LOG("SimMCIntegrator", pERROR)
            << "Maximum numerical error allowed = " << fMaxPcntErr << " %";
        LOG("SimMCIntegrator", pFATAL)
              << "Integral didn't converge to required numerical accuracy";
        LOG("SimMCIntegrator", pFATAL)
            << "Estimated Error = " << err
                  << " % - Aborting @ " << ntot << " function evaluations";
        abort();
      }
      if(err < fMaxPcntErr) {
        LOG("SimMCIntegrator", pNOTICE)
          << "Integral = " << sum << " / Estimated err = " << err << " %";
        return sum;
      }
    }
  }
  return 0;
}
//____________________________________________________________________________
void SimMCIntegrator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void SimMCIntegrator::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfigData();
}
//____________________________________________________________________________
void SimMCIntegrator::LoadConfigData(void)
{
  fNMax       = (unsigned int) fConfig->GetInt("NMaxAllowedEval");
  fNMin       = (unsigned int) fConfig->GetInt("NEvalPerStep");
  fMaxPcntErr = fConfig->GetDouble("MaxErr");
}
//____________________________________________________________________________

