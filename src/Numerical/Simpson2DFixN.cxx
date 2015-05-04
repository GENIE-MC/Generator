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

#include "Numerical/Simpson2DFixN.h"
#include "Numerical/GSFunc.h"
#include "Messenger/Messenger.h"

using namespace genie;

//____________________________________________________________________________
Simpson2DFixN::Simpson2DFixN():
IntegratorI("genie::Simpson2DFixN")
{

}
//____________________________________________________________________________
Simpson2DFixN::Simpson2DFixN(string config) :
IntegratorI("genie::Simpson2DFixN", config)
{

}
//____________________________________________________________________________
Simpson2DFixN::~Simpson2DFixN()
{

}
//____________________________________________________________________________
double Simpson2DFixN::Integrate(GSFunc & gsfunc) const
{
  unsigned int ndim = gsfunc.NParams();
  assert(ndim==2); // Simpson2DFixN requires an 2-D function

  vector<double> r(ndim);    
  vector<double> sum1d(fNBinsD0);

  Range1D_t xrange = gsfunc.ParamLimits(0);
  Range1D_t yrange = gsfunc.ParamLimits(1);

  double xmin = xrange.min;
  double xmax = xrange.max;
  double ymin = yrange.min;
  double ymax = yrange.max;

  if(fSpacingD0 == kGSpLoge) {
    xmin = TMath::Log(xmin);
    xmax = TMath::Log(xmax);
  }
  if(fSpacingD1 == kGSpLoge) {
    ymin = TMath::Log(ymin);
    ymax = TMath::Log(ymax);
  }

  double dx = (xmax - xmin)/(fNBinsD0-1);
  double dy = (ymax - ymin)/(fNBinsD1-1);

  for(unsigned int i=0; i<fNBinsD0; i++) {

    sum1d[i] = 0;
    r[0] = (fSpacingD0 == kGSpLoge) ? 
                     TMath::Exp(xmin+i*dx) : xmin + i*dx;

    r[1] = (fSpacingD1 == kGSpLoge) ? TMath::Exp(ymin) : ymin;
    sum1d[i] += 0.5 * (this->J(r) * gsfunc(r));
    r[1] = (fSpacingD1 == kGSpLoge) ? TMath::Exp(ymax) : ymax;
    sum1d[i] += 0.5 * (this->J(r) * gsfunc(r));
    
    for(unsigned int j=1; j<fNBinsD1-1; j++) {
        r[1] = (fSpacingD1 == kGSpLoge) ? 
                     TMath::Exp(ymin+j*dy) : ymin + j*dy;

        sum1d[i] += (this->J(r) * gsfunc(r) * (j%2 + 1));
    }//j
    sum1d[i] *= (2.*dy/3.);
  }

  double sum2d = (sum1d[0]+sum1d[fNBinsD0-1])/2.;
  for(unsigned int i=1; i<fNBinsD0-1; i++) {
    sum2d += (sum1d[i] * (i%2 + 1));
  }
  sum2d *= (2.*dx/3.);

  return sum2d;
}
//____________________________________________________________________________
double Simpson2DFixN::J(const vector<double> & r) const
{
  double J = 1.;

  if(fSpacingD0 == kGSpLoge) J *= r[0];
  if(fSpacingD1 == kGSpLoge) J *= r[1];

  return J;
}
//____________________________________________________________________________
void Simpson2DFixN::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void Simpson2DFixN::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfigData();
}
//____________________________________________________________________________
void Simpson2DFixN::LoadConfigData(void)
{
  bool inloge0 = fConfig->GetBoolDef("InLogeDim0", false);
  bool inloge1 = fConfig->GetBoolDef("InLogeDim1", false);

  if(inloge0) fSpacingD0 = kGSpLoge;
  else        fSpacingD0 = kGSpLinear;

  if(inloge1) fSpacingD1 = kGSpLoge;
  else        fSpacingD1 = kGSpLinear;

  fNBinsD0 = fConfig->GetIntDef("NBinsDim0", 401);
  fNBinsD1 = fConfig->GetIntDef("NBinsDim1", 401);
}
//____________________________________________________________________________

