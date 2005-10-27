//____________________________________________________________________________
/*!

\class    genie::DISPXSec

\brief    Computes single differential DIS cross section.

          It can be configured to compute either
            \li dxsec / dx  where \c x is the Bjorken scaling variable, or
            \li dxsec / dy  where \c y is the Inelasticity

          This is a model-independent algorithm. It merely integrates the
          specified double differential cross section.

          Is a concrete implementation of the XSecAlgorithmI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  July 01, 2004

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Algorithm/AlgFactory.h"
#include "Conventions/Constants.h"
#include "CrossSections/DISPXSec.h"
#include "Messenger/Messenger.h"
#include "Numerical/UnifGrid.h"
#include "Numerical/FunctionMap.h"
#include "Numerical/IntegratorI.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
DISPXSec::DISPXSec() :
XSecAlgorithmI("genie::DISPXSec")
{

}
//____________________________________________________________________________
DISPXSec::DISPXSec(string config) :
XSecAlgorithmI("genie::DISPXSec", config)
{

}
//____________________________________________________________________________
DISPXSec::~DISPXSec()
{

}
//____________________________________________________________________________
double DISPXSec::XSec(const Interaction * interaction) const
{
  LOG("DISPXSec", pDEBUG) << *fConfig;

  //-- Make sure it knows what kind of partial (dxsec/d?) xsec algorithm it is

  assert( fConfig->Exists("is-differential-over") );

  string variable = fConfig->GetString("is-differential-over");

  LOG("DISPXSec", pDEBUG) << "XSec is differential over var: " << variable;

  //-- Get a d^2xsec/dxdy xsec algorithm

  fConfig->AssertExistence("partial-xsec-alg-name", "partial-xsec-param-set");

  const Algorithm * xsec_alg_base = this->SubAlg(
                           "partial-xsec-alg-name", "partial-xsec-param-set");
  const XSecAlgorithmI * partial_xsec_alg =
                         dynamic_cast<const XSecAlgorithmI *> (xsec_alg_base);

  //----- Set default & check for user defined integration range

  int    nlogt = 61;     //------- default integration range
  double tmin  = 0.001;  // t = x or y
  double tmax  = 0.999;

  //-- Get x or y integration range from config (if exists)

  if ( variable == "y" ) {

     //-- it is dxsec/dy, get intergation limits over x

     if ( fConfig->Exists("n-log-x") ) fConfig->Get("n-log-x", nlogt );
     if ( fConfig->Exists("x-min")   ) fConfig->Get("x-min",   tmin  );
     if ( fConfig->Exists("x-max")   ) fConfig->Get("x-max",   tmax  );

  } else if ( variable == "x" ) {

     //-- it is dxsec/dx, get intergation limits over y

     if ( fConfig->Exists("n-log-y") ) fConfig->Get("n-log-y", nlogt );
     if ( fConfig->Exists("y-min")   ) fConfig->Get("y-min",   tmin  );
     if ( fConfig->Exists("y-max")   ) fConfig->Get("y-max",   tmax  );
  }

  LOG("DISPXSec", pDEBUG) << "Integration: n(log)bins = "
                     << nlogt << ", min = " << tmin << ", max = " << tmax;

  //-- Check that t (x or y) range is meaningful
  assert( nlogt > 1 );
  assert( tmax > tmin && tmax < 1 && tmin < 1 && tmax > 0 & tmin > 0 );

  //----- Define the integration area & step
  double log_tmax = TMath::Log(tmax);
  double log_tmin = TMath::Log(tmin);
  double dlogt    = (log_tmax - log_tmin) / (nlogt-1);

  //----- Define the integration grid & instantiate a FunctionMap
  UnifGrid grid;
  grid.AddDimension(nlogt, log_tmin, log_tmax); // 1-D

  FunctionMap tdxsec(grid); // t * (dxsec/dt), t = x,y

  //----- Loop over t (x or y) and compute dxsec/dt

  for(int it = 0; it < nlogt; it++) {

    double t  = TMath::Exp(log_tmin + it * dlogt);

    //-- update the "running" scattering parameter
    if      (variable == "y") interaction->GetKinematicsPtr()->Setx(t);
    else if (variable == "x") interaction->GetKinematicsPtr()->Sety(t);

    //-- compute d^2xsec/dxdy
    double pxsec = partial_xsec_alg->XSec(interaction);

    tdxsec.AddPoint(t*pxsec, it); // t * dxsec/dt (t = x or y)

  } //t

  //----- Numerical integration

  //-- get specified integration algorithm from the config. registry
  //   or use Simpson1D if no one else is defined
  string intgr = fConfig->GetStringDef("integrator", "genie::Simpson1D");

  AlgFactory * algf = AlgFactory::Instance();
  const IntegratorI * integrator =
          dynamic_cast<const IntegratorI *> (algf->GetAlgorithm(intgr));

  double xsec = integrator->Integrate(tdxsec);
  return xsec;
}
//____________________________________________________________________________

