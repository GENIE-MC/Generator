//____________________________________________________________________________
/*!

\class    genie::RESPXSec

\brief    Computes single differential RES cross section.

          It can be configured to compute either
            \li dxsec / dQ2  where \c Q2 is the momentum transfer, or
            \li dxsec / dW   where \c W is the hadronic invariant mass

          This is a model-independent algorithm. It merely integrates the
          specified double differential cross section.

          Is a concrete implementation of the XSecAlgorithmI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  March 03, 2005

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Algorithm/AlgFactory.h"
#include "Conventions/Constants.h"
#include "CrossSections/RESPXSec.h"
#include "Messenger/Messenger.h"
#include "Numerical/UnifGrid.h"
#include "Numerical/FunctionMap.h"
#include "Numerical/IntegratorI.h"
#include "Utils/KineUtils.h"
#include "Utils/Range1.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
RESPXSec::RESPXSec() :
XSecAlgorithmI()
{
  fName = "genie::RESPXSec";
}
//____________________________________________________________________________
RESPXSec::RESPXSec(const char * param_set) :
XSecAlgorithmI(param_set)
{
  fName = "genie::RESPXSec";

  FindConfig();
}
//____________________________________________________________________________
RESPXSec::~RESPXSec()
{

}
//____________________________________________________________________________
double RESPXSec::XSec(const Interaction * interaction) const
{
  //-- Make sure it knows what kind of partial (dxsec/d?) xsec algorithm it is

  fConfig->AssertExistence("is-differential-over");

  string variable = fConfig->GetString("is-differential-over");

  LOG("RESPXSec", pINFO) << "Computing dxsec_{RES}/d" << variable;

  //-- Get the specified double differential (d^2xsec/dQ2dW) xsec algorithm

  const Algorithm * xsec_alg_base = this->SubAlg(
                           "partial-xsec-alg-name", "partial-xsec-param-set");
  const XSecAlgorithmI * d2xsec_alg =
                         dynamic_cast<const XSecAlgorithmI *> (xsec_alg_base);

  //-- Get t (W,logQ2) integration range from config (if it exists
  //   or set default values). It should be OK if default range extends to
  //   unphysical values - the double differential xsec alg should return
  //   0. for these phase space points.
  //   Note: Note that for RES, W values > ~ 1.4 GeV are not unphysical.
  //   Set the external config value if you want to impose a upper W limit.

  int nt = 0; double tmin = 0, tmax = 0, dt = 0, t = 0;

  if ( variable.find("W") != string::npos ) {

     //-- it is dxsec/dW, get intergation limits over Q2

     // default is physical Q2 range for input W
     Range1D_t rQ2 = utils::kinematics::Q2Range_W(interaction);

     if(rQ2.min<0 && rQ2.max<0) return 0.;

     nt   = ( fConfig->Exists("nLogQ2") ) ? fConfig->GetInt("nLogQ2")    : 151;
     tmin = ( fConfig->Exists("Q2min")  ) ? fConfig->GetDouble("Q2min")  : rQ2.min;
     tmax = ( fConfig->Exists("Q2max")  ) ? fConfig->GetDouble("Q2max")  : rQ2.max;

  } else if ( variable.find("Q2") != string::npos ) {

     //-- it is dxsec/dQ2, get intergation limits over W

     // default is physical W range for the given energy
     Range1D_t rW = utils::kinematics::WRange(interaction);

     nt   = ( fConfig->Exists("nW") )    ? fConfig->GetInt("nW")       : 51;
     tmin = ( fConfig->Exists("Wmin")  ) ? fConfig->GetDouble("Wmin")  : rW.min;
     tmax = ( fConfig->Exists("Wmax")  ) ? fConfig->GetDouble("Wmax")  : rW.max;

  } else return 0.;

  //-- Check that t (W or Q2) range is meaningful

  LOG("RESPXSec", pINFO) << "tmin = " << tmin << ", tmax = " << tmax << ", nt = " << nt;

  assert( nt > 1 && tmax > tmin);

  //-- Define the integration grid & instantiate a FunctionMap

  UnifGrid grid;

  if ( variable == "W" ) {

    // It is a dxsec/dW -- do the Q2 integration over *dlogQ2*

    dt = (TMath::Log(tmax) - TMath::Log(tmin)) / (nt-1);
    grid.AddDimension(nt, TMath::Log(tmin), TMath::Log(tmax)); // 1-D

  } else if ( variable == "Q2" ) {

    // It is a dxsec/dQ2 -- do the W integration

    dt = (tmax - tmin) / (nt-1);
    grid.AddDimension(nt, tmin, tmax); // 1-D
  }

  FunctionMap funcmap(grid); // Q2*(d^2xsec/dlogQ2) or d^2xsec/dW

  //-- Loop over t (W or Q2) and compute dxsec/dt

  for(int it = 0; it < nt; it++) {

    if ( variable == "W" ) {

       t = TMath::Exp( TMath::Log(tmin) + it * dt);
       interaction->GetScatParamsPtr()->Set("Q2", t);
       double d2xsec = d2xsec_alg->XSec(interaction);
       funcmap.AddPoint(t*d2xsec, it);

    } else if ( variable == "Q2" ) {

       t = tmin + it * dt;
       interaction->GetScatParamsPtr()->Set("W", t);
       double d2xsec = d2xsec_alg->XSec(interaction);
       funcmap.AddPoint(d2xsec, it);
    }
  } //t

  //-- Numerical integration

  const IntegratorI * integrator = this->Integrator();

  double xsec = integrator->Integrate(funcmap);

  return xsec;
}
//____________________________________________________________________________
const IntegratorI * RESPXSec::Integrator(void) const
{
  //-- get specified integration algorithm from the config. registry
  //   or use Simpson1D if no one else is defined

  string integrator_name;

  if( fConfig->Exists("integrator") )
                          fConfig->Get("integrator", integrator_name );
  else integrator_name = "genie::Simpson1D";

  //-- Get an instance of the AlgFactory & request the integrator algorithm

  AlgFactory * algf = AlgFactory::Instance();

  const Algorithm * intg_alg_base = algf->GetAlgorithm(integrator_name);

  const IntegratorI * integrator =
                           dynamic_cast<const IntegratorI *> (intg_alg_base);

  assert(integrator);

  return integrator;
}
//____________________________________________________________________________

