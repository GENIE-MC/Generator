//____________________________________________________________________________
/*!

\class    genie::DISXSec

\brief    Computes the DIS Cross Section.

          Is a concrete implementation of the XSecAlgorithmI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 04, 2004

*/
//____________________________________________________________________________

#include <TMath.h>

#include "AlgFactory/AlgFactory.h"
#include "Conventions/Constants.h"
#include "CrossSections/DISXSec.h"
#include "Messenger/Messenger.h"
#include "Numerical/UnifGrid.h"
#include "Numerical/FunctionMap.h"
#include "Utils/MathUtils.h"
#include "Utils/KineLimits.h"
#include "Utils/Range1.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
DISXSec::DISXSec() :
XSecAlgorithmI()
{
  fName = "genie::DISXSec";
}
//____________________________________________________________________________
DISXSec::DISXSec(const char * param_set) :
XSecAlgorithmI(param_set)
{
  fName = "genie::DISXSec";

  FindConfig();
}
//____________________________________________________________________________
DISXSec::~DISXSec()
{

}
//____________________________________________________________________________
double DISXSec::XSec(const Interaction * interaction) const
{
  //-- Get the requested d^2xsec/dxdy xsec algorithm to use

  const Algorithm * xsec_alg_base = this->SubAlg(
                           "partial-xsec-alg-name", "partial-xsec-param-set");
  const XSecAlgorithmI * partial_xsec_alg =
                         dynamic_cast<const XSecAlgorithmI *> (xsec_alg_base);
                         
  LOG("DISXSec", pDEBUG) << *partial_xsec_alg;

  //-- Get neutrino energy in the struck nucleon rest frame

  const InitialState & init_state = interaction -> GetInitialState();

  TLorentzVector * p4 = init_state.GetProbeP4(kRfStruckNucAtRest);
  
  double Ev  = p4->Energy();

  delete p4;

  //-- Check the energy threshold

  double Ethr = kine_limits::EnergyThreshold(interaction);
  
  if(Ev <= Ethr) {   
     LOG("DISXSec", pINFO) << "E = " << Ev << " < Ethreshold = " << Ethr;
     return 0.;
  }
  //-- Get x,y from config (if they exist) or set defaults

  int    nlogx = this -> NLogX ();
  int    nlogy = this -> NLogY ();
  double xmin  = this -> Xmin  ();
  double xmax  = this -> Xmax  ();
  double ymin  = this -> Ymin  ();
  double ymax  = this -> Ymax  ();

  //-- Check that x,y range is meaningful

  assert( nlogx > 1 && nlogy > 1 );
  
  assert( xmax > xmin && xmax < 1 && xmin < 1 && xmax > 0 & xmin > 0 );
  assert( ymax > ymin && ymax < 1 && ymin < 1 && ymax > 0 & ymin > 0 );
                
  //-- Define the integration area
                
  double log_xmax = TMath::Log(xmax);
  double log_xmin = TMath::Log(xmin);
  double log_ymax = TMath::Log(ymax);
  double log_ymin = TMath::Log(ymin);

  double dlogx = (log_xmax - log_xmin) / (nlogx-1);
  double dlogy = (log_ymax - log_ymin) / (nlogy-1);

  //-- Define the integration grid & instantiate a FunctionMap

  UnifGrid grid;

  grid.AddDimension(nlogx, log_xmin, log_xmax);
  grid.AddDimension(nlogy, log_ymin, log_ymax);

  FunctionMap xyd2xsec(grid);

  //----- Loop over x,y & compute the differential xsec

  LOG("DISXSec", pDEBUG) 
         << "integration limits: x = (" << xmin << ", " << xmax << ") "
                                    << "y = (" << ymin << ", " << ymax << ")";

  for(int ix = 0; ix < nlogx; ix++) {
    double x  = TMath::Exp(log_xmin + ix * dlogx);

    for(int iy = 0; iy < nlogy; iy++) {
       double y = TMath::Exp(log_ymin + iy * dlogy);
                                        
       //-- update the scattering parameters
       interaction->GetScatParamsPtr()->Set("x", x);
       interaction->GetScatParamsPtr()->Set("y", y);

       double pxsec = 0.;

       if ( this->IsWithinIntegrationRange(interaction) ) {
         
          //-- compute d^2xsec/dxdy
          pxsec = partial_xsec_alg->XSec(interaction);

          LOG("DISXSec", pDEBUG)
              << "dxsec/dxdy (x = " << x << ", y = " << y
                                        << ", Ev = " << Ev << ") = " << pxsec;
       }                
       //-- push x*y*(d^2xsec/dxdy) to the FunctionMap       
       xyd2xsec.AddPoint(x*y*pxsec, ix, iy);
              
    } //y
  } //x

  //----- Perform the numerical integration

  LOG("DISXSec", pDEBUG) << "Performing numerical integration";

  //-- get specified integration algorithm from the config. registry
  //   or use Simpson2D if none is defined
  
  const IntegratorI * integrator = this->Integrator();

  double xsec = integrator->Integrate(xyd2xsec);

  LOG("DISXSec", pINFO)  << "xsec_dis (E = " << Ev << " GeV) = " << xsec;

  return xsec;
}
//____________________________________________________________________________
const IntegratorI * DISXSec::Integrator(void) const
{
// Returns the specified (in the config. registry) integration algorithm 
// If none is specified it returns a Simpson2D Integration algorithm

  string integrator_name;

  if( fConfig->Exists("integrator") )
                          fConfig->Get("integrator", integrator_name );
  else integrator_name = "genie::Simpson2D";

  //----- Get the requested algorithm from the algorithm factory

  AlgFactory * algf = AlgFactory::Instance();

  const Algorithm * alg_base = algf->GetAlgorithm(integrator_name);

  const IntegratorI * integrator =
                              dynamic_cast<const IntegratorI *> (alg_base);

  assert(integrator);

  return integrator;
}
//____________________________________________________________________________
bool DISXSec::IsWithinIntegrationRange(const Interaction * interaction) const
{
  //-- Allows the user to set integration limits on W, Q2. For each x,y
  //   pair, the corresponding W,Q2 pair is calculated. The physical
  //   W, Q2 range is calculated and the user-cuts taken into account for
  //   narrowing it down.

  // get physical integration range for W and Q2

  Range1D_t rW  = kine_limits::WRange     (interaction);
  Range1D_t rQ2 = kine_limits::Q2Range_xy (interaction);

  LOG("DISXSec", pDEBUG)
       << "\n Physical integration range: "
             << "W = [" << rW.min << ", " << rW.max << "] GeV "
                      << "Q2 = [" << rQ2.min << ", " << rQ2.max << "] GeV^2";

  // check whether the user wants to override these values

  double Wmin  = (fConfig->Exists("Wmin"))  ? fConfig->GetDouble("Wmin")  : -1;
  double Wmax  = (fConfig->Exists("Wmax"))  ? fConfig->GetDouble("Wmax")  : 1e9;
  double Q2min = (fConfig->Exists("Q2min")) ? fConfig->GetDouble("Q2min") : -1;
  double Q2max = (fConfig->Exists("Q2max")) ? fConfig->GetDouble("Q2max") : 1e9;

  // apply cuts

  kine_limits::ApplyCutsToKineLimits(rW,  Wmin,  Wmax );
  kine_limits::ApplyCutsToKineLimits(rQ2, Q2min, Q2max);

  LOG("DISXSec", pDEBUG)
       << "\n Physical && User integration range: "
             << "W = [" << rW.min << ", " << rW.max << "] GeV "
                      << "Q2 = [" << rQ2.min << ", " << rQ2.max << "] GeV^2";
  // current W, Q2

  double x  = interaction->GetScatteringParams().x();
  double y  = interaction->GetScatteringParams().y();

  const InitialState & init_state = interaction -> GetInitialState();

  TLorentzVector * p4 = init_state.GetProbeP4(kRfStruckNucAtRest);

  double Ev  = p4->Energy();

  delete p4;

  double M  = init_state.GetTarget().StruckNucleonMass();
  double M2 = M*M;
  
  double currW2 = TMath::Max(0., M2 + 2*Ev*M*y*(1-x));
  double currW  = TMath::Sqrt(currW2);
  double currQ2 = TMath::Max(0., 2*M*Ev*x*y);
  
  bool in_range = math_utils::IsWithinLimits(currQ2, rQ2)
                                     && math_utils::IsWithinLimits(currW, rW);

  if(!in_range) {
      LOG("DISXSec", pDEBUG) << "*** excluding from phase space: ";
  } else {
      LOG("DISXSec", pDEBUG) << "*** including in phase space: ";
  }
  LOG("DISXSec", pDEBUG) 
            << "x = " << x << ", y = " << y << ", Q2 = " << currQ2
                                      << ", W = " << currW << ", Ev = " << Ev;

  return in_range;
}
//____________________________________________________________________________
int DISXSec::NLogX(void) const
{
  return ( fConfig->Exists("n-log-x") ) ? fConfig->GetInt("n-log-x") : 61;
}
//____________________________________________________________________________
double DISXSec::Xmin(void) const
{
  return ( fConfig->Exists("x-min")   ) ? fConfig->GetDouble("x-min") : 0.001;
}  
//____________________________________________________________________________
double DISXSec::Xmax(void) const
{
  return ( fConfig->Exists("x-max")   ) ? fConfig->GetDouble("x-max") : 0.999;
}
//____________________________________________________________________________
int DISXSec::NLogY(void) const
{
  return ( fConfig->Exists("n-log-y") ) ? fConfig->GetInt("n-log-y")  : 61;
}
//____________________________________________________________________________
double DISXSec::Ymin(void) const
{  
  return ( fConfig->Exists("y-min")   ) ? fConfig->GetDouble("y-min") : 0.001;
}
//____________________________________________________________________________
double DISXSec::Ymax(void) const
{
  return ( fConfig->Exists("y-max")   ) ? fConfig->GetDouble("y-max") : 0.999;
}
//____________________________________________________________________________
