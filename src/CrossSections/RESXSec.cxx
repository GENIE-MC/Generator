//____________________________________________________________________________
/*!

\class    genie::RESXSec

\brief    Computed the RES Cross Section.

          Is a concrete implementation of the XSecAlgorithmI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 04, 2004

*/
//____________________________________________________________________________

#include "AlgFactory/AlgFactory.h"
#include "Conventions/Constants.h"
#include "CrossSections/RESXSec.h"
#include "Messenger/Messenger.h"
#include "Numerical/UnifGrid.h"
#include "Numerical/FunctionMap.h"
#include "Numerical/IntegratorI.h"
#include "Utils/MathUtils.h"
#include "Utils/KineLimits.h"
#include "Utils/Range1.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
RESXSec::RESXSec() :
XSecAlgorithmI()
{
  fName = "genie::RESXSec";
}
//____________________________________________________________________________
RESXSec::RESXSec(const char * param_set) :
XSecAlgorithmI(param_set)
{
  fName = "genie::RESXSec";

  FindConfig();
}
//____________________________________________________________________________
RESXSec::~RESXSec()
{

}
//____________________________________________________________________________
double RESXSec::XSec(const Interaction * interaction) const
{
  LOG("RESXSec", pINFO) << "HERE";

  //-- Get the requested d^2xsec/dxdy xsec algorithm to use

  const XSecAlgorithmI * partial_xsec_alg = this->PartialXSecAlgorithm();

  //-- Get neutrino energy in the struck nucleon rest frame

  const InitialState & init_state = interaction -> GetInitialState();

  TLorentzVector * p4 = init_state.GetProbeP4(kRfStruckNucAtRest);
  
  double Ev  = p4->Energy();

  delete p4;

  //double s   = M2 + 2*M*Ev;

  //-- Get x,y from config (if they exist) or set defaults

  int    nW      = this -> NW     ();
  int    nlogQ2  = this -> NLogQ2 ();
  double Wmin    = this -> Wmin   ();
  double Wmax    = this -> Wmax   ();
  double Q2min   = this -> Q2min  ();
  double Q2max   = this -> Q2max  ();

  //-- Check that x,y range is meaningful

  assert( nW > 1 && nlogQ2 > 1 );
  
  assert( Wmax  > Wmin  && Wmin  > 0 );
  assert( Q2max > Q2min && Q2min > 0 );
                
  //-- Define the integration area

  double logQ2min = TMath::Log( Q2min );
  double logQ2max = TMath::Log( Q2max );
  
  double dW     = ( Wmax     - Wmin    ) / (nW-1);
  double dlogQ2 = ( logQ2max - logQ2min) / (nlogQ2-1);

  //-- Define the integration grid & instantiate a FunctionMap

  UnifGrid grid;

  grid.AddDimension(nW,     Wmin,     Wmax    );
  grid.AddDimension(nlogQ2, logQ2min, logQ2max);

  FunctionMap Q2d2xsec(grid);

  //----- Loop over x,y & compute the differential xsec

  LOG("RESXSec", pINFO) 
           << "Integrating over  W: (" << Wmin  << ", " << Wmax  << ")";
  LOG("RESXSec", pINFO) 
           << "Integrating over Q2: (" << Q2min << ", " << Q2max << ")";

  for(int i = 0; i < nW; i++) {

    double W = Wmin + i*dW;

    for(int j = 0; j < nlogQ2; j++) {

       double Q2 = TMath::Exp( logQ2min + j * dlogQ2 );

       double pxsec = 0;
                                        
       //-- update the scattering parameters
       interaction->GetScatParamsPtr()->Set("W",  W );
       interaction->GetScatParamsPtr()->Set("Q2", Q2);

       if ( this->IsKinematicallyAllowed(interaction) ) {
         
          //-- compute d^2xsec/dxdy
          pxsec = partial_xsec_alg->XSec(interaction);

          LOG("RESXSec", pINFO)
              << "dxsec[RES]/dQ2dW (Q2 = " << Q2
                      << ", W = " << W << ", Ev = " << Ev << ") = " << pxsec;
       }
                
       //-- push Q2*(d^2xsec/dWdQ2) to the FunctionMap
       
       Q2d2xsec.AddPoint(Q2*pxsec, i, j);
              
    } //Q2
  } //W

  //----- Perform the numerical integration

  LOG("RESXSec", pDEBUG) << "Performing numerical integration";

  //-- get specified integration algorithm from the config. registry
  //   or use Simpson2D if none is defined
  
  const IntegratorI * integrator = this->Integrator();

  double xsec = integrator->Integrate(Q2d2xsec);

  LOG("RESXSec", pINFO)  << "xsec[RES] (Ev = " << Ev << " GeV) = " << xsec;

  return xsec;
}
//____________________________________________________________________________
const XSecAlgorithmI * RESXSec::PartialXSecAlgorithm(void) const
{
  assert(
     fConfig->Exists("partial-xsec-alg-name") &&
                                     fConfig->Exists("partial-xsec-param-set")
  );

  //-- Get the partial xsec alg-name & param-set from the config. registry

  string alg_name, param_set;

  fConfig->Get("partial-xsec-alg-name",  alg_name  );
  fConfig->Get("partial-xsec-param-set", param_set );

  //----- Get the requested algorithm from the algorithm factory

  AlgFactory * algf = AlgFactory::Instance();

  const Algorithm * algbase = algf->GetAlgorithm(alg_name, param_set);

  const XSecAlgorithmI * xsec_alg =
                            dynamic_cast<const XSecAlgorithmI *> (algbase);

  assert(xsec_alg);

  return xsec_alg;
}
//____________________________________________________________________________
const IntegratorI * RESXSec::Integrator(void) const
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
bool RESXSec::IsKinematicallyAllowed(const Interaction * interaction) const
{
// checks whether the current Q2, W is kinematically allowed

  //-- get physical integration range for W and Q2

  Range1D_t rW  = kine_limits::WRange    (interaction);
  Range1D_t rQ2 = kine_limits::Q2Range_W (interaction);

  LOG("RESXSec", pDEBUG)
       << "\n Physical W integration range: "
                               << "[" << rW.min << ", " << rW.max << "] GeV";
  LOG("RESXSec", pDEBUG)
          << "\n Physical Q2 integration range: "
                           << "[" << rQ2.min << ", " << rQ2.max << "] GeV^2";

  //-- get current W, Q2

  double W   = interaction->GetScatteringParams().W();
  double Q2  = interaction->GetScatteringParams().Q2();

  bool allowed = math_utils::IsWithinLimits(Q2, rQ2)
                                         && math_utils::IsWithinLimits(W, rW);
  if(!allowed) {
       LOG("RESXSec", pDEBUG)
            << "*** ommitted from integral: "<< ", Q2 = "<< Q2<< ", W = "<< W;
  }
  return allowed;
}
//____________________________________________________________________________
int RESXSec::NW(void) const
{
  return ( fConfig->Exists("nW") ) ? fConfig->GetInt("nW") : 61;
}
//____________________________________________________________________________
double RESXSec::Wmin(void) const
{
  return ( fConfig->Exists("Wmin")   ) ? fConfig->GetDouble("Wmin") : 0.5;
}  
//____________________________________________________________________________
double RESXSec::Wmax(void) const
{
  return ( fConfig->Exists("Wmax")   ) ? fConfig->GetDouble("Wmax") : 2.0;
}
//____________________________________________________________________________
int RESXSec::NLogQ2(void) const
{
  return ( fConfig->Exists("nLogQ2") ) ? fConfig->GetInt("nLogQ2") : 81;
}
//____________________________________________________________________________
double RESXSec::Q2min(void) const
{  
  return ( fConfig->Exists("Q2min")   ) ? fConfig->GetDouble("Q2min") : 0.1;
}
//____________________________________________________________________________
double RESXSec::Q2max(void) const
{
  return ( fConfig->Exists("Q2max")   ) ? fConfig->GetDouble("Q2max") : 5.0;
}
//____________________________________________________________________________


