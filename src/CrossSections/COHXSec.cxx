//____________________________________________________________________________
/*!

\class    genie::COHXSec

\brief    Computes the cross section for COH neutrino-nucleus scattering.

          Is a concrete implementation of the XSecAlgorithmI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 04, 2004

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Algorithm/AlgFactory.h"
#include "Conventions/Constants.h"
#include "CrossSections/COHXSec.h"
#include "Messenger/Messenger.h"
#include "Numerical/UnifGrid.h"
#include "Numerical/FunctionMap.h"
#include "Numerical/IntegratorI.h"
#include "PDG/PDGUtils.h"
#include "Utils/MathUtils.h"
#include "Utils/KineUtils.h"
#include "Utils/Range1.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
COHXSec::COHXSec() :
XSecAlgorithmI("genie::COHXSec")
{

}
//____________________________________________________________________________
COHXSec::COHXSec(string config) :
XSecAlgorithmI("genie::COHXSec", config)
{

}
//____________________________________________________________________________
COHXSec::~COHXSec()
{

}
//____________________________________________________________________________
double COHXSec::XSec(const Interaction * in) const
{
  if(! this -> ValidProcess    (in) ) return 0.;
  if(! this -> ValidKinematics (in) ) return 0.;

  Interaction * interaction = new Interaction(*in);
  interaction->SetBit(kISkipProcessChk);
  interaction->SetBit(kISkipKinematicChk);

  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  // Get the neutrino energy in LAB
  const InitialState & init_state = interaction -> GetInitialState();
  double Ev = init_state.GetProbeE(kRfLab);

  // Define the integration grid & instantiate a FunctionMap
  double Mpi     = kPionMass;
  double e       = 1e-3;
  double ymin    = Mpi/Ev + e;
  double ymax    = 1. - e;
  double xmin    = 0. + e;
  double xmax    = 1. - e;
  double logxmin = TMath::Log(xmin);
  double logxmax = TMath::Log(xmax);
  double logymin = TMath::Log(ymin);
  double logymax = TMath::Log(ymax);
  double dlogx   = (logxmax-logxmin)/(fNlogx-1);
  double dlogy   = (logymax-logymin)/(fNlogy-1);

  UnifGrid grid;
  grid.AddDimension(fNlogx, logxmin, logxmax);
  grid.AddDimension(fNlogy, logymin, logymax);

  FunctionMap xyd2xsec(grid);

  // Loop over x,y & compute the differential cross section
  for(int ix = 0; ix < fNlogx; ix++) {
    double x = TMath::Exp(logxmin + ix * dlogx);

    for(int iy = 0; iy < fNlogy; iy++) {
       double y = TMath::Exp(logymin + iy * dlogy);

       //-- update the scattering parameters
       interaction->GetKinematicsPtr()->Setx(x);
       interaction->GetKinematicsPtr()->Sety(y);
       //-- compute d^2xsec/dxdy
       double pxsec = fPartialXSecAlg->XSec(interaction);
       LOG("COHXSec", pDEBUG)
              << "dxsec/dxdy (x = " << x << ", y = " << y
                                 << ", Ev = " << Ev << ") = " << pxsec;
       //-- update max differential xsec
       //max_pxsec = TMath::Max(max_pxsec, pxsec);
       //-- push x*y*(d^2xsec/dxdy) to the FunctionMap
       xyd2xsec.AddPoint(x*y*pxsec, ix, iy);
    } //y
  } //x

  delete interaction;

  // Perform the numerical integration
  LOG("COHXSec", pDEBUG) << "Performing numerical integration";
  double xsec = fIntegrator->Integrate(xyd2xsec);

  LOG("COHXSec", pDEBUG)  << "XSec[COH] (E = " << Ev << " GeV) = " << xsec;

  return xsec;
}
//____________________________________________________________________________
bool COHXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const InitialState & init_state = interaction->GetInitialState();
  const ProcessInfo &  proc_info  = interaction->GetProcessInfo();

  if (!proc_info.IsCoherent()) return false;
  if (!proc_info.IsWeakNC())   return false;

  int  nu = init_state.GetProbePDGCode();
  if (!pdg::IsNeutrino(nu) && !pdg::IsAntiNeutrino(nu)) return false;

  bool hitnuc = init_state.GetTarget().StruckNucleonIsSet();
  if(hitnuc) return false;

  return true;
}
//____________________________________________________________________________
bool COHXSec::ValidKinematics(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipKinematicChk)) return true;

  const InitialState & init_state = interaction -> GetInitialState();
  double Ev   = init_state.GetProbeE(kRfLab);
  double Ethr = kPionMass;

  if(Ev <= Ethr) {
     LOG("COHXSec", pINFO) << "E = " << Ev << " <= Ethreshold = " << Ethr;
     return false;
  }
  return true;
}
//____________________________________________________________________________
void COHXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
  this->LoadSubAlg();
}
//____________________________________________________________________________
void COHXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
  this->LoadSubAlg();
}
//____________________________________________________________________________
void COHXSec::LoadConfigData(void)
{
  //-- Get x,y from config (if they exist) or set defaults
  fNlogx = fConfig -> GetIntDef    ("n-log-x", 61   );
  fNlogy = fConfig -> GetIntDef    ("n-log-y", 61   );
  fXmin  = fConfig -> GetDoubleDef ("x-min",   0.001);
  fXmax  = fConfig -> GetDoubleDef ("x-max",   0.999);
  fYmin  = fConfig -> GetDoubleDef ("y-min",   0.001);
  fYmax  = fConfig -> GetDoubleDef ("y-max",   0.999);

  //-- Check that x,y range is meaningful
  assert(fNlogx > 1 && fNlogy > 1);
  assert(fXmax > fXmin && fXmax < 1 && fXmin < 1 && fXmax > 0 & fXmin > 0);
  assert(fYmax > fYmin && fYmax < 1 && fYmin < 1 && fYmax > 0 & fYmin > 0);
}
//____________________________________________________________________________
void COHXSec::LoadSubAlg(void)
{
  fPartialXSecAlg = 0;
  fIntegrator     = 0;

  //-- Get the requested d^2xsec/dxdy xsec algorithm to use
  fPartialXSecAlg =
         dynamic_cast<const XSecAlgorithmI *> (this->SubAlg(
                         "partial-xsec-alg-name", "partial-xsec-param-set"));
  assert(fPartialXSecAlg);

  LOG("COHXSec", pDEBUG) << *fPartialXSecAlg;

  //-- get specified integration algorithm from the config. registry
  //   or use Simpson2D if none is defined
  string intgr = fConfig->GetStringDef("integrator", "genie::Simpson2D");

  AlgFactory * algf = AlgFactory::Instance();
  fIntegrator = dynamic_cast<const IntegratorI *> (algf->GetAlgorithm(intgr));
  assert(fIntegrator);
}
//____________________________________________________________________________
