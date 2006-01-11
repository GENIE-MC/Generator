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
#include "PDG/PDGUtils.h"
#include "Utils/KineUtils.h"

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

  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  //----- Define the integration grid & instantiate a FunctionMap
  UnifGrid grid;
  grid.AddDimension(fNLogt, fLogtmin, fLogtmax); // 1-D

  FunctionMap tdxsec(grid); // t * (dxsec/dt), t = x,y

  //----- Loop over t (x or y) and compute dxsec/dt
  for(int it = 0; it < fNLogt; it++) {
    double t  = TMath::Exp(fLogtmin + it * fdLogt);

    //-- update the "running" scattering parameter
    if      (fKineVar == "y") interaction->GetKinematicsPtr()->Setx(t);
    else if (fKineVar == "x") interaction->GetKinematicsPtr()->Sety(t);
    else    abort();

    //-- compute d^2xsec/dxdy
    double pxsec = fPartialXSecAlg->XSec(interaction);
    tdxsec.AddPoint(t*pxsec, it); // t * dxsec/dt (t = x or y)
  } //t

  //----- Numerical integration
  double xsec = fIntegrator->Integrate(tdxsec);
  return xsec;
}
//____________________________________________________________________________
bool DISPXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const InitialState & init_state = interaction->GetInitialState();
  const ProcessInfo &  proc_info  = interaction->GetProcessInfo();

  int  nuc = init_state.GetTarget().StruckNucleonPDGCode();
  int  nu  = init_state.GetProbePDGCode();

  if (!pdg::IsProton(nuc)  && !pdg::IsNeutron(nuc))     return false;
  if (!pdg::IsNeutrino(nu) && !pdg::IsAntiNeutrino(nu)) return false;
  if (!proc_info.IsDeepInelastic()) return false;
  if (!proc_info.IsWeak())          return false;

  return true;
}
//____________________________________________________________________________
bool DISPXSec::ValidKinematics(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipKinematicChk)) return true;

  //-- Get neutrino energy in the struck nucleon rest frame
  const InitialState & init_state = interaction -> GetInitialState();
  double Ev = init_state.GetProbeE(kRfStruckNucAtRest);

  //-- Check the energy threshold
  double Ethr = utils::kinematics::EnergyThreshold(interaction);
  if(Ev <= Ethr) {
     LOG("DISXSec", pINFO) << "E = " << Ev << " <= Ethreshold = " << Ethr;
     return false;
  }

  return true;
}
//____________________________________________________________________________
void DISPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
  this->LoadSubAlg();
}
//____________________________________________________________________________
void DISPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
  this->LoadSubAlg();
}
//____________________________________________________________________________
void DISPXSec::LoadConfigData(void)
{
  //-- Make sure it knows what kind of partial (dxsec/d?) xsec algorithm it is

  assert( fConfig->Exists("is-differential-over") );
  fKineVar = fConfig->GetString("is-differential-over");
  LOG("DISPXSec", pDEBUG) << "XSec is differential over var: " << fKineVar;

  //-- Get x or y integration range from config (if exists)
  int    nlogt;
  double tmin, tmax;
  double e = 1e-4;

  if ( fKineVar == "y" ) {
     //-- it is dxsec/dy, get intergation limits over x
     nlogt = fConfig->GetIntDef    ("n-log-x", 61 );
     tmin  = fConfig->GetDoubleDef ("x-min",   e  );
     tmax  = fConfig->GetDoubleDef ("x-max",   1-e);

  } else if ( fKineVar == "x" ) {
     //-- it is dxsec/dx, get intergation limits over y
     nlogt = fConfig->GetIntDef    ("n-log-y", 61 );
     tmin  = fConfig->GetDoubleDef ("y-min",   e  );
     tmax  = fConfig->GetDoubleDef ("y-max",   1-e);
  }

  LOG("DISPXSec", pDEBUG) << "Integration: n(log)bins = "
                     << nlogt << ", min = " << tmin << ", max = " << tmax;

  //-- Check that t (x or y) range is meaningful
  assert( nlogt > 1 );
  assert( tmax > tmin && tmax < 1 && tmin < 1 && tmax > 0 & tmin > 0 );

  fNLogt   = nlogt;
  fLogtmax = TMath::Log(tmax);
  fLogtmin = TMath::Log(tmin);
  fdLogt   = (fLogtmax - fLogtmin) / (fNLogt-1);
}
//____________________________________________________________________________
void DISPXSec::LoadSubAlg(void)
{
  fPartialXSecAlg = 0;
  fIntegrator     = 0;

  //-- Get the requested d^2xsec/dxdy xsec algorithm to use
  fPartialXSecAlg =
         dynamic_cast<const XSecAlgorithmI *> (this->SubAlg(
                         "partial-xsec-alg-name", "partial-xsec-param-set"));
  assert(fPartialXSecAlg);

  LOG("DISPXSec", pDEBUG) << *fPartialXSecAlg;

  //-- get specified integration algorithm from the config. registry
  //   or use Simpson2D if none is defined
  string intgr = fConfig->GetStringDef("integrator", "genie::Simpson2D");

  AlgFactory * algf = AlgFactory::Instance();
  fIntegrator = dynamic_cast<const IntegratorI *> (algf->GetAlgorithm(intgr));
  assert(fIntegrator);
}
//____________________________________________________________________________

