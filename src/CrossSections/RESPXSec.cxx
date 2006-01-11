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
#include "PDG/PDGUtils.h"
#include "Utils/KineUtils.h"
#include "Utils/MathUtils.h"
#include "Utils/Range1.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
RESPXSec::RESPXSec() :
XSecAlgorithmI("genie::RESPXSec")
{

}
//____________________________________________________________________________
RESPXSec::RESPXSec(string config) :
XSecAlgorithmI("genie::RESPXSec", config)
{

}
//____________________________________________________________________________
RESPXSec::~RESPXSec()
{

}
//____________________________________________________________________________
double RESPXSec::XSec(const Interaction * interaction) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  //-- Get the physicical t (W or Q2) integration range and apply user
  //   cuts (if they were specified)
  double tmin, tmax;
  if ( fKineVar.find("W") != string::npos ) {
     //-- it is dxsec/dW, get intergation limits over Q2

     // default is physical Q2 range for input W
     Range1D_t rQ2 = utils::kinematics::Q2Range_W(interaction);
     // apply kinematic cuts
     if ( utils::math::IsWithinLimits(fKineMinCut, rQ2) ) rQ2.min = fKineMinCut;
     if ( utils::math::IsWithinLimits(fKineMaxCut, rQ2) ) rQ2.max = fKineMaxCut;
     assert(rQ2.min > 0 && rQ2.max > rQ2.min);

     tmin = rQ2.min;
     tmax = rQ2.max;

  } else if ( fKineVar.find("Q2") != string::npos ) {
     //-- it is dxsec/dQ2, get intergation limits over W

     // default is physical W range for the given energy
     Range1D_t rW = utils::kinematics::WRange(interaction);
     // apply kinematic cuts
     if ( utils::math::IsWithinLimits(fKineMinCut, rW) ) rW.min = fKineMinCut;
     if ( utils::math::IsWithinLimits(fKineMaxCut, rW) ) rW.max = fKineMaxCut;
     assert(rW.min < rW.max);

     tmin = rW.min;
     tmax = rW.max;

  } else abort();

  LOG("RESPXSec", pDEBUG) << "Integration: n(log)bins = "
                       << fNLogt << ", min = " << tmin << ", max = " << tmax;
  assert(tmax > tmin);

  double logtmax = TMath::Log(tmax);
  double logtmin = TMath::Log(tmin);
  double dlogt   = (logtmax - logtmin) / (fNLogt-1);

  //-- Define the integration grid & instantiate a FunctionMap
  //   Do the integration in log(t)
  UnifGrid grid;
  grid.AddDimension(fNLogt, logtmin, logtmax); // 1-D

  FunctionMap funcmap(grid); // Q2*(d^2xsec/dlogQ2) or d^2xsec/dW

  //-- Loop over t (W or Q2) and compute dxsec/dt
  for(int it = 0; it < fNLogt; it++) {
    double t = TMath::Exp(logtmin + it * dlogt);

    if      (fKineVar == "W" ) interaction->GetKinematicsPtr()->SetQ2(t);
    else if (fKineVar == "Q2") interaction->GetKinematicsPtr()->SetW(t);
    else    abort();

    double d2xsec = fPartialXSecAlg->XSec(interaction);
    funcmap.AddPoint(t*d2xsec, it);
  } //t

  //-- Numerical integration
  double xsec = fIntegrator->Integrate(funcmap);
  LOG("RESPXSec", pDEBUG) << "dxsec[RES]/d" << fKineVar << " = " << xsec;
  return xsec;
}
//____________________________________________________________________________
bool RESPXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const InitialState & init_state = interaction->GetInitialState();
  const ProcessInfo &  proc_info  = interaction->GetProcessInfo();

  if(!proc_info.IsResonant()) return false;
  if(!proc_info.IsWeak())     return false;

  int  nuc = init_state.GetTarget().StruckNucleonPDGCode();
  int  nu  = init_state.GetProbePDGCode();

  if (!pdg::IsProton(nuc)  && !pdg::IsNeutron(nuc))     return false;
  if (!pdg::IsNeutrino(nu) && !pdg::IsAntiNeutrino(nu)) return false;

  return true;
}
//____________________________________________________________________________
bool RESPXSec::ValidKinematics(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipKinematicChk)) return true;

  const InitialState & init_state = interaction -> GetInitialState();
  double Ev  = init_state.GetProbeE(kRfStruckNucAtRest);

  double EvThr = utils::kinematics::EnergyThreshold(interaction);
  if(Ev <= EvThr) return false;

  return true;
}
//____________________________________________________________________________
void RESPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
  this->LoadSubAlg();
}
//____________________________________________________________________________
void RESPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
  this->LoadSubAlg();
}
//____________________________________________________________________________
void RESPXSec::LoadConfigData(void)
{
  assert( fConfig->Exists("is-differential-over") );
  fKineVar = fConfig->GetString("is-differential-over");
  LOG("RESPXSec", pDEBUG) << "XSec is differential over var: " << fKineVar;

  int    nt;
  double tmin, tmax;

  if ( fKineVar.find("W") != string::npos ) {
     //-- it is dxsec/dW, get user cuts over Q2 and integration steps
     nt    = fConfig -> GetIntDef    ("nLogQ2", 151);
     tmin  = fConfig -> GetDoubleDef ("Q2min",  -1 );
     tmax  = fConfig -> GetDoubleDef ("Q2max",  -1 );

  } else if ( fKineVar.find("Q2") != string::npos ) {
     //-- it is dxsec/dQ2, get user cuts over W and integration steps
     nt   = fConfig -> GetIntDef    ("nW",   51 );
     tmin = fConfig -> GetDoubleDef ("Wmin", -1 );
     tmax = fConfig -> GetDoubleDef ("Wmax", -1 );

  } else abort();

  assert(nt > 1);

  fNLogt      = nt;
  fKineMinCut = tmin;
  fKineMaxCut = tmax;
}
//____________________________________________________________________________
void RESPXSec::LoadSubAlg(void)
{
  fPartialXSecAlg = 0;
  fIntegrator     = 0;

  //-- Get the requested d^2xsec/dxdy xsec algorithm to use
  fPartialXSecAlg =
         dynamic_cast<const XSecAlgorithmI *> (this->SubAlg(
                         "partial-xsec-alg-name", "partial-xsec-param-set"));
  assert(fPartialXSecAlg);

  LOG("DISXSec", pDEBUG) << *fPartialXSecAlg;

  //-- get specified integration algorithm from the config. registry
  //   or use Simpson2D if none is defined
  string intgr = fConfig->GetStringDef("integrator", "genie::Simpson1D");

  AlgFactory * algf = AlgFactory::Instance();
  fIntegrator = dynamic_cast<const IntegratorI *> (algf->GetAlgorithm(intgr));
  assert(fIntegrator);
}
//____________________________________________________________________________

