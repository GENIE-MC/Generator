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

#include "Conventions/Constants.h"
#include "CrossSections/COHXSec.h"
#include "CrossSections/GXSecFunc.h"
#include "Messenger/Messenger.h"
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

  GXSecFunc * func =
       new Integrand_D2XSec_DxDy_E(fPartialXSecAlg, interaction);
  func->SetParam(0,"x",xmin,xmax);
  func->SetParam(1,"y",ymin,ymax);
  double xsec = fIntegrator->Integrate(*func);

  LOG("COHXSec", pDEBUG)  << "XSec[COH] (E = " << Ev << " GeV) = " << xsec;

  delete interaction;
  delete func;
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
  this->LoadConfig();
}
//____________________________________________________________________________
void COHXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void COHXSec::LoadConfig(void)
{
  fPartialXSecAlg = 0;
  fIntegrator     = 0;

  //-- Get the requested d^2xsec/dxdy xsec algorithm to use
  fPartialXSecAlg =
         dynamic_cast<const XSecAlgorithmI *> (this->SubAlg(
                         "partial-xsec-alg-name", "partial-xsec-param-set"));
  assert(fPartialXSecAlg);

  LOG("COHXSec", pDEBUG) << *fPartialXSecAlg;

  //-- get specified integration algorithm
  fIntegrator = dynamic_cast<const IntegratorI *> (
                 this->SubAlg("integrator-alg-name", "integrator-param-set"));
  assert(fIntegrator);
}
//____________________________________________________________________________
