//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - May 04, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Conventions/Constants.h"
#include "Conventions/KineVar.h"
#include "CrossSections/RESXSec.h"
#include "CrossSections/GXSecFunc.h"
#include "Messenger/Messenger.h"
#include "Numerical/IntegratorI.h"
#include "PDG/PDGUtils.h"
#include "Utils/MathUtils.h"
#include "Utils/KineUtils.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
RESXSec::RESXSec() :
XSecAlgorithmI("genie::RESXSec")
{

}
//____________________________________________________________________________
RESXSec::RESXSec(string config) :
XSecAlgorithmI("genie::RESXSec", config)
{

}
//____________________________________________________________________________
RESXSec::~RESXSec()
{

}
//____________________________________________________________________________
double RESXSec::XSec(const Interaction * in, KinePhaseSpace_t kps) const
{
  assert(kps==kPSfE);

  LOG("RESXSec", pDEBUG) << *fConfig;

  if(! this -> ValidProcess    (in) ) return 0.;
  if(! this -> ValidKinematics (in) ) return 0.;

  Interaction * interaction = new Interaction(*in);
  interaction->SetBit(kISkipProcessChk);
  //interaction->SetBit(kISkipKinematicChk);

  // Get neutrino energy in the struck nucleon rest frame
  const InitialState & init_state = interaction -> GetInitialState();
  double Ev  = init_state.GetProbeE(kRfStruckNucAtRest);

  // Get W integration range
  Range1D_t rW = this->WRange(interaction);
  // Get the wider possible Q2 range for the input W range
  Range1D_t rQ2 = utils::kinematics::Q2Range_W(interaction, rW);

  GXSecFunc * func = new Integrand_D2XSec_DWDQ2_E(fPartialXSecAlg, interaction);
  func->SetParam(0,"W",  rW);
  func->SetParam(1,"Q2", rQ2);
  double xsec = fIntegrator->Integrate(*func);

  SLOG("RESXSec", pINFO)  << "XSec[RES] (Ev = " << Ev << " GeV) = " << xsec;

  delete func;
  return xsec;
}
//____________________________________________________________________________
bool RESXSec::ValidProcess(const Interaction * interaction) const
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
bool RESXSec::ValidKinematics(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipKinematicChk)) return true;

  const InitialState & init_state = interaction -> GetInitialState();
  double Ev  = init_state.GetProbeE(kRfStruckNucAtRest);

  double EvThr = utils::kinematics::EnergyThreshold(interaction);
  if(Ev <= EvThr) return false;

  return true;
}
//____________________________________________________________________________
void RESXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void RESXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void RESXSec::LoadConfig(void)
{
  fPartialXSecAlg = 0;
  fIntegrator     = 0;

  //-- get the requested d^2xsec/dxdy xsec algorithm to use
  fPartialXSecAlg =
          dynamic_cast<const XSecAlgorithmI *> (SubAlg(
                          "partial-xsec-alg-name", "partial-xsec-param-set"));

  //-- get the specified integration algorithm
  fIntegrator = dynamic_cast<const IntegratorI *> (
                 this->SubAlg("integrator-alg-name", "integrator-param-set"));

  assert( fPartialXSecAlg );
  assert( fIntegrator     );
}
//____________________________________________________________________________
Range1D_t RESXSec::WRange(const Interaction * interaction) const
{
  //-- Get the physically allowed W range for this interaction and allow the
  //   user inputs (if any) to narrow it

  Range1D_t rW = utils::kinematics::KineRange(interaction, kKVW); 
  LOG("RESXSec", pDEBUG)
       << "Physical W range: " << "[" << rW.min << ", " << rW.max << "] GeV";

  // user cuts
  double Wmin = fConfig->GetDoubleDef("Wmin", -1.0);
  double Wmax = fConfig->GetDoubleDef("Wmax",  1e9);

  utils::kinematics::ApplyCutsToKineLimits(rW,  Wmin,  Wmax );

  LOG("RESXSec", pDEBUG)
       << "Physical & User W range: "
                               << "[" << rW.min << ", " << rW.max << "] GeV";
  return rW;
}
//___________________________________________________________________________
Range1D_t RESXSec::Q2Range(const Interaction * interaction) const
{
  //-- Get the physically allowed Q2 range for this interaction and allow the
  //   user inputs (if any) to narrow it

  Range1D_t rQ2 = utils::kinematics::KineRange(interaction, kKVQ2); 
  LOG("RESXSec", pDEBUG) << "Physical Q2 range: "
                 << "[" << rQ2.min << ", " << rQ2.max << "] GeV^2";

  // user cuts
  double Q2min = fConfig->GetDoubleDef("Q2min", -1.0);
  double Q2max = fConfig->GetDoubleDef("Q2max",  1e9);

  utils::kinematics::ApplyCutsToKineLimits(rQ2, Q2min, Q2max );

  LOG("RESXSec", pDEBUG)
       << "Physical & User Q2 range: "
                << "[" << rQ2.min << ", " << rQ2.max << "] GeV^2";
  return rQ2;
}
//___________________________________________________________________________


