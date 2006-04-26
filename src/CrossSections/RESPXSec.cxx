//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - March 03, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Conventions/Constants.h"
#include "CrossSections/RESPXSec.h"
#include "CrossSections/GXSecFunc.h"
#include "Messenger/Messenger.h"
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

  GXSecFunc * func = 0;
  double xsec = 0;

  if (fKineVar == "W") {

    double Q2 = interaction->GetKinematics().Q2();
    func = new Integrand_D2XSec_DWDQ2_EQ2(
                         fPartialXSecAlg, interaction, Q2);

    // default is physical W range for the given energy
    Range1D_t rW = utils::kinematics::WRange(interaction);
    // apply kinematic cuts
    if ( utils::math::IsWithinLimits(fKineMinCut, rW) ) rW.min = fKineMinCut;
    if ( utils::math::IsWithinLimits(fKineMaxCut, rW) ) rW.max = fKineMaxCut;
    assert(rW.min < rW.max);

    func->SetParam(0,"W",rW);
    xsec = fIntegrator->Integrate(*func);

  } else if (fKineVar == "Q2") {

    double W = interaction->GetKinematics().W();
    func = new Integrand_D2XSec_DWDQ2_EW(
                         fPartialXSecAlg, interaction, W);

    // default is physical Q2 range for input W
    Range1D_t rQ2 = utils::kinematics::Q2Range_W(interaction);
    // apply kinematic cuts
    if ( utils::math::IsWithinLimits(fKineMinCut, rQ2) ) rQ2.min = fKineMinCut;
    if ( utils::math::IsWithinLimits(fKineMaxCut, rQ2) ) rQ2.max = fKineMaxCut;
    assert(rQ2.min > 0 && rQ2.max > rQ2.min);

    func->SetParam(0,"Q2",rQ2);
    xsec = fIntegrator->Integrate(*func);

  } else abort();

  delete func;
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

  double tmin, tmax;

  if ( fKineVar.find("W") != string::npos ) {
     //-- it is dxsec/dW, get user cuts over Q2 and integration steps
     tmin  = fConfig -> GetDoubleDef ("Q2min",  -1 );
     tmax  = fConfig -> GetDoubleDef ("Q2max",  -1 );

  } else if ( fKineVar.find("Q2") != string::npos ) {
     //-- it is dxsec/dQ2, get user cuts over W and integration steps
     tmin = fConfig -> GetDoubleDef ("Wmin", -1 );
     tmax = fConfig -> GetDoubleDef ("Wmax", -1 );

  } else abort();

  fKineMinCut = tmin;
  fKineMaxCut = tmax;
}
//____________________________________________________________________________
void RESPXSec::LoadSubAlg(void)
{
  fPartialXSecAlg = 0;
  fIntegrator     = 0;

  //-- get the requested d^2xsec/dxdy xsec algorithm to use
  fPartialXSecAlg =
         dynamic_cast<const XSecAlgorithmI *> (this->SubAlg(
                         "partial-xsec-alg-name", "partial-xsec-param-set"));
  assert(fPartialXSecAlg);

  LOG("DISXSec", pDEBUG) << *fPartialXSecAlg;
  //-- get the specified integration algorithm
  fIntegrator = dynamic_cast<const IntegratorI *> (
                 this->SubAlg("integrator-alg-name", "integrator-param-set"));
  assert(fIntegrator);
}
//____________________________________________________________________________

