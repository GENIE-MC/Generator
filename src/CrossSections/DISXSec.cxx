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
#include "CrossSections/DISXSec.h"
#include "CrossSections/GXSecFunc.h"
#include "Messenger/Messenger.h"
#include "Numerical/IntegratorI.h"
#include "PDG/PDGUtils.h"
#include "Utils/MathUtils.h"
#include "Utils/Range1.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
DISXSec::DISXSec() :
XSecAlgorithmI("genie::DISXSec")
{

}
//____________________________________________________________________________
DISXSec::DISXSec(string config) :
XSecAlgorithmI("genie::DISXSec", config)
{

}
//____________________________________________________________________________
DISXSec::~DISXSec()
{

}
//____________________________________________________________________________
double DISXSec::XSec(const Interaction * in, KinePhaseSpace_t kps) const
{
  assert(kps==kPSfE);

  if(! this -> ValidProcess    (in) ) return 0.;
  if(! this -> ValidKinematics (in) ) return 0.;

  Interaction * interaction = new Interaction(*in);
  interaction->SetBit(kISkipProcessChk);
  interaction->SetBit(kISkipKinematicChk);

  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  // Get neutrino energy in the struck nucleon rest frame
  const InitialState & init_state = interaction -> InitState();
  double Ev = init_state.ProbeE(kRfHitNucRest);

  Range1D_t WCuts (fWmin, fWmax );
  Range1D_t Q2Cuts(fQ2min,fQ2max);

  GXSecFunc * func = new
    Integrand_D2XSec_DxDy_E_WQ2Cuts(fPartialXSecAlg, interaction,WCuts,Q2Cuts);

  func->SetParam(0,"x",fXmin,fXmax);
  func->SetParam(1,"y",fYmin,fYmax);
  double xsec = fIntegrator->Integrate(*func);

  LOG("DISXSec", pDEBUG)  << "XSec[DIS] (E = " << Ev << " GeV) = " << xsec;

  delete interaction;
  delete func;
  return xsec;
}
//____________________________________________________________________________
bool DISXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const InitialState & init_state = interaction->InitState();
  const ProcessInfo &  proc_info  = interaction->ProcInfo();

  int  nuc = init_state.Tgt().HitNucPdg();
  int  nu  = init_state.ProbePdg();

  if (!pdg::IsProton(nuc)  && !pdg::IsNeutron(nuc))     return false;
  if (!pdg::IsNeutrino(nu) && !pdg::IsAntiNeutrino(nu)) return false;
  if (!proc_info.IsDeepInelastic()) return false;
  if (!proc_info.IsWeak())          return false;

  return true;
}
//____________________________________________________________________________
bool DISXSec::ValidKinematics(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipKinematicChk)) return true;

  //-- Get neutrino energy in the struck nucleon rest frame
  const InitialState & init_state = interaction -> InitState();
  double Ev = init_state.ProbeE(kRfHitNucRest);

  //-- Check the energy threshold
  double Ethr = interaction->EnergyThreshold();
  if(Ev <= Ethr) {
     LOG("DISXSec", pINFO) << "E = " << Ev << " <= Ethreshold = " << Ethr;
     return false;
  }
  return true;
}
//____________________________________________________________________________
void DISXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DISXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DISXSec::LoadConfig(void)
{
  //-- Get the requested d^2xsec/dxdy xsec algorithm to use
  fPartialXSecAlg = 0;
  fPartialXSecAlg =
         dynamic_cast<const XSecAlgorithmI *> (this->SubAlg("DiffXSecAlg"));
  assert(fPartialXSecAlg);

  //-- get specified integration algorithm
  fIntegrator = 0;
  fIntegrator = dynamic_cast<const IntegratorI *> (this->SubAlg("Integrator"));
  assert(fIntegrator);

  double e=1E-3;
  fXmin = fConfig -> GetDoubleDef ("Kine-Xmin", e);
  fXmax = fConfig -> GetDoubleDef ("Kine-Xmax", 1-e);
  fYmin = fConfig -> GetDoubleDef ("Kine-Ymin", e);
  fYmax = fConfig -> GetDoubleDef ("Kine-Ymax", 1-e);

  //-- Check that x,y range is meaningful
  assert(fXmax > fXmin && fXmax < 1 && fXmin < 1 && fXmax > 0 & fXmin > 0);
  assert(fYmax > fYmin && fYmax < 1 && fYmin < 1 && fYmax > 0 & fYmin > 0);

  // check whether the user imposes kinematic cuts
  fWmin  = fConfig->GetDoubleDef( "Kine-Wmin",  -1  );
  fWmax  = fConfig->GetDoubleDef( "Kine-Wmax",  1e9 );
  fQ2min = fConfig->GetDoubleDef( "Kine-Q2min", -1  );
  fQ2max = fConfig->GetDoubleDef( "Kine-Q2max", 1e9 );
}
//____________________________________________________________________________

