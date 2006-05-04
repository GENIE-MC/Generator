//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - July 01, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <cstdlib>

#include <TMath.h>

#include "Conventions/Constants.h"
#include "Conventions/Controls.h"
#include "CrossSections/DISPXSec.h"
#include "CrossSections/GXSecFunc.h"
#include "Messenger/Messenger.h"
#include "Numerical/IntegratorI.h"
#include "PDG/PDGUtils.h"
#include "Utils/KineUtils.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;

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
double DISPXSec::XSec(
                  const Interaction * interaction, KinePhaseSpace_t kps) const
{
  assert(kps==kPSxfE || kps==kPSyfE);

  LOG("DISPXSec", pDEBUG) << *fConfig;

  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  GXSecFunc * func = 0;
  double xsec = 0;

  if (kps==kPSxfE) {

    double x = interaction->GetKinematics().x();
    func = new Integrand_D2XSec_DxDy_Ex(
                         fPartialXSecAlg, interaction, x);
    func->SetParam(0,"y",kMinY,kMaxY);
    xsec = fIntegrator->Integrate(*func);

  } else if (kps==kPSyfE) {

    double y = interaction->GetKinematics().y();
    func = new Integrand_D2XSec_DxDy_Ey(
                         fPartialXSecAlg, interaction, y);
    func->SetParam(0,"x",kMinX,kMaxX);
    xsec = fIntegrator->Integrate(*func);

  } else {
    LOG("DISPXSec", pFATAL) 
         << "Can not handle phase space = "  << KinePhaseSpace::AsString(kps);
    exit(1);
  }
  delete func;
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
  this->LoadConfig();
}
//____________________________________________________________________________
void DISPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DISPXSec::LoadConfig(void)
{
  fPartialXSecAlg = 0;
  fIntegrator     = 0;

  //-- get the requested d^2xsec/dxdy xsec algorithm to use
  fPartialXSecAlg =
         dynamic_cast<const XSecAlgorithmI *> (this->SubAlg(
                         "partial-xsec-alg-name", "partial-xsec-param-set"));
  assert(fPartialXSecAlg);

  LOG("DISPXSec", pDEBUG) << *fPartialXSecAlg;

  //-- get the specified integration algorithm
  fIntegrator = dynamic_cast<const IntegratorI *> (
                 this->SubAlg("integrator-alg-name", "integrator-param-set"));
  assert(fIntegrator);
}
//____________________________________________________________________________

