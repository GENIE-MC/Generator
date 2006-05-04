//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - February 14, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Conventions/Constants.h"
#include "Conventions/RefFrame.h"
#include "CrossSections/IMDXSec.h"
#include "CrossSections/GXSecFunc.h"
#include "Messenger/Messenger.h"
#include "Numerical/Simpson1D.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
IMDXSec::IMDXSec() :
XSecAlgorithmI("genie::IMDXSec")
{

}
//____________________________________________________________________________
IMDXSec::IMDXSec(string config) :
XSecAlgorithmI("genie::IMDXSec", config)
{

}
//____________________________________________________________________________
IMDXSec::~IMDXSec()
{

}
//____________________________________________________________________________
double IMDXSec::XSec(
              const Interaction * interaction, KinePhaseSpace_t kps) const
{
  assert(kps==kPSfE);

  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  const InitialState & init_state = interaction->GetInitialState();
  double E  = init_state.GetProbeE(kRfLab);

  double e = 1e-4;
  Range1D_t y(e, 1.-e);
  GXSecFunc * func = new Integrand_DXSec_Dy_E(fDiffXSecModel, interaction);
  func->SetParam(0,"y",y);
  double xsec = fIntegrator->Integrate(*func);

  LOG("IMD", pDEBUG) << "XSec[IMD] (E = " << E << ") = " << xsec;

  delete func;
  return xsec;
}
//____________________________________________________________________________
bool IMDXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;
  return true;
}
//____________________________________________________________________________
bool IMDXSec::ValidKinematics(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipKinematicChk)) return true;

  const InitialState & init_state = interaction -> GetInitialState();

  double E = init_state.GetProbeE(kRfLab);
  double s = kElectronMass2 + 2*kElectronMass*E;

  //-- check if it is kinematically allowed
  if(s < kMuonMass2) {
     LOG("IMD", pINFO)
        << "Ev = " << E << " (s = " << s << ") is below threshold (s-min = "
        << kMuonMass2 << ") for IMD";
     return false;
  }
  return true;
}
//____________________________________________________________________________
void IMDXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void IMDXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void IMDXSec::LoadConfig(void)
{
  fDiffXSecModel = 0;
  fIntegrator    = 0;

  //-- get an algorithm to calculate differential cross sections dxsec/dQ2
  fDiffXSecModel =
          dynamic_cast<const XSecAlgorithmI *> (this->SubAlg(
                         "partial-xsec-alg-name", "partial-xsec-param-set"));
  assert(fDiffXSecModel);

  //-- get specified integration algorithm
  fIntegrator = dynamic_cast<const IntegratorI *> (
                this->SubAlg("integrator-alg-name", "integrator-param-set"));
  assert(fIntegrator);
}
//____________________________________________________________________________


