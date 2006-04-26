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

#include <TMath.h>

#include "Conventions/Constants.h"
#include "CrossSections/DISPXSec.h"
#include "CrossSections/GXSecFunc.h"
#include "Messenger/Messenger.h"
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

  GXSecFunc * func = 0;
  double xsec = 0;

  if (fKineVar == "y") {
    double x = interaction->GetKinematics().x();
    func = new Integrand_D2XSec_DxDy_Ex(
                         fPartialXSecAlg, interaction, x);
    func->SetParam(0,"y",ftmin,ftmax);
    xsec = fIntegrator->Integrate(*func);

  } else if (fKineVar == "x") {

    double y = interaction->GetKinematics().y();
    func = new Integrand_D2XSec_DxDy_Ey(
                         fPartialXSecAlg, interaction, y);
    func->SetParam(0,"x",ftmin,ftmax);
    xsec = fIntegrator->Integrate(*func);

  } else abort();

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
  double e = 1e-4;

  if ( fKineVar == "y" ) {
     //-- it is dxsec/dy, get intergation limits over x
     ftmin  = fConfig->GetDoubleDef ("x-min",   e  );
     ftmax  = fConfig->GetDoubleDef ("x-max",   1-e);

  } else if ( fKineVar == "x" ) {
     //-- it is dxsec/dx, get intergation limits over y
     ftmin  = fConfig->GetDoubleDef ("y-min",   e  );
     ftmax  = fConfig->GetDoubleDef ("y-max",   1-e);
  }

  LOG("DISPXSec", pDEBUG)
           << "Integration range: " << "(" << ftmin << ", " << ftmax << ")";

  //-- Check that t (x or y) range is meaningful
  assert( ftmax > ftmin && ftmax < 1 && ftmin < 1 && ftmax > 0 & ftmin > 0 );
}
//____________________________________________________________________________
void DISPXSec::LoadSubAlg(void)
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

