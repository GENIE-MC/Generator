//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool
*/
//____________________________________________________________________________

#include <TMath.h>
#include <Math/IFunction.h>
#include <Math/Integrator.h>

#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/RefFrame.h"
#include "Physics/NuElectron/XSection/NuElectronXSec.h"
#include "Physics/XSectionIntegration/GSLXSecFunc.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/GSLUtils.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
NuElectronXSec::NuElectronXSec() :
XSecIntegratorI("genie::NuElectronXSec")
{

}
//____________________________________________________________________________
NuElectronXSec::NuElectronXSec(string config) :
XSecIntegratorI("genie::NuElectronXSec", config)
{

}
//____________________________________________________________________________
NuElectronXSec::~NuElectronXSec()
{

}
//____________________________________________________________________________
double NuElectronXSec::Integrate(
                 const XSecAlgorithmI * model, const Interaction * in) const
{
  if(! model->ValidProcess(in) ) return 0.;

  const KPhaseSpace & kps = in->PhaseSpace();
  if(!kps.IsAboveThreshold()) {
     LOG("NuEXSec", pDEBUG)  << "*** Below energy threshold";
     return 0;
  }
  Range1D_t yl = kps.Limits(kKVy);

  LOG("NuEXSec", pDEBUG) << "y = (" << yl.min << ", " << yl.max << ")";

  Interaction * interaction = new Interaction(*in);
  interaction->SetBit(kISkipProcessChk);
  //interaction->SetBit(kISkipKinematicChk);

  ROOT::Math::IBaseFunctionOneDim * func =
     new utils::gsl::dXSec_dy_E(model, interaction);
  ROOT::Math::IntegrationOneDim::Type ig_type =
     utils::gsl::Integration1DimTypeFromString(fGSLIntgType);
  ROOT::Math::Integrator ig(*func,ig_type,1,fGSLRelTol,fGSLMaxEval);
  double xsec = ig.Integral(yl.min, yl.max) * (1E-38 * units::cm2);

  //LOG("NuEXSec", pDEBUG) << "*** XSec[ve-] (E=" << E << ") = " << xsec;

  delete interaction;
  delete func;
  return xsec;
}
//____________________________________________________________________________
void NuElectronXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void NuElectronXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void NuElectronXSec::LoadConfig(void)
{
  // Get GSL integration type & relative tolerance
	GetParamDef( "gsl-integration-type", fGSLIntgType, string( "adaptive" ) ) ;
	GetParamDef( "gsl-relative-tolerance", fGSLRelTol, 0.01 ) ;
	int max;
	GetParamDef( "gsl-max-eval", max, 100000 ) ;
	fGSLMaxEval  = (unsigned int) max ;
}
//____________________________________________________________________________
