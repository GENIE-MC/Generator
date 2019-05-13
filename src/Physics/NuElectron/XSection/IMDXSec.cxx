//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 07, 2009 - CA
   Integrated with GNU Numerical Library (GSL) via ROOT's MathMore library.

*/
//____________________________________________________________________________

#include <TMath.h>
#include <Math/IFunction.h>
#include <Math/Integrator.h>

#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/RefFrame.h"
#include "Physics/NuElectron/XSection/IMDXSec.h"
#include "Physics/XSectionIntegration/GSLXSecFunc.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/GSLUtils.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
IMDXSec::IMDXSec() :
XSecIntegratorI("genie::IMDXSec")
{

}
//____________________________________________________________________________
IMDXSec::IMDXSec(string config) :
XSecIntegratorI("genie::IMDXSec", config)
{

}
//____________________________________________________________________________
IMDXSec::~IMDXSec()
{

}
//____________________________________________________________________________
double IMDXSec::Integrate(
                 const XSecAlgorithmI * model, const Interaction * in) const
{
  if(! model->ValidProcess(in) ) return 0.;

  const KPhaseSpace & kps = in->PhaseSpace();
  if(!kps.IsAboveThreshold()) {
     LOG("IMDXSec", pDEBUG)  << "*** below energy threshold";
     return 0;
  }
  Range1D_t yl = kps.Limits(kKVy);

  Interaction * interaction = new Interaction(*in);
  interaction->SetBit(kISkipProcessChk);
  //interaction->SetBit(kISkipKinematicChk);

  double abstol = 1; // We mostly care about relative tolerance
  ROOT::Math::IBaseFunctionOneDim * func = new utils::gsl::dXSec_dy_E(model, interaction);
  ROOT::Math::IntegrationOneDim::Type ig_type = utils::gsl::Integration1DimTypeFromString(fGSLIntgType);
  ROOT::Math::Integrator ig(*func,ig_type,abstol,fGSLRelTol,fGSLMaxEval);
  double xsec = ig.Integral(yl.min, yl.max) * (1E-38 * units::cm2);

  //LOG("IMDXSec", pDEBUG) << "XSec[IMD] (E = " << E << ") = " << xsec;

  delete interaction;
  delete func;

  return xsec;
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
  // Get GSL integration type & relative tolerance
	GetParamDef( "gsl-integration-type", fGSLIntgType, string( "adaptive" ) ) ;
	GetParamDef( "gsl-relative-tolerance", fGSLRelTol, 1E-4 ) ;
	int max;
	GetParamDef( "gsl-max-eval", max, 100000 ) ;
	fGSLMaxEval  = (unsigned int) max ;

}
//____________________________________________________________________________


