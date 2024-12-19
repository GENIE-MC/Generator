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

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Conventions/RefFrame.h"
#include "Physics/QuasiElastic/XSection/QELXSec.h"

#include "Physics/NuclearState/NuclearModelI.h"
#include "Physics/XSectionIntegration/GSLXSecFunc.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Physics/NuclearState/NuclearUtils.h"
#include "Framework/Utils/Range1.h"
#include "Framework/Numerical/GSLUtils.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
QELXSec::QELXSec() :
XSecIntegratorI("genie::QELXSec")
{

}
//____________________________________________________________________________
QELXSec::QELXSec(string config) :
XSecIntegratorI("genie::QELXSec", config)
{

}
//____________________________________________________________________________
QELXSec::~QELXSec()
{

}
//____________________________________________________________________________
double QELXSec::Integrate(
                  const XSecAlgorithmI * model, const Interaction * in) const
{
  LOG("QELXSec",pDEBUG) << "Beginning integrate";
  if(! model->ValidProcess(in)) return 0.;

  const KPhaseSpace & kps = in->PhaseSpace();
  if(!kps.IsAboveThreshold()) {
     LOG("QELXSec", pDEBUG)  << "*** Below energy threshold";
     return 0;
  }
  Range1D_t rQ2 = kps.Limits(kKVQ2);
  if(rQ2.min<0 || rQ2.max<0) return 0;
  LOG("QELXSec", pDEBUG)
          << "Q2 integration range = (" << rQ2.min << ", " << rQ2.max << ")";

  Interaction * interaction = new Interaction(*in);
  interaction->SetBit(kISkipProcessChk);
  interaction->SetBit(kISkipKinematicChk);

  ROOT::Math::IBaseFunctionOneDim * func = new
      utils::gsl::dXSec_dQ2_E(model, interaction);
  ROOT::Math::IntegrationOneDim::Type ig_type =
      utils::gsl::Integration1DimTypeFromString(fGSLIntgType);

  double abstol = 0; //We mostly care about relative tolerance
  ROOT::Math::Integrator ig(*func,ig_type,abstol,fGSLRelTol,fGSLMaxSizeOfSubintervals, fGSLRule);
  double xsec = ig.Integral(rQ2.min, rQ2.max) * (1E-38 * units::cm2);

  //LOG("QELXSec", pDEBUG) << "XSec[QEL] (E = " << E << ") = " << xsec;

  delete func;
  delete interaction;

  return xsec;
}
//____________________________________________________________________________
void QELXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void QELXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void QELXSec::LoadConfig(void)
{
  // Get GSL integration type & relative tolerance
	GetParamDef( "gsl-integration-type", fGSLIntgType, string("adaptive"));
	GetParamDef( "gsl-relative-tolerance", fGSLRelTol, 0.001 ) ;
	int max_size_of_subintervals;
	GetParamDef( "gsl-max-size-of-subintervals", max_size_of_subintervals, 40000);
	fGSLMaxSizeOfSubintervals = (unsigned int) max_size_of_subintervals;
	int rule;
	GetParamDef( "gsl-rule", rule, 3);
	fGSLRule = (unsigned int) rule;
    if (fGSLRule>6) fGSLRule=3;
}
//____________________________________________________________________________
