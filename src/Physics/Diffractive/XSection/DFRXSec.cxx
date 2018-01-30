//____________________________________________________________________________
/*
 Copyright (c) 2003-2018, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include <TMath.h>
#include <Math/IFunction.h>
#include <Math/IntegratorMultiDim.h>
#include "Math/AdaptiveIntegratorMultiDim.h"

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/RefFrame.h"
#include "Physics/XSectionIntegration/GSLXSecFunc.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Physics/Coherent/XSection/BergerSehgalCOHPiPXSec2015.h"
#include "Framework/Numerical/GSLUtils.h"
#include "Framework/Utils/HadXSUtils.h"
#include "Framework/Utils/KineUtils.h"

#include "Physics/Diffractive/XSection/DFRXSec.h"

using namespace genie;
using namespace genie::utils;

//____________________________________________________________________________
DFRXSec::DFRXSec ()
 : XSecIntegratorI("genie::DFRXSec")
{

}

//____________________________________________________________________________
DFRXSec::DFRXSec (std::string config)
 : XSecIntegratorI("genie::DFRXSec", config)
{

}

//____________________________________________________________________________
DFRXSec::~DFRXSec()
{

}

//____________________________________________________________________________
double DFRXSec::Integrate (const XSecAlgorithmI* model, const Interaction* in) const
{
  if(! model->ValidProcess(in) ) return 0.;

  const KPhaseSpace & kps = in->PhaseSpace();
  if(!kps.IsAboveThreshold()) {
     LOG("DFRXSec", pDEBUG)  << "*** Below energy threshold";
     return 0;
  }
  Range1D_t xl = kps.Limits(kKVx);
  Range1D_t yl = kps.Limits(kKVy);
  // the actual t lower limit depends on Q^2 and nu, or, equivalently, x and y.
  // defer to KPhaseSpace::IsAlllowed() on where to start the t integral.
  Range1D_t tl;
  tl.min = controls::kASmallNum;
  tl.max = fTMax;

  Interaction * interaction = new Interaction(*in);
  interaction->SetBit(kISkipProcessChk);   // todo: was enabled in COH model.  do I need it here?
  //interaction->SetBit(kISkipKinematicChk);

  double xsec = 0.0;
  LOG("DFRXSec", pINFO)
    << "x integration range = [" << xl.min << ", " << xl.max << "]";
  LOG("DFRXSec", pINFO)
    << "y integration range = [" << yl.min << ", " << yl.max << "]";
  LOG("DFRXSec", pINFO)
    << "t integration range = [" << tl.min << ", " << tl.max << "]";

  ROOT::Math::IBaseFunctionMultiDim * func =
    new utils::gsl::d3XSec_dxdydt_E(model, interaction);
  ROOT::Math::IntegrationMultiDim::Type ig_type =
    utils::gsl::IntegrationNDimTypeFromString(fGSLIntgType);

  double abstol = 1; //We mostly care about relative tolerance.
  ROOT::Math::IntegratorMultiDim ig(*func, ig_type, abstol, fGSLRelTol, fGSLMaxEval);
  if (ig_type == ROOT::Math::IntegrationMultiDim::kADAPTIVE) {
    ROOT::Math::AdaptiveIntegratorMultiDim * cast =
      dynamic_cast<ROOT::Math::AdaptiveIntegratorMultiDim*>( ig.GetIntegrator() );
    assert(cast);
    cast->SetMinPts(fGSLMinEval);
  }

  double kine_min[3] = { xl.min, yl.min, tl.min };
  double kine_max[3] = { xl.max, yl.max, tl.max };
  xsec = ig.Integral(kine_min, kine_max) * (1E-38 * units::cm2);
  return xsec;
}

//____________________________________________________________________________
void DFRXSec::Configure (const Registry& config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}

//____________________________________________________________________________
void DFRXSec::Configure (std::string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}

//____________________________________________________________________________
void DFRXSec::LoadConfig (void)
{

  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();

  // Get GSL integration type & relative tolerance
  fGSLIntgType = fConfig->GetStringDef("gsl-integration-type",  "adaptive");
  fGSLRelTol   = fConfig->GetDoubleDef("gsl-relative-tolerance", 1E-2);
  fGSLMaxEval  = (unsigned int) fConfig->GetIntDef ("gsl-max-eval", 500000);
  fGSLMinEval  = (unsigned int) fConfig->GetIntDef ("gsl-min-eval",  5000);

  //-- DFR model parameter t_max for t = (q - p_pi)^2
  fTMax = fConfig->GetDoubleDef("DFR-t-max", gc->GetDouble("DFR-t-max"));
}
