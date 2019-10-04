//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Mar 03, 2009 - CA
   Renamed COHPiXSec -> COHXSec. Adapt to naming changes made to the coherent 
   generator for including coherent vector meson production.
 @ Sep 07, 2009 - CA
   Integrated with GNU Numerical Library (GSL) via ROOT's MathMore library.

*/
//____________________________________________________________________________

#include <TMath.h>
#include <Math/IFunction.h>
#include <Math/IntegratorMultiDim.h>
#include "Math/AdaptiveIntegratorMultiDim.h"

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Units.h"
#include "Physics/Coherent/XSection/COHXSec.h"
#include "Physics/XSectionIntegration/GSLXSecFunc.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/Utils/Range1.h"
#include "Framework/Numerical/GSLUtils.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::utils;

//____________________________________________________________________________
COHXSec::COHXSec() :
XSecIntegratorI("genie::COHXSec")
{

}
//____________________________________________________________________________
COHXSec::COHXSec(string config) :
XSecIntegratorI("genie::COHXSec", config)
{

}
//____________________________________________________________________________
COHXSec::~COHXSec()
{

}
//____________________________________________________________________________
double COHXSec::Integrate(
      const XSecAlgorithmI * model, const Interaction * in) const
{
  if(! model->ValidProcess(in) ) return 0.;

  const KPhaseSpace & kps = in->PhaseSpace();
  if(!kps.IsAboveThreshold()) {
     LOG("COHXSec", pDEBUG)  << "*** Below energy threshold";
     return 0;
  }
  Range1D_t xl = kps.Limits(kKVx);
  Range1D_t yl = kps.Limits(kKVy);
  Range1D_t Q2l;
  Q2l.min = controls::kASmallNum;
  Q2l.max = fQ2Max;

  Interaction * interaction = new Interaction(*in);
  interaction->SetBit(kISkipProcessChk);
  //interaction->SetBit(kISkipKinematicChk);

  double xsec = 0.0;

  if (model->Id().Name() == "genie::ReinSehgalCOHPiPXSec") {
    LOG("COHXSec", pINFO)
      << "x integration range = [" << xl.min << ", " << xl.max << "]";
    LOG("COHXSec", pINFO)
      << "y integration range = [" << yl.min << ", " << yl.max << "]";

    ROOT::Math::IBaseFunctionMultiDim * func = 
      new utils::gsl::d2XSec_dxdy_E(model, interaction);
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
  
    double kine_min[2] = { xl.min, yl.min };
    double kine_max[2] = { xl.max, yl.max };
    xsec = ig.Integral(kine_min, kine_max) * (1E-38 * units::cm2);
    delete func;
  } 
  else if (model->Id().Name() == "genie::BergerSehgalCOHPiPXSec2015")
  {
    ROOT::Math::IBaseFunctionMultiDim * func = 
      new utils::gsl::d2XSec_dQ2dy_E(model, interaction);
    ROOT::Math::IntegrationMultiDim::Type ig_type = 
      utils::gsl::IntegrationNDimTypeFromString(fGSLIntgType);
    ROOT::Math::IntegratorMultiDim ig(ig_type);
    ig.SetRelTolerance(fGSLRelTol);
    ig.SetFunction(*func);
    if (ig_type == ROOT::Math::IntegrationMultiDim::kADAPTIVE) {
    ROOT::Math::AdaptiveIntegratorMultiDim * cast =
      dynamic_cast<ROOT::Math::AdaptiveIntegratorMultiDim*>( ig.GetIntegrator() );
      assert(cast);
      cast->SetMinPts(fGSLMinEval);
    }
    double kine_min[2] = { Q2l.min, yl.min };
    double kine_max[2] = { Q2l.max, yl.max };
    xsec = ig.Integral(kine_min, kine_max) * (1E-38 * units::cm2);
    delete func;
  }
  else if (model->Id().Name() == "genie::BergerSehgalFMCOHPiPXSec2015") 
  {
    Range1D_t tl;
    tl.min = controls::kASmallNum;
    tl.max = fTMax;

    ROOT::Math::IBaseFunctionMultiDim * func = 
      new utils::gsl::d2XSec_dQ2dydt_E(model, interaction);
    ROOT::Math::IntegrationMultiDim::Type ig_type = 
      utils::gsl::IntegrationNDimTypeFromString(fGSLIntgType);
    ROOT::Math::IntegratorMultiDim ig(ig_type);
    ig.SetRelTolerance(fGSLRelTol);
    ig.SetFunction(*func);
    if (ig_type == ROOT::Math::IntegrationMultiDim::kADAPTIVE) {
    ROOT::Math::AdaptiveIntegratorMultiDim * cast =
      dynamic_cast<ROOT::Math::AdaptiveIntegratorMultiDim*>( ig.GetIntegrator() );
      assert(cast);
      cast->SetMinPts(fGSLMinEval);
    }
    double kine_min[3] = { Q2l.min, yl.min, tl.min };
    double kine_max[3] = { Q2l.max, yl.max, tl.max };
    xsec = ig.Integral(kine_min, kine_max) * (1E-38 * units::cm2);
    delete func;
  }

  const InitialState & init_state = in->InitState();
  double Ev = init_state.ProbeE(kRfLab);
  LOG("COHXSec", pINFO) << "XSec[COH] (E = " << Ev << " GeV) = " << xsec;

  delete interaction;

  return xsec;
}
//____________________________________________________________________________
void COHXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void COHXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void COHXSec::LoadConfig(void)
{

  // Get GSL integration type & relative tolerance
  GetParamDef( "gsl-integration-type", fGSLIntgType, string("adaptive") ) ;
  GetParamDef( "gsl-relative-tolerance", fGSLRelTol, 1E-2 ) ;

  int max_eval, min_eval ;
  GetParamDef( "gsl-max-eval", max_eval, 500000 ) ;
  GetParamDef( "gsl-min-eval", min_eval, 5000 ) ;

  fGSLMaxEval  = (unsigned int) max_eval ;
  fGSLMinEval  = (unsigned int) min_eval ;

  //-- COH model parameter t_max for t = (q - p_pi)^2
  GetParam("COH-t-max", fTMax ) ;

  //-- COH model bounds of integration for Q^2
  GetParam( "COH-Q2-min", fQ2Min ) ;
  GetParam( "COH-Q2-max", fQ2Max ) ;

}
//____________________________________________________________________________
