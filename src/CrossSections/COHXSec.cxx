//____________________________________________________________________________
/*
 Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
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

#include "Conventions/GBuild.h"
#include "Conventions/Constants.h"
#include "Conventions/Controls.h"
#include "Conventions/Units.h"
#include "CrossSections/COHXSec.h"
#include "CrossSections/GSLXSecFunc.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGUtils.h"
#include "Utils/MathUtils.h"
#include "Utils/Range1.h"
#include "Utils/GSLUtils.h"

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
  Q2l.max = 1.0;  // TODO - No magic numbers! 

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
  } 
  else if (model->Id().Name() == "genie::BergerSehgalCOHPiPXSec" || 
           model->Id().Name() == "genie::BergerSehgalFMCOHPiPXSec") 
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
  fGSLIntgType = fConfig->GetStringDef("gsl-integration-type",  "adaptive");
  fGSLRelTol   = fConfig->GetDoubleDef("gsl-relative-tolerance", 1E-2); 
  fGSLMaxEval  = (unsigned int) fConfig->GetIntDef ("gsl-max-eval", 500000); 
  fGSLMinEval  = (unsigned int) fConfig->GetIntDef ("gsl-min-eval",  5000); 
}
//____________________________________________________________________________
