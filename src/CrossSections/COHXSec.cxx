//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory 

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

#include "Conventions/GBuild.h"
#include "Conventions/Constants.h"
#include "Conventions/Units.h"
#include "CrossSections/COHXSec.h"
#include "CrossSections/GXSecFunc.h"
#include "CrossSections/GSLXSecFunc.h"
#include "Messenger/Messenger.h"
#include "Numerical/IntegratorI.h"
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

  LOG("COHXSec", pINFO)
       << "x integration range = [" << xl.min << ", " << xl.max << "]";
  LOG("COHXSec", pINFO)
       << "y integration range = [" << yl.min << ", " << yl.max << "]";

  Interaction * interaction = new Interaction(*in);
  interaction->SetBit(kISkipProcessChk);
  //interaction->SetBit(kISkipKinematicChk);

#ifdef __GENIE_GSL_ENABLED__
  ROOT::Math::IBaseFunctionMultiDim * func = 
      new utils::gsl::wrap::d2XSec_dxdy_E(model, interaction);
  ROOT::Math::IntegrationMultiDim::Type ig_type = 
      utils::gsl::IntegrationNDimTypeFromString(fGSLIntgType);
  ROOT::Math::IntegratorMultiDim ig(ig_type);
  ig.SetRelTolerance(fGSLRelTol);
  ig.SetFunction(*func);
  double kine_min[2] = { xl.min, yl.min };
  double kine_max[2] = { xl.max, yl.max };
  double xsec = ig.Integral(kine_min, kine_max) * (1E-38 * units::cm2);

#else
  GXSecFunc * func = new Integrand_D2XSec_DxDy_E(model, interaction);
  func->SetParam(0,"x",xl);
  func->SetParam(1,"y",yl);
  double xsec = fIntegrator->Integrate(*func);

#endif

  const InitialState & init_state = in->InitState();
  double Ev = init_state.ProbeE(kRfLab);
  LOG("COHXSec", pINFO) << "XSec[COH] (E = " << Ev << " GeV) = " << xsec;

  delete interaction;
  delete func;
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
  // Get specified GENIE integration algorithm
  fIntegrator = dynamic_cast<const IntegratorI *> (this->SubAlg("Integrator"));
  assert(fIntegrator);

  // Get GSL integration type & relative tolerance
  fGSLIntgType = fConfig->GetStringDef("gsl-integration-type",  "adaptive");
  fGSLRelTol   = fConfig->GetDoubleDef("gsl-relative-tolerance", 0.01);
}
//____________________________________________________________________________
