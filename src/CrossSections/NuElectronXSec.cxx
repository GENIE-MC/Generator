//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - February 10, 2006

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 07, 2009 - CA
   Integrated with GNU Numerical Library (GSL) via ROOT's MathMore library.

*/
//____________________________________________________________________________

#include <TMath.h>
#include <Math/IFunction.h>
#include <Math/Integrator.h>

#include "Conventions/GBuild.h"
#include "Conventions/Constants.h"
#include "Conventions/Units.h"
#include "Conventions/RefFrame.h"
#include "CrossSections/NuElectronXSec.h"
#include "CrossSections/GXSecFunc.h"
#include "CrossSections/GSLXSecFunc.h"
#include "Messenger/Messenger.h"
#include "Numerical/IntegratorI.h"
#include "Utils/GSLUtils.h"
  
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

#ifdef __GENIE_GSL_ENABLED__
  ROOT::Math::IBaseFunctionOneDim * func = 
     new utils::gsl::wrap::dXSec_dy_E(model, interaction);
  ROOT::Math::IntegrationOneDim::Type ig_type = 
     utils::gsl::Integration1DimTypeFromString(fGSLIntgType);
  ROOT::Math::Integrator ig(ig_type);
  ig.SetFunction(*func);
  ig.SetRelTolerance(fGSLRelTol);
  double xsec = ig.Integral(yl.min, yl.max) * (1E-38 * units::cm2);

#else
  GXSecFunc * func = new Integrand_DXSec_Dy_E(model, interaction);
  func->SetParam(0,"y",yl);
  double xsec = fIntegrator->Integrate(*func);

#endif

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
  // Get the specified GENIE integration algorithm
  fIntegrator = dynamic_cast<const IntegratorI *>(this->SubAlg("Integrator"));
  assert(fIntegrator);

  // Get GSL integration type & relative tolerance
  fGSLIntgType = fConfig->GetStringDef("gsl-integration-type",  "adaptive");
  fGSLRelTol   = fConfig->GetDoubleDef("gsl-relative-tolerance", 0.01);
}
//____________________________________________________________________________

