//____________________________________________________________________________
/*
 Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
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
#include <Math/IntegratorMultiDim.h>

#include "Conventions/GBuild.h"
#include "Conventions/Constants.h"
#include "Conventions/Units.h"
#include "Conventions/KineVar.h"
#include "CrossSections/RESXSec.h"
#include "CrossSections/GSLXSecFunc.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGUtils.h"
#include "Utils/MathUtils.h"
#include "Utils/KineUtils.h"
#include "Utils/GSLUtils.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
RESXSec::RESXSec() :
XSecIntegratorI("genie::RESXSec")
{

}
//____________________________________________________________________________
RESXSec::RESXSec(string config) :
XSecIntegratorI("genie::RESXSec", config)
{

}
//____________________________________________________________________________
RESXSec::~RESXSec()
{

}
//____________________________________________________________________________
double RESXSec::Integrate(
                 const XSecAlgorithmI * model, const Interaction * in) const
{
  if(! model->ValidProcess(in) ) return 0.;

  const KPhaseSpace & kps = in->PhaseSpace();
  if(!kps.IsAboveThreshold()) {
     LOG("COHXSec", pDEBUG)  << "*** Below energy threshold";
     return 0;
  }
  Range1D_t Q2l = kps.Limits(kKVQ2);
  Range1D_t Wl  = kps.Limits(kKVW);

  Interaction * interaction = new Interaction(*in);
  interaction->SetBit(kISkipProcessChk);
  //interaction->SetBit(kISkipKinematicChk);

  ROOT::Math::IBaseFunctionMultiDim * func = 
      new utils::gsl::d2XSec_dWdQ2_E(model, interaction);
  
  ROOT::Math::IntegrationMultiDim::Type ig_type = 
      utils::gsl::IntegrationNDimTypeFromString(fGSLIntgType);
  double abstol = 1E-16; //We mostly care about relative tolerance.
  
  ROOT::Math::IntegratorMultiDim ig(*func, ig_type, abstol, fGSLRelTol, fGSLMaxEval);

  double kine_min[2] = { Wl.min, Q2l.min };
  double kine_max[2] = { Wl.max, Q2l.max };
  double xsec = ig.Integral(kine_min, kine_max) * (1E-38 * units::cm2);

  LOG("RESXSec", pERROR)  << "Integrator opt / Integrator = " <<  ig.Options().Integrator();

  if(xsec < 0) {
    LOG("RESXSec", pERROR)  << "Algorithm " << *model << " returns a negative cross-section (xsec = " << xsec << " 1E-38 * cm2)";
    LOG("RESXSec", pERROR)  << "for process" << *interaction;
    LOG("RESXSec", pERROR)  << "Integrator status code = " << ig.Status();
    LOG("RESXSec", pERROR)  << "Integrator error code = " << ig.Error();
  }
         
  //LOG("RESXSec", pINFO)  << "XSec[RES] (Ev = " << Ev << " GeV) = " << xsec;

  delete interaction;
  delete func;
  return xsec;
}
//____________________________________________________________________________
void RESXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void RESXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void RESXSec::LoadConfig(void)
{
  // Get GSL integration type & relative tolerance
  fGSLIntgType = fConfig->GetStringDef("gsl-integration-type"  ,  "adaptive");
  fGSLRelTol   = fConfig->GetDoubleDef("gsl-relative-tolerance",   1E-2);
  fGSLMaxEval  = (unsigned int) fConfig->GetIntDef("gsl-max-eval", 500000);
  fGSLMinEval  = (unsigned int) fConfig->GetIntDef("gsl-min-eval", 5000);
}
//____________________________________________________________________________
