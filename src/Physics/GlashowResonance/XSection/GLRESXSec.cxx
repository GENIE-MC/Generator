//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Alfonso Garcia <alfonsog \at nikhef.nl>
         NIKHEF (Amsterdam)

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TMath.h>
#include <Math/IFunction.h>
#include <Math/Integrator.h>

#include "Physics/GlashowResonance/XSection/GLRESXSec.h"
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Physics/XSectionIntegration/GSLXSecFunc.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/RunOpt.h"
#include "Framework/Utils/Range1.h"
#include "Framework/Utils/Cache.h"
#include "Framework/Utils/CacheBranchFx.h"
#include "Framework/Utils/XSecSplineList.h"
#include "Framework/Numerical/GSLUtils.h"

using namespace genie;
using namespace genie::controls;
using namespace genie::constants;

//____________________________________________________________________________
GLRESXSec::GLRESXSec() :
XSecIntegratorI("genie::GLRESXSec")
{

}
//____________________________________________________________________________
GLRESXSec::GLRESXSec(string config) :
XSecIntegratorI("genie::GLRESXSec", config)
{

}
//____________________________________________________________________________
GLRESXSec::~GLRESXSec()
{

}
//____________________________________________________________________________
double GLRESXSec::Integrate(
                 const XSecAlgorithmI * model, const Interaction * in) const
{
  if(! model->ValidProcess(in) ) return 0.;

  const KPhaseSpace & kps = in->PhaseSpace();
  if(!kps.IsAboveThreshold()) {
     LOG("GLRESXSec", pDEBUG)  << "*** Below energy threshold";
     return 0;
  }


  const InitialState & init_state = in->InitState();
  double Ev = init_state.ProbeE(kRfLab);

  int NNucl   = init_state.Tgt().Z();
  
  // If the input interaction is off a nuclear target, then chek whether 
  // the corresponding free nucleon cross section already exists at the 
  // cross section spline list. 
  // If yes, calculate the nuclear cross section based on that value.
  //
  XSecSplineList * xsl = XSecSplineList::Instance();
  if( !xsl->IsEmpty() ) {
    Interaction * interaction = new Interaction(*in);
    Target * target = interaction->InitStatePtr()->TgtPtr();
    target->SetId(kPdgTgtFreeP);
    if(xsl->SplineExists(model,interaction)) {
      const Spline * spl = xsl->GetSpline(model, interaction);
      double xsec = spl->Evaluate(Ev);
      LOG("GLRESXSec", pINFO)  << "From XSecSplineList: XSec[ve-,free nucleon] (E = " << Ev << " GeV) = " << xsec;
      if( !interaction->TestBit(kIAssumeFreeNucleon) ) { 
      	xsec *= NNucl; 
        LOG("GLRESXSec", pINFO)  << "XSec[ve-] (E = " << Ev << " GeV) = " << xsec;
      }
      delete interaction;
      return xsec;
    }
    delete interaction;
  }


  Range1D_t yl = kps.Limits(kKVy);

  LOG("GLRESXSec", pDEBUG) << "y = (" << yl.min << ", " << yl.max << ")";

  Interaction * interaction = new Interaction(*in);
  interaction->SetBit(kISkipProcessChk);
  //interaction->SetBit(kISkipKinematicChk);

  ROOT::Math::IBaseFunctionOneDim * func = 
     new utils::gsl::dXSec_dy_E(model, interaction);
  ROOT::Math::IntegrationOneDim::Type ig_type = 
     utils::gsl::Integration1DimTypeFromString(fGSLIntgType);
  ROOT::Math::Integrator ig(*func,ig_type,1,fGSLRelTol,fGSLMaxEval);
  double xsec = ig.Integral(yl.min, yl.max) * (1E-38 * units::cm2);

  LOG("GLRESXSec", pDEBUG) << "*** XSec[ve-] (E=" << interaction->InitState().ProbeE(kRfLab) << ") = " << xsec;

  delete interaction;
  delete func;
  return xsec;
}
//____________________________________________________________________________
void GLRESXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void GLRESXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void GLRESXSec::LoadConfig(void)
{
  // Get GSL integration type & relative tolerance
  GetParamDef("gsl-integration-type", fGSLIntgType, string("adaptive") ) ;
  GetParamDef("gsl-relative-tolerance", fGSLRelTol, 1E-2 ) ;
  int max_eval ;
  GetParamDef( "gsl-max-eval", max_eval, 500000 ) ;
  fGSLMaxEval  = (unsigned int) max_eval ;
}
//____________________________________________________________________________

