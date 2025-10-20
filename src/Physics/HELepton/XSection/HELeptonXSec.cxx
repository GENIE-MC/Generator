//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Alfonso Garcia <aagarciasoto \at km3net.de>
 IFIC & Harvard University
*/
//____________________________________________________________________________

#include <TMath.h>
#include <Math/IFunction.h>
#include <Math/Integrator.h>

#include "Physics/HELepton/XSection/HELeptonXSec.h"
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
HELeptonXSec::HELeptonXSec() :
XSecIntegratorI("genie::HELeptonXSec")
{

}
//____________________________________________________________________________
HELeptonXSec::HELeptonXSec(string config) :
XSecIntegratorI("genie::HELeptonXSec", config)
{

}
//____________________________________________________________________________
HELeptonXSec::~HELeptonXSec()
{

}
//____________________________________________________________________________
double HELeptonXSec::Integrate(
                 const XSecAlgorithmI * model, const Interaction * in) const
{
  if(! model->ValidProcess(in) ) return 0.;

  const KPhaseSpace & kps = in->PhaseSpace();
  if(!kps.IsAboveThreshold()) {
     LOG("HELeptonXSec", pDEBUG)  << "*** Below energy threshold";
     return 0;
  }


  const ProcessInfo &  proc_info  = in->ProcInfo();
  const InitialState & init_state = in->InitState();
  double Ev = init_state.ProbeE(kRfLab);
  
  // If the input interaction is off a nuclear target, then chek whether 
  // the corresponding free nucleon cross section already exists at the 
  // cross section spline list. 
  // If yes, calculate the nuclear cross section based on that value.
  //
  XSecSplineList * xsl = XSecSplineList::Instance();
  if( !xsl->IsEmpty() && !proc_info.IsPhotonCoherent() ) {

    Interaction * interaction = new Interaction(*in);
    Target * target = interaction->InitStatePtr()->TgtPtr();

    int NNucl = 0;
    bool skip = false;
    if ( proc_info.IsGlashowResonance() ) {
      NNucl = init_state.Tgt().Z();
      target->SetId(kPdgTgtFreeP);
    }
    else if ( proc_info.IsPhotonResonance() ) {
      int nucpdgc = init_state.Tgt().HitNucPdg();
      if (pdg::IsProton(nucpdgc)) {
        NNucl = init_state.Tgt().Z();
        target->SetId(kPdgTgtFreeP);
      }
      else {
        NNucl = init_state.Tgt().N();
        target->SetId(kPdgTgtFreeN);
      }
    }

    if( xsl->SplineExists(model,interaction) ) {
      const Spline * spl = xsl->GetSpline(model, interaction);
      double xsec = spl->Evaluate(Ev);
      LOG("HELeptonXSec", pINFO)  << "From XSecSplineList: XSec["<<init_state.ProbePdg()<<","<<proc_info.ScatteringTypeAsString()<<","<<interaction->ExclTag().FinalLeptonPdg()<< ",free nucleon] (E = " << Ev << " GeV) = " << xsec;
      if( !interaction->TestBit(kIAssumeFreeNucleon) ) { 
        xsec *= NNucl; 
        LOG("HELeptonXSec", pINFO)  << "XSec["<<init_state.ProbePdg()<<","<<proc_info.ScatteringTypeAsString()<<","<<interaction->ExclTag().FinalLeptonPdg()<< "] (E = " << Ev << " GeV) = " << xsec;
      }
      delete interaction;
      return xsec;
    }
    delete interaction;
    
  }
   
  Interaction * interaction = new Interaction(*in);
  interaction->SetBit(kISkipProcessChk);

  ROOT::Math::IntegrationMultiDim::Type ig_type = utils::gsl::IntegrationNDimTypeFromString(fGSLIntgType);

  double xsec = 0.;
  ROOT::Math::IBaseFunctionMultiDim * func;
  if ( proc_info.IsPhotonCoherent() ) {
    double kine_min[3] = { -1.,  0.,  0. }; 
    double kine_max[3] = {  1.,  1.,  1. }; 
    func = new utils::gsl::d2Xsec_dn1dn2dn3_E(model, interaction);
    ROOT::Math::IntegratorMultiDim ig(*func,ig_type,0,0,fGSLMaxEval);
    xsec = ig.Integral(kine_min, kine_max);
  }
  else {
    double kine_min[2] = { -1.,  0. }; 
    double kine_max[2] = {  1.,  1. }; 
    func = new utils::gsl::d2Xsec_dn1dn2_E(model, interaction);
    ROOT::Math::IntegratorMultiDim ig(*func,ig_type,0,0,fGSLMaxEval);
    xsec = ig.Integral(kine_min, kine_max);
  }

  delete func;

  xsec *= (1E-38 * units::cm2);

  LOG("HELeptonXSec", pDEBUG) << "*** XSec["<<init_state.ProbePdg()<<","<<proc_info.ScatteringTypeAsString()<<","<<interaction->ExclTag().FinalLeptonPdg()<< "] (E=" << interaction->InitState().ProbeE(kRfLab) << ") = " << xsec / (1E-38 * units::cm2);

  delete interaction;
  return xsec;
}
//____________________________________________________________________________
void HELeptonXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void HELeptonXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void HELeptonXSec::LoadConfig(void)
{
  // Get GSL integration type & relative tolerance
  GetParamDef("gsl-integration-type", fGSLIntgType, string("vegas") ) ;

  int max_eval ;
  GetParamDef( "gsl-max-eval", max_eval, 500000 ) ;
  fGSLMaxEval  = (unsigned int) max_eval ;
}
//____________________________________________________________________________

