//____________________________________________________________________________
/*
 Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Chris Marshall and Martti Nirkko

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include <TMath.h>
#include <Math/IFunction.h>
#include <Math/IntegratorMultiDim.h>
#include "Math/AdaptiveIntegratorMultiDim.h"

#include "SingleKaon/AlamSimoAtharVacasSKXSec.h"
#include "Conventions/GBuild.h"
#include "Conventions/Constants.h"
#include "Conventions/Controls.h"
#include "Conventions/Units.h"
#include "Conventions/KineVar.h"
#include "Conventions/KinePhaseSpace.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "Numerical/Spline.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGLibrary.h"
#include "Utils/MathUtils.h"
#include "Utils/KineUtils.h"
#include "Utils/Range1.h"
#include "Utils/GSLUtils.h"
#include "Utils/XSecSplineList.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;
using namespace genie::utils;

//____________________________________________________________________________
AlamSimoAtharVacasSKXSec::AlamSimoAtharVacasSKXSec() :
XSecIntegratorI("genie::AlamSimoAtharVacasSKXSec")
{

}
//____________________________________________________________________________
AlamSimoAtharVacasSKXSec::AlamSimoAtharVacasSKXSec(string config) :
XSecIntegratorI("genie::AlamSimoAtharVacasSKXSec", config)
{

}
//____________________________________________________________________________
AlamSimoAtharVacasSKXSec::~AlamSimoAtharVacasSKXSec()
{

}
//____________________________________________________________________________
double AlamSimoAtharVacasSKXSec::Integrate(
      const XSecAlgorithmI * model, const Interaction * in) const
{
  LOG("SKXSec", pDEBUG) << "Integrating the Alam Simo Athar Vacas model";

  const InitialState & init_state = in -> InitState();
  
  if(! model->ValidProcess(in) ) return 0.;
  
  const KPhaseSpace & kps = in->PhaseSpace(); // only OK phase space for this
  if(!kps.IsAboveThreshold()) {
     LOG("SKXSec", pDEBUG)  << "*** Below energy threshold";
     return 0;
  }

  // If the input interaction is off a nuclear target, then chek whether 
  // the corresponding free nucleon cross section already exists at the 
  // cross section spline list. 
  // Cross section for PP scales with number of protons, NP and NN scale 
  // with number of neutrons
  int nucpdgc = init_state.Tgt().HitNucPdg();
  int NNucl   = (pdg::IsProton(nucpdgc)) ? 
                   init_state.Tgt().Z() : init_state.Tgt().N();
  double Ev = init_state.ProbeE(kRfHitNucRest);

  XSecSplineList * xsl = XSecSplineList::Instance();
  if(init_state.Tgt().IsNucleus() && !xsl->IsEmpty() ) {
    Interaction * interaction = new Interaction(*in);
    Target * target = interaction->InitStatePtr()->TgtPtr();
    if(pdg::IsProton(nucpdgc)) { target->SetId(kPdgTgtFreeP); }
    else                       { target->SetId(kPdgTgtFreeN); }
    if(xsl->SplineExists(model,interaction)) {
      const Spline * spl = xsl->GetSpline(model, interaction);
      double xsec = spl->Evaluate(Ev);
      LOG("SKXSec", pINFO)  
        << "From XSecSplineList: XSec[SK,free nucleon] (E = " << Ev << " GeV) = " << xsec;
      if(! interaction->TestBit(kIAssumeFreeNucleon) ) { 
          xsec *= NNucl; 
          LOG("SKXSec", pINFO)  << "XSec[SK] (E = " << Ev << " GeV) = " << xsec;
      }
      delete interaction;
      return xsec;
    }
    delete interaction;
  }

  // no free nucelon spline exists -- do the integration
 
  // Check this
  double Enu = init_state.ProbeE(kRfLab);
  int kpdg = in->ExclTag().StrangeHadronPdg();
  double mk   = PDGLibrary::Instance()->Find(kpdg)->Mass();
  double ml   = PDGLibrary::Instance()->Find(in->FSPrimLeptonPdg())->Mass();

  // integration bounds for T (kinetic energy)
  double zero    = 0.0;
  double tmax = Enu - mk - ml;

  LOG("SKXSec", pDEBUG)
       << "Lepton/Kaon KE integration range = [" << 0.0 << ", " << tmax << "]";

  Interaction * interaction = new Interaction(*in);
  interaction->SetBit(kISkipProcessChk);
  interaction->SetBit(kISkipKinematicChk);
  
  double xsec = 0;

  // do the integration over log(1-costheta) so it's not so sharply peaked
    
  ROOT::Math::IBaseFunctionMultiDim * func = 
        new utils::gsl::d3Xsec_dTldTkdCosThetal(model, interaction);
  double kine_min[3] = { zero, zero, -20 }; // Tlep, Tkaon, cosine theta lep
  double kine_max[3] = { tmax, tmax,  0.69314718056 }; // Tlep, Tkaon, cosine theta lep
    
  ROOT::Math::IntegrationMultiDim::Type ig_type = 
    utils::gsl::IntegrationNDimTypeFromString(fGSLIntgType);
        
  double abstol = 1; //We mostly care about relative tolerance.
  ROOT::Math::IntegratorMultiDim ig(*func, ig_type, abstol, fGSLRelTol, fGSLMaxEval);
  
  xsec = ig.Integral(kine_min, kine_max) * (1E-38 * units::cm2);
  delete func;

  delete interaction;

  return xsec;
}
//____________________________________________________________________________
void AlamSimoAtharVacasSKXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void AlamSimoAtharVacasSKXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void AlamSimoAtharVacasSKXSec::LoadConfig(void)
{
  // Get GSL integration type & relative tolerance
  fGSLIntgType   = fConfig->GetStringDef("gsl-integration-type" ,    "vegas");
  fGSLMaxEval    = (unsigned int) fConfig->GetIntDef("gsl-max-evals", 20000);
  fGSLRelTol     = fConfig->GetDoubleDef("gsl-relative-tolerance",    0.01);
  fSplitIntegral = fConfig->GetBoolDef("split-integral",              true);
}
//_____________________________________________________________________________
// GSL wrappers
//____________________________________________________________________________

genie::utils::gsl::d3Xsec_dTldTkdCosThetal::d3Xsec_dTldTkdCosThetal(
     const XSecAlgorithmI * m, const Interaction * i) :
ROOT::Math::IBaseFunctionMultiDim(),
fModel(m),
fInteraction(i)
{
  
}

genie::utils::gsl::d3Xsec_dTldTkdCosThetal::~d3Xsec_dTldTkdCosThetal()
{
  
}   

unsigned int genie::utils::gsl::d3Xsec_dTldTkdCosThetal::NDim(void) const
{
  // phi_kq is important for kinematics generation
  // But dependence is weak so we will not use it in the integration
  return 3;
}

double genie::utils::gsl::d3Xsec_dTldTkdCosThetal::DoEval(const double * xin) const
{
// inputs:
//    Tl [GeV]
//    Tk [GeV]
//    cosine theta l
//    * calculate phi_kq based on neutrino energy -- this is for the integral only
// outputs:
//   differential cross section [10^-38 cm^2]
//
 
  double Enu = fInteraction->InitState().ProbeE(kRfLab);
    
  double phikq = 0.5*constants::kPi;
  if( Enu > 3.0 ) phikq = 0.55*constants::kPi;
  else if( Enu > 1.0 ) phikq = constants::kPi*(0.5 + 0.025*(Enu-1.0));

  Kinematics * kinematics = fInteraction->KinePtr();

  double T_l         = xin[0];
  double T_k         = xin[1];
  //double cos_theta_l = xin[2];
  double log_oneminuscostheta = xin[2];
  double cos_theta_l = 1.0 - TMath::Exp(log_oneminuscostheta);
  double J = 1.0 - cos_theta_l; // Jacobian for transformation
    
  kinematics->SetKV(kKVTl, T_l);
  kinematics->SetKV(kKVTk, T_k);
  kinematics->SetKV(kKVctl, cos_theta_l);
  kinematics->SetKV(kKVphikq, phikq);
  
  double xsec = fModel->XSec(fInteraction);
  LOG( "GXSecFunc", pDEBUG ) 
     << "t_l = " << T_l << " t_k = " << T_k << " costhetal = " << cos_theta_l << " phikq = " << phikq 
     << " enu = " << Enu << " Xsec = " << xsec;

  // integrate out phi_kq by multiplying by 2pi
  
  xsec *= 2.0 * genie::constants::kPi * J;
 
  return xsec/(1E-38 * units::cm2);
}

ROOT::Math::IBaseFunctionMultiDim *
   genie::utils::gsl::d3Xsec_dTldTkdCosThetal::Clone() const
{
  return
    new genie::utils::gsl::d3Xsec_dTldTkdCosThetal(fModel,fInteraction);
}
//____________________________________________________________________________


