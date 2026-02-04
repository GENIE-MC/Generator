//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Alfonso Garcia <aagarciasoto \at km3net.de>
 IFIC & Harvard University
*/
//____________________________________________________________________________

#include <TMath.h>

#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/KineUtils.h"
#include "Physics/HELepton/XSection/GLRESPXSec.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
GLRESPXSec::GLRESPXSec() :
XSecAlgorithmI("genie::GLRESPXSec")
{
  born = new Born();

}
//____________________________________________________________________________
GLRESPXSec::GLRESPXSec(string config) :
XSecAlgorithmI("genie::GLRESPXSec", config)
{

}
//____________________________________________________________________________
GLRESPXSec::~GLRESPXSec()
{

}
//____________________________________________________________________________
double GLRESPXSec::XSec(
   const Interaction * interaction, KinePhaseSpace_t kps) const
{

  if(! this -> ValidProcess    (interaction) ) return 0.;

  const InitialState & init_state = interaction -> InitState();
  const Kinematics &   kinematics = interaction -> Kine();
  const XclsTag &      xclstag    = interaction -> ExclTag();

  int loutpdg  = xclstag.FinalLeptonPdg();

  double mlout = interaction->FSPrimLepton()->Mass(); //mass of charged lepton
  double mlin  = kElectronMass;                       //mass of incoming charged lepton

  double Enuin = init_state.ProbeE(kRfLab);
  double s = born->GetS(mlin,Enuin);

  double n1 = kinematics.GetKV(kKVn1);
  double n2 = kinematics.GetKV(kKVn2);
  double t  = born->GetT(mlin,mlout,s,n1);
  if (t>0) return 0.;

  //nlo correction
  double zeta     = born->GetReAlpha()/kPi*(2.*TMath::Log(TMath::Sqrt(-t)/kElectronMass)-1.);
  double omx      = TMath::Power(n2, 1./zeta );
  double pdf_soft = TMath::Exp(zeta*(3./4.-TMath::EulerGamma()))/TMath::Gamma(1.+zeta) + omx*(omx-2.)/2./n2;
  if ( omx<0. || omx>1. ) return 0.;
  double s_r = s*(1. - omx);
  double t_r = t*(1. - omx);

  //http://users.jyu.fi/~tulappi/fysh300sl11/l2.pdf [slide 22]
  //remember we always define nuout as p4
  double Enuout = (mlin*mlin-t_r)/2./mlin;
  if ( !born->IsInPhaseSpace(mlin,mlout,Enuin,Enuout) ) return 0.;

  double xsec = kPi/4./(s-mlin*mlin) * pdf_soft ;
  
  if ( pdg::IsPion(loutpdg) ) {
    if ( TMath::Sqrt(s_r)<fWmin ) return 0.; // The W limit is because hadronization might have issues at low W (as in PYTHIA6).
    xsec *= 64.41/10.63;    
  }

  double ME = 0.;
  if ( loutpdg == kPdgElectron ) ME = born->PXSecCCRNC(s_r,t_r,mlin,mlout);
  else                           ME = born->PXSecCCR  (s_r,t_r,mlin,mlout); 
  xsec *= TMath::Max(0.,ME);

  //----- If requested return the free electron xsec even for nuclear target
  if( interaction->TestBit(kIAssumeFreeElectron) ) return xsec;
   
  //----- Scale for the number of scattering centers at the target
  int Ne = init_state.Tgt().Z(); // num of scattering centers
  xsec *= Ne;

  if(kps!=kPSn1n2fE) {
      LOG("GLRESPXSec", pWARN)
          << "Doesn't support transformation from "
          << KinePhaseSpace::AsString(kPSn1n2fE) << " to "
          << KinePhaseSpace::AsString(kps);
      xsec = 0;
  }

  LOG("GLRESPXSec", pINFO) << "dxsec/dn1dn2 (E= " << Enuin << ", n1= " << n1 << ", n2=" << n2 << ") = " << xsec;

  return xsec;

}
//____________________________________________________________________________
double GLRESPXSec::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);

  return xsec;
}
//____________________________________________________________________________
bool GLRESPXSec::ValidProcess(const Interaction* interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const ProcessInfo & proc_info  = interaction->ProcInfo();
  if(!proc_info.IsGlashowResonance()) return false;
  if(!proc_info.IsWeakCC()) return false;

  const InitialState & init_state = interaction -> InitState();
  if(!pdg::IsAntiNuE(init_state.ProbePdg())) return false;

  if(init_state.Tgt().HitNucIsSet()) return false;
 
  return true;
}
//____________________________________________________________________________
void GLRESPXSec::Configure(const Registry & config)
{

  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void GLRESPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void GLRESPXSec::LoadConfig(void)
{

  //-- load the differential cross section integrator
  fXSecIntegrator = dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);

  GetParam( "Xsec-Wmin", fWmin ) ;


}