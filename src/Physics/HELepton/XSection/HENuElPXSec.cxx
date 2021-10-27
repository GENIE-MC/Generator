//____________________________________________________________________________
/*
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
#include "Physics/HELepton/XSection/HENuElPXSec.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
HENuElPXSec::HENuElPXSec() :
XSecAlgorithmI("genie::HENuElPXSec")
{
  born = new Born();

}
//____________________________________________________________________________
HENuElPXSec::HENuElPXSec(string config) :
XSecAlgorithmI("genie::HENuElPXSec", config)
{

}
//____________________________________________________________________________
HENuElPXSec::~HENuElPXSec()
{

}
//____________________________________________________________________________
double HENuElPXSec::XSec(
   const Interaction * interaction, KinePhaseSpace_t kps) const
{

  if(! this -> ValidProcess    (interaction) ) return 0.;

  const ProcessInfo &  proc_info  = interaction->ProcInfo();
  const InitialState & init_state = interaction -> InitState();
  const Kinematics &   kinematics = interaction -> Kine();
  const XclsTag &      xclstag    = interaction -> ExclTag();

  bool isCC = proc_info.isCC();

  double mlout = interaction->FSPrimLepton()->Mass(); //mass of charged lepton
  
  double E = init_state.ProbeE(kRfLab);
  double s = 2 * kElectronMass * E + kElectronMass2;

  double n1 = kinematics.GetKV(kKVn1);
  double n2 = kinematics.GetKV(kKVn2);
  double t  = 1;
  if ( isCC ) t = born->GetT( 0., kElectronMass, mlout, 0.   , s, n1 );
  else        t = born->GetT( 0., kElectronMass, 0.   , mlout, s, n1 );
  if (t>0) return 0.;

  //nlo correction
  double zeta = born->GetReAlpha()/kPi*(2.*TMath::Log(TMath::Sqrt(-t)/kElectronMass)-1.);
  double omx  = TMath::Power(n2, 1./zeta );
  if ( omx<0. || omx>1. ) return 0.;

  double s_r = s*(1. - omx);
  double t_r = t*(1. - omx);

  double pdf_soft = TMath::Exp(zeta*(3./4.-TMath::EulerGamma()))/TMath::Gamma(1.+zeta) + omx*(omx-2.)/2./n2;
  double xsec = kPi/4./(s-kElectronMass2) * pdf_soft ;
  
  double ME = 0;
  if ( pdg::IsNuE(nu) ) ME = PXSecCCVNC(s,t,kElectronMass2,mlout*mlout);
  else {
    if (isCC) ME = born->PXSecCCV(s_r,t_r,kElectronMass2,mlout*mlout);
    else      ME = born->PXSecNCV(s_r,t_r,kElectronMass2,mlout*mlout);
  }  
  xsec *= TMath::Max(0.,ME);

  //----- If requested return the free electron xsec even for nuclear target
  if( interaction->TestBit(kIAssumeFreeElectron) ) return xsec;
   
  //----- Scale for the number of scattering centers at the target
  int Ne = init_state.Tgt().Z(); // num of scattering centers
  xsec *= Ne;

  if(kps!=kPSn1n2fE) {
      LOG("HENuElPXSec", pWARN)
          << "Doesn't support transformation from "
          << KinePhaseSpace::AsString(kPSn1n2fE) << " to "
          << KinePhaseSpace::AsString(kps);
      xsec = 0;
  }

  LOG("HENuElPXSec", pINFO) << "dxsec/dn1dn2 (E= " << E << ", n1= " << n1 << ", n2=" << n2 << ") = " << xsec;

  return xsec;

}
//____________________________________________________________________________
double HENuElPXSec::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);

  return xsec;
}
//____________________________________________________________________________
bool HENuElPXSec::ValidProcess(const Interaction* interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const ProcessInfo & proc_info  = interaction->ProcInfo();
  if(!proc_info.IsGlashowResonance()) return false;

  const InitialState & init_state = interaction -> InitState();
  if(pdg::IsAntiNuE(init_state.ProbePdg())) return false;

  if(init_state.Tgt().HitNucIsSet()) return false;
 
  return true;
}
//____________________________________________________________________________
void HENuElPXSec::Configure(const Registry & config)
{

  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void HENuElPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void HENuElPXSec::LoadConfig(void)
{

  //-- load the differential cross section integrator
  fXSecIntegrator = dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);

}