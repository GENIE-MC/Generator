//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Alfonso Garcia <aagarciasoto \at km3net.de>
 IFIC & Harvard University
*/
//____________________________________________________________________________

#include "Physics/HELepton/XSection/PhotonCOHPXSec.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/KineUtils.h"

using namespace genie;
using namespace genie::constants;

const double a  = 0.523 * (units::fermi);
const double r0 = 1.126 * (units::fermi);

//____________________________________________________________________________
PhotonCOHPXSec::PhotonCOHPXSec() :
XSecAlgorithmI("genie::PhotonCOHPXSec")
{
  born = new Born();
}
//____________________________________________________________________________
PhotonCOHPXSec::PhotonCOHPXSec(string config) :
XSecAlgorithmI("genie::PhotonCOHPXSec", config)
{

}
//____________________________________________________________________________
PhotonCOHPXSec::~PhotonCOHPXSec()
{

}
//____________________________________________________________________________
double PhotonCOHPXSec::XSec(
   const Interaction * interaction, KinePhaseSpace_t kps) const
{

  if(! this -> ValidProcess    (interaction) ) return 0.;

  const InitialState & init_state = interaction -> InitState();
  const Kinematics &   kinematics = interaction -> Kine();

  int probepdg = init_state.ProbePdg();
  double E     = init_state.ProbeE(kRfLab);

  double mlout = 0;
  if      (pdg::IsNuE  (TMath::Abs(probepdg))) mlout = kElectronMass;
  else if (pdg::IsNuMu (TMath::Abs(probepdg))) mlout = kMuonMass;
  else if (pdg::IsNuTau(TMath::Abs(probepdg))) mlout = kTauMass;
  double mlout2 = mlout*mlout;

  int A        = init_state.Tgt().A(); 
  int Z        = init_state.Tgt().Z(); 
  double Mtgt = Z*kProtonMass + (A-Z)*kNeutronMass;

  double n1 = kinematics.GetKV(kKVn1);
  double n2 = kinematics.GetKV(kKVn2);
  double n3 = kinematics.GetKV(kKVn3);

  double mL = mlout+kMw;
  
  double Delta = TMath::Power(2.*E*Mtgt-mL*mL,2)-4.*TMath::Power(Mtgt*mL,2);
  if (Delta<0) return 0.;

  double s12_min = E/(2.*E+Mtgt)*(mL*mL+2.*E*Mtgt-TMath::Sqrt(Delta));
  double s12_max = E/(2.*E+Mtgt)*(mL*mL+2.*E*Mtgt+TMath::Sqrt(Delta));
  double s12 = TMath::Exp( TMath::Log(s12_min)+(TMath::Log(s12_max)-TMath::Log(s12_min))*n2);

  double Q2_min = TMath::Power(s12,2)*Mtgt/2./E/(2.*E*Mtgt-s12);
  double Q2_max = s12 - mL*mL;
  double Q2 = TMath::Exp( TMath::Log(Q2_min) + (TMath::Log(Q2_max)-TMath::Log(Q2_min))*n3 );

  double s = s12 - Q2;
  double s13 = s12/2.*((1.+kMw2/s-mlout2/s)-TMath::Sqrt(born->Lambda(1.,kMw2/s,mlout2/s))*n1);

  double ME_T = born->PXSecPhoton_T(s12,s13,Q2,mlout2) * (1.-s12/2./E/Mtgt-TMath::Power(s12/2./E,2)/Q2);
  double ME_L = born->PXSecPhoton_L(s12,s13,Q2,mlout2) * TMath::Power(1.-s12/4./E/Mtgt,2);

  double ME   = ME_T+ME_L;
  double dps2 = 1/32./kPi2/s12 * TMath::Sqrt( born->Lambda(1.,kMw2/s,mlout2/s) ) * (TMath::Log(Q2_max)-TMath::Log(Q2_min)) * (TMath::Log(s12_max)-TMath::Log(s12_min));    
  double dP   = born->GetReAlpha()*TMath::Power(Z,2)*F2_Q( TMath::Sqrt(Q2), r0*TMath::Power(A,1./3.) );

  double xsec = ME * dps2 * dP;
   
  if(kps!=kPSn1n2n3fE) {
      LOG("PhotonCOHPXSec", pWARN)
          << "Doesn't support transformation from "
          << KinePhaseSpace::AsString(kPSn1n2n3fE) << " to "
          << KinePhaseSpace::AsString(kps);
      xsec = 0;
  }

  // If requested return the free nucleon xsec even for input nuclear tgt 
  if( interaction->TestBit(kIAssumeFreeNucleon) ) return xsec;

  LOG("PhotonCOHPXSec", pINFO) << "dxsec/dn1dn2dn3 (E= " << E << ", n1= " << n1 << ", n2=" << n2 << ", n3=" << n3 << ") = " << xsec;

  return xsec;
}
//____________________________________________________________________________
double PhotonCOHPXSec::F2_Q(double Q, double r0) const
{
  // Analytic Woods-Saxon, A.3 of https://arxiv.org/pdf/1807.10973.pdf
  double FF = 0.0;
  double coth = 1./TMath::TanH(kPi*Q*a);
  FF = 3.*kPi*a/(TMath::Power(r0,2)+TMath::Power(kPi*a,2)) / (Q*r0* TMath::SinH(kPi*Q*a));
  FF *= (kPi*a*coth*TMath::Sin(Q*r0) - r0 * TMath::Cos(Q*r0));        
  return TMath::Power(FF,2);
}
//____________________________________________________________________________
double PhotonCOHPXSec::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool PhotonCOHPXSec::ValidProcess(const Interaction* interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const ProcessInfo & proc_info  = interaction->ProcInfo();
  if(!proc_info.IsPhotonCoherent()) return false;

  const InitialState & init_state = interaction -> InitState();
  if(!pdg::IsLepton(init_state.ProbePdg())) return false;

  if(init_state.Tgt().HitNucIsSet()) return false;
  
  return true;
}
//____________________________________________________________________________


//____________________________________________________________________________
void PhotonCOHPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void PhotonCOHPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void PhotonCOHPXSec::LoadConfig(void)
{

  //-- load the differential cross section integrator
  fXSecIntegrator = dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);

}
