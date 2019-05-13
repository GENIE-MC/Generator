//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Feb 12, 2013 - CA (code from Rosen Matev)
   Tweak kinematics. The final state primary lepton is always the electron.
   The kinematical variable y has the definition used in Marciano's paper.
*/
//____________________________________________________________________________

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Controls.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/RefFrame.h"
#include "Physics/NuElectron/XSection/NuElectronPXSec.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;

//____________________________________________________________________________
NuElectronPXSec::NuElectronPXSec() :
XSecAlgorithmI("genie::NuElectronPXSec")
{

}
//____________________________________________________________________________
NuElectronPXSec::NuElectronPXSec(string config) :
XSecAlgorithmI("genie::NuElectronPXSec", config)
{

}
//____________________________________________________________________________
NuElectronPXSec::~NuElectronPXSec()
{

}
//____________________________________________________________________________
double NuElectronPXSec::XSec(
                 const Interaction * interaction, KinePhaseSpace_t kps) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  //----- get initial state & kinematics
  const InitialState & init_state = interaction -> InitState();
  const Kinematics &   kinematics = interaction -> Kine();
  const ProcessInfo &  proc_info  = interaction -> ProcInfo();

  double Ev = init_state.ProbeE(kRfLab);
  double me = kElectronMass;
  double y  = kinematics.y();
  double A  = kGF2*2*me*Ev/kPi;

  y = 1 - me/Ev - y; // FSPL = electron. XSec below are expressed in Marciano's y!
  if(y > 1/(1+0.5*me/Ev)) return 0;
  if(y < 0) return 0;

  double xsec = 0; // <-- dxsec/dy

  int inu = init_state.ProbePdg();

  // nue + e- -> nue + e- [CC + NC + interference]
  if(pdg::IsNuE(inu)) 
  {
    double em = -0.5 - fSin28w;
    double ep = -fSin28w;
    xsec = TMath::Power(em,2) + TMath::Power(ep*(1-y),2) - ep*em*me*y/Ev;
    xsec *= A;
  }

  // nuebar + e- -> nue + e- [CC + NC + interference]
  if(pdg::IsAntiNuE(inu)) 
  {
    double em = -0.5 - fSin28w;
    double ep = -fSin28w;
    xsec = TMath::Power(ep,2) + TMath::Power(em*(1-y),2) - ep*em*me*y/Ev;
    xsec *= A;
  }

  // numu/nutau + e- -> numu/nutau + e- [NC]
  if( (pdg::IsNuMu(inu)||pdg::IsNuTau(inu)) && proc_info.IsWeakNC() ) 
  {
    double em = 0.5 - fSin28w;
    double ep = -fSin28w;
    xsec = TMath::Power(em,2) + TMath::Power(ep*(1-y),2) - ep*em*me*y/Ev;
    xsec *= A;
  }

  // numubar/nutaubar + e- -> numubar/nutaubar + e- [NC]
  if( (pdg::IsAntiNuMu(inu)||pdg::IsAntiNuTau(inu)) && proc_info.IsWeakNC() )
  {
    double em = 0.5 - fSin28w;
    double ep = -fSin28w;
    xsec = TMath::Power(ep,2) + TMath::Power(em*(1-y),2) - ep*em*me*y/Ev;
    xsec *= A;
  }

  // numu/nutau + e- -> l- + nu_e [CC}
  if( (pdg::IsNuMu(inu)||pdg::IsNuTau(inu)) && proc_info.IsWeakCC() ) xsec=0;
/*
    double ml  = (pdg::IsNuMu(inu)) ? kMuonMass : kTauMass;
    double ml2 = TMath::Power(ml,2);
    xsec = (kGF2*s/kPi)*(1-ml2/s);
    xsec = TMath::Max(0.,xsec); // if s<ml2 => xsec<0 : force to xsec=0
*/

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("Elastic", pDEBUG)
    << "*** dxsec(ve-)/dy [free e-](Ev="<< Ev << ", y= "<< y<< ") = "<< xsec;
#endif

  //----- The algorithm computes dxsec/dy
  //      Check whether variable tranformation is needed
  if(kps!=kPSyfE) {
    double J = utils::kinematics::Jacobian(interaction,kPSyfE,kps);
    xsec *= J;
  }

  //----- If requested return the free electron xsec even for nuclear target
  if( interaction->TestBit(kIAssumeFreeElectron) ) return xsec;

  //----- Scale for the number of scattering centers at the target
  int Ne = init_state.Tgt().Z(); // num of scattering centers
  xsec *= Ne;

  return xsec;
}
//____________________________________________________________________________
double NuElectronPXSec::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool NuElectronPXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;
  return true;
}
//____________________________________________________________________________
bool NuElectronPXSec::ValidKinematics(const Interaction* interaction) const
{
  if(interaction->TestBit(kISkipKinematicChk)) return true;
  return true;
}
//____________________________________________________________________________
void NuElectronPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void NuElectronPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void NuElectronPXSec::LoadConfig(void)
{
  // weinberg angle
  double thw ;
  GetParam( "WeinbergAngle", thw ) ;
  fSin28w = TMath::Power(TMath::Sin(thw), 2);
  fSin48w = TMath::Power(TMath::Sin(thw), 4);

  // load XSec Integrator
  fXSecIntegrator =
      dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);
}
//____________________________________________________________________________

