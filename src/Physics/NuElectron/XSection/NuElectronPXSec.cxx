//____________________________________________________________________________
/*
 Copyright (c) 2003-2022, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory

 Changes required to implement the Electron Velocity module
 were installed by Brinden Carlson (Univ. of Florida)
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
#include "Physics/NuclearState/ElectronVelocity.h"

#include "TFile.h"
#include "TGraph.h"
#include <fstream>
#include <iterator>

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
  // std::ofstream output_file("./sigmas_errors.txt");
  // output_file <<"x\tsigma\terr\tavg\n";
  // vector<double> x;
  // for (unsigned int i = 1; i<=fSigmas.size();i++){
  //   x.push_back(i);
  //   output_file <<i<<"\t"<<fSigmas[i]<<"\t"<<
  //   fErrors[i]<<"\t"<<fAvg[i]<<"\n"; //Put beam eng in here
  // }
  //TGraph sigma(fSigmas.size(),x.data(),fSigmas.data());
  //TGraph errors(fErrors.size(),x.data(),fErrors.data());

  //sigma.Write("sigma");
  //errors.Write("errors");
  //f.Close();

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

  double Ev = init_state.ProbeE(kRfHitElRest); //Electron rest frame

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
  double xsec_sum = 0; 
  double xsec_sum2 = 0; 
  int NInt = 0; //Count number of integrations
  do{
    NInt++;
    Interaction in_curr(*interaction); //Copy interaction object
    fElectronVelocity->InitializeVelocity(in_curr); //Modify interaction to give electron random velocity from selected distribution
    double xsec = fXSecIntegrator->Integrate(this,&in_curr); // In ele at rest
    //get beta - comps orthogonal to beam x,y
    //scale = sqrt(1-b_t^2)
    TVector3 beta = in_curr.InitState().Tgt().HitEleP4().BoostVector(); // beta
    double beta_tangent = sqrt(TMath::Power(beta[0],2)+TMath::Power(beta[1],2)); //Component tangential to beam
    xsec *= sqrt(1-TMath::Power(beta_tangent,2)); //Correct for lorentz factor
    xsec_sum+=xsec;
    xsec_sum2+=TMath::Power(xsec,2);
    double xsec_mean = xsec_sum/NInt;
    //var = (sum(xi^2)/N-xsec_mean^2)
    //rel_err = sigma/sqrt(n)*mean
    //double xsec_sigma = sqrt(xsec_sum2/NInt-TMath::Power(xsec_mean,2));
    double xsec_err = sqrt((xsec_sum2/NInt-TMath::Power(xsec_mean,2))/NInt)/xsec_mean;
    // fSigmas.push_back(xsec_sigma);
    // fErrors.push_back(xsec_err);
    // fAvg.push_back(xsec_mean);
    if (NInt > 1 && xsec_err < fErrTolerance){break;} //Break condition for dipping below set tolerance
  }
  while ( NInt < fNIntegration); 
  double xsec_avg = xsec_sum/NInt;
  return xsec_avg;
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
  GetParam( "N-Integration-Samples", fNIntegration ) ; //
  GetParam( "NuE-XSecRelError" , fErrTolerance ) ; //
  fSin28w = TMath::Power(TMath::Sin(thw), 2);
  fSin48w = TMath::Power(TMath::Sin(thw), 4);

  //XSecOnElectron::LoadConfig(); //Gets fNInt, fErr, fXSec, fEle - only keep fSin in lower level modules

  // load XSec Integrator
  fXSecIntegrator =
      dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator")); //
  fElectronVelocity =
      dynamic_cast<const ElectronVelocity *> (this->SubAlg("Electron-Velocity")); //
  assert(fXSecIntegrator);
}
//____________________________________________________________________________
