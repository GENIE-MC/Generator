//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - February 10, 2006

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "Algorithm/AlgConfigPool.h"
#include "Conventions/Constants.h"
#include "Conventions/RefFrame.h"
#include "Elastic/NuElectronPXSec.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGUtils.h"
#include "Utils/KineUtils.h"

using namespace genie;
using namespace genie::constants;

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

  double E     = init_state.ProbeE(kRfLab);
  double s     = 2*kElectronMass*E;
  double y     = kinematics.y();
  double ydep  = TMath::Power(1.-y,2.);
  double C     = 0.25 * (kGF2*s/kPi);
  double m2    = kElectronMass2;
  double mterm = -0.5*s*kGF2*m2*y*(fCv*fCv-fCa*fCa)/kPi; // small if m2/s<<1

  double xsec = 0; // <-- dxsec/dy

  int inu = init_state.ProbePdg();

  // nue + e- -> nue + e- [CC + NC + interference]
  if(pdg::IsNuE(inu))
    xsec = C*(TMath::Power(fCv+fCa+2,2) + 
                                   TMath::Power(fCv-fCa,2) * ydep) + mterm;

  // nuebar + e- -> nue + e- [CC + NC + interference]
  if(pdg::IsAntiNuE(inu))
    xsec = C*(TMath::Power(fCv-fCa,2) + 
                                 TMath::Power(fCv+fCa+2,2) * ydep) + mterm;

  // numu/nutau + e- -> numu/nutau + e- [NC]
  if( (pdg::IsNuMu(inu)||pdg::IsNuTau(inu)) && proc_info.IsWeakNC() )
    xsec = C*(TMath::Power(fCv+fCa,2) + 
                                   TMath::Power(fCv-fCa,2) * ydep) + mterm;

  // numubar/nutaubar + e- -> numubar/nutaubar + e- [NC]
  if( (pdg::IsAntiNuMu(inu)||pdg::IsAntiNuTau(inu)) && proc_info.IsWeakNC() )
    xsec = C*(TMath::Power(fCv-fCa,2) + 
                                   TMath::Power(fCv+fCa,2) * ydep) + mterm;

  // numu/nutau + e- -> l- + nu_e [CC}
  if( (pdg::IsNuMu(inu)||pdg::IsNuTau(inu)) && proc_info.IsWeakCC() ) {
    double ml  = (pdg::IsNuMu(inu)) ? kMuonMass : kTauMass;
    double ml2 = TMath::Power(ml,2);
    xsec = (kGF2*s/kPi)*(1-ml2/s);
    xsec = TMath::Max(0.,xsec); // if s<ml2 => xsec<0 : force to xsec=0
  }

  LOG("Elastic", pDEBUG)
     << "*** dxsec(ve-)/dy [free e-](E="<< E << ", y= "<< y<< ") = "<< xsec;

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
  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();

  fCv = fConfig->GetDoubleDef("cv", gc->GetDouble("NuElecEL-CV"));
  fCa = fConfig->GetDoubleDef("ca", gc->GetDouble("NuElecEL-CA"));
}
//____________________________________________________________________________

