//____________________________________________________________________________
/*!

\class    genie::NuElectronPXSec

\brief    nu/nubar + e- scattering differential cross section (dxsec/dy) \n

          The cross section algorithm handles:
             - nue/nuebar + e- -> nue/nuebar + e- [CC + NC + interference]
             - numu/nutau + e- -> numu/nutau + e- [NC]
             - numubar/nutaubar + e- -> numubar/nutaubar + e- [NC]
             - numu/nutau + e- -> l- + nu_e [CC]

          NuElectronPXSec is a concrete implementation of the XSecAlgorithmI
          interface. \n

\ref      W.Greiner and B.Muller, Gauge Theory of Weak  Interactions, Springer
          F.Halzen and A.Martin, Quarks and Leptons, J.Wiley

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  February 10, 2006

*/
//____________________________________________________________________________

#include "Conventions/Constants.h"
#include "Conventions/RefFrame.h"
#include "Elastic/NuElectronPXSec.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGUtils.h"

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
double NuElectronPXSec::XSec(const Interaction * interaction) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  //----- get initial state & kinematics
  const InitialState & init_state = interaction -> GetInitialState ();
  const Kinematics &   kinematics = interaction -> GetKinematics   ();
  const ProcessInfo &  proc_info  = interaction -> GetProcessInfo  ();

  double E     = init_state.GetProbeE(kRfLab);
  double s     = 2*kElectronMass*E;
  double cv    = fCv;
  double ca    = fCa;
  double y     = kinematics.y();
  double ydep  = TMath::Power(1.-y,2.);
  double C     = 0.25 * (kGF_2*s/kPi);
  double m2    = kElectronMass_2;
  double mterm = -0.5*s*kGF_2*m2*y*(cv*cv-ca*ca)/kPi; // small if m2/s<<1

  double xsec = 0; // <-- dxsec/dy

  int inu = init_state.GetProbePDGCode();

  // nue + e- -> nue + e- [CC + NC + interference]
  if(pdg::IsNuE(inu))
    xsec = C*(TMath::Power(cv+ca+2,2) + TMath::Power(cv-ca,2) * ydep) + mterm;

  // nuebar + e- -> nue + e- [CC + NC + interference]
  if(pdg::IsAntiNuE(inu))
    xsec = C*(TMath::Power(cv-ca,2) + TMath::Power(cv+ca+2,2) * ydep) + mterm;

  // numu/nutau + e- -> numu/nutau + e- [NC]
  if( (pdg::IsNuMu(inu)||pdg::IsNuTau(inu)) && proc_info.IsWeakNC() )
    xsec = C*(TMath::Power(cv+ca,2) + TMath::Power(cv-ca,2) * ydep) + mterm;

  // numubar/nutaubar + e- -> numubar/nutaubar + e- [NC]
  if( (pdg::IsAntiNuMu(inu)||pdg::IsAntiNuTau(inu)) && proc_info.IsWeakNC() )
    xsec = C*(TMath::Power(cv-ca,2) + TMath::Power(cv+ca,2) * ydep) + mterm;

  // numu/nutau + e- -> l- + nu_e [CC}
  if( (pdg::IsNuMu(inu)||pdg::IsNuTau(inu)) && proc_info.IsWeakCC() ) {
    double ml  = (pdg::IsNuMu(inu)) ? kMuonMass : kTauMass;
    double ml2 = TMath::Power(ml,2);
    xsec = (kGF_2*s/kPi)*(1-ml2/s);
    xsec = TMath::Max(0.,xsec); // if s<ml2 => xsec<0 : force to xsec=0
  }

  LOG("Elastic", pDEBUG)
         << "*** dxsec[ve-]/dy (E=" << E << ", y = " << y << ") = " << xsec;
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
// Reads its configuration from its Registry

  fCv = fConfig->GetDoubleDef("cv", -0.52);
  fCa = fConfig->GetDoubleDef("ca",  0.06);
}
//____________________________________________________________________________

