//____________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory
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
#include "Physics/GlashowResonance/XSection/GLRESPXSec.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
GLRESPXSec::GLRESPXSec() :
XSecAlgorithmI("genie::GLRESPXSec")
{

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
  if(! this -> ValidKinematics (interaction) ) return 0.;

  const InitialState & init_state = interaction -> InitState();

  double E     = init_state.ProbeE(kRfLab);
  double s     = 2 * kElectronMass * E + kElectronMass2;

  // The W limit is because hadronization might have issues at low W (as in PYTHIA6).
  // To be consistent the cross section must be computed in the W range where the events are generated.
  if ( TMath::Sqrt(s)<fWmin ) return 0.;

  const Kinematics &   kinematics = interaction -> Kine();
  const XclsTag &      xclstag    = interaction -> ExclTag();
  int  final_lepton = xclstag.FinalLeptonPdg();

  double y     = kinematics.y();
  double s0    = 2 * kElectronMass * E;

  double Gw      = PDGLibrary::Instance()->Find(kPdgWM)->Width(); //PDGLibrary::Instance()->Find(kPdgWM)->Width() //genhen=2.124
  double gf      = kGF2/kPi;                                      //kGF2/kPi                                      //genhen=1.16637e-05*1.16637e-05/3.141592653589793
  double m_w     = kMw;                                           //kMw                                           //genhen=80.425
  double m_z     = kMz;                                           //kMz                                           //genhen=91.1876
  double Sin2thw = 1 - kMw2 / kMz2;                               //1 - kMw2 / kMz2                               //genhen=0.2312
  double branch  = 64.41/10.63;                                   //branching ratio hadrons/muons                 //genhen=67.96/10.57 

  double Gw2     = TMath::Power(Gw,  2);
  double m_w2    = TMath::Power(m_w,2); 
  double m_z2    = TMath::Power(m_z,2); 
  double Sin4thw = TMath::Power(Sin2thw,2);

  double prop =  TMath::Power(1-s/m_w2, 2) + Gw2/m_w2;

  double xsec  = 0;
  if      ( pdg::IsMuon(final_lepton) )     {
    double ml2 = kMuonMass2;
    xsec = s0*TMath::Power(1-y, 2) + (3*kElectronMass2+ml2)*(1-y) + kElectronMass*(kElectronMass2+ml2)/E;
    xsec *= gf/prop; 
  }
  else if ( pdg::IsTau(final_lepton) )      {
    double ml2 = kTauMass2;
    xsec = s0*TMath::Power(1-y, 2) + (3*kElectronMass2+ml2)*(1-y) + kElectronMass*(kElectronMass2+ml2)/E;
    xsec *= gf/prop; 
  }
  else if ( pdg::IsPion(final_lepton) )     {
    double ml2 = kMuonMass2;
    xsec = ( s0*TMath::Power(1-y, 2) + (3*kElectronMass2+ml2)*(1-y) + kElectronMass*(kElectronMass2+ml2)/E ) * branch;
    xsec *= gf/prop; 
  }
  else if ( pdg::IsElectron(final_lepton) ) {
    double u  = 2 * kElectronMass2 - 2 * kElectronMass * E * y;
    double t  = 2 * kElectronMass2 - s - u;
    double L  = Sin2thw - 1./2.;

    double x1 = L*m_z2*(s-m_w2) + m_w2*(u-m_z2);
    double y1 = L * m_z2 * Gw * m_w;
    double x2 = (u - m_z2) * (s - m_w2);
    double y2 = (u - m_z2) * Gw * m_w;
    double y3 = ( x1*x2 + y1*y2 ) / ( TMath::Power(x2,2) + TMath::Power(y2,2) );
    double x3 = ( x2*y1 - x1*y2 ) / ( TMath::Power(x2,2) + TMath::Power(y2,2) );

    xsec = gf * s0 * ( Sin4thw/TMath::Power(u/m_z2-1,2) + (TMath::Power(x3,2)+TMath::Power(y3,2))*TMath::Power(t-kElectronMass2,2)/TMath::Power(s0,2) );
  }

  
  LOG("GLRESPXSec", pINFO) << "dxsec/dy (E= " << E << ", y= " << y << ") = " << xsec;

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
double GLRESPXSec::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);

  return xsec;
}
//____________________________________________________________________________
bool GLRESPXSec::ValidProcess(const Interaction* interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const InitialState & init_state = interaction->InitState();
  const ProcessInfo &  proc_info  = interaction->ProcInfo();

  bool nuok    = pdg::IsAntiNuE(init_state.ProbePdg());
  bool nucok   = !(init_state.Tgt().HitNucIsSet());
  bool ccprcok = proc_info.IsWeakCC();

  if ( !nuok    ) return false;
  if ( !nucok   ) return false;
  if ( !ccprcok ) return false;

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