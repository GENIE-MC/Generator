//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Rosen Matev (r.matev@gmail.com)

 For the class documentation see the corresponding header file.

 Important revisions after version 2.8.0 :

*/
//____________________________________________________________________________

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Controls.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/RefFrame.h"
#include "Physics/NuElectron/XSection/IMDAnnihilationPXSec.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;

//____________________________________________________________________________
IMDAnnihilationPXSec::IMDAnnihilationPXSec() :
XSecAlgorithmI("genie::IMDAnnihilationPXSec")
{

}
//____________________________________________________________________________
IMDAnnihilationPXSec::IMDAnnihilationPXSec(string config) :
XSecAlgorithmI("genie::IMDAnnihilationPXSec", config)
{

}
//____________________________________________________________________________
IMDAnnihilationPXSec::~IMDAnnihilationPXSec()
{

}
//____________________________________________________________________________
double IMDAnnihilationPXSec::XSec(
                 const Interaction * interaction, KinePhaseSpace_t kps) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  //----- get initial state & kinematics
  const InitialState & init_state = interaction -> InitState();
  const Kinematics &   kinematics = interaction -> Kine();

  double Ev = init_state.ProbeE(kRfLab);
  double twoMeEv = 2*kElectronMass*Ev;
  double ymax = 1 - kMuonMass2/(twoMeEv + kElectronMass2);
  double y  = kinematics.y();
  double A  = kGF2/kPi;

  LOG("IMDAnnihilation", pDEBUG) << "Ev = " << Ev << ", y = " << y << ", ymax = " << ymax;
  //Note: y = (Ev-El)/Ev but in Marciano's paper y=(El-(m_l^2+m_e^2)/2m_e)/Ev.
  y = 1 - y - (kMuonMass2 + kElectronMass2)/twoMeEv;

  if(y > ymax) return 0;
  if(y < 0) return 0;

  double xsec = A*(twoMeEv*TMath::Power(1-y,2) - (kMuonMass2 - kElectronMass2)*(1-y)); // <-- dxsec/dy

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("IMDAnnihilation", pDEBUG)
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
double IMDAnnihilationPXSec::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool IMDAnnihilationPXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;
  return true;
}
//____________________________________________________________________________
bool IMDAnnihilationPXSec::ValidKinematics(const Interaction* interaction) const
{
  if(interaction->TestBit(kISkipKinematicChk)) return true;
  return true;
}
//____________________________________________________________________________
void IMDAnnihilationPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void IMDAnnihilationPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void IMDAnnihilationPXSec::LoadConfig(void)
{
  // load XSec Integrator
  fXSecIntegrator =
      dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);
}
//____________________________________________________________________________

