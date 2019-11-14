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
#include "Physics/BoostedDarkMatter/XSection/DMElectronPXSec.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/KineUtils.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;

//____________________________________________________________________________
DMElectronPXSec::DMElectronPXSec() :
XSecAlgorithmI("genie::DMElectronPXSec")
{

}
//____________________________________________________________________________
DMElectronPXSec::DMElectronPXSec(string config) :
XSecAlgorithmI("genie::DMElectronPXSec", config)
{

}
//____________________________________________________________________________
DMElectronPXSec::~DMElectronPXSec()
{

}
//____________________________________________________________________________
double DMElectronPXSec::XSec(
                 const Interaction * interaction, KinePhaseSpace_t kps) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  //----- get initial state & kinematics
  const InitialState & init_state = interaction -> InitState();
  const Kinematics &   kinematics = interaction -> Kine();
  const ProcessInfo &  proc_info  = interaction -> ProcInfo();

  double Ev = init_state.ProbeE(kRfLab);
  double Ev2 = TMath::Power(Ev,2);
  double me = kElectronMass;
  double me2 = TMath::Power(me,2);
  double ml = interaction->FSPrimLepton()->Mass();
  double ml2 = TMath::Power(ml,2);
  double y  = kinematics.y();
  double MZ2 = TMath::Power(fMedMass,2);
  double A  = fgZp4 * TMath::Power(Ev,3) * me / (4.0 * kPi * (Ev2 - ml2) * TMath::Power(MZ2 + 2.0*Ev*me*(1.0 - y),2));

  //  y = 1 - me/Ev - y; // FSPL = electron. XSec below are expressed in Marciano's y!
  //  if(y > 1/(1+0.5*me/Ev)) return 0;
  //  if(y < 0) return 0;

  double QeV2 = TMath::Power(0.5*(fQeL+fQeR),2);
  double QeA2 = TMath::Power(0.5*(-fQeL+fQeR),2);
  double QeVA = 0.25*(TMath::Power(fQeR,2) - TMath::Power(fQeL,2));
    
  double xsec = 0; // <-- dxsec/dy

  int inu = init_state.ProbePdg();

  if( proc_info.IsDarkMatter() && fVelMode == 0 ) 
  {
    double QdmV2 = TMath::Power(0.5*(fQdmL+fQdmR),2);
    double QdmA2 = TMath::Power(0.5*(-fQdmL+fQdmR),2);
    double QdmVA = 0.25*(TMath::Power(fQdmR,2) - TMath::Power(fQdmL,2));
    double T1 = 1 + TMath::Power(y,2); 
    double T2 = (1.0 - y)* me / Ev;
    double T3 = (1.0 - y)*ml2 / Ev / me;
    double T4 = 2.0 * ml2 / Ev2;
    double TL = 2.0*TMath::Power(2.0*(1.0 - y) * me * ml / MZ2 + ml/Ev,2);
    xsec  = (T1 - T2 - T3)*QdmV2*QeV2;
    xsec += (T1 - T2 + T3 - T4)*QdmA2*QeV2;
    xsec += (T1 + T2 - T3 - T4)*QdmV2*QeA2;
    xsec += (T1 + T2 + T3 + T4 + TL)*QdmA2*QeA2;
    if ( pdg::IsDarkMatter(inu) ) {
      xsec += 4.0 * (1.0 - TMath::Power(y,2))*QdmVA*QeVA;
    }
    else {
      xsec -= 4.0 * (1.0 - TMath::Power(y,2))*QdmVA*QeVA;
    }
    xsec *= A;
  }

  if( proc_info.IsDarkMatter() && fVelMode == 2 ) 
  {
    double QdmS2 = TMath::Power(fQdmS,2);
    double T1 = 2.0 * y;
    double T2 = (1.0 - y)*me2 / Ev / me;
    double T3 = (1.0 - y)*ml2 / Ev / me;
    double T4 = 2.0 * ml2 / Ev2;
    xsec  = (T1 - T3)*QeV2*QdmS2;
    xsec += (T1 - T2 - T3 - T4)*QeA2*QdmS2;
    xsec *= A;
  }

  #ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("Elastic", pDEBUG)
    << "*** dxsec(ve-)/dy [free e-](Ev="<< Ev << ", y= "<< y<< ") = "<< xsec;
  #endif

  //----- The algorithm computes dxsec/dy
  //      Check whether variable tranformation is needed
  if(kps!=kPSyfE) {
    double J = utils::kinematics::Jacobian(interaction,kPSyfE,kps);
    LOG("Elastic", pDEBUG) << "Multiplying by jacobian " << J;
    xsec *= J;
  }

  //----- If requested return the free electron xsec even for nuclear target
  if( interaction->TestBit(kIAssumeFreeElectron) ) return xsec;

  //----- Scale for the number of scattering centers at the target
  int Ne = init_state.Tgt().Z(); // num of scattering centers
  xsec *= Ne;
  LOG("Elastic", pDEBUG) << "Multiplying by Ne " << Ne;  

  return xsec;
}
//____________________________________________________________________________
double DMElectronPXSec::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool DMElectronPXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;
  return true;
}
//____________________________________________________________________________
bool DMElectronPXSec::ValidKinematics(const Interaction* interaction) const
{
  if(interaction->TestBit(kISkipKinematicChk)) return true;
  return true;
}
//____________________________________________________________________________
void DMElectronPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DMElectronPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DMElectronPXSec::LoadConfig(void)
{
  // Dark matter couplings
  double gZp;
  GetParam("ZpCoupling", gZp);
  GetParam("DarkLeftCharge", fQdmL);
  GetParam("DarkRightCharge", fQdmR);
  GetParam("DarkScalarCharge", fQdmS);
  GetParam("ElectronLeftCharge", fQeL);
  GetParam("ElectronRightCharge", fQeR);

  fgZp4 = TMath::Power(gZp, 4);
  
  // Mediator mass
  fMedMass = PDGLibrary::Instance()->Find(kPdgMediator)->Mass();

  // velocity dependence of interaction
  GetParamDef("velocity-mode", fVelMode, 0 );

  // load XSec Integrator
  fXSecIntegrator =
      dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);
}
//____________________________________________________________________________

