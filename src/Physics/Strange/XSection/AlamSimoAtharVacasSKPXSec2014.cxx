//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Authors: Chris Marshall and Martti Nirkko

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Numerical/MathUtils.h"
#include "Physics/Strange/XSection/AlamSimoAtharVacasSKPXSec2014.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Physics/NuclearState/NuclearUtils.h"

using namespace genie;
using namespace genie::utils;
using namespace genie::constants;

//____________________________________________________________________________
AlamSimoAtharVacasSKPXSec2014::AlamSimoAtharVacasSKPXSec2014() :
XSecAlgorithmI("genie::AlamSimoAtharVacasSKPXSec2014")
{

}
//____________________________________________________________________________
AlamSimoAtharVacasSKPXSec2014::AlamSimoAtharVacasSKPXSec2014(string config) :
XSecAlgorithmI("genie::AlamSimoAtharVacasSKPXSec2014", config)
{

}
//____________________________________________________________________________
AlamSimoAtharVacasSKPXSec2014::~AlamSimoAtharVacasSKPXSec2014()
{

}
//____________________________________________________________________________
double AlamSimoAtharVacasSKPXSec2014::XSec(
              const Interaction * interaction, KinePhaseSpace_t /*kps*/) const
{
  // Check whether interaction is valid
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;
  
  const InitialState & init_state = interaction -> InitState();
  const Kinematics &   kinematics = interaction -> Kine();
  Target * target  = interaction->InitStatePtr()->TgtPtr();
  int nTargetProtons = target->Z();
  int nTargetNeutrons = target->N();
  
  // Initialisation begins here
  // --------------------------
  // mutable variables used in the matrix element calculation
  Enu = init_state.ProbeE(kRfHitNucRest); // use the struck nucleon rest frame
  leptonPDG = interaction->FSPrimLeptonPdg();
  
  // Fill reaction type (NN=1, NP=2, PP=3)
  reactionType = 0; // mutable
  int numProtons  = interaction->ExclTagPtr()->NProtons();
  int numNeutrons = interaction->ExclTagPtr()->NNeutrons();
  int kaonPDG     = interaction->ExclTagPtr()->StrangeHadronPdg();
  if      (numProtons==0 && numNeutrons==1 && kaonPDG==kPdgKP) reactionType=1;
  else if (numProtons==1 && numNeutrons==0 && kaonPDG==kPdgK0) reactionType=2;
  else if (numProtons==1 && numNeutrons==0 && kaonPDG==kPdgKP) reactionType=3;
  else {  
    LOG("SKXSec", pERROR) 
      << "Reaction not defined! This should NEVER happen!"
      << "\n - numProtons  = " << numProtons
      << "\n - numNeutrons = " << numNeutrons
      << "\n - kaonPDG     = " << kaonPDG;
    return 0.;
  }
  
  double Tlep = kinematics.GetKV(kKVTl);
  double Tkaon = kinematics.GetKV(kKVTk);
  double costheta = kinematics.GetKV(kKVctl);
  double phikq = kinematics.GetKV(kKVphikq);

  // Set lepton mass
  aml = PDGLibrary::Instance()->Find(leptonPDG)->Mass(); // mutable

  double theta = TMath::ACos(costheta);
  
  // Set reaction parameters, which are mutables used in the matrix element calculations
  if (reactionType == 1) {
    amSig = PDGLibrary::Instance()->Find(kPdgSigmaM)->Mass();
    amk   = PDGLibrary::Instance()->Find(kPdgKP)->Mass();
    ampi  = kPi0Mass;
    am    = kNeutronMass;
  }
  else if (reactionType == 2) {
    amSig = PDGLibrary::Instance()->Find(kPdgSigma0)->Mass();
    amk   = PDGLibrary::Instance()->Find(kPdgK0)->Mass();
    ampi  = kPionMass;
    am    = kNeutronMass;
  }
  else if (reactionType == 3) {
    amSig = PDGLibrary::Instance()->Find(kPdgSigma0)->Mass();
    amk   = PDGLibrary::Instance()->Find(kPdgKP)->Mass();
    ampi  = kPi0Mass;
    am    = kProtonMass;
  }
  else {
    return 0.;
  }
  
  // Cross section calculation begins here
  // -------------------------------------
  double tkmax = Enu - amk - aml - Tlep;        // maximal allowed kaon energy
  if (Tkaon > tkmax) return 0.;

  Ekaon = Tkaon+amk; // mutable
  pkvec = sqrt(Ekaon*Ekaon-amk*amk); // mutable
  
  double xsec = 0.;
  
  Elep = Tlep + aml; // mutable
  alepvec = sqrt(Elep*Elep - aml*aml); // mutable
  aq0 = Enu-Elep; // mutable
  double a1 = aq0+am-Ekaon; // mutable

  // aqvec is the magnitude of the three-momentum transfer to hadron system
  aqvec = sqrt(alepvec*alepvec+Enu*Enu-2.0*Enu*alepvec*costheta);

  // this check is basically the longitudinal component of the kaon momentum (in the momentum transfer frame) divided by the total
  // if it is larger than 1, the kinematics are non-physical so we should return zero
  double check = (aqvec*aqvec+pkvec*pkvec+am*am-a1*a1)/(2.0*aqvec*pkvec);
  
  double amat2; // the matrix element
  if (fabs(check) <= 1.0) { // so it has to be smaller than 1, but it could be negative if the kaon backscatters in com frame
    angkq = check;
    if      (reactionType == 1) amat2 = this->Amatrix_NN(theta, phikq);
    else if (reactionType == 2) amat2 = this->Amatrix_NP(theta, phikq);
    else if (reactionType == 3) amat2 = this->Amatrix_PP(theta, phikq);
    else    return 0.;
    xsec = alepvec*alepvec*amat2/(32.0*pow(2.0*kPi,4)*am*Enu*Elep*aqvec); 
  }
  else {
    xsec = 0.;
  }

  // the matrix element calculation is for d4sigma/dtk dpl dcosthetal dphi_kq
  // we have T_l instead of p_l so we multiply by dp/dT = E/p
  xsec *= Elep / alepvec;
  
  // xsec is now the nucleon-level cross section
  // There are no fancy nuclear effects for this model, so the nucleus cross section is just 
  // the nucleon XS times the number of nucleons of the appropriate isospin for the process selected
  if( reactionType == 1 || reactionType == 2 ) xsec *= nTargetNeutrons; // NN or NP
  else xsec *= nTargetProtons; // PP

  return xsec;
}
//____________________________________________________________________________
double AlamSimoAtharVacasSKPXSec2014::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool AlamSimoAtharVacasSKPXSec2014::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const ProcessInfo & proc_info = interaction->ProcInfo();  
  if(!proc_info.IsSingleKaon() || !proc_info.IsWeakCC()) return false;

  return true;
}
//____________________________________________________________________________
void AlamSimoAtharVacasSKPXSec2014::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void AlamSimoAtharVacasSKPXSec2014::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void AlamSimoAtharVacasSKPXSec2014::LoadConfig(void)
{
  fXSecIntegrator =
      dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);
  
  amLam = PDGLibrary::Instance()->Find(kPdgLambda)->Mass();
  am = kNeutronMass; // this will be nucleon mass, set event by event
  amEta = PDGLibrary::Instance()->Find(kPdgEta)->Mass();

  GetParam( "CKM-Vus", Vus ) ;
  // fpi is 0.0924 in Athar's code, use the same one that is already in UserPhysicsOptions
  GetParam( "PionDecayConstant",fpi ) ;

  GetParam( "SU3-D", d ) ;
  GetParam( "SU3-F", f ) ;
  g = kGF;                              // Fermi coupling
  // we really want the anomolous moment, but the one in UserPhysicsOptions is the full moment, despite the name
  GetParam( "AnomMagnMoment-P", amup ) ;
  amup -= 1 ;
  GetParam( "AnomMagnMoment-N", amun ) ;

  Fm1 = -(amup+2.0*amun)/(2.0*am);
  Fm2 = -3.0*amup/(2.0*am);

}
//____________________________________________________________________________





























    // NOTE: Please to not continue reading.





























    // WARNING: Go back, or you may be eaten by raptors.





























    // DANGER: HERE BE DRAGONS




























    
// *****************************************************************************
// The following are lengthy matrix calculations that must not be tampered with.    
// *****************************************************************************

//____________________________________________________________________________
double AlamSimoAtharVacasSKPXSec2014::Amatrix_NN(double theta, double phikq) const
{
  
  double sol = 0.;
  
  double akk1=0., zdotq=0., qdotpk=0., akcrosk1=0., qcrospk=0.;
  double zdotpk=0., azpk=0., aqkaon=0., akpk=0., apkk1=0.;
  double C1=0., C2=0., C3=0., C4=0., /*C5=0., C6=0.,*/ C7=0., C8=0., C9=0.;
  double aq2=0., gform=0., con=0., t1=0., t2=0., t3=0., /*t4=0.,*/ t5=0., t6=0.;
  
  akk1=Enu*Elep-Enu*alepvec*cos(theta);
  zdotq=(Enu*Enu-alepvec*alepvec)/2.0;
  qdotpk=aqvec*pkvec*angkq;
  akcrosk1 = Enu*alepvec*sin(theta);
  qcrospk=aqvec*pkvec*sqrt(1.0-angkq*angkq);
  zdotpk=(akcrosk1*qcrospk*cos(phikq)+zdotq*qdotpk)/(aqvec*aqvec);
  azpk=Ekaon*(Enu+Elep)/2.0-zdotpk;
  aqkaon=aq0*Ekaon-aqvec*pkvec*angkq;
  akpk=azpk + aqkaon/2.0;
  apkk1=azpk - aqkaon/2.0;
                     
  C1=1.0/(am*am+amk*amk-2.0*am*Ekaon-amSig*amSig);
  C2=d+f;
  C3=1./(aml*aml-2.0*akk1-amk*amk);
  C4=d-f;
  //C5=d+3.*f;
  //C6=1.0/(am*am+amk*amk-2.0*am*Ekaon-amLam*amLam);       
  C7=2.0*am/(aml*aml-2.0*akk1+amk*amk-2.*aqkaon-ampi*ampi);
  C8=2.0*am/(aml*aml-2.0*akk1+amk*amk-2.*aqkaon-amEta*amEta);
  C9=d - 3.*f;
  aq2=aml*aml-2.0*akk1;
  gform=1.0/pow(1.0-aq2/(1.0*1.0),4);
  
  con=g*g*Vus*Vus/(4.0*fpi*fpi);

  t1=1.0;
  t2=1.0;
  t3=1.0;
  //t4=1.0;
  t5=1.0;
  t6=1.0;
  
  sol = con*gform*(4.*am*(-2.*(akk1*akk1)*
           ((aml*aml)*(C3*C3)*
           (16.*(am*am)*(C1*C1)*(C4*C4*C4*C4)*(Ekaon*Ekaon*Ekaon)*(t3*t3) + 
           4.*(am + amSig)*C1*(C4*C4)*t3*
           ((-1.*akpk - 2.*(amk*amk) + apkk1)*t2 + 
           2.*(amk*amk)*
           (akpk + (amk*amk) + am*(am + amSig) - 
           1.*apkk1)*C1*(C4*C4)*t3) - 
           8.*am*C1*(C4*C4)*(Ekaon*Ekaon)*t3*
           (-2.*t2 + 
           C1*(C4*C4)*
           ((am*am) + 2.*(amk*amk) + (amSig*amSig) + 
           2.*am*(amSig + Elep - 1.*Enu))*t3) + 
           Enu*((t2*t2) - 4.*(amk*amk)*C1*(C4*C4)*t2*t3 - 
           4.*(amk*amk)*(am - 1.*amk + amSig)*
           (am + amk + amSig)*(C1*C1)*(C4*C4*C4*C4)*(t3*t3)) + 
           Elep*(-1.*(t2*t2) + 4.*(amk*amk)*C1*(C4*C4)*t2*t3 + 
           4.*(amk*amk)*(am - 1.*amk + amSig)*
           (am + amk + amSig)*(C1*C1)*(C4*C4*C4*C4)*(t3*t3)) + 
           Ekaon*(3.*(t2*t2) - 
           8.*C1*(C4*C4)*((amk*amk) + am*(Elep - 1.*Enu))*
           t2*t3 + 
           4.*(C1*C1)*(C4*C4*C4*C4)*
           ((amk*amk*amk*amk) + (amk*amk)*(amSig*amSig) + 
           2.*akpk*(-1.*(am*am) + (amSig*amSig)) - 
           2.*(amSig*amSig)*apkk1 + 
           (am*am)*(-3.*(amk*amk) + 2.*apkk1) - 
           2.*am*(amk*amk)*(amSig - 2.*Elep + 2.*Enu))
           *(t3*t3))) - 
           8.*C1*C4*Fm1*t3*
           (-1.*C4*((amk*amk*amk*amk)*C1*
           (-2. + (am - 1.*Elep + Enu)*Fm1)*t3 + 
           (amk*amk)*
           (2.*t1 + 
           C1*
           ((am*am*am)*Fm1 + 
           am*
           (-4.*amSig + 4.*Ekaon + 
           2.*akpk*Fm1 + 
           (-2.*apkk1 + 
           (amSig + 2.*Ekaon)*
           (amSig + 2.*Elep - 2.*Enu))*Fm1) + 
           (am*am)*
           (-2. + 
           (2.*amSig - 2.*Ekaon + Elep - 
           1.*Enu)*Fm1) + 
           amSig*
           (-4.*Ekaon + 
           2.*(akpk - 1.*apkk1)*Fm1 + 
           amSig*(-2. + Elep*Fm1 - 1.*Enu*Fm1)
           ))*t3) + 
           2.*Ekaon*
           ((am*am)*C1*
           (-1.*akpk + apkk1 - 
           2.*Ekaon*(amSig + Elep - 1.*Enu))*
           Fm1*t3 - 
           1.*am*(t1 - 4.*amSig*C1*Ekaon*t3) + 
           amSig*
           (t1 + 
           amSig*(akpk - 1.*apkk1)*C1*Fm1*t3)))
           + ((amk*amk)*(amSig + Ekaon) + 
           am*((amk*amk) - 2.*(Ekaon*Ekaon)))*
           (C2*C7*t5 + C8*C9*t6))) + 
           4.*(Enu*(-2.*am*Elep*
           (-2.*(1. + (C4*C4))*(t1*t1) + (aml*aml)*C3*t1*t2 - 
           2.*C1*(-4. + (aml*aml)*C3)*(C4*C4)*
           ((amk*amk) - 2.*am*Ekaon)*t1*t3 + 
           C1*(C4*C4)*t3*
           (2.*C1*(1. + (C4*C4))*
           (amk*(am + amk + amSig) - 2.*am*Ekaon)*
           (amk*(am - 1.*amk + amSig) + 
           2.*am*Ekaon)*t3 + 
           (aml*aml)*C3*
           (2.*(amk*amk*amk*amk)*C1*(C4*C4)*t3 + 
           2.*am*Ekaon*
           (t2 + 4.*am*C1*(C4*C4)*Ekaon*t3) - 
           1.*(amk*amk)*
           (t2 + 
           2.*C1*(C4*C4)*
           (((am+amSig)*(am+amSig)) + 4.*am*Ekaon)*t3)
           ))) + 
           (aml*aml)*(-2.*((-1.+C4)*(-1.+C4))*(t1*t1) + 
           t1*(C3*((aml*aml) + 2.*am*(Ekaon + Enu))*t2 + 
           4.*am*C1*C4*Ekaon*
           (2. + 2.*(C4*C4) - 1.*(am + amSig)*Fm1 + 
           C4*
           (-4. + (aml*aml)*C3 + 
           2.*am*(C3*Enu + Fm1)))*t3 + 
           (amk*amk)*
           (-1.*C3*t2 - 
           2.*C1*C4*
           (2. + 
           C4*
           (-4. + (aml*aml)*C3 + 2.*C4 + 
           2.*am*C3*Enu + 3.*am*Fm1 + 
           amSig*Fm1))*t3)) + 
           C1*C4*t3*
           (-2.*C1*C4*
           ((amk*amk)*
           (-1.*(am - 1.*amk + amSig)*
           (am + amk + amSig)*((-1.+C4)*(-1.+C4)) - 
           2.*am*
           ((amk*amk) + 2.*((am+amSig)*(am+amSig))*C4)*
           Fm1 - 
           2.*am*(am - 1.*amk + amSig)*
           (am + amk + amSig)*Enu*(Fm1*Fm1)) + 
           4.*(am*am)*(Ekaon*Ekaon)*
           (1. + (C4*C4) + amSig*Fm1 + 
           2.*C4*(-1. + (am + amSig)*Fm1) + 
           am*Fm1*(-1. + 2.*Enu*Fm1)) - 
           2.*am*(amk*amk)*Ekaon*
           (2. + 2.*(C4*C4) + amSig*Fm1 + 
           2.*C4*(-2. + (am + amSig)*Fm1) + 
           am*Fm1*(-3. + 4.*Enu*Fm1)))*t3 + 
           C3*C4*
           (((amk*amk) - 2.*am*Ekaon)*
           ((amk*amk) - 1.*(aml*aml) - 
           2.*am*(Ekaon + Enu))*t2 - 
           2.*C1*(C4*C4)*
           (amk*(am + amk + amSig) - 
           2.*am*Ekaon)*
           (amk*(am - 1.*amk + amSig) + 
           2.*am*Ekaon)*((aml*aml) + 2.*am*Enu)*
           t3) - 
           1.*
           ((amk*amk)*(am + amSig)*(-1. + 2.*C4) - 
           1.*(aml*aml)*((amk*amk) - 2.*am*Ekaon)*Fm1)
           *(C2*C7*t5 + C8*C9*t6))) + 
           apkk1*(-2.*((-1.+C4)*(-1.+C4))*(t1*t1) + 
           2.*C1*C4*
           (-1.*(amk*amk)*
           (2. + C4*(-4. + (aml*aml)*C3 + 2.*C4)) + 
           2.*am*
           ((am + amSig)*C4*
           ((aml*aml)*C3 - 2.*(1. + C4)) + 
           (2. + C4*(-4. + (aml*aml)*C3 + 2.*C4))*
           Ekaon) - 
           1.*(am + amSig)*
           (-1.*(aml*aml)*(-2. + C4) + 
           2.*am*
           ((-1. + C4)*Elep + Enu + C4*Enu))*
           Fm1)*t1*t3 + 
           C1*C4*t3*
           (2.*C1*C4*
           ((amk*amk)*
           (3.*(am*am) - 1.*(amk*amk) + 
           4.*am*amSig + (amSig*amSig) + 
           2.*((am*am) + (amk*amk) - 1.*(amSig*amSig))*
           C4 - 
           1.*
           (-1.*(amk*amk) + 
           (am + amSig)*(3.*am + amSig))*
           (-1. + (aml*aml)*C3)*(C4*C4)) + 
           2.*am*
           (-1.*(am*am) + 2.*(amk*amk) + (amSig*amSig) - 
           2.*
           ((am*am) + 2.*(amk*amk) - 1.*(amSig*amSig))*
           C4 + 
           ((am*am) - 2.*(amk*amk) - 1.*(amSig*amSig))*
           (-1. + (aml*aml)*C3)*(C4*C4))*Ekaon + 
           4.*(am*am)*
           (-1. + 
           C4*(2. + (-1. + (aml*aml)*C3)*C4))*
           (Ekaon*Ekaon))*t3 + 
           ((2.*(amk*amk) + (aml*aml))*(am + amSig) + 
           3.*(aml*aml)*((amk*amk) - 2.*am*Ekaon)*Fm1)
           *(C2*C7*t5 + C8*C9*t6))) + 
           2.*(apkk1*apkk1)*C1*C4*t3*
           ((am + amSig)*(1. + C4)*
           (C2*C7*t5 + C8*C9*t6) + 
           Fm1*(2.*amSig*(-1. + C4)*t1 + 
           (amk*amk)*(C2*C7*t5 + C8*C9*t6) - 
           2.*am*
           (t1 - 1.*C4*t1 + C2*C7*Ekaon*t5 + 
           C8*C9*Ekaon*t6)))) + 
           2.*(akpk*akpk)*C1*C4*t3*
           ((aml*aml)*(am + amSig)*
           (C3*C4*(2.*t1 + t2 - 
           4.*C1*(C4*C4)*
           ((amk*amk) + (-1.*am + amSig)*Ekaon)*t3)
           + Fm1*
           (4.*C1*C4*
           ((amk*amk) + (-1.*am + amSig)*Ekaon)*Fm1*
           t3 + C2*C7*t5 + C8*C9*t6)) - 
           1.*Elep*(-1.*(am + amSig)*(-1. + C4)*
           (C2*C7*t5 + C8*C9*t6) + 
           Fm1*(2.*amSig*(1. + C4)*t1 + 
           (amk*amk)*(C2*C7*t5 + C8*C9*t6) + 
           2.*am*
           ((1. + C4)*t1 - 
           1.*Ekaon*(C2*C7*t5 + C8*C9*t6))))) + 
           akpk*(4.*am*(am + amSig)*C1*(-1. + C4)*C4*
           (Elep*Elep)*Fm1*t1*t3 - 
           1.*(aml*aml*aml*aml)*(am + amSig)*C1*C4*t3*
           (C3*C4*(2.*t1 + t2 - 
           4.*C1*(C4*C4)*
           ((amk*amk) + (-1.*am + amSig)*Ekaon)*t3)
           - 1.*Fm1*(C2*C7*t5 + C8*C9*t6)) + 
           (aml*aml)*(-2.*(amk*amk*amk*amk)*(C1*C1)*(C4*C4)*Fm1*
           (1. + 2.*C4 - 2.*Enu*Fm1)*(t3*t3) + 
           4.*(am*am*am)*(C1*C1)*(C4*C4)*Ekaon*Fm1*
           (-1. + 2.*C4 + 2.*Enu*Fm1)*(t3*t3) - 
           2.*(amSig*amSig)*(C1*C1)*(C4*C4)*
           (2.*((-1.+C4)*(-1.+C4))*Ekaon + 
           (amk*amk)*Fm1*(1. + 2.*C4 + 2.*Enu*Fm1))*
           (t3*t3) + 
           2.*(am*am)*C1*C4*t3*
           (2.*C1*((-1.+C4)*(-1.+C4))*C4*Ekaon*t3 + 
           2.*C1*C4*(-3.*(amk*amk) + 4.*(Ekaon*Ekaon))*Enu*
           (Fm1*Fm1)*t3 + 
           Fm1*
           (-2.*(-1. + C4)*t1 + 
           (amk*amk)*C1*(1. - 6.*C4)*C4*t3)) - 
           1.*C3*
           (2.*Enu*
           (t1*t2 + 
           C1*(C4*C4)*
           (-1.*(amk*amk) + 
           am*(am + amSig + 2.*Ekaon))*
           (2.*t1 + t2)*t3 + 
           2.*(C1*C1)*(C4*C4*C4*C4)*
           ((amk*amk*amk*amk) - 1.*(amk*amk)*(amSig*amSig) + 
           2.*(am*am*am)*Ekaon + 
           (am*am)*(-3.*(amk*amk) + 4.*(Ekaon*Ekaon)) - 
           2.*am*
           ((amSig*amSig)*Ekaon + 
           2.*(amk*amk)*(amSig + Ekaon)))*(t3*t3))
           + (am + amSig)*C1*(C4*C4)*t3*
           (-1.*((amk*amk) - 2.*am*Ekaon)*t2 + 
           2.*apkk1*
           (2.*t1 + t2 - 
           4.*C1*(C4*C4)*
           ((amk*amk) + (-1.*am + amSig)*Ekaon)*
           t3))) - 
           1.*(C2*C7*t5 + C8*C9*t6)*
           (-3.*C4*t1 + 
           (Ekaon - 1.*Enu)*(C2*C7*t5 + C8*C9*t6))
           + 2.*amSig*C1*C4*t3*
           (2.*(-1.*((-1.+C4)*(-1.+C4)) + 
           (1. + C4)*Ekaon*Fm1)*t1 - 
           2.*(amk*amk)*C1*C4*
           (((-1.+C4)*(-1.+C4)) + 
           (1. + 2.*C4)*Ekaon*Fm1)*t3 + 
           ((-2. + C4)*Ekaon + apkk1*Fm1)*
           (C2*C7*t5 + C8*C9*t6)) + 
           2.*am*C1*C4*t3*
           (-2.*
           (((-1.+C4)*(-1.+C4)) + 
           (amSig*(-1. + C4) + Ekaon)*Fm1)*t1 + 
           4.*amSig*C1*C4*(1. + 2.*C4)*(Ekaon*Ekaon)*
           Fm1*t3 - 
           2.*(amSig*amSig)*C1*C4*Ekaon*Fm1*
           (-1. + 2.*C4 + 2.*Enu*Fm1)*t3 - 
           2.*(amk*amk)*C1*C4*
           (((-1.+C4)*(-1.+C4)) + 
           (4.*amSig*C4 - 
           1.*(1. + 2.*C4)*Ekaon)*Fm1 + 
           4.*(amSig + Ekaon)*Enu*(Fm1*Fm1))*t3 + 
           (apkk1*Fm1 + 
           Ekaon*(1. + C4 - 2.*Enu*Fm1))*
           (C2*C7*t5 + C8*C9*t6)) + 
           (amk*amk)*C1*C4*t3*
           (-3.*(C2*C7*t5 + C8*C9*t6) + 
           2.*Fm1*
           ((2. + C4)*t1 + 
           Enu*(C2*C7*t5 + C8*C9*t6)))) - 
           1.*Elep*(2.*((1.+C4)*(1.+C4))*(t1*t1) + 
           2.*(amk*amk*amk*amk)*(C1*C1)*(C4*C4)*
           (1. + C4*(2. + C4 - 1.*(aml*aml)*C3*C4))*
           (t3*t3) + 
           4.*(am*am*am)*(C1*C1)*(C4*C4)*
           (1. + C4*(-2. + C4 - 1.*(aml*aml)*C3*C4))*
           Ekaon*(t3*t3) + 
           2.*(am*am)*C1*(C4*C4)*t3*
           (C1*((amk*amk)*(-3. + (2. - 3.*C4)*C4) + 
           4.*((1.+C4)*(1.+C4))*(Ekaon*Ekaon))*t3 - 
           1.*(aml*aml)*C3*
           (t2 + 
           C1*(C4*C4)*(-3.*(amk*amk) + 4.*(Ekaon*Ekaon))*
           t3)) - 
           2.*t1*
           (2.*C1*C4*
           (2.*am*(am + amSig)*(-1. + C4)*C4 + 
           (amk*amk)*((1.+C4)*(1.+C4)) - 
           2.*am*((1.+C4)*(1.+C4))*Ekaon + 
           (am + amSig)*
           (apkk1 - 1.*apkk1*C4 + 
           am*(1. + C4)*Enu)*Fm1)*t3 + 
           (aml*aml)*
           (-1.*(am + amSig)*C1*(-2. + C4)*C4*Fm1*
           t3 + 
           C3*
           (t2 - 
           1.*C1*(C4*C4)*((amk*amk) - 2.*am*Ekaon)*
           t3))) + 
           (C2*C7*t5 + C8*C9*t6)*
           ((aml*aml)*
           (amSig*C1*C4*t3 + C2*C7*t5 + C8*C9*t6)
           + 2.*apkk1*
           (amSig*C1*C4*(1. + C4)*t3 + C2*C7*t5 + 
           C8*C9*t6)) + 
           (amk*amk)*C1*C4*t3*
           ((aml*aml)*
           (2.*C3*C4*
           (t2 + (amSig*amSig)*C1*(C4*C4)*t3) + 
           Fm1*(C2*C7*t5 + C8*C9*t6)) - 
           2.*
           ((amSig*amSig)*C1*C4*((1.+C4)*(1.+C4))*t3 + 
           amSig*(C2*C7*t5 + C8*C9*t6) - 
           1.*apkk1*Fm1*(C2*C7*t5 + C8*C9*t6)))
           + am*C1*C4*t3*
           ((aml*aml)*
           (4.*(amSig*amSig)*C1*C3*(C4*C4*C4)*Ekaon*t3 + 
           4.*C3*C4*Ekaon*
           (-1.*t2 + 2.*(amk*amk)*C1*(C4*C4)*t3) + 
           amSig*C3*
           (-2.*C4*t2 + 8.*(amk*amk)*C1*(C4*C4*C4)*t3)
           - 1.*(-1. + 2.*Ekaon*Fm1)*
           (C2*C7*t5 + C8*C9*t6)) - 
           2.*
           (2.*(amSig*amSig)*C1*((-1.+C4)*(-1.+C4))*C4*Ekaon*
           t3 - 
           1.*apkk1*(1. + C4 - 2.*Ekaon*Fm1)*
           (C2*C7*t5 + C8*C9*t6) + 
           (amk*amk)*
           (4.*C1*C4*
           (amSig*(1. + (C4*C4)) + 
           ((1.+C4)*(1.+C4))*Ekaon)*t3 + 
           C2*C7*t5 + C8*C9*t6)))) + 
           2.*apkk1*
           (4.*(am*am)*(C1*C1)*(C4*C4)*(1. + (C4*C4))*Ekaon*
           (t3*t3) - 
           4.*(amSig*amSig)*(C1*C1)*(C4*C4)*(1. + (C4*C4))*Ekaon*
           (t3*t3) + 
           (C2*C7*t5 + C8*C9*t6)*
           (C4*(2.*t1 + 
           (amk*amk)*C1*(-2. + Enu*Fm1)*t3) - 
           1.*(Ekaon - 1.*Enu)*
           (C2*C7*t5 + C8*C9*t6)) + 
           amSig*C1*C4*t3*
           (2.*Enu*Fm1*t1 - 4.*(amk*amk)*C1*(C4*C4*C4)*t3 - 
           4.*C2*C7*Ekaon*t5 + C2*C7*Enu*t5 + 
           C8*C9*(-4.*Ekaon + Enu)*t6 + 
           C4*
           (2.*(4. + Enu*Fm1)*t1 - 
           4.*(amk*amk)*C1*t3 - 
           1.*Enu*(C2*C7*t5 + C8*C9*t6))) + 
           am*C1*C4*t3*
           (-4.*(amk*amk)*C1*(C4*C4*C4)*t3 + 
           C4*
           (2.*(4. + Enu*Fm1)*t1 - 
           4.*(amk*amk)*C1*t3 - 
           1.*Enu*(C2*C7*t5 + C8*C9*t6)) + 
           Enu*
           (C2*C7*t5 + C8*C9*t6 + 
           2.*Fm1*
           (t1 - 
           1.*Ekaon*(C2*C7*t5 + C8*C9*t6))))))
           ) + akk1*(8.*(am*am*am)*(C1*C1)*(C4*C4)*
           ((amk*amk)*(1. + 
           (C4*C4)*
           (-1. + 
           (aml*aml)*C3*
           (-2. + (aml*aml)*C3 - 
           1.*C3*((Elep - 1.*Enu)*(Elep - 1.*Enu)))) - 
           4.*C4*(Elep + Enu)*Fm1 + 
           2.*((aml*aml) - 1.*(Elep*Elep) - 1.*(Enu*Enu))*(Fm1*Fm1)
           ) - 1.*Ekaon*
           ((aml*aml*aml*aml)*(C3*C3)*(C4*C4)*Ekaon + 
           (aml*aml)*
           (-2.*C3*(C4*C4)*Ekaon - 
           2.*(C3*C3)*(C4*C4)*(Elep - 1.*Enu)*
           (akpk - 1.*apkk1 + 2.*Ekaon*Elep - 
           2.*Ekaon*Enu) + Ekaon*(Fm1*Fm1)) + 
           4.*Fm1*
           (apkk1*(1. + C4 + Elep*Fm1) + 
           akpk*(-1. + C4 + Enu*Fm1) - 
           2.*Ekaon*
           ((1. + C4)*Elep + (Elep*Elep)*Fm1 + 
           Enu*(-1. + C4 + Enu*Fm1)))))*(t3*t3) + 
           (aml*aml*aml*aml)*(C3*C3)*
           (3.*Ekaon*(t2*t2) - 
           4.*C1*(C4*C4)*
           (akpk*amSig - 1.*amSig*apkk1 + 
           2.*(amk*amk)*(amSig + Ekaon))*t2*t3 + 
           4.*(C1*C1)*(C4*C4*C4*C4)*
           (-2.*(amSig*amSig)*apkk1*Ekaon + 
           (amk*amk*amk*amk)*(2.*amSig + Ekaon) + 
           2.*akpk*amSig*((amk*amk) + amSig*Ekaon) + 
           (amk*amk)*amSig*(-2.*apkk1 + amSig*Ekaon))*
           (t3*t3)) + 
           4.*(am*am)*C1*C4*t3*
           (4.*(aml*aml)*C1*C4*(Ekaon*Ekaon*Ekaon)*
           (C3*(-2. + (aml*aml)*C3)*(C4*C4) + (Fm1*Fm1))*t3 + 
           4.*(akpk*akpk)*C1*C4*Ekaon*
           ((aml*aml)*(C3*C3)*(C4*C4) + 2.*(Fm1*Fm1))*t3 + 
           4.*(apkk1*apkk1)*C1*C4*Ekaon*
           ((aml*aml)*(C3*C3)*(C4*C4) + 2.*(Fm1*Fm1))*t3 + 
           (amk*amk)*C1*C4*
           (2.*(Elep - 1.*Enu) + 4.*C4*(Elep + Enu) + 
           (C4*C4)*
           ((2. + (aml*aml*aml*aml)*(C3*C3))*Elep - 
           1.*(2. + (aml*aml)*C3*(-4. + (aml*aml)*C3))*
           Enu) + 
           3.*(aml*aml)*Fm1*
           (-2. + Elep*Fm1 - 1.*Enu*Fm1) + 
           4.*amSig*
           (1. + 
           (C4*C4)*
           (-1. + 
           (aml*aml)*C3*
           (-2. + (aml*aml)*C3 - 
           1.*C3*((Elep - 1.*Enu)*(Elep - 1.*Enu)))) - 
           4.*C4*(Elep + Enu)*Fm1 + 
           2.*((aml*aml) - 1.*(Elep*Elep) - 1.*(Enu*Enu))*
           (Fm1*Fm1)))*t3 - 
           4.*C4*(Ekaon*Ekaon)*
           ((aml*aml*aml*aml)*C1*(C3*C3)*(C4*C4)*
           (amSig + Elep - 1.*Enu)*t3 + 
           2.*C1*
           (((1.+C4)*(1.+C4))*Elep - 
           1.*((-1.+C4)*(-1.+C4))*Enu - 
           1.*amSig*
           (-1. + (C4*C4) + 
           2.*((-1. + C4)*Elep + Enu + C4*Enu)*
           Fm1))*t3 + 
           (aml*aml)*
           ((C3*C3)*(Elep - 1.*Enu)*t2 - 
           2.*C1*C3*(C4*C4)*(amSig - 2.*Enu)*t3 + 
           3.*C1*(amSig + Elep - 1.*Enu)*(Fm1*Fm1)*t3)
           ) - 1.*Ekaon*
           (-4.*Enu*
           ((1. - 3.*C4)*Fm1*t1 + 
           (aml*aml)*(C3*C3)*C4*Enu*t2) + 
           (amk*amk)*C1*C4*
           (4. + 
           (C4*C4)*
           (-4. + 
           (aml*aml)*C3*
           (-6. + 3.*(aml*aml)*C3 + 8.*C3*(Enu*Enu)))
           + 8.*C4*Enu*Fm1 + 
           Fm1*
           (7.*(aml*aml)*Fm1 + 
           8.*Enu*(-3. + 2.*Enu*Fm1)))*t3 + 
           4.*C4*(Elep*Elep)*
           (4.*(amk*amk)*C1*(Fm1*Fm1)*t3 - 
           1.*(aml*aml)*(C3*C3)*
           (t2 - 2.*(amk*amk)*C1*(C4*C4)*t3)) - 
           4.*Elep*
           (2.*(aml*aml)*(C3*C3)*C4*Enu*
           (-1.*t2 + 2.*(amk*amk)*C1*(C4*C4)*t3) + 
           Fm1*
           (t1 + 3.*C4*t1 - 
           2.*(amk*amk)*C1*C4*(3. + C4)*t3))) + 
           2.*apkk1*
           (-2.*(1. + C4)*Fm1*t1 + 
           2.*(amk*amk)*C1*C4*(1. + 3.*C4)*Fm1*t3 + 
           C1*C4*
           (3.*(aml*aml)*Ekaon + 6.*(amk*amk)*Elep - 
           8.*(Ekaon*Ekaon)*Elep)*(Fm1*Fm1)*t3 + 
           C4*((aml*aml*aml*aml)*C1*(C3*C3)*(C4*C4)*Ekaon*t3 + 
           2.*C1*((1.+C4)*(1.+C4))*Ekaon*t3 + 
           (aml*aml)*(C3*C3)*
           ((Ekaon - 1.*Elep + Enu)*t2 + 
           C1*(C4*C4)*(3.*(amk*amk) - 4.*(Ekaon*Ekaon))*
           (Elep - 1.*Enu)*t3))) - 
           2.*akpk*(C1*C4*
           (3.*(aml*aml)*Ekaon - 6.*(amk*amk)*Enu + 
           8.*(Ekaon*Ekaon)*Enu)*(Fm1*Fm1)*t3 + 
           2.*Fm1*
           (t1 - 1.*C4*t1 + 
           (amk*amk)*C1*(1. - 3.*C4)*C4*t3) + 
           C4*(2.*C1*((-1.+C4)*(-1.+C4))*Ekaon*t3 + 
           (aml*aml*aml*aml)*C1*(C3*C3)*(C4*C4)*Ekaon*t3 + 
           (aml*aml)*C3*
           (C3*(Ekaon - 1.*Elep + Enu)*t2 + 
           C1*(C4*C4)*
           (4.*(-1. + apkk1*C3)*Ekaon + 
           3.*(amk*amk)*C3*(Elep - 1.*Enu) + 
           4.*C3*(Ekaon*Ekaon)*(-1.*Elep + Enu))*t3
           )))) + 
           (aml*aml)*(24.*(amSig*amSig)*(akpk - 1.*apkk1)*(C1*C1)*
           (C4*C4)*Ekaon*(Fm1*Fm1)*(t3*t3) + 
           4.*(amk*amk*amk*amk)*(C1*C1)*(C4*C4)*Fm1*
           (-6. + 2.*amSig*Fm1 + Ekaon*Fm1)*(t3*t3) + 
           (C3*C3)*((2.*akpk - 1.*(amk*amk) - 2.*apkk1)*Ekaon*
           (t2*t2) + 
           4.*amSig*(2.*akpk + (amk*amk) - 2.*apkk1)*
           (akpk - 1.*apkk1)*C1*(C4*C4)*t2*t3 - 
           16.*amSig*((akpk - 1.*apkk1)*(akpk - 1.*apkk1))*(C1*C1)*(C4*C4*C4*C4)*
           ((amk*amk) + amSig*Ekaon)*(t3*t3)) + 
           8.*C3*(-1.*Ekaon*t1*t2 + 
           C1*(C4*C4)*
           ((amk*amk)*(amSig + Ekaon)*(t1 + t2) + 
           akpk*amSig*(2.*t1 + t2))*t3 - 
           1.*(C1*C1)*(C4*C4*C4*C4)*
           ((amk*amk)*(amSig*amSig)*Ekaon + 
           (amk*amk*amk*amk)*(2.*amSig + Ekaon) + 
           4.*akpk*amSig*((amk*amk) + amSig*Ekaon))*
           (t3*t3)) - 
           1.*(C2*C7*t5 + C8*C9*t6)*
           (-4.*C4*t1 + C2*C7*Ekaon*t5 + C8*C9*Ekaon*t6)
           - 4.*amSig*C1*C4*t3*
           (-6.*C4*Ekaon*Fm1*t1 + 
           (Ekaon + (akpk + apkk1)*Fm1)*
           (C2*C7*t5 + C8*C9*t6)) + 
           4.*(amk*amk)*C1*C4*t3*
           (C4*Fm1*
           (6.*t1 + 
           amSig*C1*
           (-6.*(amSig + 2.*Ekaon) + 
           (6.*akpk - 6.*apkk1 + amSig*Ekaon)*
           Fm1)*t3) - 
           1.*(1. + 2.*(amSig + Ekaon)*Fm1)*
           (C2*C7*t5 + C8*C9*t6))) - 
           1.*Elep*(8.*((1.+C4)*(1.+C4))*(t1*t1) + 
           (aml*aml)*(-2.*akpk - 3.*(amk*amk) + (aml*aml) + 
           2.*apkk1)*(C3*C3)*(t2*t2) - 
           16.*C1*C4*
           ((amk*amk)*((1.+C4)*(1.+C4)) - 2.*akpk*amSig*Fm1)*t1*
           t3 + (aml*aml)*((C2*C7*t5 + C8*C9*t6)*(C2*C7*t5 + C8*C9*t6)) + 
           8.*akpk*(amk*amk)*C1*C4*t3*
           ((aml*aml)*(C3*C3)*C4*
           (t2 + (-1.*(amk*amk) + (amSig*amSig))*C1*(C4*C4)*t3)
           + 2.*Fm1*(C2*C7*t5 + C8*C9*t6)) - 
           4.*(amk*amk)*
           (2.*(amSig*amSig)*(C1*C1)*(C4*C4)*
           (((1.+C4)*(1.+C4)) + 2.*apkk1*(Fm1*Fm1))*(t3*t3) + 
           (aml*aml*aml*aml)*C1*(C3*C3)*(C4*C4)*t3*
           (t2 + (amSig*amSig)*C1*(C4*C4)*t3) + 
           2.*amSig*C1*C4*(1. + C4)*t3*
           (C2*C7*t5 + C8*C9*t6) + 
           ((C2*C7*t5 + C8*C9*t6)*(C2*C7*t5 + C8*C9*t6)) + 
           (aml*aml)*C1*C4*t3*
           (2.*apkk1*(C3*C3)*C4*
           (t2 + (amSig*amSig)*C1*(C4*C4)*t3) + 
           Fm1*
           (3.*(amSig*amSig)*C1*C4*Fm1*t3 - 
           1.*C2*C7*t5 - 1.*C8*C9*t6))) + 
           4.*(amk*amk*amk*amk)*C1*C4*t3*
           ((aml*aml*aml*aml)*C1*(C3*C3)*(C4*C4*C4)*t3 + 
           (aml*aml)*C4*
           (3.*C1*(Fm1*Fm1)*t3 + 
           (C3*C3)*(t2 + 2.*apkk1*C1*(C4*C4)*t3)) + 
           2.*(C1*C4*(((1.+C4)*(1.+C4)) + 2.*apkk1*(Fm1*Fm1))*
           t3 - 1.*Fm1*(C2*C7*t5 + C8*C9*t6)))) - 
           4.*(8.*(akpk*akpk)*amSig*(C1*C1)*(C4*C4)*
           ((amk*amk) + amSig*Ekaon)*(Fm1*Fm1)*(t3*t3) + 
           4.*(amSig*amSig)*apkk1*(C1*C1)*(C4*C4)*
           (((1.+C4)*(1.+C4))*Ekaon - 
           1.*(amk*amk)*(-1. + C4)*Fm1 + 
           2.*apkk1*Ekaon*(Fm1*Fm1))*(t3*t3) - 
           1.*(amk*amk)*
           ((C2*C7*t5 + C8*C9*t6)*
           (-2.*C4*t1 + 2.*(amk*amk)*C1*C4*t3 + 
           C2*C7*Ekaon*t5 + C8*C9*Ekaon*t6) + 
           2.*apkk1*C1*C4*Fm1*t3*
           (2.*(-1. + C4)*(t1 + (amk*amk)*C1*C4*t3) - 
           1.*Ekaon*(C2*C7*t5 + C8*C9*t6))) + 
           2.*amSig*C1*C4*t3*
           (4.*(amk*amk)*(apkk1*apkk1)*C1*C4*(Fm1*Fm1)*t3 - 
           1.*(amk*amk)*Ekaon*(C2*C7*t5 + C8*C9*t6) + 
           apkk1*
           (-2.*
           (((1.+C4)*(1.+C4)) + (-1. + C4)*Ekaon*Fm1)*
           t1 + 
           2.*(amk*amk)*C1*C4*
           (((1.+C4)*(1.+C4)) - 
           2.*(-1. + C4)*Ekaon*Fm1)*t3 + 
           (1. + C4)*Ekaon*(C2*C7*t5 + C8*C9*t6)))
           - 2.*akpk*C1*C4*t3*
           (2.*(amSig*amSig)*C1*C4*
           (((-1.+C4)*(-1.+C4))*Ekaon + 
           (amk*amk)*(1. + C4)*Fm1)*t3 + 
           (amk*amk)*Fm1*
           (2.*(1. + C4)*
           (-1.*t1 + (amk*amk)*C1*C4*t3) + 
           Ekaon*(C2*C7*t5 + C8*C9*t6)) + 
           amSig*
           (2.*(((-1.+C4)*(-1.+C4)) - 
           1.*(1. + C4)*Ekaon*Fm1)*t1 + 
           2.*(amk*amk)*C1*C4*
           (((-1.+C4)*(-1.+C4)) + 
           2.*(1. + C4)*Ekaon*Fm1)*t3 - 
           1.*((-1. + C4)*Ekaon + 4.*apkk1*Fm1)*
           (C2*C7*t5 + C8*C9*t6)))) + 
           Enu*(8.*((-1.+C4)*(-1.+C4))*(t1*t1) - 
           2.*akpk*(aml*aml)*(C3*C3)*(t2*t2) + 
           8.*akpk*(amk*amk)*(aml*aml)*C1*(C3*C3)*(C4*C4)*t2*t3 - 
           8.*akpk*(amk*amk)*((amk*amk) - 1.*(amSig*amSig))*(C1*C1)*
           (C4*C4)*((aml*aml)*(C3*C3)*(C4*C4) + 2.*(Fm1*Fm1))*(t3*t3)
           + 8.*t1*
           (2.*C1*C4*
           ((amk*amk)*((-1.+C4)*(-1.+C4)) + 
           2.*amSig*apkk1*Fm1)*t3 - 
           1.*(aml*aml)*C3*(t2 - 2.*(amk*amk)*C1*(C4*C4)*t3))
           + (aml*aml)*
           (((aml*aml) + 2.*apkk1)*(C3*C3)*(t2*t2) + 
           ((C2*C7*t5 + C8*C9*t6)*(C2*C7*t5 + C8*C9*t6))) + 
           4.*(amk*amk*amk*amk)*C1*C4*t3*
           ((aml*aml*aml*aml)*C1*(C3*C3)*(C4*C4*C4)*t3 + 
           (aml*aml)*C4*
           (-4.*C1*C3*(C4*C4)*t3 + 3.*C1*(Fm1*Fm1)*t3 + 
           (C3*C3)*(t2 + 2.*apkk1*C1*(C4*C4)*t3)) + 
           2.*(C1*((-1.+C4)*(-1.+C4))*C4*t3 - 
           1.*Fm1*(C2*C7*t5 + C8*C9*t6))) - 
           1.*(amk*amk)*
           (4.*(aml*aml*aml*aml)*C1*(C3*C3)*(C4*C4)*t3*
           (t2 + (amSig*amSig)*C1*(C4*C4)*t3) + 
           (aml*aml)*
           (-8.*C1*C3*(C4*C4)*t3*
           (t2 + 2.*(amSig*amSig)*C1*(C4*C4)*t3) + 
           (C3*C3)*
           (3.*(t2*t2) + 8.*apkk1*C1*(C4*C4)*t2*t3 + 
           8.*(amSig*amSig)*apkk1*(C1*C1)*(C4*C4*C4*C4)*(t3*t3))
           + 4.*C1*C4*Fm1*t3*
           (3.*(amSig*amSig)*C1*C4*Fm1*t3 + C2*C7*t5 + 
           C8*C9*t6)) + 
           4.*(2.*(amSig*amSig)*(C1*C1)*((-1.+C4)*(-1.+C4))*(C4*C4)*
           (t3*t3) - 
           2.*amSig*C1*(-1. + C4)*C4*t3*
           (C2*C7*t5 + C8*C9*t6) + 
           (C2*C7*t5 + C8*C9*t6)*
           (4.*apkk1*C1*C4*Fm1*t3 + C2*C7*t5 + 
           C8*C9*t6)))) + 
           2.*am*((aml*aml)*(C3*C3)*((Ekaon - 1.*Elep + Enu)*(Ekaon - 1.*Elep + Enu))*
           (t2*t2) - 
           2.*(-2.*(-1. + (C4*C4))*(t1*t1) - 
           2.*C1*C4*
           (2.*akpk*
           (1. + C4*(-2. + (aml*aml)*C3 + C4) - 
           1.*amSig*Fm1 + 
           (Ekaon + C4*(amSig + Ekaon) - 
           2.*Elep)*Fm1) + 
           2.*apkk1*
           (((1.+C4)*(1.+C4)) - 
           1.*
           (amSig*(1. + C4) + 
           (-1. + C4)*Ekaon - 2.*Enu)*Fm1) + 
           (amk*amk)*C4*
           ((aml*aml)*C3 + 4.*(-1.*Elep + Enu)*Fm1)
           - 1.*Ekaon*
           (4.*((1.+C4)*(1.+C4))*Elep + 
           4.*((-1.+C4)*(-1.+C4))*Enu + 
           2.*amSig*(-1. + C4)*Elep*Fm1 - 
           2.*amSig*(1. + C4)*Enu*Fm1 + 
           (aml*aml)*C4*
           (2.*C3*(Ekaon + 2.*Enu) + 3.*Fm1)))*
           t1*t3 + 
           C1*C4*t3*
           (-2.*(amk*amk*amk*amk)*C1*C4*
           (1. + 
           (C4*C4)*
           (-1. + 
           (aml*aml)*C3*
           (-2. + 
           C3*((aml*aml) + ((Elep - 1.*Enu)*(Elep - 1.*Enu)))))
           + 4.*(Elep - 1.*Enu)*Fm1 + 
           2.*((aml*aml) + (Elep*Elep) + (Enu*Enu))*(Fm1*Fm1)
           )*t3 + 
           2.*(akpk*akpk)*C4*
           (4.*(amk*amk)*C1*(Fm1*Fm1)*t3 - 
           1.*(aml*aml)*(C3*C3)*
           (t2 - 2.*(amk*amk)*C1*(C4*C4)*t3)) + 
           (aml*aml*aml*aml)*(C3*C3)*C4*
           (-1.*apkk1*t2 + 
           2.*Ekaon*
           (-2.*Ekaon*t2 + Elep*t2 - 
           1.*Enu*t2 + 
           (amSig*amSig)*C1*(C4*C4)*Ekaon*t3)) + 
           2.*apkk1*Ekaon*
           (8.*amSig*C1*(-1. + C4)*C4*Ekaon*Fm1*
           t3 - 
           4.*(amSig*amSig)*C1*C4*Fm1*
           (1. + C4 + Elep*Fm1)*t3 + 
           (1. + C4 - 2.*(Ekaon + 2.*Enu)*Fm1)*
           (C2*C7*t5 + C8*C9*t6)) - 
           1.*(aml*aml)*
           (2.*(apkk1*apkk1)*(C3*C3)*C4*t2 + 
           apkk1*
           (2.*(C3*C3)*C4*
           (2.*Ekaon*(-1.*Elep + Enu)*t2 + 
           amSig*(Ekaon - 1.*Elep + Enu)*t2 + 
           2.*(amSig*amSig)*C1*(C4*C4)*Ekaon*
           (Elep - 1.*Enu)*t3) - 
           1.*Fm1*(C2*C7*t5 + C8*C9*t6)) + 
           Ekaon*
           (-4.*C3*C4*(Ekaon + Enu)*t2 + 
           4.*(amSig*amSig)*C1*C3*(C4*C4*C4)*Ekaon*t3 + 
           24.*amSig*C1*C4*Ekaon*Fm1*t3 - 
           2.*(amSig*amSig)*C1*C4*Ekaon*(Fm1*Fm1)*
           t3 + 
           (1. + 
           2.*(2.*Ekaon + Elep + Enu)*Fm1)*
           (C2*C7*t5 + C8*C9*t6))) + 
           akpk*
           ((aml*aml*aml*aml)*(C3*C3)*C4*
           (t2 - 2.*(amk*amk)*C1*(C4*C4)*t3) - 
           2.*
           (-8.*amSig*C1*C4*(1. + C4)*(Ekaon*Ekaon)*
           Fm1*t3 + 
           4.*(amSig*amSig)*C1*C4*Ekaon*Fm1*
           (-1. + C4 + Enu*Fm1)*t3 + 
           2.*(amk*amk)*C1*C4*
           (((-1.+C4)*(-1.+C4)) + 
           2.*
           (2.*amSig*C4 - 1.*(1. + C4)*Ekaon)*
           Fm1 + 4.*(amSig + Ekaon)*Enu*(Fm1*Fm1)
           )*t3 - 
           1.*
           (4.*apkk1*Fm1 + 
           Ekaon*
           (-1. + C4 + 2.*Ekaon*Fm1 - 
           4.*Elep*Fm1))*(C2*C7*t5 + C8*C9*t6)
           ) + 
           (aml*aml)*
           (C3*
           (-2.*C4*t2 + 8.*(amk*amk)*C1*(C4*C4*C4)*t3)
           + (C3*C3)*C4*
           (-1.*(amk*amk)*
           (t2 + 
           8.*C1*(C4*C4)*
           (apkk1 - 
           1.*(amSig + Ekaon)*(Elep - 1.*Enu))
           *t3) + 
           2.*
           (2.*apkk1*t2 + 
           2.*Ekaon*(-1.*Elep + Enu)*t2 + 
           amSig*(Ekaon - 1.*Elep + Enu)*t2 + 
           2.*(amSig*amSig)*C1*(C4*C4)*Ekaon*
           (Elep - 1.*Enu)*t3)) + 
           Fm1*
           (-6.*(amk*amk)*C1*C4*Fm1*t3 + 
           C2*C7*t5 + C8*C9*t6))) + 
           (amk*amk)*
           (2.*(aml*aml*aml*aml)*(C3*C3)*C4*
           (t2 + 
           C1*(C4*C4)*
           (-1.*(amSig*amSig) + apkk1 + 
           amSig*(Ekaon - 1.*Elep + Enu) + 
           2.*Ekaon*(Ekaon - 1.*Elep + Enu))*
           t3) + 
           2.*
           (-4.*C1*C4*Ekaon*
           (((1.+C4)*(1.+C4))*Elep - 
           1.*((-1.+C4)*(-1.+C4))*Enu)*t3 + 
           4.*(apkk1*apkk1)*C1*C4*(Fm1*Fm1)*t3 + 
           2.*apkk1*C1*C4*
           (((1.+C4)*(1.+C4)) + 
           2.*(-1. + C4)*Ekaon*Fm1 - 
           4.*Ekaon*Elep*(Fm1*Fm1))*t3 + 
           (amSig*amSig)*C1*C4*
           (-1. + (C4*C4) + 
           4.*C4*(Elep + Enu)*Fm1 + 
           2.*((Elep*Elep) + (Enu*Enu))*(Fm1*Fm1))*t3 - 
           2.*amSig*C1*C4*
           (-1.*((-1.+C4)*(-1.+C4))*Enu + 
           4.*apkk1*C4*Fm1 - 
           1.*Ekaon*
           (-1. + (C4*C4) + 
           2.*
           ((-1. + C4)*Elep + Enu + C4*Enu)*
           Fm1) + 
           Elep*
           (((1.+C4)*(1.+C4)) + 4.*apkk1*(Fm1*Fm1)))*
           t3 + 
           (Ekaon - 1.*Elep - 1.*C4*Elep + 
           Enu - 1.*C4*Enu + 
           2.*Ekaon*(Elep - 1.*Enu)*Fm1)*
           (C2*C7*t5 + C8*C9*t6)) + 
           (aml*aml)*
           (2.*C3*C4*
           (-1.*t2 + 
           2.*C1*(C4*C4)*
           ((amSig*amSig) - 
           1.*amSig*(Ekaon + 2.*Enu) - 
           2.*Ekaon*(Ekaon + 2.*Enu))*t3) + 
           (C3*C3)*C4*
           (4.*(apkk1*apkk1)*C1*(C4*C4)*t3 + 
           apkk1*
           (t2 - 
           8.*C1*(C4*C4)*(amSig + Ekaon)*
           (Elep - 1.*Enu)*t3) + 
           2.*(Elep - 1.*Enu)*
           (-2.*Ekaon*t2 + 
           (Elep - 1.*Enu)*
           (t2 + (amSig*amSig)*C1*(C4*C4)*t3))) + 
           2.*Fm1*
           (C1*C4*
           (-2.*(amSig*amSig)*Fm1 + 3.*apkk1*Fm1 + 
           amSig*
           (6. - 
           1.*(Ekaon + 3.*Elep - 3.*Enu)*Fm1)
           + 2.*Ekaon*
           (-3. + 
           (Ekaon - 3.*Elep + 3.*Enu)*Fm1))*t3
           + C2*C7*t5 + C8*C9*t6)))))))));
  
  if (sol <= 0.0) {
    LOG( "SKXSec", pERROR )
      << "Matrix element is 0 or negative! This is bad and really should throw an exception!";
  }  
  
  return sol;
}

//____________________________________________________________________________
double AlamSimoAtharVacasSKPXSec2014::Amatrix_NP(double theta, double phikq) const
{

  double sol = 0.;
  
  double akk1=0., zdotq=0., qdotpk=0., akcrosk1=0., qcrospk=0.;
  double zdotpk=0., azpk=0., aqkaon=0., akpk=0., apkk1=0.;
  double C1=0., C2=0., C3=0., C4=0., C5=0., C6=0., C7=0./*, C8=0., C9=0.*/;
  double aq2=0., gform=0., con=0., t1=0., t2=0., t3=0., t4=0., t5=0./*, t6=0.*/;
  
  akk1=Enu*Elep-Enu*alepvec*cos(theta);
  zdotq=(Enu*Enu-alepvec*alepvec)/2.0;
  qdotpk=aqvec*pkvec*angkq;
  akcrosk1 = Enu*alepvec*sin(theta);
  qcrospk=aqvec*pkvec*sqrt(1.0-angkq*angkq);
  zdotpk=(akcrosk1*qcrospk*cos(phikq)+zdotq*qdotpk)/(aqvec*aqvec);
  azpk=Ekaon*(Enu+Elep)/2.0-zdotpk;
  aqkaon=aq0*Ekaon-aqvec*pkvec*angkq;
  akpk=azpk + aqkaon/2.0;
  apkk1=azpk - aqkaon/2.0;                                          
  C1=1.0/(am*am+amk*amk-2.0*am*Ekaon-amSig*amSig);
  C2=d+f;
  C3=1./(aml*aml-2.0*akk1-amk*amk);
  C4=d-f;
  C5=d+3.*f;
  C6=1.0/(am*am+amk*amk-2.0*am*Ekaon-amLam*amLam);
  C7=2.0*am/(aml*aml-2.0*akk1+amk*amk-2.*aqkaon-ampi*ampi);
  //C8=am/(aml*aml-2.0*akk1+amk*amk-2.*aqkaon-amEta*amEta);
  //C9=d - 3.*f;
  aq2=aml*aml-2.0*akk1;             
  gform=1.0/pow(1.0-aq2/(1.1*1.1),4);

  con=g*g*Vus*Vus/(4.0*fpi*fpi);

  t1=1.0;
  t2=1.0;               // !Full Term
  t3=1.0;
  t4=1.0;
  t5=1.0;
  //t6=1.0;
  
  sol = con*gform*(0.4444444444444444*am*
           (akk1*(-72.*Elep*(t1*t1) + 144.*C2*Elep*(t1*t1) - 
           72.*(C2*C2)*Elep*(t1*t1) + 72.*Enu*(t1*t1) + 
           144.*C2*Enu*(t1*t1) + 72.*(C2*C2)*Enu*(t1*t1) - 
           72.*(aml*aml)*C3*Ekaon*t1*t2 - 
           72.*(aml*aml)*C3*Enu*t1*t2 + 
           18.*akpk*(aml*aml)*(C3*C3)*Ekaon*(t2*t2) - 
           9.*(amk*amk)*(aml*aml)*(C3*C3)*Ekaon*(t2*t2) + 
           27.*(aml*aml*aml*aml)*(C3*C3)*Ekaon*(t2*t2) - 
           18.*(aml*aml)*apkk1*(C3*C3)*Ekaon*(t2*t2) + 
           18.*akpk*(aml*aml)*(C3*C3)*Elep*(t2*t2) + 
           27.*(amk*amk)*(aml*aml)*(C3*C3)*Elep*(t2*t2) - 
           9.*(aml*aml*aml*aml)*(C3*C3)*Elep*(t2*t2) - 
           18.*(aml*aml)*apkk1*(C3*C3)*Elep*(t2*t2) - 
           18.*akpk*(aml*aml)*(C3*C3)*Enu*(t2*t2) - 
           27.*(amk*amk)*(aml*aml)*(C3*C3)*Enu*(t2*t2) + 
           9.*(aml*aml*aml*aml)*(C3*C3)*Enu*(t2*t2) + 
           18.*(aml*aml)*apkk1*(C3*C3)*Enu*(t2*t2) - 
           72.*akpk*amSig*C1*C4*t1*t3 - 
           72.*amSig*apkk1*C1*C4*t1*t3 - 
           72.*akpk*amSig*C1*C2*C4*t1*t3 + 
           72.*amSig*apkk1*C1*C2*C4*t1*t3 + 
           72.*akpk*amSig*C1*(C4*C4)*t1*t3 - 
           72.*amSig*apkk1*C1*(C4*C4)*t1*t3 + 
           72.*akpk*amSig*C1*C2*(C4*C4)*t1*t3 + 
           72.*amSig*apkk1*C1*C2*(C4*C4)*t1*t3 - 
           72.*akpk*(aml*aml)*amSig*C1*C3*(C4*C4)*t1*t3 - 
           36.*(amk*amk)*(aml*aml)*amSig*C1*C3*(C4*C4)*t1*t3 - 
           36.*(amk*amk)*(aml*aml)*C1*C3*(C4*C4)*Ekaon*t1*t3 - 
           72.*(amk*amk)*C1*C4*Elep*t1*t3 + 
           72.*(amk*amk)*C1*C2*C4*Elep*t1*t3 - 
           72.*(amk*amk)*C1*(C4*C4)*Elep*t1*t3 + 
           72.*(amk*amk)*C1*C2*(C4*C4)*Elep*t1*t3 - 
           72.*(amk*amk)*C1*C4*Enu*t1*t3 - 
           72.*(amk*amk)*C1*C2*C4*Enu*t1*t3 + 
           72.*(amk*amk)*C1*(C4*C4)*Enu*t1*t3 + 
           72.*(amk*amk)*C1*C2*(C4*C4)*Enu*t1*t3 - 
           72.*(amk*amk)*(aml*aml)*C1*C3*(C4*C4)*Enu*t1*t3 + 
           72.*akpk*(amk*amk)*C1*C4*Fm1*t1*t3 + 
           72.*(amk*amk)*apkk1*C1*C4*Fm1*t1*t3 - 
           72.*akpk*(amk*amk)*C1*C2*C4*Fm1*t1*t3 + 
           108.*(amk*amk)*(aml*aml)*C1*C2*C4*Fm1*t1*t3 + 
           72.*(amk*amk)*apkk1*C1*C2*C4*Fm1*t1*t3 + 
           72.*akpk*amSig*C1*C4*Ekaon*Fm1*t1*t3 + 
           72.*amSig*apkk1*C1*C4*Ekaon*Fm1*t1*t3 - 
           72.*akpk*amSig*C1*C2*C4*Ekaon*Fm1*t1*t3 + 
           108.*(aml*aml)*amSig*C1*C2*C4*Ekaon*Fm1*t1*t3 + 
           72.*amSig*apkk1*C1*C2*C4*Ekaon*Fm1*t1*t3 + 
           144.*akpk*amSig*C1*C4*Elep*Fm1*t1*t3 - 
           144.*amSig*apkk1*C1*C4*Enu*Fm1*t1*t3 - 
           36.*akpk*(aml*aml)*amSig*C1*C3*(C4*C4)*t2*t3 - 
           36.*(amk*amk)*(aml*aml)*amSig*C1*C3*(C4*C4)*t2*t3 - 
           36.*(akpk*akpk)*(aml*aml)*amSig*C1*(C3*C3)*(C4*C4)*t2*t3 - 
           18.*akpk*(amk*amk)*(aml*aml)*amSig*C1*(C3*C3)*(C4*C4)*t2*
           t3 + 18.*akpk*(aml*aml*aml*aml)*amSig*C1*(C3*C3)*(C4*C4)*t2*
           t3 + 36.*(amk*amk)*(aml*aml*aml*aml)*amSig*C1*(C3*C3)*(C4*C4)*t2*
           t3 + 72.*akpk*(aml*aml)*amSig*apkk1*C1*(C3*C3)*(C4*C4)*
           t2*t3 + 18.*(amk*amk)*(aml*aml)*amSig*apkk1*C1*(C3*C3)*
           (C4*C4)*t2*t3 - 
           18.*(aml*aml*aml*aml)*amSig*apkk1*C1*(C3*C3)*(C4*C4)*t2*t3 - 
           36.*(aml*aml)*amSig*(apkk1*apkk1)*C1*(C3*C3)*(C4*C4)*t2*t3 - 
           36.*(amk*amk)*(aml*aml)*C1*C3*(C4*C4)*Ekaon*t2*t3 + 
           36.*(amk*amk)*(aml*aml*aml*aml)*C1*(C3*C3)*(C4*C4)*Ekaon*t2*t3 + 
           36.*akpk*(amk*amk)*(aml*aml)*C1*(C3*C3)*(C4*C4)*Elep*t2*
           t3 + 18.*(amk*amk*amk*amk)*(aml*aml)*C1*(C3*C3)*(C4*C4)*Elep*t2*
           t3 - 18.*(amk*amk)*(aml*aml*aml*aml)*C1*(C3*C3)*(C4*C4)*Elep*t2*
           t3 - 36.*(amk*amk)*(aml*aml)*apkk1*C1*(C3*C3)*(C4*C4)*
           Elep*t2*t3 - 
           36.*(amk*amk)*(aml*aml)*C1*C3*(C4*C4)*Enu*t2*t3 - 
           36.*akpk*(amk*amk)*(aml*aml)*C1*(C3*C3)*(C4*C4)*Enu*t2*
           t3 - 18.*(amk*amk*amk*amk)*(aml*aml)*C1*(C3*C3)*(C4*C4)*Enu*t2*
           t3 + 18.*(amk*amk)*(aml*aml*aml*aml)*C1*(C3*C3)*(C4*C4)*Enu*t2*
           t3 + 36.*(amk*amk)*(aml*aml)*apkk1*C1*(C3*C3)*(C4*C4)*Enu*
           t2*t3 + 36.*akpk*(amk*amk)*amSig*(C1*C1)*(C4*C4)*
           (t3*t3) - 36.*(amk*amk)*amSig*apkk1*(C1*C1)*(C4*C4)*
           (t3*t3) - 72.*akpk*(amk*amk)*amSig*(C1*C1)*(C4*C4*C4)*
           (t3*t3) - 72.*(amk*amk)*amSig*apkk1*(C1*C1)*(C4*C4*C4)*
           (t3*t3) + 36.*akpk*(amk*amk)*amSig*(C1*C1)*(C4*C4*C4*C4)*
           (t3*t3) - 36.*(amk*amk)*amSig*apkk1*(C1*C1)*(C4*C4*C4*C4)*
           (t3*t3) - 72.*akpk*(amk*amk)*(aml*aml)*amSig*(C1*C1)*C3*
           (C4*C4*C4*C4)*(t3*t3) - 
           36.*(amk*amk*amk*amk)*(aml*aml)*amSig*(C1*C1)*C3*(C4*C4*C4*C4)*(t3*t3) - 
           36.*(akpk*akpk)*(amk*amk)*(aml*aml)*amSig*(C1*C1)*(C3*C3)*
           (C4*C4*C4*C4)*(t3*t3) + 
           18.*akpk*(amk*amk)*(aml*aml*aml*aml)*amSig*(C1*C1)*(C3*C3)*(C4*C4*C4*C4)*
           (t3*t3) + 18.*(amk*amk*amk*amk)*(aml*aml*aml*aml)*amSig*(C1*C1)*(C3*C3)*
           (C4*C4*C4*C4)*(t3*t3) + 
           72.*akpk*(amk*amk)*(aml*aml)*amSig*apkk1*(C1*C1)*(C3*C3)*
           (C4*C4*C4*C4)*(t3*t3) - 
           18.*(amk*amk)*(aml*aml*aml*aml)*amSig*apkk1*(C1*C1)*(C3*C3)*(C4*C4*C4*C4)*
           (t3*t3) - 36.*(amk*amk)*(aml*aml)*amSig*(apkk1*apkk1)*(C1*C1)*
           (C3*C3)*(C4*C4*C4*C4)*(t3*t3) + 
           36.*akpk*(amSig*amSig)*(C1*C1)*(C4*C4)*Ekaon*(t3*t3) - 
           36.*(amSig*amSig)*apkk1*(C1*C1)*(C4*C4)*Ekaon*(t3*t3) - 
           72.*akpk*(amSig*amSig)*(C1*C1)*(C4*C4*C4)*Ekaon*(t3*t3) - 
           72.*(amSig*amSig)*apkk1*(C1*C1)*(C4*C4*C4)*Ekaon*(t3*t3) + 
           36.*akpk*(amSig*amSig)*(C1*C1)*(C4*C4*C4*C4)*Ekaon*(t3*t3) - 
           36.*(amSig*amSig)*apkk1*(C1*C1)*(C4*C4*C4*C4)*Ekaon*(t3*t3) - 
           18.*(amk*amk*amk*amk)*(aml*aml)*(C1*C1)*C3*(C4*C4*C4*C4)*Ekaon*(t3*t3) - 
           72.*akpk*(aml*aml)*(amSig*amSig)*(C1*C1)*C3*(C4*C4*C4*C4)*Ekaon*
           (t3*t3) - 18.*(amk*amk)*(aml*aml)*(amSig*amSig)*(C1*C1)*C3*
           (C4*C4*C4*C4)*Ekaon*(t3*t3) + 
           9.*(amk*amk*amk*amk)*(aml*aml*aml*aml)*(C1*C1)*(C3*C3)*(C4*C4*C4*C4)*Ekaon*(t3*t3) - 
           36.*(akpk*akpk)*(aml*aml)*(amSig*amSig)*(C1*C1)*(C3*C3)*(C4*C4*C4*C4)*
           Ekaon*(t3*t3) + 
           18.*akpk*(aml*aml*aml*aml)*(amSig*amSig)*(C1*C1)*(C3*C3)*(C4*C4*C4*C4)*Ekaon*
           (t3*t3) + 9.*(amk*amk)*(aml*aml*aml*aml)*(amSig*amSig)*(C1*C1)*(C3*C3)*
           (C4*C4*C4*C4)*Ekaon*(t3*t3) + 
           72.*akpk*(aml*aml)*(amSig*amSig)*apkk1*(C1*C1)*(C3*C3)*(C4*C4*C4*C4)*
           Ekaon*(t3*t3) - 
           18.*(aml*aml*aml*aml)*(amSig*amSig)*apkk1*(C1*C1)*(C3*C3)*(C4*C4*C4*C4)*
           Ekaon*(t3*t3) - 
           36.*(aml*aml)*(amSig*amSig)*(apkk1*apkk1)*(C1*C1)*(C3*C3)*(C4*C4*C4*C4)*
           Ekaon*(t3*t3) - 
           18.*(amk*amk*amk*amk)*(C1*C1)*(C4*C4)*Elep*(t3*t3) + 
           18.*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4)*Elep*(t3*t3) - 
           36.*(amk*amk*amk*amk)*(C1*C1)*(C4*C4*C4)*Elep*(t3*t3) + 
           36.*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4*C4)*Elep*(t3*t3) - 
           18.*(amk*amk*amk*amk)*(C1*C1)*(C4*C4*C4*C4)*Elep*(t3*t3) + 
           18.*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4*C4*C4)*Elep*(t3*t3) + 
           18.*akpk*(amk*amk*amk*amk)*(aml*aml)*(C1*C1)*(C3*C3)*(C4*C4*C4*C4)*Elep*
           (t3*t3) - 9.*(amk*amk*amk*amk)*(aml*aml*aml*aml)*(C1*C1)*(C3*C3)*(C4*C4*C4*C4)*Elep*
           (t3*t3) - 18.*akpk*(amk*amk)*(aml*aml)*(amSig*amSig)*(C1*C1)*
           (C3*C3)*(C4*C4*C4*C4)*Elep*(t3*t3) + 
           9.*(amk*amk)*(aml*aml*aml*aml)*(amSig*amSig)*(C1*C1)*(C3*C3)*(C4*C4*C4*C4)*Elep*
           (t3*t3) - 18.*(amk*amk*amk*amk)*(aml*aml)*apkk1*(C1*C1)*(C3*C3)*
           (C4*C4*C4*C4)*Elep*(t3*t3) + 
           18.*(amk*amk)*(aml*aml)*(amSig*amSig)*apkk1*(C1*C1)*(C3*C3)*
           (C4*C4*C4*C4)*Elep*(t3*t3) + 
           18.*(amk*amk*amk*amk)*(C1*C1)*(C4*C4)*Enu*(t3*t3) - 
           18.*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4)*Enu*(t3*t3) - 
           36.*(amk*amk*amk*amk)*(C1*C1)*(C4*C4*C4)*Enu*(t3*t3) + 
           36.*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4*C4)*Enu*(t3*t3) + 
           18.*(amk*amk*amk*amk)*(C1*C1)*(C4*C4*C4*C4)*Enu*(t3*t3) - 
           18.*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4*C4*C4)*Enu*(t3*t3) - 
           36.*(amk*amk*amk*amk)*(aml*aml)*(C1*C1)*C3*(C4*C4*C4*C4)*Enu*(t3*t3) + 
           36.*(amk*amk)*(aml*aml)*(amSig*amSig)*(C1*C1)*C3*(C4*C4*C4*C4)*Enu*
           (t3*t3) - 18.*akpk*(amk*amk*amk*amk)*(aml*aml)*(C1*C1)*(C3*C3)*
           (C4*C4*C4*C4)*Enu*(t3*t3) + 
           9.*(amk*amk*amk*amk)*(aml*aml*aml*aml)*(C1*C1)*(C3*C3)*(C4*C4*C4*C4)*Enu*(t3*t3) + 
           18.*akpk*(amk*amk)*(aml*aml)*(amSig*amSig)*(C1*C1)*(C3*C3)*
           (C4*C4*C4*C4)*Enu*(t3*t3) - 
           9.*(amk*amk)*(aml*aml*aml*aml)*(amSig*amSig)*(C1*C1)*(C3*C3)*(C4*C4*C4*C4)*Enu*
           (t3*t3) + 18.*(amk*amk*amk*amk)*(aml*aml)*apkk1*(C1*C1)*(C3*C3)*
           (C4*C4*C4*C4)*Enu*(t3*t3) - 
           18.*(amk*amk)*(aml*aml)*(amSig*amSig)*apkk1*(C1*C1)*(C3*C3)*
           (C4*C4*C4*C4)*Enu*(t3*t3) + 
           36.*akpk*(amk*amk*amk*amk)*(C1*C1)*(C4*C4)*Fm1*(t3*t3) - 
           54.*(amk*amk*amk*amk)*(aml*aml)*(C1*C1)*(C4*C4)*Fm1*(t3*t3) + 
           36.*akpk*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4)*Fm1*(t3*t3) - 
           54.*(amk*amk)*(aml*aml)*(amSig*amSig)*(C1*C1)*(C4*C4)*Fm1*
           (t3*t3) - 36.*(amk*amk*amk*amk)*apkk1*(C1*C1)*(C4*C4)*Fm1*
           (t3*t3) - 36.*(amk*amk)*(amSig*amSig)*apkk1*(C1*C1)*(C4*C4)*
           Fm1*(t3*t3) + 
           36.*akpk*(amk*amk*amk*amk)*(C1*C1)*(C4*C4*C4)*Fm1*(t3*t3) + 
           36.*akpk*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4*C4)*Fm1*(t3*t3) + 
           36.*(amk*amk*amk*amk)*apkk1*(C1*C1)*(C4*C4*C4)*Fm1*(t3*t3) + 
           36.*(amk*amk)*(amSig*amSig)*apkk1*(C1*C1)*(C4*C4*C4)*Fm1*
           (t3*t3) + 72.*akpk*(amk*amk)*amSig*(C1*C1)*(C4*C4)*Ekaon*
           Fm1*(t3*t3) - 
           108.*(amk*amk)*(aml*aml)*amSig*(C1*C1)*(C4*C4)*Ekaon*Fm1*
           (t3*t3) - 72.*(amk*amk)*amSig*apkk1*(C1*C1)*(C4*C4)*
           Ekaon*Fm1*(t3*t3) + 
           72.*akpk*(amk*amk)*amSig*(C1*C1)*(C4*C4*C4)*Ekaon*Fm1*
           (t3*t3) + 72.*(amk*amk)*amSig*apkk1*(C1*C1)*(C4*C4*C4)*
           Ekaon*Fm1*(t3*t3) - 
           72.*(akpk*akpk)*(amk*amk)*amSig*(C1*C1)*(C4*C4)*(Fm1*Fm1)*
           (t3*t3) + 54.*akpk*(amk*amk)*(aml*aml)*amSig*(C1*C1)*
           (C4*C4)*(Fm1*Fm1)*(t3*t3) + 
           18.*(amk*amk*amk*amk)*(aml*aml)*amSig*(C1*C1)*(C4*C4)*(Fm1*Fm1)*
           (t3*t3) - 54.*(amk*amk)*(aml*aml)*amSig*apkk1*(C1*C1)*
           (C4*C4)*(Fm1*Fm1)*(t3*t3) - 
           72.*(amk*amk)*amSig*(apkk1*apkk1)*(C1*C1)*(C4*C4)*(Fm1*Fm1)*
           (t3*t3) + 9.*(amk*amk*amk*amk)*(aml*aml)*(C1*C1)*(C4*C4)*Ekaon*
           (Fm1*Fm1)*(t3*t3) - 
           72.*(akpk*akpk)*(amSig*amSig)*(C1*C1)*(C4*C4)*Ekaon*(Fm1*Fm1)*
           (t3*t3) + 54.*akpk*(aml*aml)*(amSig*amSig)*(C1*C1)*(C4*C4)*
           Ekaon*(Fm1*Fm1)*(t3*t3) + 
           9.*(amk*amk)*(aml*aml)*(amSig*amSig)*(C1*C1)*(C4*C4)*Ekaon*
           (Fm1*Fm1)*(t3*t3) - 
           54.*(aml*aml)*(amSig*amSig)*apkk1*(C1*C1)*(C4*C4)*Ekaon*
           (Fm1*Fm1)*(t3*t3) - 
           72.*(amSig*amSig)*(apkk1*apkk1)*(C1*C1)*(C4*C4)*Ekaon*(Fm1*Fm1)*
           (t3*t3) - 27.*(amk*amk*amk*amk)*(aml*aml)*(C1*C1)*(C4*C4)*Elep*
           (Fm1*Fm1)*(t3*t3) + 
           27.*(amk*amk)*(aml*aml)*(amSig*amSig)*(C1*C1)*(C4*C4)*Elep*
           (Fm1*Fm1)*(t3*t3) - 
           36.*(amk*amk*amk*amk)*apkk1*(C1*C1)*(C4*C4)*Elep*(Fm1*Fm1)*(t3*t3) + 
           36.*(amk*amk)*(amSig*amSig)*apkk1*(C1*C1)*(C4*C4)*Elep*
           (Fm1*Fm1)*(t3*t3) - 
           36.*akpk*(amk*amk*amk*amk)*(C1*C1)*(C4*C4)*Enu*(Fm1*Fm1)*(t3*t3) + 
           27.*(amk*amk*amk*amk)*(aml*aml)*(C1*C1)*(C4*C4)*Enu*(Fm1*Fm1)*(t3*t3) + 
           36.*akpk*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4)*Enu*(Fm1*Fm1)*
           (t3*t3) - 27.*(amk*amk)*(aml*aml)*(amSig*amSig)*(C1*C1)*(C4*C4)*
           Enu*(Fm1*Fm1)*(t3*t3) - 72.*akpk*amLam*C5*C6*t1*t4 - 
           72.*amLam*apkk1*C5*C6*t1*t4 - 
           72.*akpk*amLam*C2*C5*C6*t1*t4 + 
           72.*amLam*apkk1*C2*C5*C6*t1*t4 - 
           24.*akpk*amLam*(C5*C5)*C6*t1*t4 + 
           24.*amLam*apkk1*(C5*C5)*C6*t1*t4 - 
           24.*akpk*amLam*C2*(C5*C5)*C6*t1*t4 - 
           24.*amLam*apkk1*C2*(C5*C5)*C6*t1*t4 + 
           24.*akpk*(aml*aml)*amLam*(C5*C5*C5)*C6*t1*t4 + 
           12.*(amk*amk)*(aml*aml)*amLam*(C5*C5*C5)*C6*t1*t4 + 
           12.*(amk*amk)*(aml*aml)*(C5*C5*C5)*C6*Ekaon*t1*t4 - 
           72.*(amk*amk)*C5*C6*Elep*t1*t4 + 
           72.*(amk*amk)*C2*C5*C6*Elep*t1*t4 + 
           24.*(amk*amk)*(C5*C5)*C6*Elep*t1*t4 - 
           24.*(amk*amk)*C2*(C5*C5)*C6*Elep*t1*t4 - 
           72.*(amk*amk)*C5*C6*Enu*t1*t4 - 
           72.*(amk*amk)*C2*C5*C6*Enu*t1*t4 - 
           24.*(amk*amk)*(C5*C5)*C6*Enu*t1*t4 - 
           24.*(amk*amk)*C2*(C5*C5)*C6*Enu*t1*t4 + 
           24.*(amk*amk)*(aml*aml)*(C5*C5*C5)*C6*Enu*t1*t4 + 
           24.*akpk*(amk*amk)*C5*C6*Fm2*t1*t4 + 
           24.*(amk*amk)*apkk1*C5*C6*Fm2*t1*t4 - 
           24.*akpk*(amk*amk)*C2*C5*C6*Fm2*t1*t4 + 
           36.*(amk*amk)*(aml*aml)*C2*C5*C6*Fm2*t1*t4 + 
           24.*(amk*amk)*apkk1*C2*C5*C6*Fm2*t1*t4 + 
           24.*akpk*amLam*C5*C6*Ekaon*Fm2*t1*t4 + 
           24.*amLam*apkk1*C5*C6*Ekaon*Fm2*t1*t4 - 
           24.*akpk*amLam*C2*C5*C6*Ekaon*Fm2*t1*t4 + 
           36.*(aml*aml)*amLam*C2*C5*C6*Ekaon*Fm2*t1*t4 + 
           24.*amLam*apkk1*C2*C5*C6*Ekaon*Fm2*t1*t4 + 
           48.*akpk*amLam*C5*C6*Elep*Fm2*t1*t4 - 
           48.*amLam*apkk1*C5*C6*Enu*Fm2*t1*t4 + 
           12.*akpk*(aml*aml)*amLam*C3*(C5*C5)*C6*t2*t4 + 
           12.*(amk*amk)*(aml*aml)*amLam*C3*(C5*C5)*C6*t2*t4 + 
           12.*(akpk*akpk)*(aml*aml)*amLam*C3*(C5*C5*C5)*C6*t2*t4 + 
           6.*akpk*(amk*amk)*(aml*aml)*amLam*C3*(C5*C5*C5)*C6*t2*t4 - 
           6.*akpk*(aml*aml*aml*aml)*amLam*C3*(C5*C5*C5)*C6*t2*t4 - 
           12.*(amk*amk)*(aml*aml*aml*aml)*amLam*C3*(C5*C5*C5)*C6*t2*t4 - 
           24.*akpk*(aml*aml)*amLam*apkk1*C3*(C5*C5*C5)*C6*t2*t4 - 
           6.*(amk*amk)*(aml*aml)*amLam*apkk1*C3*(C5*C5*C5)*C6*t2*t4 + 
           6.*(aml*aml*aml*aml)*amLam*apkk1*C3*(C5*C5*C5)*C6*t2*t4 + 
           12.*(aml*aml)*amLam*(apkk1*apkk1)*C3*(C5*C5*C5)*C6*t2*t4 + 
           12.*(amk*amk)*(aml*aml)*C3*(C5*C5)*C6*Ekaon*t2*t4 - 
           12.*(amk*amk)*(aml*aml*aml*aml)*C3*(C5*C5*C5)*C6*Ekaon*t2*t4 - 
           12.*akpk*(amk*amk)*(aml*aml)*C3*(C5*C5*C5)*C6*Elep*t2*t4 - 
           6.*(amk*amk*amk*amk)*(aml*aml)*C3*(C5*C5*C5)*C6*Elep*t2*t4 + 
           6.*(amk*amk)*(aml*aml*aml*aml)*C3*(C5*C5*C5)*C6*Elep*t2*t4 + 
           12.*(amk*amk)*(aml*aml)*apkk1*C3*(C5*C5*C5)*C6*Elep*t2*t4 + 
           12.*(amk*amk)*(aml*aml)*C3*(C5*C5)*C6*Enu*t2*t4 + 
           12.*akpk*(amk*amk)*(aml*aml)*C3*(C5*C5*C5)*C6*Enu*t2*t4 + 
           6.*(amk*amk*amk*amk)*(aml*aml)*C3*(C5*C5*C5)*C6*Enu*t2*t4 - 
           6.*(amk*amk)*(aml*aml*aml*aml)*C3*(C5*C5*C5)*C6*Enu*t2*t4 - 
           12.*(amk*amk)*(aml*aml)*apkk1*C3*(C5*C5*C5)*C6*Enu*t2*t4 + 
           36.*akpk*(amk*amk)*amLam*C1*C4*C5*C6*t3*t4 + 
           36.*akpk*(amk*amk)*amSig*C1*C4*C5*C6*t3*t4 - 
           36.*(amk*amk)*amLam*apkk1*C1*C4*C5*C6*t3*t4 - 
           36.*(amk*amk)*amSig*apkk1*C1*C4*C5*C6*t3*t4 - 
           36.*akpk*(amk*amk)*amLam*C1*(C4*C4)*C5*C6*t3*t4 - 
           36.*akpk*(amk*amk)*amSig*C1*(C4*C4)*C5*C6*t3*t4 - 
           36.*(amk*amk)*amLam*apkk1*C1*(C4*C4)*C5*C6*t3*t4 - 
           36.*(amk*amk)*amSig*apkk1*C1*(C4*C4)*C5*C6*t3*t4 + 
           12.*akpk*(amk*amk)*amLam*C1*C4*(C5*C5)*C6*t3*t4 + 
           12.*akpk*(amk*amk)*amSig*C1*C4*(C5*C5)*C6*t3*t4 + 
           12.*(amk*amk)*amLam*apkk1*C1*C4*(C5*C5)*C6*t3*t4 + 
           12.*(amk*amk)*amSig*apkk1*C1*C4*(C5*C5)*C6*t3*t4 - 
           12.*akpk*(amk*amk)*amLam*C1*(C4*C4)*(C5*C5)*C6*t3*t4 - 
           12.*akpk*(amk*amk)*amSig*C1*(C4*C4)*(C5*C5)*C6*t3*t4 + 
           12.*(amk*amk)*amLam*apkk1*C1*(C4*C4)*(C5*C5)*C6*t3*t4 + 
           12.*(amk*amk)*amSig*apkk1*C1*(C4*C4)*(C5*C5)*C6*t3*t4 + 
           12.*akpk*(amk*amk)*(aml*aml)*amLam*C1*C3*(C4*C4)*(C5*C5)*
           C6*t3*t4 + 
           6.*(amk*amk*amk*amk)*(aml*aml)*amLam*C1*C3*(C4*C4)*(C5*C5)*C6*t3*
           t4 + 12.*akpk*(amk*amk)*(aml*aml)*amSig*C1*C3*(C4*C4)*
           (C5*C5)*C6*t3*t4 + 
           6.*(amk*amk*amk*amk)*(aml*aml)*amSig*C1*C3*(C4*C4)*(C5*C5)*C6*t3*
           t4 + 12.*akpk*(amk*amk)*(aml*aml)*amLam*C1*(C4*C4)*
           (C5*C5*C5)*C6*t3*t4 + 
           6.*(amk*amk*amk*amk)*(aml*aml)*amLam*C1*(C4*C4)*(C5*C5*C5)*C6*t3*t4 + 
           12.*akpk*(amk*amk)*(aml*aml)*amSig*C1*(C4*C4)*(C5*C5*C5)*C6*
           t3*t4 + 6.*(amk*amk*amk*amk)*(aml*aml)*amSig*C1*(C4*C4)*(C5*C5*C5)*
           C6*t3*t4 + 
           12.*(akpk*akpk)*(amk*amk)*(aml*aml)*amLam*C1*C3*(C4*C4)*
           (C5*C5*C5)*C6*t3*t4 - 
           6.*akpk*(amk*amk)*(aml*aml*aml*aml)*amLam*C1*C3*(C4*C4)*(C5*C5*C5)*C6*
           t3*t4 - 6.*(amk*amk*amk*amk)*(aml*aml*aml*aml)*amLam*C1*C3*(C4*C4)*
           (C5*C5*C5)*C6*t3*t4 + 
           12.*(akpk*akpk)*(amk*amk)*(aml*aml)*amSig*C1*C3*(C4*C4)*
           (C5*C5*C5)*C6*t3*t4 - 
           6.*akpk*(amk*amk)*(aml*aml*aml*aml)*amSig*C1*C3*(C4*C4)*(C5*C5*C5)*C6*
           t3*t4 - 6.*(amk*amk*amk*amk)*(aml*aml*aml*aml)*amSig*C1*C3*(C4*C4)*
           (C5*C5*C5)*C6*t3*t4 - 
           24.*akpk*(amk*amk)*(aml*aml)*amLam*apkk1*C1*C3*(C4*C4)*
           (C5*C5*C5)*C6*t3*t4 + 
           6.*(amk*amk)*(aml*aml*aml*aml)*amLam*apkk1*C1*C3*(C4*C4)*(C5*C5*C5)*
           C6*t3*t4 - 
           24.*akpk*(amk*amk)*(aml*aml)*amSig*apkk1*C1*C3*(C4*C4)*
           (C5*C5*C5)*C6*t3*t4 + 
           6.*(amk*amk)*(aml*aml*aml*aml)*amSig*apkk1*C1*C3*(C4*C4)*(C5*C5*C5)*
           C6*t3*t4 + 
           12.*(amk*amk)*(aml*aml)*amLam*(apkk1*apkk1)*C1*C3*(C4*C4)*
           (C5*C5*C5)*C6*t3*t4 + 
           12.*(amk*amk)*(aml*aml)*amSig*(apkk1*apkk1)*C1*C3*(C4*C4)*
           (C5*C5*C5)*C6*t3*t4 + 
           72.*akpk*amLam*amSig*C1*C4*C5*C6*Ekaon*t3*t4 - 
           72.*amLam*amSig*apkk1*C1*C4*C5*C6*Ekaon*t3*t4 - 
           72.*akpk*amLam*amSig*C1*(C4*C4)*C5*C6*Ekaon*t3*
           t4 - 72.*amLam*amSig*apkk1*C1*(C4*C4)*C5*C6*Ekaon*
           t3*t4 + 24.*akpk*amLam*amSig*C1*C4*(C5*C5)*C6*
           Ekaon*t3*t4 + 
           24.*amLam*amSig*apkk1*C1*C4*(C5*C5)*C6*Ekaon*t3*
           t4 - 24.*akpk*amLam*amSig*C1*(C4*C4)*(C5*C5)*C6*
           Ekaon*t3*t4 + 
           24.*amLam*amSig*apkk1*C1*(C4*C4)*(C5*C5)*C6*Ekaon*t3*
           t4 + 6.*(amk*amk*amk*amk)*(aml*aml)*C1*C3*(C4*C4)*(C5*C5)*C6*
           Ekaon*t3*t4 + 
           24.*akpk*(aml*aml)*amLam*amSig*C1*C3*(C4*C4)*(C5*C5)*C6*
           Ekaon*t3*t4 + 
           6.*(amk*amk)*(aml*aml)*amLam*amSig*C1*C3*(C4*C4)*(C5*C5)*
           C6*Ekaon*t3*t4 + 
           6.*(amk*amk*amk*amk)*(aml*aml)*C1*(C4*C4)*(C5*C5*C5)*C6*Ekaon*t3*t4 + 
           24.*akpk*(aml*aml)*amLam*amSig*C1*(C4*C4)*(C5*C5*C5)*C6*
           Ekaon*t3*t4 + 
           6.*(amk*amk)*(aml*aml)*amLam*amSig*C1*(C4*C4)*(C5*C5*C5)*C6*
           Ekaon*t3*t4 - 
           6.*(amk*amk*amk*amk)*(aml*aml*aml*aml)*C1*C3*(C4*C4)*(C5*C5*C5)*C6*Ekaon*t3*
           t4 + 24.*(akpk*akpk)*(aml*aml)*amLam*amSig*C1*C3*(C4*C4)*
           (C5*C5*C5)*C6*Ekaon*t3*t4 - 
           12.*akpk*(aml*aml*aml*aml)*amLam*amSig*C1*C3*(C4*C4)*(C5*C5*C5)*C6*
           Ekaon*t3*t4 - 
           6.*(amk*amk)*(aml*aml*aml*aml)*amLam*amSig*C1*C3*(C4*C4)*(C5*C5*C5)*
           C6*Ekaon*t3*t4 - 
           48.*akpk*(aml*aml)*amLam*amSig*apkk1*C1*C3*(C4*C4)*
           (C5*C5*C5)*C6*Ekaon*t3*t4 + 
           12.*(aml*aml*aml*aml)*amLam*amSig*apkk1*C1*C3*(C4*C4)*(C5*C5*C5)*
           C6*Ekaon*t3*t4 + 
           24.*(aml*aml)*amLam*amSig*(apkk1*apkk1)*C1*C3*(C4*C4)*
           (C5*C5*C5)*C6*Ekaon*t3*t4 - 
           36.*(amk*amk*amk*amk)*C1*C4*C5*C6*Elep*t3*t4 + 
           36.*(amk*amk)*amLam*amSig*C1*C4*C5*C6*Elep*t3*t4 - 
           36.*(amk*amk*amk*amk)*C1*(C4*C4)*C5*C6*Elep*t3*t4 + 
           36.*(amk*amk)*amLam*amSig*C1*(C4*C4)*C5*C6*Elep*t3*
           t4 + 12.*(amk*amk*amk*amk)*C1*C4*(C5*C5)*C6*Elep*t3*t4 - 
           12.*(amk*amk)*amLam*amSig*C1*C4*(C5*C5)*C6*Elep*t3*
           t4 + 12.*(amk*amk*amk*amk)*C1*(C4*C4)*(C5*C5)*C6*Elep*t3*t4 - 
           12.*(amk*amk)*amLam*amSig*C1*(C4*C4)*(C5*C5)*C6*Elep*t3*
           t4 - 12.*akpk*(amk*amk*amk*amk)*(aml*aml)*C1*C3*(C4*C4)*(C5*C5*C5)*
           C6*Elep*t3*t4 + 
           6.*(amk*amk*amk*amk)*(aml*aml*aml*aml)*C1*C3*(C4*C4)*(C5*C5*C5)*C6*Elep*t3*
           t4 + 12.*akpk*(amk*amk)*(aml*aml)*amLam*amSig*C1*C3*
           (C4*C4)*(C5*C5*C5)*C6*Elep*t3*t4 - 
           6.*(amk*amk)*(aml*aml*aml*aml)*amLam*amSig*C1*C3*(C4*C4)*(C5*C5*C5)*
           C6*Elep*t3*t4 + 
           12.*(amk*amk*amk*amk)*(aml*aml)*apkk1*C1*C3*(C4*C4)*(C5*C5*C5)*C6*
           Elep*t3*t4 - 
           12.*(amk*amk)*(aml*aml)*amLam*amSig*apkk1*C1*C3*(C4*C4)*
           (C5*C5*C5)*C6*Elep*t3*t4 + 
           36.*(amk*amk*amk*amk)*C1*C4*C5*C6*Enu*t3*t4 - 
           36.*(amk*amk)*amLam*amSig*C1*C4*C5*C6*Enu*t3*t4 - 
           36.*(amk*amk*amk*amk)*C1*(C4*C4)*C5*C6*Enu*t3*t4 + 
           36.*(amk*amk)*amLam*amSig*C1*(C4*C4)*C5*C6*Enu*t3*
           t4 + 12.*(amk*amk*amk*amk)*C1*C4*(C5*C5)*C6*Enu*t3*t4 - 
           12.*(amk*amk)*amLam*amSig*C1*C4*(C5*C5)*C6*Enu*t3*
           t4 - 12.*(amk*amk*amk*amk)*C1*(C4*C4)*(C5*C5)*C6*Enu*t3*t4 + 
           12.*(amk*amk)*amLam*amSig*C1*(C4*C4)*(C5*C5)*C6*Enu*t3*
           t4 + 12.*(amk*amk*amk*amk)*(aml*aml)*C1*C3*(C4*C4)*(C5*C5)*C6*Enu*
           t3*t4 - 12.*(amk*amk)*(aml*aml)*amLam*amSig*C1*C3*
           (C4*C4)*(C5*C5)*C6*Enu*t3*t4 + 
           12.*(amk*amk*amk*amk)*(aml*aml)*C1*(C4*C4)*(C5*C5*C5)*C6*Enu*t3*t4 - 
           12.*(amk*amk)*(aml*aml)*amLam*amSig*C1*(C4*C4)*(C5*C5*C5)*C6*
           Enu*t3*t4 + 
           12.*akpk*(amk*amk*amk*amk)*(aml*aml)*C1*C3*(C4*C4)*(C5*C5*C5)*C6*Enu*
           t3*t4 - 6.*(amk*amk*amk*amk)*(aml*aml*aml*aml)*C1*C3*(C4*C4)*(C5*C5*C5)*C6*
           Enu*t3*t4 - 
           12.*akpk*(amk*amk)*(aml*aml)*amLam*amSig*C1*C3*(C4*C4)*
           (C5*C5*C5)*C6*Enu*t3*t4 + 
           6.*(amk*amk)*(aml*aml*aml*aml)*amLam*amSig*C1*C3*(C4*C4)*(C5*C5*C5)*
           C6*Enu*t3*t4 - 
           12.*(amk*amk*amk*amk)*(aml*aml)*apkk1*C1*C3*(C4*C4)*(C5*C5*C5)*C6*Enu*
           t3*t4 + 12.*(amk*amk)*(aml*aml)*amLam*amSig*apkk1*C1*
           C3*(C4*C4)*(C5*C5*C5)*C6*Enu*t3*t4 + 
           36.*akpk*(amk*amk*amk*amk)*C1*C4*C5*C6*Fm1*t3*t4 - 
           54.*(amk*amk*amk*amk)*(aml*aml)*C1*C4*C5*C6*Fm1*t3*t4 + 
           36.*akpk*(amk*amk)*amLam*amSig*C1*C4*C5*C6*Fm1*t3*
           t4 - 54.*(amk*amk)*(aml*aml)*amLam*amSig*C1*C4*C5*C6*
           Fm1*t3*t4 - 
           36.*(amk*amk*amk*amk)*apkk1*C1*C4*C5*C6*Fm1*t3*t4 - 
           36.*(amk*amk)*amLam*amSig*apkk1*C1*C4*C5*C6*Fm1*t3*
           t4 - 12.*akpk*(amk*amk*amk*amk)*C1*C4*(C5*C5)*C6*Fm1*t3*t4 - 
           12.*akpk*(amk*amk)*amLam*amSig*C1*C4*(C5*C5)*C6*Fm1*
           t3*t4 - 12.*(amk*amk*amk*amk)*apkk1*C1*C4*(C5*C5)*C6*Fm1*t3*
           t4 - 12.*(amk*amk)*amLam*amSig*apkk1*C1*C4*(C5*C5)*
           C6*Fm1*t3*t4 + 
           36.*akpk*(amk*amk)*amLam*C1*C4*C5*C6*Ekaon*Fm1*t3*
           t4 - 54.*(amk*amk)*(aml*aml)*amLam*C1*C4*C5*C6*Ekaon*
           Fm1*t3*t4 + 
           36.*akpk*(amk*amk)*amSig*C1*C4*C5*C6*Ekaon*Fm1*t3*
           t4 - 54.*(amk*amk)*(aml*aml)*amSig*C1*C4*C5*C6*Ekaon*
           Fm1*t3*t4 - 
           36.*(amk*amk)*amLam*apkk1*C1*C4*C5*C6*Ekaon*Fm1*t3*
           t4 - 36.*(amk*amk)*amSig*apkk1*C1*C4*C5*C6*Ekaon*
           Fm1*t3*t4 - 
           12.*akpk*(amk*amk)*amLam*C1*C4*(C5*C5)*C6*Ekaon*Fm1*
           t3*t4 - 12.*akpk*(amk*amk)*amSig*C1*C4*(C5*C5)*C6*
           Ekaon*Fm1*t3*t4 - 
           12.*(amk*amk)*amLam*apkk1*C1*C4*(C5*C5)*C6*Ekaon*Fm1*
           t3*t4 - 12.*(amk*amk)*amSig*apkk1*C1*C4*(C5*C5)*C6*
           Ekaon*Fm1*t3*t4 + 
           24.*akpk*(amk*amk)*amLam*C1*C4*(C5*C5)*C6*Elep*Fm1*t3*
           t4 - 24.*akpk*(amk*amk)*amSig*C1*C4*(C5*C5)*C6*Elep*
           Fm1*t3*t4 - 
           24.*(amk*amk)*amLam*apkk1*C1*C4*(C5*C5)*C6*Enu*Fm1*t3*
           t4 + 24.*(amk*amk)*amSig*apkk1*C1*C4*(C5*C5)*C6*Enu*
           Fm1*t3*t4 + 
           12.*akpk*(amk*amk*amk*amk)*C1*C4*C5*C6*Fm2*t3*t4 - 
           18.*(amk*amk*amk*amk)*(aml*aml)*C1*C4*C5*C6*Fm2*t3*t4 + 
           12.*akpk*(amk*amk)*amLam*amSig*C1*C4*C5*C6*Fm2*t3*
           t4 - 18.*(amk*amk)*(aml*aml)*amLam*amSig*C1*C4*C5*C6*
           Fm2*t3*t4 - 
           12.*(amk*amk*amk*amk)*apkk1*C1*C4*C5*C6*Fm2*t3*t4 - 
           12.*(amk*amk)*amLam*amSig*apkk1*C1*C4*C5*C6*Fm2*t3*
           t4 + 12.*akpk*(amk*amk*amk*amk)*C1*(C4*C4)*C5*C6*Fm2*t3*t4 + 
           12.*akpk*(amk*amk)*amLam*amSig*C1*(C4*C4)*C5*C6*Fm2*
           t3*t4 + 12.*(amk*amk*amk*amk)*apkk1*C1*(C4*C4)*C5*C6*Fm2*t3*
           t4 + 12.*(amk*amk)*amLam*amSig*apkk1*C1*(C4*C4)*C5*
           C6*Fm2*t3*t4 + 
           12.*akpk*(amk*amk)*amLam*C1*C4*C5*C6*Ekaon*Fm2*t3*
           t4 - 18.*(amk*amk)*(aml*aml)*amLam*C1*C4*C5*C6*Ekaon*
           Fm2*t3*t4 + 
           12.*akpk*(amk*amk)*amSig*C1*C4*C5*C6*Ekaon*Fm2*t3*
           t4 - 18.*(amk*amk)*(aml*aml)*amSig*C1*C4*C5*C6*Ekaon*
           Fm2*t3*t4 - 
           12.*(amk*amk)*amLam*apkk1*C1*C4*C5*C6*Ekaon*Fm2*t3*
           t4 - 12.*(amk*amk)*amSig*apkk1*C1*C4*C5*C6*Ekaon*
           Fm2*t3*t4 + 
           12.*akpk*(amk*amk)*amLam*C1*(C4*C4)*C5*C6*Ekaon*Fm2*
           t3*t4 + 12.*akpk*(amk*amk)*amSig*C1*(C4*C4)*C5*C6*
           Ekaon*Fm2*t3*t4 + 
           12.*(amk*amk)*amLam*apkk1*C1*(C4*C4)*C5*C6*Ekaon*Fm2*
           t3*t4 + 12.*(amk*amk)*amSig*apkk1*C1*(C4*C4)*C5*C6*
           Ekaon*Fm2*t3*t4 + 
           24.*akpk*(amk*amk)*amLam*C1*(C4*C4)*C5*C6*Elep*Fm2*t3*
           t4 - 24.*akpk*(amk*amk)*amSig*C1*(C4*C4)*C5*C6*Elep*
           Fm2*t3*t4 - 
           24.*(amk*amk)*amLam*apkk1*C1*(C4*C4)*C5*C6*Enu*Fm2*t3*
           t4 + 24.*(amk*amk)*amSig*apkk1*C1*(C4*C4)*C5*C6*Enu*
           Fm2*t3*t4 - 
           24.*(akpk*akpk)*(amk*amk)*amLam*C1*C4*C5*C6*Fm1*Fm2*t3*
           t4 + 18.*akpk*(amk*amk)*(aml*aml)*amLam*C1*C4*C5*C6*
           Fm1*Fm2*t3*t4 + 
           6.*(amk*amk*amk*amk)*(aml*aml)*amLam*C1*C4*C5*C6*Fm1*Fm2*t3*
           t4 - 24.*(akpk*akpk)*(amk*amk)*amSig*C1*C4*C5*C6*Fm1*
           Fm2*t3*t4 + 
           18.*akpk*(amk*amk)*(aml*aml)*amSig*C1*C4*C5*C6*Fm1*Fm2*
           t3*t4 + 6.*(amk*amk*amk*amk)*(aml*aml)*amSig*C1*C4*C5*C6*Fm1*
           Fm2*t3*t4 - 
           18.*(amk*amk)*(aml*aml)*amLam*apkk1*C1*C4*C5*C6*Fm1*
           Fm2*t3*t4 - 
           18.*(amk*amk)*(aml*aml)*amSig*apkk1*C1*C4*C5*C6*Fm1*
           Fm2*t3*t4 - 
           24.*(amk*amk)*amLam*(apkk1*apkk1)*C1*C4*C5*C6*Fm1*Fm2*t3*
           t4 - 24.*(amk*amk)*amSig*(apkk1*apkk1)*C1*C4*C5*C6*Fm1*
           Fm2*t3*t4 + 
           6.*(amk*amk*amk*amk)*(aml*aml)*C1*C4*C5*C6*Ekaon*Fm1*Fm2*t3*
           t4 - 48.*(akpk*akpk)*amLam*amSig*C1*C4*C5*C6*Ekaon*
           Fm1*Fm2*t3*t4 + 
           36.*akpk*(aml*aml)*amLam*amSig*C1*C4*C5*C6*Ekaon*
           Fm1*Fm2*t3*t4 + 
           6.*(amk*amk)*(aml*aml)*amLam*amSig*C1*C4*C5*C6*Ekaon*
           Fm1*Fm2*t3*t4 - 
           36.*(aml*aml)*amLam*amSig*apkk1*C1*C4*C5*C6*Ekaon*
           Fm1*Fm2*t3*t4 - 
           48.*amLam*amSig*(apkk1*apkk1)*C1*C4*C5*C6*Ekaon*Fm1*
           Fm2*t3*t4 - 
           18.*(amk*amk*amk*amk)*(aml*aml)*C1*C4*C5*C6*Elep*Fm1*Fm2*t3*
           t4 + 18.*(amk*amk)*(aml*aml)*amLam*amSig*C1*C4*C5*C6*
           Elep*Fm1*Fm2*t3*t4 - 
           24.*(amk*amk*amk*amk)*apkk1*C1*C4*C5*C6*Elep*Fm1*Fm2*t3*
           t4 + 24.*(amk*amk)*amLam*amSig*apkk1*C1*C4*C5*C6*
           Elep*Fm1*Fm2*t3*t4 - 
           24.*akpk*(amk*amk*amk*amk)*C1*C4*C5*C6*Enu*Fm1*Fm2*t3*t4 + 
           18.*(amk*amk*amk*amk)*(aml*aml)*C1*C4*C5*C6*Enu*Fm1*Fm2*t3*
           t4 + 24.*akpk*(amk*amk)*amLam*amSig*C1*C4*C5*C6*
           Enu*Fm1*Fm2*t3*t4 - 
           18.*(amk*amk)*(aml*aml)*amLam*amSig*C1*C4*C5*C6*Enu*
           Fm1*Fm2*t3*t4 + 
           36.*akpk*(amk*amk)*amLam*(C5*C5)*(C6*C6)*(t4*t4) - 
           36.*(amk*amk)*amLam*apkk1*(C5*C5)*(C6*C6)*(t4*t4) + 
           24.*akpk*(amk*amk)*amLam*(C5*C5*C5)*(C6*C6)*(t4*t4) + 
           24.*(amk*amk)*amLam*apkk1*(C5*C5*C5)*(C6*C6)*(t4*t4) + 
           4.*akpk*(amk*amk)*amLam*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) - 
           4.*(amk*amk)*amLam*apkk1*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) - 
           8.*akpk*(amk*amk)*(aml*aml)*amLam*(C5*C5*C5*C5*C5)*(C6*C6)*(t4*t4) - 
           4.*(amk*amk*amk*amk)*(aml*aml)*amLam*(C5*C5*C5*C5*C5)*(C6*C6)*(t4*t4) - 
           4.*(akpk*akpk)*(amk*amk)*(aml*aml)*amLam*(C5*C5*C5*C5*C5*C5)*(C6*C6)*
           (t4*t4) + 2.*akpk*(amk*amk)*(aml*aml*aml*aml)*amLam*(C5*C5*C5*C5*C5*C5)*(C6*C6)*
           (t4*t4) + 2.*(amk*amk*amk*amk)*(aml*aml*aml*aml)*amLam*(C5*C5*C5*C5*C5*C5)*(C6*C6)*
           (t4*t4) + 8.*akpk*(amk*amk)*(aml*aml)*amLam*apkk1*(C5*C5*C5*C5*C5*C5)*
           (C6*C6)*(t4*t4) - 
           2.*(amk*amk)*(aml*aml*aml*aml)*amLam*apkk1*(C5*C5*C5*C5*C5*C5)*(C6*C6)*(t4*t4) - 
           4.*(amk*amk)*(aml*aml)*amLam*(apkk1*apkk1)*(C5*C5*C5*C5*C5*C5)*(C6*C6)*
           (t4*t4) + 36.*akpk*(amLam*amLam)*(C5*C5)*(C6*C6)*Ekaon*
           (t4*t4) - 36.*(amLam*amLam)*apkk1*(C5*C5)*(C6*C6)*Ekaon*
           (t4*t4) + 24.*akpk*(amLam*amLam)*(C5*C5*C5)*(C6*C6)*Ekaon*
           (t4*t4) + 24.*(amLam*amLam)*apkk1*(C5*C5*C5)*(C6*C6)*Ekaon*
           (t4*t4) + 4.*akpk*(amLam*amLam)*(C5*C5*C5*C5)*(C6*C6)*Ekaon*
           (t4*t4) - 4.*(amLam*amLam)*apkk1*(C5*C5*C5*C5)*(C6*C6)*Ekaon*
           (t4*t4) - 2.*(amk*amk*amk*amk)*(aml*aml)*(C5*C5*C5*C5*C5)*(C6*C6)*Ekaon*
           (t4*t4) - 8.*akpk*(aml*aml)*(amLam*amLam)*(C5*C5*C5*C5*C5)*(C6*C6)*
           Ekaon*(t4*t4) - 
           2.*(amk*amk)*(aml*aml)*(amLam*amLam)*(C5*C5*C5*C5*C5)*(C6*C6)*Ekaon*
           (t4*t4) + (amk*amk*amk*amk)*(aml*aml*aml*aml)*(C5*C5*C5*C5*C5*C5)*(C6*C6)*Ekaon*(t4*t4) - 
           4.*(akpk*akpk)*(aml*aml)*(amLam*amLam)*(C5*C5*C5*C5*C5*C5)*(C6*C6)*Ekaon*
           (t4*t4) + 2.*akpk*(aml*aml*aml*aml)*(amLam*amLam)*(C5*C5*C5*C5*C5*C5)*(C6*C6)*
           Ekaon*(t4*t4) + 
           (amk*amk)*(aml*aml*aml*aml)*(amLam*amLam)*(C5*C5*C5*C5*C5*C5)*(C6*C6)*Ekaon*(t4*t4) + 
           8.*akpk*(aml*aml)*(amLam*amLam)*apkk1*(C5*C5*C5*C5*C5*C5)*(C6*C6)*Ekaon*
           (t4*t4) - 2.*(aml*aml*aml*aml)*(amLam*amLam)*apkk1*(C5*C5*C5*C5*C5*C5)*(C6*C6)*
           Ekaon*(t4*t4) - 
           4.*(aml*aml)*(amLam*amLam)*(apkk1*apkk1)*(C5*C5*C5*C5*C5*C5)*(C6*C6)*Ekaon*
           (t4*t4) - 18.*(amk*amk*amk*amk)*(C5*C5)*(C6*C6)*Elep*(t4*t4) + 
           18.*(amk*amk)*(amLam*amLam)*(C5*C5)*(C6*C6)*Elep*(t4*t4) + 
           12.*(amk*amk*amk*amk)*(C5*C5*C5)*(C6*C6)*Elep*(t4*t4) - 
           12.*(amk*amk)*(amLam*amLam)*(C5*C5*C5)*(C6*C6)*Elep*(t4*t4) - 
           2.*(amk*amk*amk*amk)*(C5*C5*C5*C5)*(C6*C6)*Elep*(t4*t4) + 
           2.*(amk*amk)*(amLam*amLam)*(C5*C5*C5*C5)*(C6*C6)*Elep*(t4*t4) + 
           2.*akpk*(amk*amk*amk*amk)*(aml*aml)*(C5*C5*C5*C5*C5*C5)*(C6*C6)*Elep*(t4*t4) - 
           1.*(amk*amk*amk*amk)*(aml*aml*aml*aml)*(C5*C5*C5*C5*C5*C5)*(C6*C6)*Elep*(t4*t4) - 
           2.*akpk*(amk*amk)*(aml*aml)*(amLam*amLam)*(C5*C5*C5*C5*C5*C5)*(C6*C6)*Elep*
           (t4*t4) + (amk*amk)*(aml*aml*aml*aml)*(amLam*amLam)*(C5*C5*C5*C5*C5*C5)*(C6*C6)*Elep*
           (t4*t4) - 2.*(amk*amk*amk*amk)*(aml*aml)*apkk1*(C5*C5*C5*C5*C5*C5)*(C6*C6)*Elep*
           (t4*t4) + 2.*(amk*amk)*(aml*aml)*(amLam*amLam)*apkk1*(C5*C5*C5*C5*C5*C5)*
           (C6*C6)*Elep*(t4*t4) + 
           18.*(amk*amk*amk*amk)*(C5*C5)*(C6*C6)*Enu*(t4*t4) - 
           18.*(amk*amk)*(amLam*amLam)*(C5*C5)*(C6*C6)*Enu*(t4*t4) + 
           12.*(amk*amk*amk*amk)*(C5*C5*C5)*(C6*C6)*Enu*(t4*t4) - 
           12.*(amk*amk)*(amLam*amLam)*(C5*C5*C5)*(C6*C6)*Enu*(t4*t4) + 
           2.*(amk*amk*amk*amk)*(C5*C5*C5*C5)*(C6*C6)*Enu*(t4*t4) - 
           2.*(amk*amk)*(amLam*amLam)*(C5*C5*C5*C5)*(C6*C6)*Enu*(t4*t4) - 
           4.*(amk*amk*amk*amk)*(aml*aml)*(C5*C5*C5*C5*C5)*(C6*C6)*Enu*(t4*t4) + 
           4.*(amk*amk)*(aml*aml)*(amLam*amLam)*(C5*C5*C5*C5*C5)*(C6*C6)*Enu*
           (t4*t4) - 2.*akpk*(amk*amk*amk*amk)*(aml*aml)*(C5*C5*C5*C5*C5*C5)*(C6*C6)*Enu*
           (t4*t4) + (amk*amk*amk*amk)*(aml*aml*aml*aml)*(C5*C5*C5*C5*C5*C5)*(C6*C6)*Enu*(t4*t4) + 
           2.*akpk*(amk*amk)*(aml*aml)*(amLam*amLam)*(C5*C5*C5*C5*C5*C5)*(C6*C6)*Enu*
           (t4*t4) - 1.*(amk*amk)*(aml*aml*aml*aml)*(amLam*amLam)*(C5*C5*C5*C5*C5*C5)*(C6*C6)*
           Enu*(t4*t4) + 
           2.*(amk*amk*amk*amk)*(aml*aml)*apkk1*(C5*C5*C5*C5*C5*C5)*(C6*C6)*Enu*(t4*t4) - 
           2.*(amk*amk)*(aml*aml)*(amLam*amLam)*apkk1*(C5*C5*C5*C5*C5*C5)*(C6*C6)*Enu*
           (t4*t4) + 12.*akpk*(amk*amk*amk*amk)*(C5*C5)*(C6*C6)*Fm2*(t4*t4) - 
           18.*(amk*amk*amk*amk)*(aml*aml)*(C5*C5)*(C6*C6)*Fm2*(t4*t4) + 
           12.*akpk*(amk*amk)*(amLam*amLam)*(C5*C5)*(C6*C6)*Fm2*(t4*t4) - 
           18.*(amk*amk)*(aml*aml)*(amLam*amLam)*(C5*C5)*(C6*C6)*Fm2*
           (t4*t4) - 12.*(amk*amk*amk*amk)*apkk1*(C5*C5)*(C6*C6)*Fm2*
           (t4*t4) - 12.*(amk*amk)*(amLam*amLam)*apkk1*(C5*C5)*(C6*C6)*
           Fm2*(t4*t4) - 
           4.*akpk*(amk*amk*amk*amk)*(C5*C5*C5)*(C6*C6)*Fm2*(t4*t4) - 
           4.*akpk*(amk*amk)*(amLam*amLam)*(C5*C5*C5)*(C6*C6)*Fm2*(t4*t4) - 
           4.*(amk*amk*amk*amk)*apkk1*(C5*C5*C5)*(C6*C6)*Fm2*(t4*t4) - 
           4.*(amk*amk)*(amLam*amLam)*apkk1*(C5*C5*C5)*(C6*C6)*Fm2*(t4*t4) + 
           24.*akpk*(amk*amk)*amLam*(C5*C5)*(C6*C6)*Ekaon*Fm2*
           (t4*t4) - 36.*(amk*amk)*(aml*aml)*amLam*(C5*C5)*(C6*C6)*
           Ekaon*Fm2*(t4*t4) - 
           24.*(amk*amk)*amLam*apkk1*(C5*C5)*(C6*C6)*Ekaon*Fm2*
           (t4*t4) - 8.*akpk*(amk*amk)*amLam*(C5*C5*C5)*(C6*C6)*Ekaon*
           Fm2*(t4*t4) - 
           8.*(amk*amk)*amLam*apkk1*(C5*C5*C5)*(C6*C6)*Ekaon*Fm2*
           (t4*t4) - 8.*(akpk*akpk)*(amk*amk)*amLam*(C5*C5)*(C6*C6)*
           (Fm2*Fm2)*(t4*t4) + 
           6.*akpk*(amk*amk)*(aml*aml)*amLam*(C5*C5)*(C6*C6)*(Fm2*Fm2)*
           (t4*t4) + 2.*(amk*amk*amk*amk)*(aml*aml)*amLam*(C5*C5)*(C6*C6)*
           (Fm2*Fm2)*(t4*t4) - 
           6.*(amk*amk)*(aml*aml)*amLam*apkk1*(C5*C5)*(C6*C6)*(Fm2*Fm2)*
           (t4*t4) - 8.*(amk*amk)*amLam*(apkk1*apkk1)*(C5*C5)*(C6*C6)*
           (Fm2*Fm2)*(t4*t4) + 
           (amk*amk*amk*amk)*(aml*aml)*(C5*C5)*(C6*C6)*Ekaon*(Fm2*Fm2)*(t4*t4) - 
           8.*(akpk*akpk)*(amLam*amLam)*(C5*C5)*(C6*C6)*Ekaon*(Fm2*Fm2)*
           (t4*t4) + 6.*akpk*(aml*aml)*(amLam*amLam)*(C5*C5)*(C6*C6)*
           Ekaon*(Fm2*Fm2)*(t4*t4) + 
           (amk*amk)*(aml*aml)*(amLam*amLam)*(C5*C5)*(C6*C6)*Ekaon*(Fm2*Fm2)*
           (t4*t4) - 6.*(aml*aml)*(amLam*amLam)*apkk1*(C5*C5)*(C6*C6)*
           Ekaon*(Fm2*Fm2)*(t4*t4) - 
           8.*(amLam*amLam)*(apkk1*apkk1)*(C5*C5)*(C6*C6)*Ekaon*(Fm2*Fm2)*
           (t4*t4) - 3.*(amk*amk*amk*amk)*(aml*aml)*(C5*C5)*(C6*C6)*Elep*
           (Fm2*Fm2)*(t4*t4) + 
           3.*(amk*amk)*(aml*aml)*(amLam*amLam)*(C5*C5)*(C6*C6)*Elep*
           (Fm2*Fm2)*(t4*t4) - 
           4.*(amk*amk*amk*amk)*apkk1*(C5*C5)*(C6*C6)*Elep*(Fm2*Fm2)*(t4*t4) + 
           4.*(amk*amk)*(amLam*amLam)*apkk1*(C5*C5)*(C6*C6)*Elep*(Fm2*Fm2)*
           (t4*t4) - 4.*akpk*(amk*amk*amk*amk)*(C5*C5)*(C6*C6)*Enu*(Fm2*Fm2)*
           (t4*t4) + 3.*(amk*amk*amk*amk)*(aml*aml)*(C5*C5)*(C6*C6)*Enu*(Fm2*Fm2)*
           (t4*t4) + 4.*akpk*(amk*amk)*(amLam*amLam)*(C5*C5)*(C6*C6)*Enu*
           (Fm2*Fm2)*(t4*t4) - 
           3.*(amk*amk)*(aml*aml)*(amLam*amLam)*(C5*C5)*(C6*C6)*Enu*(Fm2*Fm2)*
           (t4*t4) + 2.*(am*am*am)*
           ((amk*amk)*(9.*(C1*C1)*(C4*C4)*
           (1. + 
           (C4*C4)*
           (-1. + (aml*aml*aml*aml)*(C3*C3) - 
           1.*(aml*aml)*C3*
           (2. + C3*((Elep - 1.*Enu)*(Elep - 1.*Enu)))) - 
           4.*C4*(Elep + Enu)*Fm1 + 
           2.*(aml*aml)*(Fm1*Fm1) - 2.*(Elep*Elep)*(Fm1*Fm1) - 
           2.*(Enu*Enu)*(Fm1*Fm1))*(t3*t3) - 
           6.*C1*C4*C5*C6*
           (-3. - 2.*C5*Elep*Fm1 - 2.*C5*Enu*Fm1 - 
           2.*(aml*aml)*Fm1*Fm2 + 
           2.*(Elep*Elep)*Fm1*Fm2 + 
           2.*(Enu*Enu)*Fm1*Fm2 + 
           C4*
           (-1.*(1. + (aml*aml)*C3)*C5 + 
           (aml*aml)*(C5*C5)*
           (-1. + (aml*aml)*C3 - 
           1.*C3*((Elep - 1.*Enu)*(Elep - 1.*Enu))) + 
           2.*(Elep + Enu)*Fm2))*t3*t4 + 
           (C5*C5)*(C6*C6)*
           (9. - 1.*(C5*C5) - 2.*(aml*aml)*(C5*C5*C5) + 
           (aml*aml)*(C5*C5*C5*C5)*
           ((aml*aml) - 1.*((Elep - 1.*Enu)*(Elep - 1.*Enu))) + 
           4.*C5*(Elep + Enu)*Fm2 + 
           2.*(aml*aml)*(Fm2*Fm2) - 2.*(Elep*Elep)*(Fm2*Fm2) - 
           2.*(Enu*Enu)*(Fm2*Fm2))*(t4*t4)) - 
           1.*Ekaon*
           ((aml*aml*aml*aml)*Ekaon*
           ((-3.*C1*C3*(C4*C4)*t3 + (C5*C5*C5)*C6*t4)*(-3.*C1*C3*(C4*C4)*t3 + (C5*C5*C5)*C6*t4)) + 
           (aml*aml)*
           (-9.*(C1*C1)*(C4*C4)*
           (2.*C3*(C4*C4)*Ekaon + 
           2.*(C3*C3)*(C4*C4)*(Elep - 1.*Enu)*
           (akpk - 1.*apkk1 + 2.*Ekaon*Elep - 
           2.*Ekaon*Enu) - 1.*Ekaon*(Fm1*Fm1))*
           (t3*t3) + 
           6.*C1*C4*C5*C6*
           (C3*C4*C5*
           (Ekaon*
           (1. + 4.*C5*((Elep - 1.*Enu)*(Elep - 1.*Enu))) + 
           2.*(akpk - 1.*apkk1)*C5*
           (Elep - 1.*Enu)) + 
           Ekaon*(C4*(C5*C5) + Fm1*Fm2))*t3*t4 + 
           (C5*C5)*(C6*C6)*
           (-2.*(C5*C5*C5)*Ekaon - 
           2.*(C5*C5*C5*C5)*(Elep - 1.*Enu)*
           (akpk - 1.*apkk1 + 2.*Ekaon*Elep - 
           2.*Ekaon*Enu) + Ekaon*(Fm2*Fm2))*(t4*t4)
           ) + 
           4.*(3.*C1*C4*Fm1*t3 + C5*C6*Fm2*t4)*
           (apkk1*
           (3.*C1*C4*(1. + C4 + Elep*Fm1)*t3 + 
           C5*C6*(3. - 1.*C5 + Elep*Fm2)*t4) + 
           akpk*
           (3.*C1*C4*(-1. + C4 + Enu*Fm1)*t3 - 
           1.*C5*C6*(3. + C5 - 1.*Enu*Fm2)*t4)
           - 2.*Ekaon*
           (3.*C1*C4*
           ((1. + C4)*Elep + (Elep*Elep)*Fm1 + 
           Enu*(-1. + C4 + Enu*Fm1))*t3 + 
           C5*C6*
           (-1.*(-3. + C5)*Elep + 
           (Elep*Elep)*Fm2 + 
           Enu*(-3. - 1.*C5 + Enu*Fm2))*t4))))
           + (am*am)*(9.*(C1*C1)*(C4*C4)*
           ((amk*amk)*
           (2.*Elep + 4.*C4*Elep + 2.*(C4*C4)*Elep - 
           6.*akpk*(aml*aml)*(C3*C3)*(C4*C4)*Elep + 
           (aml*aml*aml*aml)*(C3*C3)*(C4*C4)*Elep + 
           6.*(aml*aml)*apkk1*(C3*C3)*(C4*C4)*Elep - 
           2.*Enu + 4.*C4*Enu - 2.*(C4*C4)*Enu + 
           4.*(aml*aml)*C3*(C4*C4)*Enu + 
           6.*akpk*(aml*aml)*(C3*C3)*(C4*C4)*Enu - 
           1.*(aml*aml*aml*aml)*(C3*C3)*(C4*C4)*Enu - 
           6.*(aml*aml)*apkk1*(C3*C3)*(C4*C4)*Enu - 
           4.*akpk*Fm1 - 6.*(aml*aml)*Fm1 + 
           4.*apkk1*Fm1 + 12.*akpk*C4*Fm1 + 
           12.*apkk1*C4*Fm1 + 
           3.*(aml*aml)*Elep*(Fm1*Fm1) + 
           12.*apkk1*Elep*(Fm1*Fm1) + 
           12.*akpk*Enu*(Fm1*Fm1) - 
           3.*(aml*aml)*Enu*(Fm1*Fm1) + 
           4.*amSig*
           (1. + 
           (C4*C4)*
           (-1. + (aml*aml*aml*aml)*(C3*C3) - 
           1.*(aml*aml)*C3*
           (2. + C3*((Elep - 1.*Enu)*(Elep - 1.*Enu)))) - 
           4.*C4*(Elep + Enu)*Fm1 + 
           2.*(aml*aml)*(Fm1*Fm1) - 
           2.*(Elep*Elep)*(Fm1*Fm1) - 2.*(Enu*Enu)*(Fm1*Fm1))
           - 1.*Ekaon*
           (4. + 
           (C4*C4)*
           (-4. + 3.*(aml*aml*aml*aml)*(C3*C3) + 
           2.*(aml*aml)*C3*
           (-3. + 4.*C3*((Elep - 1.*Enu)*(Elep - 1.*Enu))))
           + 24.*Elep*Fm1 - 24.*Enu*Fm1 + 
           8.*C4*(Elep + Enu)*Fm1 + 
           7.*(aml*aml)*(Fm1*Fm1) + 
           16.*(Elep*Elep)*(Fm1*Fm1) + 
           16.*(Enu*Enu)*(Fm1*Fm1))) + 
           2.*Ekaon*
           (2.*(akpk*akpk)*
           ((aml*aml)*(C3*C3)*(C4*C4) + 2.*(Fm1*Fm1)) + 
           2.*(apkk1*apkk1)*
           ((aml*aml)*(C3*C3)*(C4*C4) + 2.*(Fm1*Fm1)) + 
           apkk1*
           (2. + 4.*C4 + 
           (C4*C4)*
           (2. + (aml*aml*aml*aml)*(C3*C3) + 
           4.*(aml*aml)*(C3*C3)*Ekaon*
           (-1.*Elep + Enu)) + 
           3.*(aml*aml)*(Fm1*Fm1) - 
           8.*Ekaon*Elep*(Fm1*Fm1)) - 
           1.*akpk*
           (2. - 4.*C4 + 
           (C4*C4)*
           (2. + (aml*aml*aml*aml)*(C3*C3) + 
           4.*(aml*aml)*C3*
           (-1. + apkk1*C3 + 
           C3*Ekaon*(-1.*Elep + Enu))) + 
           3.*(aml*aml)*(Fm1*Fm1) + 
           8.*Ekaon*Enu*(Fm1*Fm1)) + 
           2.*Ekaon*
           ((aml*aml*aml*aml)*(C3*C3)*(C4*C4)*
           (Ekaon - 1.*Elep + Enu) - 
           2.*
           (((1.+C4)*(1.+C4))*Elep - 
           1.*((-1.+C4)*(-1.+C4))*Enu) + 
           amSig*
           (-2. + 
           (2. + 2.*(aml*aml)*C3 - 
           1.*(aml*aml*aml*aml)*(C3*C3))*(C4*C4) - 
           4.*Elep*Fm1 + 4.*Enu*Fm1 + 
           4.*C4*(Elep + Enu)*Fm1 - 
           3.*(aml*aml)*(Fm1*Fm1)) + 
           (aml*aml)*
           (-2.*C3*(C4*C4)*(Ekaon + 2.*Enu) + 
           (Ekaon - 3.*Elep + 3.*Enu)*(Fm1*Fm1)))
           ))*(t3*t3) + 
           C5*C6*t4*
           (24.*(akpk + apkk1 - 1.*Ekaon*(Elep + Enu))*
           Fm2*t1 + 
           24.*C2*
           ((akpk - 1.*apkk1)*Fm2 + 
           3.*Ekaon*(1. + Elep*Fm2 - 1.*Enu*Fm2))*
           t1 + 
           2.*(aml*aml)*(C5*C5*C5*C5)*C6*
           (4.*Ekaon*
           (akpk + 
           Ekaon*(amLam - 1.*Ekaon - 2.*Enu)) + 
           (amk*amk)*(-4.*amLam + 3.*Ekaon + 2.*Enu))*
           t4 + 
           2.*(C5*C5*C5)*C6*
           ((amk*amk)*
           (-2.*amLam + 2.*Ekaon + Elep - 1.*Enu)
           + 2.*Ekaon*
           (-1.*akpk + apkk1 + 
           2.*Ekaon*(amLam - 1.*Elep + Enu)))*t4
           + (aml*aml)*(C5*C5*C5*C5*C5)*C6*
           ((amk*amk)*
           (-2.*
           (3.*akpk - 3.*apkk1 + 
           2.*(amLam + 2.*Ekaon)*
           (Elep - 1.*Enu))*(Elep - 1.*Enu) + 
           (aml*aml)*
           (4.*amLam - 3.*Ekaon + Elep - 
           1.*Enu)) + 
           2.*Ekaon*
           (2.*(akpk*akpk) + 
           2.*apkk1*
           (apkk1 + 2.*Ekaon*(-1.*Elep + Enu))
           + (aml*aml)*
           (apkk1 + 
           2.*Ekaon*
           (-1.*amLam + Ekaon - 1.*Elep + Enu)
           ) - 
           1.*akpk*
           ((aml*aml) + 
           4.*
           (apkk1 - 1.*Ekaon*Elep + Ekaon*Enu)
           )))*t4 - 
           4.*(C5*C5)*
           (-3.*(aml*aml)*C3*(Ekaon - 1.*Elep + Enu)*
           (apkk1 + 2.*Ekaon*(-1.*Elep + Enu))*t2
           + C6*
           (3.*apkk1*(2.*Ekaon + (amk*amk)*Fm2) - 
           1.*(Elep + Enu)*
           (-4.*(Ekaon*Ekaon)*(-3. + amLam*Fm2) + 
           (amk*amk)*
           (-3. + 4.*amLam*Fm2 + 2.*Ekaon*Fm2)
           ))*t4 + 
           3.*akpk*
           ((aml*aml)*C3*(Ekaon - 1.*Elep + Enu)*
           t2 + C6*(2.*Ekaon + (amk*amk)*Fm2)*t4))
           + C5*
           (4.*(aml*aml)*C6*(Ekaon*Ekaon*Ekaon)*(Fm2*Fm2)*t4 - 
           4.*C6*(Ekaon*Ekaon)*
           (Enu*
           (-18. + 4.*akpk*(Fm2*Fm2) - 
           3.*(aml*aml)*(Fm2*Fm2)) + 
           3.*amLam*
           (6. + 4.*Elep*Fm2 - 4.*Enu*Fm2 + 
           (aml*aml)*(Fm2*Fm2)) + 
           Elep*
           (18. + 3.*(aml*aml)*(Fm2*Fm2) + 
           4.*apkk1*(Fm2*Fm2)))*t4 + 
           (amk*amk)*C6*
           (4.*amLam*
           (9. + 2.*(aml*aml)*(Fm2*Fm2) - 
           2.*(Elep*Elep)*(Fm2*Fm2) - 
           2.*(Enu*Enu)*(Fm2*Fm2)) + 
           3.*
           (-6.*Enu - 4.*akpk*Fm2 - 
           6.*(aml*aml)*Fm2 + 4.*apkk1*Fm2 + 
           4.*akpk*Enu*(Fm2*Fm2) - 
           1.*(aml*aml)*Enu*(Fm2*Fm2) + 
           Elep*
           (6. + (aml*aml)*(Fm2*Fm2) + 
           4.*apkk1*(Fm2*Fm2))))*t4 - 
           1.*Ekaon*
           (24.*t1 + 
           C6*
           (-8.*(akpk*akpk)*(Fm2*Fm2) + 
           6.*akpk*(6. + (aml*aml)*(Fm2*Fm2)) - 
           2.*apkk1*
           (18. + 3.*(aml*aml)*(Fm2*Fm2) + 
           4.*apkk1*(Fm2*Fm2)) + 
           (amk*amk)*
           (36. + 72.*Elep*Fm2 - 
           72.*Enu*Fm2 + 7.*(aml*aml)*(Fm2*Fm2) + 
           16.*(Elep*Elep)*(Fm2*Fm2) + 
           16.*(Enu*Enu)*(Fm2*Fm2)))*t4))) - 
           6.*C1*C4*t3*
           (-12.*akpk*Fm1*t1 - 12.*apkk1*Fm1*t1 + 
           12.*Ekaon*Elep*Fm1*t1 + 
           12.*Ekaon*Enu*Fm1*t1 - 
           12.*C2*
           (Ekaon + akpk*Fm1 - 1.*apkk1*Fm1 + 
           3.*Ekaon*Elep*Fm1 - 3.*Ekaon*Enu*Fm1)*t1
           - 6.*(amk*amk)*amLam*C5*C6*t4 - 
           6.*(amk*amk)*amSig*C5*C6*t4 + 
           12.*akpk*C5*C6*Ekaon*t4 + 
           12.*(amk*amk)*C5*C6*Ekaon*t4 - 
           12.*apkk1*C5*C6*Ekaon*t4 + 
           4.*akpk*(C5*C5)*C6*Ekaon*t4 + 
           4.*apkk1*(C5*C5)*C6*Ekaon*t4 + 
           12.*amLam*C5*C6*(Ekaon*Ekaon)*t4 + 
           12.*amSig*C5*C6*(Ekaon*Ekaon)*t4 - 
           6.*(amk*amk)*C5*C6*Elep*t4 + 
           2.*(amk*amk)*(C5*C5)*C6*Elep*t4 + 
           24.*C5*C6*(Ekaon*Ekaon)*Elep*t4 - 
           8.*(C5*C5)*C6*(Ekaon*Ekaon)*Elep*t4 + 
           6.*(amk*amk)*C5*C6*Enu*t4 + 
           2.*(amk*amk)*(C5*C5)*C6*Enu*t4 - 
           24.*C5*C6*(Ekaon*Ekaon)*Enu*t4 - 
           8.*(C5*C5)*C6*(Ekaon*Ekaon)*Enu*t4 + 
           6.*akpk*(amk*amk)*C5*C6*Fm1*t4 + 
           9.*(amk*amk)*(aml*aml)*C5*C6*Fm1*t4 - 
           6.*(amk*amk)*apkk1*C5*C6*Fm1*t4 + 
           6.*akpk*(amk*amk)*(C5*C5)*C6*Fm1*t4 + 
           6.*(amk*amk)*apkk1*(C5*C5)*C6*Fm1*t4 - 
           4.*(amk*amk)*amLam*(C5*C5)*C6*Elep*Fm1*t4 - 
           4.*(amk*amk)*amSig*(C5*C5)*C6*Elep*Fm1*t4 + 
           36.*(amk*amk)*C5*C6*Ekaon*Elep*Fm1*t4 - 
           4.*(amk*amk)*(C5*C5)*C6*Ekaon*Elep*Fm1*t4 + 
           12.*amLam*C5*C6*(Ekaon*Ekaon)*Elep*Fm1*t4 + 
           12.*amSig*C5*C6*(Ekaon*Ekaon)*Elep*Fm1*t4 + 
           4.*amLam*(C5*C5)*C6*(Ekaon*Ekaon)*Elep*Fm1*t4 + 
           4.*amSig*(C5*C5)*C6*(Ekaon*Ekaon)*Elep*Fm1*t4 - 
           4.*(amk*amk)*amLam*(C5*C5)*C6*Enu*Fm1*t4 - 
           4.*(amk*amk)*amSig*(C5*C5)*C6*Enu*Fm1*t4 - 
           36.*(amk*amk)*C5*C6*Ekaon*Enu*Fm1*t4 - 
           4.*(amk*amk)*(C5*C5)*C6*Ekaon*Enu*Fm1*t4 - 
           12.*amLam*C5*C6*(Ekaon*Ekaon)*Enu*Fm1*t4 - 
           12.*amSig*C5*C6*(Ekaon*Ekaon)*Enu*Fm1*t4 + 
           4.*amLam*(C5*C5)*C6*(Ekaon*Ekaon)*Enu*Fm1*t4 + 
           4.*amSig*(C5*C5)*C6*(Ekaon*Ekaon)*Enu*Fm1*t4 + 
           2.*akpk*(amk*amk)*C5*C6*Fm2*t4 + 
           3.*(amk*amk)*(aml*aml)*C5*C6*Fm2*t4 - 
           2.*(amk*amk)*apkk1*C5*C6*Fm2*t4 + 
           12.*(amk*amk)*C5*C6*Ekaon*Elep*Fm2*t4 + 
           4.*amLam*C5*C6*(Ekaon*Ekaon)*Elep*Fm2*t4 + 
           4.*amSig*C5*C6*(Ekaon*Ekaon)*Elep*Fm2*t4 - 
           12.*(amk*amk)*C5*C6*Ekaon*Enu*Fm2*t4 - 
           4.*amLam*C5*C6*(Ekaon*Ekaon)*Enu*Fm2*t4 - 
           4.*amSig*C5*C6*(Ekaon*Ekaon)*Enu*Fm2*t4 - 
           4.*(amk*amk)*(aml*aml)*amLam*C5*C6*Fm1*Fm2*t4 - 
           4.*(amk*amk)*(aml*aml)*amSig*C5*C6*Fm1*Fm2*t4 - 
           8.*(akpk*akpk)*C5*C6*Ekaon*Fm1*Fm2*t4 + 
           6.*akpk*(aml*aml)*C5*C6*Ekaon*Fm1*Fm2*t4 + 
           7.*(amk*amk)*(aml*aml)*C5*C6*Ekaon*Fm1*Fm2*t4 - 
           6.*(aml*aml)*apkk1*C5*C6*Ekaon*Fm1*Fm2*t4 - 
           8.*(apkk1*apkk1)*C5*C6*Ekaon*Fm1*Fm2*t4 + 
           6.*(aml*aml)*amLam*C5*C6*(Ekaon*Ekaon)*Fm1*Fm2*
           t4 + 
           6.*(aml*aml)*amSig*C5*C6*(Ekaon*Ekaon)*Fm1*Fm2*
           t4 - 
           4.*(aml*aml)*C5*C6*(Ekaon*Ekaon*Ekaon)*Fm1*Fm2*t4 - 
           3.*(amk*amk)*(aml*aml)*C5*C6*Elep*Fm1*Fm2*t4 - 
           12.*(amk*amk)*apkk1*C5*C6*Elep*Fm1*Fm2*t4 + 
           12.*(aml*aml)*C5*C6*(Ekaon*Ekaon)*Elep*Fm1*Fm2*
           t4 + 
           16.*apkk1*C5*C6*(Ekaon*Ekaon)*Elep*Fm1*Fm2*t4 + 
           4.*(amk*amk)*amLam*C5*C6*(Elep*Elep)*Fm1*Fm2*t4 + 
           4.*(amk*amk)*amSig*C5*C6*(Elep*Elep)*Fm1*Fm2*t4 + 
           16.*(amk*amk)*C5*C6*Ekaon*(Elep*Elep)*Fm1*Fm2*
           t4 - 
           12.*akpk*(amk*amk)*C5*C6*Enu*Fm1*Fm2*t4 + 
           3.*(amk*amk)*(aml*aml)*C5*C6*Enu*Fm1*Fm2*t4 + 
           16.*akpk*C5*C6*(Ekaon*Ekaon)*Enu*Fm1*Fm2*t4 - 
           12.*(aml*aml)*C5*C6*(Ekaon*Ekaon)*Enu*Fm1*Fm2*t4 + 
           4.*(amk*amk)*amLam*C5*C6*(Enu*Enu)*Fm1*Fm2*t4 + 
           4.*(amk*amk)*amSig*C5*C6*(Enu*Enu)*Fm1*Fm2*t4 + 
           16.*(amk*amk)*C5*C6*Ekaon*(Enu*Enu)*Fm1*Fm2*t4 + 
           C4*(-6.*(aml*aml)*apkk1*(C3*C3)*Elep*t2 + 
           6.*(aml*aml)*apkk1*(C3*C3)*Enu*t2 - 
           2.*(amk*amk)*amLam*(C5*C5)*C6*t4 - 
           2.*(amk*amk)*amSig*(C5*C5)*C6*t4 - 
           2.*(amk*amk)*(aml*aml)*amLam*C3*(C5*C5)*C6*t4 - 
           2.*(amk*amk)*(aml*aml)*amSig*C3*(C5*C5)*C6*t4 - 
           2.*(amk*amk)*(aml*aml)*amLam*(C5*C5*C5)*C6*t4 - 
           2.*(amk*amk)*(aml*aml)*amSig*(C5*C5*C5)*C6*t4 + 
           2.*(amk*amk)*(aml*aml*aml*aml)*amLam*C3*(C5*C5*C5)*C6*t4 + 
           2.*(amk*amk)*(aml*aml*aml*aml)*amSig*C3*(C5*C5*C5)*C6*t4 + 
           4.*(aml*aml)*(C5*C5)*
           (-1.*C5 + C3*(-1. + (aml*aml)*C5))*C6*
           (Ekaon*Ekaon*Ekaon)*t4 - 
           6.*(amk*amk)*C5*C6*Elep*t4 + 
           2.*(amk*amk)*(C5*C5)*C6*Elep*t4 + 
           (amk*amk)*(aml*aml*aml*aml)*C3*(C5*C5*C5)*C6*Elep*t4 + 
           6.*(amk*amk)*(aml*aml)*apkk1*C3*(C5*C5*C5)*C6*Elep*
           t4 - 
           2.*(amk*amk)*(aml*aml)*amLam*C3*(C5*C5*C5)*C6*
           (Elep*Elep)*t4 - 
           2.*(amk*amk)*(aml*aml)*amSig*C3*(C5*C5*C5)*C6*
           (Elep*Elep)*t4 - 6.*(amk*amk)*C5*C6*Enu*t4 - 
           2.*(amk*amk)*(C5*C5)*C6*Enu*t4 + 
           2.*(amk*amk)*(aml*aml)*C3*(C5*C5)*C6*Enu*t4 + 
           2.*(amk*amk)*(aml*aml)*(C5*C5*C5)*C6*Enu*t4 - 
           1.*(amk*amk)*(aml*aml*aml*aml)*C3*(C5*C5*C5)*C6*Enu*t4 - 
           6.*(amk*amk)*(aml*aml)*apkk1*C3*(C5*C5*C5)*C6*Enu*
           t4 + 
           4.*(amk*amk)*(aml*aml)*amLam*C3*(C5*C5*C5)*C6*Elep*
           Enu*t4 + 
           4.*(amk*amk)*(aml*aml)*amSig*C3*(C5*C5*C5)*C6*Elep*
           Enu*t4 - 
           2.*(amk*amk)*(aml*aml)*amLam*C3*(C5*C5*C5)*C6*
           (Enu*Enu)*t4 - 
           2.*(amk*amk)*(aml*aml)*amSig*C3*(C5*C5*C5)*C6*
           (Enu*Enu)*t4 - 
           6.*(amk*amk)*apkk1*C5*C6*Fm2*t4 + 
           4.*(amk*amk)*amLam*C5*C6*Elep*Fm2*t4 + 
           4.*(amk*amk)*amSig*C5*C6*Elep*Fm2*t4 + 
           4.*(amk*amk)*amLam*C5*C6*Enu*Fm2*t4 + 
           4.*(amk*amk)*amSig*C5*C6*Enu*Fm2*t4 + 
           6.*akpk*
           (-1.*(amk*amk)*C5*C6*Fm2*t4 + 
           (aml*aml)*C3*(Elep - 1.*Enu)*
           (C3*t2 - 1.*(amk*amk)*(C5*C5*C5)*C6*t4)) - 
           2.*(Ekaon*Ekaon)*
           ((aml*aml*aml*aml)*C3*(C5*C5*C5)*C6*
           (amLam + amSig + 2.*Elep - 2.*Enu)*
           t4 - 
           2.*C5*C6*
           (6.*Elep - 2.*C5*Elep + 6.*Enu + 
           2.*C5*Enu + 
           amLam*(C5 - 1.*(Elep + Enu)*Fm2) + 
           amSig*(C5 - 1.*(Elep + Enu)*Fm2))*
           t4 + 
           (aml*aml)*
           (6.*(C3*C3)*(Elep - 1.*Enu)*t2 - 
           1.*(C5*C5*C5)*C6*
           (amLam + amSig - 4.*Enu)*t4 - 
           1.*C3*(C5*C5)*C6*
           (amLam + amSig + 4.*akpk*C5*Elep - 
           4.*apkk1*C5*Elep - 4.*Enu - 
           4.*akpk*C5*Enu + 4.*apkk1*C5*Enu)*
           t4)) + 
           Ekaon*
           (-12.*t1 + 6.*(aml*aml)*apkk1*(C3*C3)*t2 + 
           12.*(aml*aml)*(C3*C3)*(Elep*Elep)*t2 - 
           24.*(aml*aml)*(C3*C3)*Elep*Enu*t2 + 
           12.*(aml*aml)*(C3*C3)*(Enu*Enu)*t2 - 
           12.*apkk1*C5*C6*t4 + 
           4.*(amk*amk)*(C5*C5)*C6*t4 + 
           4.*apkk1*(C5*C5)*C6*t4 + 
           3.*(amk*amk)*(aml*aml)*C3*(C5*C5)*C6*t4 + 
           3.*(amk*amk)*(aml*aml)*(C5*C5*C5)*C6*t4 + 
           4.*(akpk*akpk)*(aml*aml)*C3*(C5*C5*C5)*C6*t4 - 
           3.*(amk*amk)*(aml*aml*aml*aml)*C3*(C5*C5*C5)*C6*t4 + 
           2.*(aml*aml*aml*aml)*apkk1*C3*(C5*C5*C5)*C6*t4 + 
           4.*(aml*aml)*(apkk1*apkk1)*C3*(C5*C5*C5)*C6*t4 - 
           8.*(amk*amk)*(aml*aml)*C3*(C5*C5*C5)*C6*(Elep*Elep)*
           t4 + 
           16.*(amk*amk)*(aml*aml)*C3*(C5*C5*C5)*C6*Elep*
           Enu*t4 - 
           8.*(amk*amk)*(aml*aml)*C3*(C5*C5*C5)*C6*(Enu*Enu)*
           t4 + 4.*(amk*amk)*C5*C6*Elep*Fm2*t4 + 
           4.*(amk*amk)*C5*C6*Enu*Fm2*t4 - 
           2.*akpk*
           ((aml*aml*aml*aml)*C3*(C5*C5*C5)*C6*t4 + 
           2.*C5*(3. + C5)*C6*t4 + 
           (aml*aml)*
           (3.*(C3*C3)*t2 - 2.*(C5*C5*C5)*C6*t4 + 
           2.*C3*(C5*C5)*(-1. + 2.*apkk1*C5)*C6*
           t4)))))) - 
           144.*(amk*amk)*(C2*C2)*C7*t1*t5 + 
           72.*(aml*aml)*(C2*C2)*C7*t1*t5 + 
           72.*(amk*amk*amk*amk)*C1*C2*C4*C7*t3*t5 - 
           36.*(amk*amk)*(aml*aml)*C1*C2*C4*C7*t3*t5 + 
           72.*akpk*amSig*C1*C2*C4*C7*Ekaon*t3*t5 + 
           72.*(amk*amk)*amSig*C1*C2*C4*C7*Ekaon*t3*t5 - 
           36.*(aml*aml)*amSig*C1*C2*C4*C7*Ekaon*t3*t5 - 
           72.*amSig*apkk1*C1*C2*C4*C7*Ekaon*t3*t5 - 
           72.*akpk*amSig*C1*C2*(C4*C4)*C7*Ekaon*t3*t5 - 
           72.*amSig*apkk1*C1*C2*(C4*C4)*C7*Ekaon*t3*t5 + 
           72.*(amk*amk)*amSig*C1*C2*C4*C7*Elep*t3*t5 + 
           72.*(amk*amk)*amSig*C1*C2*(C4*C4)*C7*Elep*t3*t5 - 
           72.*(amk*amk)*amSig*C1*C2*C4*C7*Enu*t3*t5 + 
           72.*(amk*amk)*amSig*C1*C2*(C4*C4)*C7*Enu*t3*t5 - 
           36.*akpk*(aml*aml)*amSig*C1*C2*C4*C7*Fm1*t3*t5 - 
           72.*(amk*amk)*(aml*aml)*amSig*C1*C2*C4*C7*Fm1*t3*t5 - 
           288.*akpk*amSig*apkk1*C1*C2*C4*C7*Fm1*t3*t5 - 
           36.*(aml*aml)*amSig*apkk1*C1*C2*C4*C7*Fm1*t3*t5 + 
           72.*akpk*(amk*amk)*C1*C2*C4*C7*Ekaon*Fm1*t3*t5 - 
           72.*(amk*amk)*(aml*aml)*C1*C2*C4*C7*Ekaon*Fm1*t3*t5 - 
           72.*(amk*amk)*apkk1*C1*C2*C4*C7*Ekaon*Fm1*t3*t5 - 
           144.*akpk*(amk*amk)*C1*C2*C4*C7*Elep*Fm1*t3*t5 + 
           72.*(amk*amk*amk*amk)*C1*C2*C4*C7*Elep*Fm1*t3*t5 - 
           36.*(amk*amk)*(aml*aml)*C1*C2*C4*C7*Elep*Fm1*t3*t5 - 
           72.*(amk*amk*amk*amk)*C1*C2*C4*C7*Enu*Fm1*t3*t5 - 
           36.*(amk*amk)*(aml*aml)*C1*C2*C4*C7*Enu*Fm1*t3*t5 - 
           144.*(amk*amk)*apkk1*C1*C2*C4*C7*Enu*Fm1*t3*t5 + 
           72.*(amk*amk*amk*amk)*C2*C5*C6*C7*t4*t5 - 
           36.*(amk*amk)*(aml*aml)*C2*C5*C6*C7*t4*t5 + 
           72.*akpk*amLam*C2*C5*C6*C7*Ekaon*t4*t5 + 
           72.*(amk*amk)*amLam*C2*C5*C6*C7*Ekaon*t4*t5 - 
           36.*(aml*aml)*amLam*C2*C5*C6*C7*Ekaon*t4*t5 - 
           72.*amLam*apkk1*C2*C5*C6*C7*Ekaon*t4*t5 + 
           24.*akpk*amLam*C2*(C5*C5)*C6*C7*Ekaon*t4*t5 + 
           24.*amLam*apkk1*C2*(C5*C5)*C6*C7*Ekaon*t4*t5 + 
           72.*(amk*amk)*amLam*C2*C5*C6*C7*Elep*t4*t5 - 
           24.*(amk*amk)*amLam*C2*(C5*C5)*C6*C7*Elep*t4*t5 - 
           72.*(amk*amk)*amLam*C2*C5*C6*C7*Enu*t4*t5 - 
           24.*(amk*amk)*amLam*C2*(C5*C5)*C6*C7*Enu*t4*t5 - 
           12.*akpk*(aml*aml)*amLam*C2*C5*C6*C7*Fm2*t4*t5 - 
           24.*(amk*amk)*(aml*aml)*amLam*C2*C5*C6*C7*Fm2*t4*t5 - 
           96.*akpk*amLam*apkk1*C2*C5*C6*C7*Fm2*t4*t5 - 
           12.*(aml*aml)*amLam*apkk1*C2*C5*C6*C7*Fm2*t4*t5 + 
           24.*akpk*(amk*amk)*C2*C5*C6*C7*Ekaon*Fm2*t4*t5 - 
           24.*(amk*amk)*(aml*aml)*C2*C5*C6*C7*Ekaon*Fm2*t4*t5 - 
           24.*(amk*amk)*apkk1*C2*C5*C6*C7*Ekaon*Fm2*t4*t5 - 
           48.*akpk*(amk*amk)*C2*C5*C6*C7*Elep*Fm2*t4*t5 + 
           24.*(amk*amk*amk*amk)*C2*C5*C6*C7*Elep*Fm2*t4*t5 - 
           12.*(amk*amk)*(aml*aml)*C2*C5*C6*C7*Elep*Fm2*t4*t5 - 
           24.*(amk*amk*amk*amk)*C2*C5*C6*C7*Enu*Fm2*t4*t5 - 
           12.*(amk*amk)*(aml*aml)*C2*C5*C6*C7*Enu*Fm2*t4*t5 - 
           48.*(amk*amk)*apkk1*C2*C5*C6*C7*Enu*Fm2*t4*t5 + 
           144.*(amk*amk)*(C2*C2)*(C7*C7)*Ekaon*(t5*t5) - 
           36.*(aml*aml)*(C2*C2)*(C7*C7)*Ekaon*(t5*t5) + 
           144.*(amk*amk)*(C2*C2)*(C7*C7)*Elep*(t5*t5) - 
           36.*(aml*aml)*(C2*C2)*(C7*C7)*Elep*(t5*t5) - 
           144.*(amk*amk)*(C2*C2)*(C7*C7)*Enu*(t5*t5) + 
           36.*(aml*aml)*(C2*C2)*(C7*C7)*Enu*(t5*t5) + 
           2.*am*(36.*(-1. + (C2*C2))*(t1*t1) + 
           18.*akpk*(amk*amk)*(C1*C1)*(C4*C4)*(t3*t3) + 
           9.*(amk*amk*amk*amk)*(C1*C1)*(C4*C4)*(t3*t3) + 
           9.*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4)*(t3*t3) - 
           18.*(amk*amk)*apkk1*(C1*C1)*(C4*C4)*(t3*t3) - 
           36.*akpk*(amk*amk)*(C1*C1)*(C4*C4*C4)*(t3*t3) - 
           36.*(amk*amk)*apkk1*(C1*C1)*(C4*C4*C4)*(t3*t3) + 
           18.*akpk*(amk*amk)*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) - 
           9.*(amk*amk*amk*amk)*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) - 
           9.*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) - 
           18.*(amk*amk)*apkk1*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) + 
           18.*(amk*amk)*amSig*(C1*C1)*(C4*C4)*Ekaon*(t3*t3) - 
           18.*(amk*amk)*amSig*(C1*C1)*(C4*C4*C4*C4)*Ekaon*(t3*t3) + 
           18.*(amk*amk)*amSig*(C1*C1)*(C4*C4)*Elep*(t3*t3) + 
           36.*(amk*amk)*amSig*(C1*C1)*(C4*C4*C4)*Elep*(t3*t3) + 
           18.*(amk*amk)*amSig*(C1*C1)*(C4*C4*C4*C4)*Elep*(t3*t3) + 
           36.*(amk*amk)*(C1*C1)*(C4*C4)*Ekaon*Elep*(t3*t3) + 
           72.*(amk*amk)*(C1*C1)*(C4*C4*C4)*Ekaon*Elep*(t3*t3) + 
           36.*(amk*amk)*(C1*C1)*(C4*C4*C4*C4)*Ekaon*Elep*(t3*t3) - 
           18.*(amk*amk)*amSig*(C1*C1)*(C4*C4)*Enu*(t3*t3) + 
           36.*(amk*amk)*amSig*(C1*C1)*(C4*C4*C4)*Enu*(t3*t3) - 
           18.*(amk*amk)*amSig*(C1*C1)*(C4*C4*C4*C4)*Enu*(t3*t3) - 
           36.*(amk*amk)*(C1*C1)*(C4*C4)*Ekaon*Enu*(t3*t3) + 
           72.*(amk*amk)*(C1*C1)*(C4*C4*C4)*Ekaon*Enu*(t3*t3) - 
           36.*(amk*amk)*(C1*C1)*(C4*C4*C4*C4)*Ekaon*Enu*(t3*t3) + 
           72.*akpk*(amk*amk)*amSig*(C1*C1)*(C4*C4*C4)*Fm1*(t3*t3) + 
           72.*(amk*amk)*amSig*apkk1*(C1*C1)*(C4*C4*C4)*Fm1*
           (t3*t3) - 
           36.*akpk*(amk*amk)*(C1*C1)*(C4*C4)*Ekaon*Fm1*(t3*t3) - 
           36.*akpk*(amSig*amSig)*(C1*C1)*(C4*C4)*Ekaon*Fm1*
           (t3*t3) + 
           36.*(amk*amk)*apkk1*(C1*C1)*(C4*C4)*Ekaon*Fm1*
           (t3*t3) + 
           36.*(amSig*amSig)*apkk1*(C1*C1)*(C4*C4)*Ekaon*Fm1*
           (t3*t3) - 
           36.*akpk*(amk*amk)*(C1*C1)*(C4*C4*C4)*Ekaon*Fm1*(t3*t3) + 
           36.*akpk*(amSig*amSig)*(C1*C1)*(C4*C4*C4)*Ekaon*Fm1*
           (t3*t3) - 
           36.*(amk*amk)*apkk1*(C1*C1)*(C4*C4*C4)*Ekaon*Fm1*
           (t3*t3) + 
           36.*(amSig*amSig)*apkk1*(C1*C1)*(C4*C4*C4)*Ekaon*Fm1*
           (t3*t3) - 
           72.*akpk*amSig*(C1*C1)*(C4*C4)*(Ekaon*Ekaon)*Fm1*
           (t3*t3) + 
           72.*amSig*apkk1*(C1*C1)*(C4*C4)*(Ekaon*Ekaon)*Fm1*
           (t3*t3) - 
           72.*akpk*amSig*(C1*C1)*(C4*C4*C4)*(Ekaon*Ekaon)*Fm1*
           (t3*t3) - 
           72.*amSig*apkk1*(C1*C1)*(C4*C4*C4)*(Ekaon*Ekaon)*Fm1*
           (t3*t3) + 
           36.*(amk*amk*amk*amk)*(C1*C1)*(C4*C4)*Elep*Fm1*(t3*t3) - 
           36.*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4*C4)*Elep*Fm1*
           (t3*t3) + 
           36.*(amk*amk)*amSig*(C1*C1)*(C4*C4)*Ekaon*Elep*Fm1*
           (t3*t3) - 
           36.*(amk*amk)*amSig*(C1*C1)*(C4*C4*C4)*Ekaon*Elep*Fm1*
           (t3*t3) - 
           36.*(amk*amk*amk*amk)*(C1*C1)*(C4*C4)*Enu*Fm1*(t3*t3) - 
           36.*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4*C4)*Enu*Fm1*
           (t3*t3) - 
           36.*(amk*amk)*amSig*(C1*C1)*(C4*C4)*Ekaon*Enu*Fm1*
           (t3*t3) - 
           36.*(amk*amk)*amSig*(C1*C1)*(C4*C4*C4)*Ekaon*Enu*Fm1*
           (t3*t3) - 
           36.*(akpk*akpk)*(amk*amk)*(C1*C1)*(C4*C4)*(Fm1*Fm1)*(t3*t3) - 
           36.*(amk*amk)*(apkk1*apkk1)*(C1*C1)*(C4*C4)*(Fm1*Fm1)*
           (t3*t3) + 
           72.*(amk*amk)*amSig*apkk1*(C1*C1)*(C4*C4)*Elep*
           (Fm1*Fm1)*(t3*t3) + 
           72.*(amk*amk)*apkk1*(C1*C1)*(C4*C4)*Ekaon*Elep*
           (Fm1*Fm1)*(t3*t3) + 
           36.*(amSig*amSig)*apkk1*(C1*C1)*(C4*C4)*Ekaon*Elep*
           (Fm1*Fm1)*(t3*t3) + 
           18.*(amk*amk*amk*amk)*(C1*C1)*(C4*C4)*(Elep*Elep)*(Fm1*Fm1)*(t3*t3) - 
           18.*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4)*(Elep*Elep)*
           (Fm1*Fm1)*(t3*t3) + 
           72.*akpk*(amk*amk)*amSig*(C1*C1)*(C4*C4)*Enu*(Fm1*Fm1)*
           (t3*t3) + 
           72.*akpk*(amk*amk)*(C1*C1)*(C4*C4)*Ekaon*Enu*(Fm1*Fm1)*
           (t3*t3) + 
           36.*akpk*(amSig*amSig)*(C1*C1)*(C4*C4)*Ekaon*Enu*
           (Fm1*Fm1)*(t3*t3) + 
           18.*(amk*amk*amk*amk)*(C1*C1)*(C4*C4)*(Enu*Enu)*(Fm1*Fm1)*(t3*t3) - 
           18.*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4)*(Enu*Enu)*(Fm1*Fm1)*
           (t3*t3) + 36.*akpk*(amk*amk)*C1*C4*C5*C6*t3*t4 + 
           18.*(amk*amk*amk*amk)*C1*C4*C5*C6*t3*t4 + 
           18.*(amk*amk)*amLam*amSig*C1*C4*C5*C6*t3*t4 - 
           36.*(amk*amk)*apkk1*C1*C4*C5*C6*t3*t4 - 
           36.*akpk*(amk*amk)*C1*(C4*C4)*C5*C6*t3*t4 - 
           36.*(amk*amk)*apkk1*C1*(C4*C4)*C5*C6*t3*t4 + 
           12.*akpk*(amk*amk)*C1*C4*(C5*C5)*C6*t3*t4 + 
           12.*(amk*amk)*apkk1*C1*C4*(C5*C5)*C6*t3*t4 - 
           12.*akpk*(amk*amk)*C1*(C4*C4)*(C5*C5)*C6*t3*t4 + 
           6.*(amk*amk*amk*amk)*C1*(C4*C4)*(C5*C5)*C6*t3*t4 + 
           6.*(amk*amk)*amLam*amSig*C1*(C4*C4)*(C5*C5)*C6*t3*
           t4 + 12.*(amk*amk)*apkk1*C1*(C4*C4)*(C5*C5)*C6*t3*
           t4 + 18.*(amk*amk)*amLam*C1*C4*C5*C6*Ekaon*t3*
           t4 + 18.*(amk*amk)*amSig*C1*C4*C5*C6*Ekaon*t3*
           t4 + 6.*(amk*amk)*amLam*C1*(C4*C4)*(C5*C5)*C6*Ekaon*
           t3*t4 + 
           6.*(amk*amk)*amSig*C1*(C4*C4)*(C5*C5)*C6*Ekaon*t3*
           t4 + 18.*(amk*amk)*amLam*C1*C4*C5*C6*Elep*t3*
           t4 + 18.*(amk*amk)*amSig*C1*C4*C5*C6*Elep*t3*
           t4 + 18.*(amk*amk)*amLam*C1*(C4*C4)*C5*C6*Elep*t3*
           t4 + 18.*(amk*amk)*amSig*C1*(C4*C4)*C5*C6*Elep*t3*
           t4 - 6.*(amk*amk)*amLam*C1*C4*(C5*C5)*C6*Elep*t3*
           t4 - 6.*(amk*amk)*amSig*C1*C4*(C5*C5)*C6*Elep*t3*
           t4 - 6.*(amk*amk)*amLam*C1*(C4*C4)*(C5*C5)*C6*Elep*
           t3*t4 - 
           6.*(amk*amk)*amSig*C1*(C4*C4)*(C5*C5)*C6*Elep*t3*
           t4 + 72.*(amk*amk)*C1*C4*C5*C6*Ekaon*Elep*t3*
           t4 + 72.*(amk*amk)*C1*(C4*C4)*C5*C6*Ekaon*Elep*t3*
           t4 - 24.*(amk*amk)*C1*C4*(C5*C5)*C6*Ekaon*Elep*t3*
           t4 - 24.*(amk*amk)*C1*(C4*C4)*(C5*C5)*C6*Ekaon*Elep*
           t3*t4 - 
           18.*(amk*amk)*amLam*C1*C4*C5*C6*Enu*t3*t4 - 
           18.*(amk*amk)*amSig*C1*C4*C5*C6*Enu*t3*t4 + 
           18.*(amk*amk)*amLam*C1*(C4*C4)*C5*C6*Enu*t3*t4 + 
           18.*(amk*amk)*amSig*C1*(C4*C4)*C5*C6*Enu*t3*t4 - 
           6.*(amk*amk)*amLam*C1*C4*(C5*C5)*C6*Enu*t3*t4 - 
           6.*(amk*amk)*amSig*C1*C4*(C5*C5)*C6*Enu*t3*t4 + 
           6.*(amk*amk)*amLam*C1*(C4*C4)*(C5*C5)*C6*Enu*t3*t4 + 
           6.*(amk*amk)*amSig*C1*(C4*C4)*(C5*C5)*C6*Enu*t3*t4 - 
           72.*(amk*amk)*C1*C4*C5*C6*Ekaon*Enu*t3*t4 + 
           72.*(amk*amk)*C1*(C4*C4)*C5*C6*Ekaon*Enu*t3*t4 - 
           24.*(amk*amk)*C1*C4*(C5*C5)*C6*Ekaon*Enu*t3*t4 + 
           24.*(amk*amk)*C1*(C4*C4)*(C5*C5)*C6*Ekaon*Enu*t3*
           t4 - 12.*akpk*(amk*amk)*amLam*C1*C4*(C5*C5)*C6*
           Fm1*t3*t4 - 
           12.*akpk*(amk*amk)*amSig*C1*C4*(C5*C5)*C6*Fm1*t3*
           t4 - 12.*(amk*amk)*amLam*apkk1*C1*C4*(C5*C5)*C6*
           Fm1*t3*t4 - 
           12.*(amk*amk)*amSig*apkk1*C1*C4*(C5*C5)*C6*Fm1*t3*
           t4 - 36.*akpk*(amk*amk)*C1*C4*C5*C6*Ekaon*Fm1*
           t3*t4 - 
           36.*akpk*amLam*amSig*C1*C4*C5*C6*Ekaon*Fm1*t3*
           t4 + 36.*(amk*amk)*apkk1*C1*C4*C5*C6*Ekaon*Fm1*
           t3*t4 + 
           36.*amLam*amSig*apkk1*C1*C4*C5*C6*Ekaon*Fm1*
           t3*t4 + 
           12.*akpk*(amk*amk)*C1*C4*(C5*C5)*C6*Ekaon*Fm1*t3*
           t4 - 12.*akpk*amLam*amSig*C1*C4*(C5*C5)*C6*
           Ekaon*Fm1*t3*t4 + 
           12.*(amk*amk)*apkk1*C1*C4*(C5*C5)*C6*Ekaon*Fm1*t3*
           t4 - 12.*amLam*amSig*apkk1*C1*C4*(C5*C5)*C6*
           Ekaon*Fm1*t3*t4 - 
           36.*akpk*amLam*C1*C4*C5*C6*(Ekaon*Ekaon)*Fm1*t3*
           t4 - 36.*akpk*amSig*C1*C4*C5*C6*(Ekaon*Ekaon)*Fm1*
           t3*t4 + 
           36.*amLam*apkk1*C1*C4*C5*C6*(Ekaon*Ekaon)*Fm1*t3*
           t4 + 36.*amSig*apkk1*C1*C4*C5*C6*(Ekaon*Ekaon)*
           Fm1*t3*t4 + 
           12.*akpk*amLam*C1*C4*(C5*C5)*C6*(Ekaon*Ekaon)*Fm1*t3*
           t4 + 12.*akpk*amSig*C1*C4*(C5*C5)*C6*(Ekaon*Ekaon)*
           Fm1*t3*t4 + 
           12.*amLam*apkk1*C1*C4*(C5*C5)*C6*(Ekaon*Ekaon)*Fm1*
           t3*t4 + 
           12.*amSig*apkk1*C1*C4*(C5*C5)*C6*(Ekaon*Ekaon)*Fm1*
           t3*t4 + 
           36.*(amk*amk*amk*amk)*C1*C4*C5*C6*Elep*Fm1*t3*t4 + 
           12.*(amk*amk)*amLam*amSig*C1*C4*(C5*C5)*C6*Elep*
           Fm1*t3*t4 + 
           18.*(amk*amk)*amLam*C1*C4*C5*C6*Ekaon*Elep*Fm1*
           t3*t4 + 
           18.*(amk*amk)*amSig*C1*C4*C5*C6*Ekaon*Elep*Fm1*
           t3*t4 - 
           24.*akpk*amLam*C1*C4*(C5*C5)*C6*Ekaon*Elep*Fm1*
           t3*t4 + 
           6.*(amk*amk)*amLam*C1*C4*(C5*C5)*C6*Ekaon*Elep*Fm1*
           t3*t4 + 
           24.*akpk*amSig*C1*C4*(C5*C5)*C6*Ekaon*Elep*Fm1*
           t3*t4 + 
           6.*(amk*amk)*amSig*C1*C4*(C5*C5)*C6*Ekaon*Elep*Fm1*
           t3*t4 - 
           36.*(amk*amk*amk*amk)*C1*C4*C5*C6*Enu*Fm1*t3*t4 + 
           12.*(amk*amk)*amLam*amSig*C1*C4*(C5*C5)*C6*Enu*Fm1*
           t3*t4 - 
           18.*(amk*amk)*amLam*C1*C4*C5*C6*Ekaon*Enu*Fm1*t3*
           t4 - 18.*(amk*amk)*amSig*C1*C4*C5*C6*Ekaon*Enu*
           Fm1*t3*t4 + 
           6.*(amk*amk)*amLam*C1*C4*(C5*C5)*C6*Ekaon*Enu*Fm1*
           t3*t4 + 
           6.*(amk*amk)*amSig*C1*C4*(C5*C5)*C6*Ekaon*Enu*Fm1*
           t3*t4 + 
           24.*amLam*apkk1*C1*C4*(C5*C5)*C6*Ekaon*Enu*Fm1*
           t3*t4 - 
           24.*amSig*apkk1*C1*C4*(C5*C5)*C6*Ekaon*Enu*Fm1*
           t3*t4 + 
           12.*akpk*(amk*amk)*amLam*C1*(C4*C4)*C5*C6*Fm2*t3*
           t4 + 12.*akpk*(amk*amk)*amSig*C1*(C4*C4)*C5*C6*
           Fm2*t3*t4 + 
           12.*(amk*amk)*amLam*apkk1*C1*(C4*C4)*C5*C6*Fm2*t3*
           t4 + 12.*(amk*amk)*amSig*apkk1*C1*(C4*C4)*C5*C6*
           Fm2*t3*t4 - 
           12.*akpk*(amk*amk)*C1*C4*C5*C6*Ekaon*Fm2*t3*t4 - 
           12.*akpk*amLam*amSig*C1*C4*C5*C6*Ekaon*Fm2*t3*
           t4 + 12.*(amk*amk)*apkk1*C1*C4*C5*C6*Ekaon*Fm2*
           t3*t4 + 
           12.*amLam*amSig*apkk1*C1*C4*C5*C6*Ekaon*Fm2*
           t3*t4 - 
           12.*akpk*(amk*amk)*C1*(C4*C4)*C5*C6*Ekaon*Fm2*t3*
           t4 + 12.*akpk*amLam*amSig*C1*(C4*C4)*C5*C6*
           Ekaon*Fm2*t3*t4 - 
           12.*(amk*amk)*apkk1*C1*(C4*C4)*C5*C6*Ekaon*Fm2*t3*
           t4 + 12.*amLam*amSig*apkk1*C1*(C4*C4)*C5*C6*
           Ekaon*Fm2*t3*t4 - 
           12.*akpk*amLam*C1*C4*C5*C6*(Ekaon*Ekaon)*Fm2*t3*
           t4 - 12.*akpk*amSig*C1*C4*C5*C6*(Ekaon*Ekaon)*Fm2*
           t3*t4 + 
           12.*amLam*apkk1*C1*C4*C5*C6*(Ekaon*Ekaon)*Fm2*t3*
           t4 + 12.*amSig*apkk1*C1*C4*C5*C6*(Ekaon*Ekaon)*
           Fm2*t3*t4 - 
           12.*akpk*amLam*C1*(C4*C4)*C5*C6*(Ekaon*Ekaon)*Fm2*t3*
           t4 - 12.*akpk*amSig*C1*(C4*C4)*C5*C6*(Ekaon*Ekaon)*
           Fm2*t3*t4 - 
           12.*amLam*apkk1*C1*(C4*C4)*C5*C6*(Ekaon*Ekaon)*Fm2*
           t3*t4 - 
           12.*amSig*apkk1*C1*(C4*C4)*C5*C6*(Ekaon*Ekaon)*Fm2*
           t3*t4 + 
           12.*(amk*amk*amk*amk)*C1*C4*C5*C6*Elep*Fm2*t3*t4 - 
           12.*(amk*amk)*amLam*amSig*C1*(C4*C4)*C5*C6*Elep*
           Fm2*t3*t4 + 
           6.*(amk*amk)*amLam*C1*C4*C5*C6*Ekaon*Elep*Fm2*t3*
           t4 + 6.*(amk*amk)*amSig*C1*C4*C5*C6*Ekaon*Elep*
           Fm2*t3*t4 - 
           24.*akpk*amLam*C1*(C4*C4)*C5*C6*Ekaon*Elep*Fm2*
           t3*t4 - 
           6.*(amk*amk)*amLam*C1*(C4*C4)*C5*C6*Ekaon*Elep*Fm2*
           t3*t4 + 
           24.*akpk*amSig*C1*(C4*C4)*C5*C6*Ekaon*Elep*Fm2*
           t3*t4 - 
           6.*(amk*amk)*amSig*C1*(C4*C4)*C5*C6*Ekaon*Elep*Fm2*
           t3*t4 - 
           12.*(amk*amk*amk*amk)*C1*C4*C5*C6*Enu*Fm2*t3*t4 - 
           12.*(amk*amk)*amLam*amSig*C1*(C4*C4)*C5*C6*Enu*Fm2*
           t3*t4 - 
           6.*(amk*amk)*amLam*C1*C4*C5*C6*Ekaon*Enu*Fm2*t3*
           t4 - 6.*(amk*amk)*amSig*C1*C4*C5*C6*Ekaon*Enu*
           Fm2*t3*t4 - 
           6.*(amk*amk)*amLam*C1*(C4*C4)*C5*C6*Ekaon*Enu*Fm2*
           t3*t4 - 
           6.*(amk*amk)*amSig*C1*(C4*C4)*C5*C6*Ekaon*Enu*Fm2*
           t3*t4 + 
           24.*amLam*apkk1*C1*(C4*C4)*C5*C6*Ekaon*Enu*Fm2*
           t3*t4 - 
           24.*amSig*apkk1*C1*(C4*C4)*C5*C6*Ekaon*Enu*Fm2*
           t3*t4 - 
           24.*(akpk*akpk)*(amk*amk)*C1*C4*C5*C6*Fm1*Fm2*t3*
           t4 - 24.*(amk*amk)*(apkk1*apkk1)*C1*C4*C5*C6*Fm1*Fm2*
           t3*t4 + 
           24.*(amk*amk)*amLam*apkk1*C1*C4*C5*C6*Elep*Fm1*
           Fm2*t3*t4 + 
           24.*(amk*amk)*amSig*apkk1*C1*C4*C5*C6*Elep*Fm1*
           Fm2*t3*t4 + 
           48.*(amk*amk)*apkk1*C1*C4*C5*C6*Ekaon*Elep*Fm1*
           Fm2*t3*t4 + 
           24.*amLam*amSig*apkk1*C1*C4*C5*C6*Ekaon*Elep*
           Fm1*Fm2*t3*t4 + 
           12.*(amk*amk*amk*amk)*C1*C4*C5*C6*(Elep*Elep)*Fm1*Fm2*t3*
           t4 - 12.*(amk*amk)*amLam*amSig*C1*C4*C5*C6*
           (Elep*Elep)*Fm1*Fm2*t3*t4 + 
           24.*akpk*(amk*amk)*amLam*C1*C4*C5*C6*Enu*Fm1*Fm2*
           t3*t4 + 
           24.*akpk*(amk*amk)*amSig*C1*C4*C5*C6*Enu*Fm1*Fm2*
           t3*t4 + 
           48.*akpk*(amk*amk)*C1*C4*C5*C6*Ekaon*Enu*Fm1*Fm2*
           t3*t4 + 
           24.*akpk*amLam*amSig*C1*C4*C5*C6*Ekaon*Enu*
           Fm1*Fm2*t3*t4 + 
           12.*(amk*amk*amk*amk)*C1*C4*C5*C6*(Enu*Enu)*Fm1*Fm2*t3*t4 - 
           12.*(amk*amk)*amLam*amSig*C1*C4*C5*C6*(Enu*Enu)*Fm1*
           Fm2*t3*t4 + 
           18.*akpk*(amk*amk)*(C5*C5)*(C6*C6)*(t4*t4) + 
           9.*(amk*amk*amk*amk)*(C5*C5)*(C6*C6)*(t4*t4) + 
           9.*(amk*amk)*(amLam*amLam)*(C5*C5)*(C6*C6)*(t4*t4) - 
           18.*(amk*amk)*apkk1*(C5*C5)*(C6*C6)*(t4*t4) + 
           12.*akpk*(amk*amk)*(C5*C5*C5)*(C6*C6)*(t4*t4) + 
           12.*(amk*amk)*apkk1*(C5*C5*C5)*(C6*C6)*(t4*t4) + 
           2.*akpk*(amk*amk)*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) - 
           1.*(amk*amk*amk*amk)*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) - 
           1.*(amk*amk)*(amLam*amLam)*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) - 
           2.*(amk*amk)*apkk1*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) + 
           18.*(amk*amk)*amLam*(C5*C5)*(C6*C6)*Ekaon*(t4*t4) - 
           2.*(amk*amk)*amLam*(C5*C5*C5*C5)*(C6*C6)*Ekaon*(t4*t4) + 
           18.*(amk*amk)*amLam*(C5*C5)*(C6*C6)*Elep*(t4*t4) - 
           12.*(amk*amk)*amLam*(C5*C5*C5)*(C6*C6)*Elep*(t4*t4) + 
           2.*(amk*amk)*amLam*(C5*C5*C5*C5)*(C6*C6)*Elep*(t4*t4) + 
           36.*(amk*amk)*(C5*C5)*(C6*C6)*Ekaon*Elep*(t4*t4) - 
           24.*(amk*amk)*(C5*C5*C5)*(C6*C6)*Ekaon*Elep*(t4*t4) + 
           4.*(amk*amk)*(C5*C5*C5*C5)*(C6*C6)*Ekaon*Elep*(t4*t4) - 
           18.*(amk*amk)*amLam*(C5*C5)*(C6*C6)*Enu*(t4*t4) - 
           12.*(amk*amk)*amLam*(C5*C5*C5)*(C6*C6)*Enu*(t4*t4) - 
           2.*(amk*amk)*amLam*(C5*C5*C5*C5)*(C6*C6)*Enu*(t4*t4) - 
           36.*(amk*amk)*(C5*C5)*(C6*C6)*Ekaon*Enu*(t4*t4) - 
           24.*(amk*amk)*(C5*C5*C5)*(C6*C6)*Ekaon*Enu*(t4*t4) - 
           4.*(amk*amk)*(C5*C5*C5*C5)*(C6*C6)*Ekaon*Enu*(t4*t4) - 
           8.*akpk*(amk*amk)*amLam*(C5*C5*C5)*(C6*C6)*Fm2*(t4*t4) - 
           8.*(amk*amk)*amLam*apkk1*(C5*C5*C5)*(C6*C6)*Fm2*(t4*t4) - 
           12.*akpk*(amk*amk)*(C5*C5)*(C6*C6)*Ekaon*Fm2*(t4*t4) - 
           12.*akpk*(amLam*amLam)*(C5*C5)*(C6*C6)*Ekaon*Fm2*
           (t4*t4) + 
           12.*(amk*amk)*apkk1*(C5*C5)*(C6*C6)*Ekaon*Fm2*
           (t4*t4) + 
           12.*(amLam*amLam)*apkk1*(C5*C5)*(C6*C6)*Ekaon*Fm2*
           (t4*t4) + 
           4.*akpk*(amk*amk)*(C5*C5*C5)*(C6*C6)*Ekaon*Fm2*(t4*t4) - 
           4.*akpk*(amLam*amLam)*(C5*C5*C5)*(C6*C6)*Ekaon*Fm2*
           (t4*t4) + 
           4.*(amk*amk)*apkk1*(C5*C5*C5)*(C6*C6)*Ekaon*Fm2*(t4*t4) - 
           4.*(amLam*amLam)*apkk1*(C5*C5*C5)*(C6*C6)*Ekaon*Fm2*
           (t4*t4) - 
           24.*akpk*amLam*(C5*C5)*(C6*C6)*(Ekaon*Ekaon)*Fm2*
           (t4*t4) + 
           24.*amLam*apkk1*(C5*C5)*(C6*C6)*(Ekaon*Ekaon)*Fm2*
           (t4*t4) + 
           8.*akpk*amLam*(C5*C5*C5)*(C6*C6)*(Ekaon*Ekaon)*Fm2*
           (t4*t4) + 
           8.*amLam*apkk1*(C5*C5*C5)*(C6*C6)*(Ekaon*Ekaon)*Fm2*
           (t4*t4) + 
           12.*(amk*amk*amk*amk)*(C5*C5)*(C6*C6)*Elep*Fm2*(t4*t4) + 
           4.*(amk*amk)*(amLam*amLam)*(C5*C5*C5)*(C6*C6)*Elep*Fm2*
           (t4*t4) + 
           12.*(amk*amk)*amLam*(C5*C5)*(C6*C6)*Ekaon*Elep*Fm2*
           (t4*t4) + 
           4.*(amk*amk)*amLam*(C5*C5*C5)*(C6*C6)*Ekaon*Elep*Fm2*
           (t4*t4) - 
           12.*(amk*amk*amk*amk)*(C5*C5)*(C6*C6)*Enu*Fm2*(t4*t4) + 
           4.*(amk*amk)*(amLam*amLam)*(C5*C5*C5)*(C6*C6)*Enu*Fm2*
           (t4*t4) - 
           12.*(amk*amk)*amLam*(C5*C5)*(C6*C6)*Ekaon*Enu*Fm2*
           (t4*t4) + 
           4.*(amk*amk)*amLam*(C5*C5*C5)*(C6*C6)*Ekaon*Enu*Fm2*
           (t4*t4) - 
           4.*(akpk*akpk)*(amk*amk)*(C5*C5)*(C6*C6)*(Fm2*Fm2)*(t4*t4) - 
           4.*(amk*amk)*(apkk1*apkk1)*(C5*C5)*(C6*C6)*(Fm2*Fm2)*(t4*t4) + 
           8.*(amk*amk)*amLam*apkk1*(C5*C5)*(C6*C6)*Elep*(Fm2*Fm2)*
           (t4*t4) + 
           8.*(amk*amk)*apkk1*(C5*C5)*(C6*C6)*Ekaon*Elep*(Fm2*Fm2)*
           (t4*t4) + 
           4.*(amLam*amLam)*apkk1*(C5*C5)*(C6*C6)*Ekaon*Elep*
           (Fm2*Fm2)*(t4*t4) + 
           2.*(amk*amk*amk*amk)*(C5*C5)*(C6*C6)*(Elep*Elep)*(Fm2*Fm2)*(t4*t4) - 
           2.*(amk*amk)*(amLam*amLam)*(C5*C5)*(C6*C6)*(Elep*Elep)*(Fm2*Fm2)*
           (t4*t4) + 
           8.*akpk*(amk*amk)*amLam*(C5*C5)*(C6*C6)*Enu*(Fm2*Fm2)*
           (t4*t4) + 
           8.*akpk*(amk*amk)*(C5*C5)*(C6*C6)*Ekaon*Enu*(Fm2*Fm2)*
           (t4*t4) + 
           4.*akpk*(amLam*amLam)*(C5*C5)*(C6*C6)*Ekaon*Enu*(Fm2*Fm2)*
           (t4*t4) + 
           2.*(amk*amk*amk*amk)*(C5*C5)*(C6*C6)*(Enu*Enu)*(Fm2*Fm2)*(t4*t4) - 
           2.*(amk*amk)*(amLam*amLam)*(C5*C5)*(C6*C6)*(Enu*Enu)*(Fm2*Fm2)*
           (t4*t4) + 
           6.*t1*(-6.*(amk*amk)*C1*C2*C4*t3 - 
           6.*(amk*amk)*C1*(C4*C4)*t3 - 
           3.*(amk*amk)*(aml*aml)*C1*C3*(C4*C4)*t3 - 
           6.*amSig*C1*C2*C4*Ekaon*t3 - 
           6.*amSig*C1*(C4*C4)*Ekaon*t3 + 
           6.*(aml*aml)*C1*C3*(C4*C4)*(Ekaon*Ekaon)*t3 + 
           12.*C1*C4*Ekaon*Elep*t3 - 
           12.*C1*C2*C4*Ekaon*Elep*t3 + 
           12.*C1*(C4*C4)*Ekaon*Elep*t3 - 
           12.*C1*C2*(C4*C4)*Ekaon*Elep*t3 + 
           12.*C1*C4*Ekaon*Enu*t3 + 
           12.*C1*C2*C4*Ekaon*Enu*t3 - 
           12.*C1*(C4*C4)*Ekaon*Enu*t3 - 
           12.*C1*C2*(C4*C4)*Ekaon*Enu*t3 + 
           12.*(aml*aml)*C1*C3*(C4*C4)*Ekaon*Enu*t3 - 
           9.*(aml*aml)*C1*C2*C4*Ekaon*Fm1*t3 - 
           12.*(amk*amk)*C1*C2*C4*Elep*Fm1*t3 - 
           6.*amSig*C1*C4*Ekaon*Elep*Fm1*t3 - 
           6.*amSig*C1*C2*C4*Ekaon*Elep*Fm1*t3 + 
           12.*(amk*amk)*C1*C2*C4*Enu*Fm1*t3 - 
           6.*amSig*C1*C4*Ekaon*Enu*Fm1*t3 + 
           6.*amSig*C1*C2*C4*Ekaon*Enu*Fm1*t3 - 
           6.*(amk*amk)*C2*C5*C6*t4 + 
           2.*(amk*amk)*(C5*C5)*C6*t4 + 
           (amk*amk)*(aml*aml)*(C5*C5*C5)*C6*t4 - 
           6.*amLam*C2*C5*C6*Ekaon*t4 + 
           2.*amLam*(C5*C5)*C6*Ekaon*t4 - 
           2.*(aml*aml)*(C5*C5*C5)*C6*(Ekaon*Ekaon)*t4 + 
           12.*C5*C6*Ekaon*Elep*t4 - 
           12.*C2*C5*C6*Ekaon*Elep*t4 - 
           4.*(C5*C5)*C6*Ekaon*Elep*t4 + 
           4.*C2*(C5*C5)*C6*Ekaon*Elep*t4 + 
           12.*C5*C6*Ekaon*Enu*t4 + 
           12.*C2*C5*C6*Ekaon*Enu*t4 + 
           4.*(C5*C5)*C6*Ekaon*Enu*t4 + 
           4.*C2*(C5*C5)*C6*Ekaon*Enu*t4 - 
           4.*(aml*aml)*(C5*C5*C5)*C6*Ekaon*Enu*t4 - 
           3.*(aml*aml)*C2*C5*C6*Ekaon*Fm2*t4 - 
           4.*(amk*amk)*C2*C5*C6*Elep*Fm2*t4 - 
           2.*amLam*C5*C6*Ekaon*Elep*Fm2*t4 - 
           2.*amLam*C2*C5*C6*Ekaon*Elep*Fm2*t4 + 
           4.*(amk*amk)*C2*C5*C6*Enu*Fm2*t4 - 
           2.*amLam*C5*C6*Ekaon*Enu*Fm2*t4 + 
           2.*amLam*C2*C5*C6*Ekaon*Enu*Fm2*t4 + 
           2.*akpk*
           (3.*C1*C4*
           (-1. + C4 - 1.*(aml*aml)*C3*C4 + 
           amSig*Fm1 - 1.*Ekaon*Fm1 + 
           2.*Elep*Fm1 + 
           C2*(-1. + C4 + amSig*Fm1 + Ekaon*Fm1)
           )*t3 + 
           C5*C6*
           (-3. - 1.*C5 + (aml*aml)*(C5*C5) + 
           amLam*Fm2 - 1.*Ekaon*Fm2 + 
           2.*Elep*Fm2 + 
           C2*
           (-3. - 1.*C5 + amLam*Fm2 + 
           Ekaon*Fm2))*t4) + 
           2.*apkk1*
           (3.*C1*C4*
           (-1. - 1.*C4 + amSig*Fm1 - 
           1.*Ekaon*Fm1 - 2.*Enu*Fm1 + 
           C2*
           (1. + C4 - 1.*amSig*Fm1 - 
           1.*Ekaon*Fm1))*t3 - 
           1.*C5*C6*
           (3. - 1.*C5 - 1.*amLam*Fm2 + 
           Ekaon*Fm2 + 2.*Enu*Fm2 + 
           C2*(-3. + C5 + amLam*Fm2 + Ekaon*Fm2)
           )*t4)) + 
           (aml*aml*aml*aml)*(-9.*apkk1*C1*(C3*C3)*(C4*C4)*t2*t3 - 
           36.*C1*(C3*C3)*(C4*C4)*(Ekaon*Ekaon)*t2*t3 + 
           18.*C1*(C3*C3)*(C4*C4)*Ekaon*Elep*t2*t3 - 
           18.*C1*(C3*C3)*(C4*C4)*Ekaon*Enu*t2*t3 - 
           9.*(amSig*amSig)*(C1*C1)*(C3*C3)*(C4*C4*C4*C4)*(Ekaon*Ekaon)*
           (t3*t3) + 3.*apkk1*C3*(C5*C5*C5)*C6*t2*t4 + 
           12.*C3*(C5*C5*C5)*C6*(Ekaon*Ekaon)*t2*t4 - 
           6.*C3*(C5*C5*C5)*C6*Ekaon*Elep*t2*t4 + 
           6.*C3*(C5*C5*C5)*C6*Ekaon*Enu*t2*t4 + 
           6.*amLam*amSig*C1*C3*(C4*C4)*(C5*C5*C5)*C6*
           (Ekaon*Ekaon)*t3*t4 - 
           1.*(amLam*amLam)*(C5*C5*C5*C5*C5*C5)*(C6*C6)*(Ekaon*Ekaon)*(t4*t4) + 
           (amk*amk*amk*amk)*
           ((-3.*C1*C3*(C4*C4)*t3 + (C5*C5*C5)*C6*t4)*(-3.*C1*C3*(C4*C4)*t3 + (C5*C5*C5)*C6*t4)) + 
           akpk*(3.*C1*C3*(C4*C4)*t3 - 1.*(C5*C5*C5)*C6*t4)*
           (3.*C3*(t2 + (amk*amk)*C1*(C4*C4)*t3) - 
           1.*(amk*amk)*(C5*C5*C5)*C6*t4) + 
           (amk*amk)*
           (9.*(C1*C1)*(C3*C3)*(C4*C4*C4*C4)*
           ((amSig*amSig) - 1.*apkk1 - 
           1.*amSig*(Ekaon - 1.*Elep + Enu) - 
           2.*Ekaon*(Ekaon - 1.*Elep + Enu))*
           (t3*t3) + 
           (C5*C5*C5)*C6*t4*
           (-6.*C3*t2 - 
           1.*(C5*C5*C5)*C6*
           (-1.*(amLam*amLam) + apkk1 + 
           amLam*(Ekaon - 1.*Elep + Enu) + 
           2.*Ekaon*(Ekaon - 1.*Elep + Enu))*
           t4) + 
           3.*C1*C3*(C4*C4)*t3*
           (6.*C3*t2 + 
           (C5*C5*C5)*C6*
           (2.*apkk1 + 
           (amSig + 4.*Ekaon)*
           (Ekaon - 1.*Elep + Enu) + 
           amLam*
           (-2.*amSig + Ekaon - 1.*Elep + Enu)
           )*t4))) + 
           36.*akpk*C1*C2*C4*C7*Ekaon*t3*t5 - 
           36.*(amk*amk)*C1*C2*C4*C7*Ekaon*t3*t5 - 
           36.*apkk1*C1*C2*C4*C7*Ekaon*t3*t5 - 
           36.*akpk*C1*C2*(C4*C4)*C7*Ekaon*t3*t5 - 
           36.*apkk1*C1*C2*(C4*C4)*C7*Ekaon*t3*t5 + 
           36.*(amk*amk)*C1*C2*C4*C7*Elep*t3*t5 + 
           36.*(amk*amk)*C1*C2*(C4*C4)*C7*Elep*t3*t5 - 
           36.*(amk*amk)*C1*C2*C4*C7*Enu*t3*t5 + 
           36.*(amk*amk)*C1*C2*(C4*C4)*C7*Enu*t3*t5 - 
           144.*akpk*apkk1*C1*C2*C4*C7*Fm1*t3*t5 - 
           72.*akpk*C1*C2*C4*C7*(Ekaon*Ekaon)*Fm1*t3*t5 + 
           72.*apkk1*C1*C2*C4*C7*(Ekaon*Ekaon)*Fm1*t3*t5 + 
           144.*akpk*C1*C2*C4*C7*Ekaon*Elep*Fm1*t3*t5 - 
           72.*(amk*amk)*C1*C2*C4*C7*Ekaon*Elep*Fm1*t3*t5 + 
           72.*(amk*amk)*C1*C2*C4*C7*Ekaon*Enu*Fm1*t3*t5 + 
           144.*apkk1*C1*C2*C4*C7*Ekaon*Enu*Fm1*t3*t5 + 
           36.*akpk*C2*C5*C6*C7*Ekaon*t4*t5 - 
           36.*(amk*amk)*C2*C5*C6*C7*Ekaon*t4*t5 - 
           36.*apkk1*C2*C5*C6*C7*Ekaon*t4*t5 + 
           12.*akpk*C2*(C5*C5)*C6*C7*Ekaon*t4*t5 + 
           12.*apkk1*C2*(C5*C5)*C6*C7*Ekaon*t4*t5 + 
           36.*(amk*amk)*C2*C5*C6*C7*Elep*t4*t5 - 
           12.*(amk*amk)*C2*(C5*C5)*C6*C7*Elep*t4*t5 - 
           36.*(amk*amk)*C2*C5*C6*C7*Enu*t4*t5 - 
           12.*(amk*amk)*C2*(C5*C5)*C6*C7*Enu*t4*t5 - 
           48.*akpk*apkk1*C2*C5*C6*C7*Fm2*t4*t5 - 
           24.*akpk*C2*C5*C6*C7*(Ekaon*Ekaon)*Fm2*t4*t5 + 
           24.*apkk1*C2*C5*C6*C7*(Ekaon*Ekaon)*Fm2*t4*t5 + 
           48.*akpk*C2*C5*C6*C7*Ekaon*Elep*Fm2*t4*t5 - 
           24.*(amk*amk)*C2*C5*C6*C7*Ekaon*Elep*Fm2*t4*t5 + 
           24.*(amk*amk)*C2*C5*C6*C7*Ekaon*Enu*Fm2*t4*t5 + 
           48.*apkk1*C2*C5*C6*C7*Ekaon*Enu*Fm2*t4*t5 + 
           (aml*aml)*(-54.*(amk*amk)*amSig*(C1*C1)*(C4*C4)*Fm1*
           (t3*t3) + 
           54.*(amk*amk)*(C1*C1)*(C4*C4)*Ekaon*Fm1*(t3*t3) + 
           108.*amSig*(C1*C1)*(C4*C4)*(Ekaon*Ekaon)*Fm1*
           (t3*t3) + 
           27.*akpk*(amk*amk)*(C1*C1)*(C4*C4)*(Fm1*Fm1)*(t3*t3) + 
           18.*(amk*amk*amk*amk)*(C1*C1)*(C4*C4)*(Fm1*Fm1)*(t3*t3) + 
           18.*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4)*(Fm1*Fm1)*
           (t3*t3) - 
           27.*(amk*amk)*apkk1*(C1*C1)*(C4*C4)*(Fm1*Fm1)*
           (t3*t3) + 
           9.*(amk*amk)*amSig*(C1*C1)*(C4*C4)*Ekaon*(Fm1*Fm1)*
           (t3*t3) - 
           18.*(amk*amk)*(C1*C1)*(C4*C4)*(Ekaon*Ekaon)*(Fm1*Fm1)*
           (t3*t3) - 
           9.*(amSig*amSig)*(C1*C1)*(C4*C4)*(Ekaon*Ekaon)*(Fm1*Fm1)*
           (t3*t3) + 
           27.*(amk*amk)*amSig*(C1*C1)*(C4*C4)*Elep*(Fm1*Fm1)*
           (t3*t3) + 
           54.*(amk*amk)*(C1*C1)*(C4*C4)*Ekaon*Elep*(Fm1*Fm1)*
           (t3*t3) - 
           27.*(amk*amk)*amSig*(C1*C1)*(C4*C4)*Enu*(Fm1*Fm1)*
           (t3*t3) - 
           54.*(amk*amk)*(C1*C1)*(C4*C4)*Ekaon*Enu*(Fm1*Fm1)*
           (t3*t3) + 
           9.*(C3*C3)*
           ((Ekaon*Ekaon)*(t2*t2) + (Enu*Enu)*(t2*t2) - 
           2.*(akpk*akpk)*C1*(C4*C4)*t2*t3 - 
           1.*akpk*(amk*amk)*C1*(C4*C4)*t2*t3 + 
           4.*akpk*apkk1*C1*(C4*C4)*t2*t3 + 
           (amk*amk)*apkk1*C1*(C4*C4)*t2*t3 - 
           2.*(apkk1*apkk1)*C1*(C4*C4)*t2*t3 + 
           2.*akpk*amSig*C1*(C4*C4)*Enu*t2*t3 - 
           2.*amSig*apkk1*C1*(C4*C4)*Enu*t2*t3 + 
           2.*(amk*amk)*C1*(C4*C4)*(Enu*Enu)*t2*t3 - 
           2.*(akpk*akpk)*(amk*amk)*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) + 
           4.*akpk*(amk*amk)*apkk1*(C1*C1)*(C4*C4*C4*C4)*
           (t3*t3) - 
           2.*(amk*amk)*(apkk1*apkk1)*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) + 
           4.*akpk*(amk*amk)*amSig*(C1*C1)*(C4*C4*C4*C4)*Enu*
           (t3*t3) - 
           4.*(amk*amk)*amSig*apkk1*(C1*C1)*(C4*C4*C4*C4)*Enu*
           (t3*t3) + 
           (amk*amk*amk*amk)*(C1*C1)*(C4*C4*C4*C4)*(Enu*Enu)*(t3*t3) - 
           1.*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4*C4*C4)*(Enu*Enu)*
           (t3*t3) + 
           (Elep*Elep)*
           ((t2*t2) + 2.*(amk*amk)*C1*(C4*C4)*t2*t3 + 
           (amk*amk)*((amk*amk) - 1.*(amSig*amSig))*(C1*C1)*
           (C4*C4*C4*C4)*(t3*t3)) - 
           2.*Elep*
           (amSig*(akpk - 1.*apkk1)*C1*(C4*C4)*t3*
           (t2 + 2.*(amk*amk)*C1*(C4*C4)*t3) + 
           Enu*
           ((t2*t2) + 2.*(amk*amk)*C1*(C4*C4)*t2*t3 + 
           (amk*amk)*((amk*amk) - 1.*(amSig*amSig))*
           (C1*C1)*(C4*C4*C4*C4)*(t3*t3))) - 
           2.*Ekaon*
           (-1.*amSig*(akpk - 1.*apkk1)*C1*(C4*C4)*
           t2*t3 + 
           Elep*
           ((t2*t2) + 
           2.*(akpk + (amk*amk) - 1.*apkk1)*C1*
           (C4*C4)*t2*t3 + 
           (2.*(amk*amk) + (amSig*amSig))*
           (akpk - 1.*apkk1)*(C1*C1)*(C4*C4*C4*C4)*(t3*t3)
           ) - 
           1.*Enu*
           ((t2*t2) + 
           2.*(akpk + (amk*amk) - 1.*apkk1)*C1*
           (C4*C4)*t2*t3 + 
           (2.*(amk*amk) + (amSig*amSig))*
           (akpk - 1.*apkk1)*(C1*C1)*(C4*C4*C4*C4)*(t3*t3)
           ))) + 
           12.*akpk*(amk*amk)*C1*(C4*C4)*(C5*C5*C5)*C6*t3*t4 + 
           6.*(amk*amk*amk*amk)*C1*(C4*C4)*(C5*C5*C5)*C6*t3*t4 + 
           6.*(amk*amk)*amLam*amSig*C1*(C4*C4)*(C5*C5*C5)*C6*t3*
           t4 - 
           3.*(amk*amk)*amLam*C1*(C4*C4)*(C5*C5*C5)*C6*Ekaon*t3*
           t4 - 
           3.*(amk*amk)*amSig*C1*(C4*C4)*(C5*C5*C5)*C6*Ekaon*t3*
           t4 - 
           12.*(amk*amk)*C1*(C4*C4)*(C5*C5*C5)*C6*(Ekaon*Ekaon)*t3*
           t4 - 
           6.*amLam*amSig*C1*(C4*C4)*(C5*C5*C5)*C6*(Ekaon*Ekaon)*
           t3*t4 - 
           6.*(amk*amk)*amLam*C1*(C4*C4)*(C5*C5*C5)*C6*Enu*t3*
           t4 - 
           6.*(amk*amk)*amSig*C1*(C4*C4)*(C5*C5*C5)*C6*Enu*t3*
           t4 - 
           24.*(amk*amk)*C1*(C4*C4)*(C5*C5*C5)*C6*Ekaon*Enu*t3*
           t4 - 
           27.*(amk*amk)*amLam*C1*C4*C5*C6*Fm1*t3*t4 - 
           27.*(amk*amk)*amSig*C1*C4*C5*C6*Fm1*t3*t4 + 
           54.*(amk*amk)*C1*C4*C5*C6*Ekaon*Fm1*t3*t4 + 
           54.*amLam*C1*C4*C5*C6*(Ekaon*Ekaon)*Fm1*t3*t4 + 
           54.*amSig*C1*C4*C5*C6*(Ekaon*Ekaon)*Fm1*t3*t4 - 
           9.*(amk*amk)*amLam*C1*C4*C5*C6*Fm2*t3*t4 - 
           9.*(amk*amk)*amSig*C1*C4*C5*C6*Fm2*t3*t4 + 
           18.*(amk*amk)*C1*C4*C5*C6*Ekaon*Fm2*t3*t4 + 
           18.*amLam*C1*C4*C5*C6*(Ekaon*Ekaon)*Fm2*t3*t4 + 
           18.*amSig*C1*C4*C5*C6*(Ekaon*Ekaon)*Fm2*t3*t4 + 
           18.*akpk*(amk*amk)*C1*C4*C5*C6*Fm1*Fm2*t3*
           t4 + 
           12.*(amk*amk*amk*amk)*C1*C4*C5*C6*Fm1*Fm2*t3*t4 + 
           12.*(amk*amk)*amLam*amSig*C1*C4*C5*C6*Fm1*Fm2*
           t3*t4 - 
           18.*(amk*amk)*apkk1*C1*C4*C5*C6*Fm1*Fm2*t3*
           t4 + 
           3.*(amk*amk)*amLam*C1*C4*C5*C6*Ekaon*Fm1*Fm2*
           t3*t4 + 
           3.*(amk*amk)*amSig*C1*C4*C5*C6*Ekaon*Fm1*Fm2*
           t3*t4 - 
           12.*(amk*amk)*C1*C4*C5*C6*(Ekaon*Ekaon)*Fm1*Fm2*t3*
           t4 - 
           6.*amLam*amSig*C1*C4*C5*C6*(Ekaon*Ekaon)*Fm1*
           Fm2*t3*t4 + 
           9.*(amk*amk)*amLam*C1*C4*C5*C6*Elep*Fm1*Fm2*
           t3*t4 + 
           9.*(amk*amk)*amSig*C1*C4*C5*C6*Elep*Fm1*Fm2*
           t3*t4 + 
           36.*(amk*amk)*C1*C4*C5*C6*Ekaon*Elep*Fm1*Fm2*
           t3*t4 - 
           9.*(amk*amk)*amLam*C1*C4*C5*C6*Enu*Fm1*Fm2*t3*
           t4 - 
           9.*(amk*amk)*amSig*C1*C4*C5*C6*Enu*Fm1*Fm2*t3*
           t4 - 
           36.*(amk*amk)*C1*C4*C5*C6*Ekaon*Enu*Fm1*Fm2*
           t3*t4 - 
           4.*akpk*(amk*amk)*(C5*C5*C5*C5*C5)*(C6*C6)*(t4*t4) - 
           2.*(amk*amk*amk*amk)*(C5*C5*C5*C5*C5)*(C6*C6)*(t4*t4) - 
           2.*(amk*amk)*(amLam*amLam)*(C5*C5*C5*C5*C5)*(C6*C6)*(t4*t4) - 
           2.*(akpk*akpk)*(amk*amk)*(C5*C5*C5*C5*C5*C5)*(C6*C6)*(t4*t4) + 
           4.*akpk*(amk*amk)*apkk1*(C5*C5*C5*C5*C5*C5)*(C6*C6)*(t4*t4) - 
           2.*(amk*amk)*(apkk1*apkk1)*(C5*C5*C5*C5*C5*C5)*(C6*C6)*(t4*t4) + 
           2.*(amk*amk)*amLam*(C5*C5*C5*C5*C5)*(C6*C6)*Ekaon*(t4*t4) + 
           4.*(amk*amk)*(C5*C5*C5*C5*C5)*(C6*C6)*(Ekaon*Ekaon)*(t4*t4) + 
           2.*(amLam*amLam)*(C5*C5*C5*C5*C5)*(C6*C6)*(Ekaon*Ekaon)*(t4*t4) - 
           4.*akpk*(amk*amk)*amLam*(C5*C5*C5*C5*C5*C5)*(C6*C6)*Elep*
           (t4*t4) + 
           4.*(amk*amk)*amLam*apkk1*(C5*C5*C5*C5*C5*C5)*(C6*C6)*Elep*
           (t4*t4) - 
           4.*akpk*(amk*amk)*(C5*C5*C5*C5*C5*C5)*(C6*C6)*Ekaon*Elep*
           (t4*t4) - 
           2.*akpk*(amLam*amLam)*(C5*C5*C5*C5*C5*C5)*(C6*C6)*Ekaon*Elep*
           (t4*t4) + 
           4.*(amk*amk)*apkk1*(C5*C5*C5*C5*C5*C5)*(C6*C6)*Ekaon*Elep*
           (t4*t4) + 
           2.*(amLam*amLam)*apkk1*(C5*C5*C5*C5*C5*C5)*(C6*C6)*Ekaon*Elep*
           (t4*t4) + 
           (amk*amk*amk*amk)*(C5*C5*C5*C5*C5*C5)*(C6*C6)*(Elep*Elep)*(t4*t4) - 
           1.*(amk*amk)*(amLam*amLam)*(C5*C5*C5*C5*C5*C5)*(C6*C6)*(Elep*Elep)*
           (t4*t4) + 
           4.*(amk*amk)*amLam*(C5*C5*C5*C5*C5)*(C6*C6)*Enu*(t4*t4) + 
           4.*akpk*(amk*amk)*amLam*(C5*C5*C5*C5*C5*C5)*(C6*C6)*Enu*
           (t4*t4) - 
           4.*(amk*amk)*amLam*apkk1*(C5*C5*C5*C5*C5*C5)*(C6*C6)*Enu*
           (t4*t4) + 
           8.*(amk*amk)*(C5*C5*C5*C5*C5)*(C6*C6)*Ekaon*Enu*(t4*t4) + 
           4.*akpk*(amk*amk)*(C5*C5*C5*C5*C5*C5)*(C6*C6)*Ekaon*Enu*
           (t4*t4) + 
           2.*akpk*(amLam*amLam)*(C5*C5*C5*C5*C5*C5)*(C6*C6)*Ekaon*Enu*
           (t4*t4) - 
           4.*(amk*amk)*apkk1*(C5*C5*C5*C5*C5*C5)*(C6*C6)*Ekaon*Enu*
           (t4*t4) - 
           2.*(amLam*amLam)*apkk1*(C5*C5*C5*C5*C5*C5)*(C6*C6)*Ekaon*Enu*
           (t4*t4) - 
           2.*(amk*amk*amk*amk)*(C5*C5*C5*C5*C5*C5)*(C6*C6)*Elep*Enu*(t4*t4) + 
           2.*(amk*amk)*(amLam*amLam)*(C5*C5*C5*C5*C5*C5)*(C6*C6)*Elep*Enu*
           (t4*t4) + (amk*amk*amk*amk)*(C5*C5*C5*C5*C5*C5)*(C6*C6)*(Enu*Enu)*(t4*t4) - 
           1.*(amk*amk)*(amLam*amLam)*(C5*C5*C5*C5*C5*C5)*(C6*C6)*(Enu*Enu)*
           (t4*t4) - 
           18.*(amk*amk)*amLam*(C5*C5)*(C6*C6)*Fm2*(t4*t4) + 
           18.*(amk*amk)*(C5*C5)*(C6*C6)*Ekaon*Fm2*(t4*t4) + 
           36.*amLam*(C5*C5)*(C6*C6)*(Ekaon*Ekaon)*Fm2*(t4*t4) + 
           3.*akpk*(amk*amk)*(C5*C5)*(C6*C6)*(Fm2*Fm2)*(t4*t4) + 
           2.*(amk*amk*amk*amk)*(C5*C5)*(C6*C6)*(Fm2*Fm2)*(t4*t4) + 
           2.*(amk*amk)*(amLam*amLam)*(C5*C5)*(C6*C6)*(Fm2*Fm2)*
           (t4*t4) - 
           3.*(amk*amk)*apkk1*(C5*C5)*(C6*C6)*(Fm2*Fm2)*(t4*t4) + 
           (amk*amk)*amLam*(C5*C5)*(C6*C6)*Ekaon*(Fm2*Fm2)*
           (t4*t4) - 
           2.*(amk*amk)*(C5*C5)*(C6*C6)*(Ekaon*Ekaon)*(Fm2*Fm2)*
           (t4*t4) - 
           1.*(amLam*amLam)*(C5*C5)*(C6*C6)*(Ekaon*Ekaon)*(Fm2*Fm2)*
           (t4*t4) + 
           3.*(amk*amk)*amLam*(C5*C5)*(C6*C6)*Elep*(Fm2*Fm2)*
           (t4*t4) + 
           6.*(amk*amk)*(C5*C5)*(C6*C6)*Ekaon*Elep*(Fm2*Fm2)*
           (t4*t4) - 
           3.*(amk*amk)*amLam*(C5*C5)*(C6*C6)*Enu*(Fm2*Fm2)*
           (t4*t4) - 
           6.*(amk*amk)*(C5*C5)*(C6*C6)*Ekaon*Enu*(Fm2*Fm2)*
           (t4*t4) + 
           3.*C3*
           (2.*(akpk*akpk)*(C5*C5*C5)*C6*
           (t2 + 2.*(amk*amk)*C1*(C4*C4)*t3)*t4 - 
           2.*(amk*amk*amk*amk)*C1*(C4*C4)*t3*
           (3.*C1*(C4*C4)*t3 + 
           (C5*C5)*C6*
           (-1. + C5*((Elep - 1.*Enu)*(Elep - 1.*Enu)))*t4) + 
           2.*
           (3.*(amSig*amSig)*(C1*C1)*(C4*C4*C4*C4)*(Ekaon*Ekaon)*
           (t3*t3) + 
           (C5*C5)*C6*
           ((apkk1*apkk1)*C5 - 
           2.*Ekaon*(Ekaon + Enu) + 
           apkk1*C5*
           (2.*Ekaon*(-1.*Elep + Enu) + 
           amLam*(Ekaon - 1.*Elep + Enu)))*t2*
           t4 + 
           C1*(C4*C4)*Ekaon*t3*
           (6.*Ekaon*t2 + 6.*Enu*t2 - 
           1.*amLam*amSig*(C5*C5)*C6*Ekaon*t4 - 
           2.*amLam*amSig*apkk1*(C5*C5*C5)*C6*Elep*
           t4 + 
           2.*amLam*amSig*apkk1*(C5*C5*C5)*C6*Enu*
           t4)) - 
           1.*(amk*amk)*
           (6.*(C1*C1)*(C4*C4*C4*C4)*
           ((amSig*amSig) - 
           1.*amSig*(Ekaon + 2.*Enu) - 
           2.*Ekaon*(Ekaon + 2.*Enu))*(t3*t3) + 
           (C5*C5)*C6*
           (-2. + apkk1*C5 + 
           2.*C5*(Elep - 1.*Enu)*
           (-2.*Ekaon + Elep - 1.*Enu))*t2*t4
           + C1*(C4*C4)*t3*
           (6.*t2 + 
           (C5*C5)*C6*
           (-4.*(apkk1*apkk1)*C5 + 
           4.*apkk1*C5*(amSig + 2.*Ekaon)*
           (Elep - 1.*Enu) + 
           (amSig + 4.*Ekaon)*
           (Ekaon + 2.*Enu) + 
           amLam*
           (Ekaon + 4.*apkk1*C5*Elep - 
           2.*amSig*
           (1. + C5*((Elep - 1.*Enu)*(Elep - 1.*Enu))) + 
           2.*Enu - 4.*apkk1*C5*Enu))*t4)) + 
           akpk*
           (-12.*(amk*amk)*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) + 
           (C5*C5)*C6*
           (2. + (amk*amk)*C5 - 4.*apkk1*C5 - 
           2.*amLam*C5*Ekaon + 
           2.*amLam*C5*Elep + 
           4.*C5*Ekaon*Elep - 
           2.*amLam*C5*Enu - 4.*C5*Ekaon*Enu)*
           t2*t4 - 
           2.*C1*(C4*C4)*t3*
           (3.*t2 - 
           2.*(C5*C5)*C6*
           (amLam*amSig*C5*Ekaon*
           (Elep - 1.*Enu) + 
           (amk*amk)*
           (1. - 2.*apkk1*C5 + 
           amSig*C5*Elep + 2.*C5*Ekaon*Elep + 
           amLam*C5*(Elep - 1.*Enu) - 
           1.*amSig*C5*Enu - 2.*C5*Ekaon*Enu))
           *t4))) + 
           18.*C1*C2*C4*C7*Ekaon*t3*t5 - 
           18.*akpk*C1*C2*C4*C7*Fm1*t3*t5 - 
           36.*(amk*amk)*C1*C2*C4*C7*Fm1*t3*t5 - 
           18.*apkk1*C1*C2*C4*C7*Fm1*t3*t5 + 
           72.*C1*C2*C4*C7*(Ekaon*Ekaon)*Fm1*t3*t5 + 
           36.*C1*C2*C4*C7*Ekaon*Elep*Fm1*t3*t5 + 
           36.*C1*C2*C4*C7*Ekaon*Enu*Fm1*t3*t5 + 
           18.*C2*C5*C6*C7*Ekaon*t4*t5 - 
           6.*akpk*C2*C5*C6*C7*Fm2*t4*t5 - 
           12.*(amk*amk)*C2*C5*C6*C7*Fm2*t4*t5 - 
           6.*apkk1*C2*C5*C6*C7*Fm2*t4*t5 + 
           24.*C2*C5*C6*C7*(Ekaon*Ekaon)*Fm2*t4*t5 + 
           12.*C2*C5*C6*C7*Ekaon*Elep*Fm2*t4*t5 + 
           12.*C2*C5*C6*C7*Ekaon*Enu*Fm2*t4*t5))) - 
           2.*(akk1*akk1)*((aml*aml)*
           (9.*(C3*C3)*(Enu*(t2*t2) + 
           2.*akpk*am*C1*(C4*C4)*t2*t3 + 
           4.*am*(amk*amk)*C1*(C4*C4)*t2*t3 + 
           2.*akpk*amSig*C1*(C4*C4)*t2*t3 + 
           4.*(amk*amk)*amSig*C1*(C4*C4)*t2*t3 - 
           2.*am*apkk1*C1*(C4*C4)*t2*t3 - 
           2.*amSig*apkk1*C1*(C4*C4)*t2*t3 + 
           2.*(amk*amk)*C1*(C4*C4)*Enu*t2*t3 + 
           2.*akpk*am*(amk*amk)*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) + 
           2.*(am*am*am)*(amk*amk)*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) + 
           2.*am*(amk*amk*amk*amk)*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) + 
           2.*akpk*(amk*amk)*amSig*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) + 
           4.*(am*am)*(amk*amk)*amSig*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) + 
           2.*(amk*amk*amk*amk)*amSig*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) + 
           2.*am*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) - 
           2.*am*(amk*amk)*apkk1*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) - 
           2.*(amk*amk)*amSig*apkk1*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) + 
           4.*(am*am)*(C1*C1)*(C4*C4*C4*C4)*(Ekaon*Ekaon*Ekaon)*(t3*t3) - 
           1.*(am*am)*(amk*amk)*(C1*C1)*(C4*C4*C4*C4)*Enu*(t3*t3) + 
           (amk*amk*amk*amk)*(C1*C1)*(C4*C4*C4*C4)*Enu*(t3*t3) - 
           2.*am*(amk*amk)*amSig*(C1*C1)*(C4*C4*C4*C4)*Enu*(t3*t3) - 
           1.*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4*C4*C4)*Enu*(t3*t3) - 
           2.*am*C1*(C4*C4)*(Ekaon*Ekaon)*t3*
           (4.*t2 + 
           C1*(C4*C4)*
           ((am*am) + 2.*(amk*amk) + (amSig*amSig) + 
           2.*am*(amSig + Elep - 1.*Enu))*t3) + 
           Elep*(-1.*(t2*t2) - 
           2.*(amk*amk)*C1*(C4*C4)*t2*t3 + 
           (amk*amk)*
           ((am*am) - 1.*(amk*amk) + 2.*am*amSig + 
           (amSig*amSig))*(C1*C1)*(C4*C4*C4*C4)*(t3*t3)) + 
           Ekaon*
           (3.*(t2*t2) + 
           4.*C1*(C4*C4)*
           ((amk*amk) + am*(Elep - 1.*Enu))*t2*t3 + 
           (C1*C1)*(C4*C4*C4*C4)*
           ((amk*amk*amk*amk) + (amk*amk)*(amSig*amSig) - 
           2.*akpk*((am*am) - 1.*(amSig*amSig)) - 
           2.*(amSig*amSig)*apkk1 + 
           (am*am)*(-3.*(amk*amk) + 2.*apkk1) - 
           2.*am*(amk*amk)*
           (amSig - 2.*Elep + 2.*Enu))*(t3*t3)))
           - 6.*C3*(C5*C5*C5)*C6*
           (2.*(amk*amk)*amLam*t2 - 1.*amLam*apkk1*t2 + 
           2.*(amk*amk)*Ekaon*t2 - 1.*(amk*amk)*Elep*t2 + 
           (amk*amk)*Enu*t2 + (amk*amk*amk*amk)*amLam*C1*(C4*C4)*t3 + 
           (amk*amk*amk*amk)*amSig*C1*(C4*C4)*t3 - 
           1.*(amk*amk)*amLam*apkk1*C1*(C4*C4)*t3 - 
           1.*(amk*amk)*amSig*apkk1*C1*(C4*C4)*t3 + 
           (amk*amk*amk*amk)*C1*(C4*C4)*Ekaon*t3 + 
           (amk*amk)*amLam*amSig*C1*(C4*C4)*Ekaon*t3 - 
           2.*amLam*amSig*apkk1*C1*(C4*C4)*Ekaon*t3 + 
           2.*(am*am*am)*C1*(C4*C4)*((amk*amk) - 1.*(Ekaon*Ekaon))*
           t3 - 1.*(amk*amk*amk*amk)*C1*(C4*C4)*Elep*t3 + 
           (amk*amk)*amLam*amSig*C1*(C4*C4)*Elep*t3 + 
           (amk*amk*amk*amk)*C1*(C4*C4)*Enu*t3 - 
           1.*(amk*amk)*amLam*amSig*C1*(C4*C4)*Enu*t3 + 
           (am*am)*C1*(C4*C4)*
           ((amk*amk)*
           (2.*amLam + 2.*amSig - 3.*Ekaon + 
           Elep - 1.*Enu) + 
           2.*Ekaon*
           (apkk1 + 
           Ekaon*
           (-1.*amLam - 1.*amSig + 2.*Ekaon - 
           2.*Elep + 2.*Enu)))*t3 + 
           akpk*((amk*amk)*amSig*C1*(C4*C4)*t3 - 
           2.*(am*am)*C1*(C4*C4)*Ekaon*t3 + 
           am*(t2 + 2.*(amk*amk)*C1*(C4*C4)*t3) + 
           amLam*
           (t2 + 
           C1*(C4*C4)*((amk*amk) + 2.*amSig*Ekaon)*t3
           )) + 
           am*(-1.*apkk1*t2 - 4.*(Ekaon*Ekaon)*t2 + 
           2.*Ekaon*Elep*t2 - 2.*Ekaon*Enu*t2 + 
           2.*(amk*amk*amk*amk)*C1*(C4*C4)*t3 - 
           2.*amLam*amSig*C1*(C4*C4)*(Ekaon*Ekaon)*t3 + 
           (amk*amk)*
           (2.*t2 - 
           1.*C1*(C4*C4)*
           (2.*apkk1 + 
           (amSig + 4.*Ekaon)*
           (Ekaon - 1.*Elep + Enu) + 
           amLam*
           (-2.*amSig + Ekaon - 1.*Elep + Enu)
           )*t3)))*t4 + 
           (C5*C5*C5*C5*C5*C5)*(C6*C6)*
           (2.*(amk*amk*amk*amk)*amLam - 2.*(amk*amk)*amLam*apkk1 + 
           (amk*amk*amk*amk)*Ekaon + (amk*amk)*(amLam*amLam)*Ekaon - 
           2.*(amLam*amLam)*apkk1*Ekaon - 
           2.*akpk*(am + amLam)*
           (-1.*(amk*amk) + am*Ekaon - 1.*amLam*Ekaon)
           + 2.*(am*am*am)*((amk*amk) - 1.*(Ekaon*Ekaon)) - 
           1.*(amk*amk*amk*amk)*Elep + (amk*amk)*(amLam*amLam)*Elep + 
           (amk*amk*amk*amk)*Enu - 1.*(amk*amk)*(amLam*amLam)*Enu + 
           2.*am*
           ((amk*amk*amk*amk) - 1.*(amLam*amLam)*(Ekaon*Ekaon) + 
           (amk*amk)*
           ((amLam*amLam) - 1.*apkk1 - 
           1.*amLam*(Ekaon - 1.*Elep + Enu) - 
           2.*Ekaon*(Ekaon - 1.*Elep + Enu))) + 
           (am*am)*
           ((amk*amk)*
           (4.*amLam - 3.*Ekaon + Elep - 1.*Enu)
           + 2.*Ekaon*
           (apkk1 + 
           2.*Ekaon*
           (-1.*amLam + Ekaon - 1.*Elep + Enu))
           ))*(t4*t4)) - 
           2.*(-1.*(amk*amk*amk*amk)*(3.*C1*C4*Fm1*t3 + C5*C6*Fm2*t4)*
           (3.*C1*C4*
           (-2. + am*Fm1 - 1.*Elep*Fm1 + Enu*Fm1)*t3
           + C5*C6*
           (-6. + am*Fm2 - 1.*Elep*Fm2 + Enu*Fm2)*t4)
           + (amk*amk)*
           (-9.*(C1*C1)*(C4*C4)*Fm1*
           ((am*am*am)*Fm1 + 
           (am*am)*
           (-2. + 2.*amSig*Fm1 - 2.*Ekaon*Fm1 + 
           Elep*Fm1 - 1.*Enu*Fm1) + 
           amSig*
           (-4.*Ekaon + 
           2.*(akpk - 1.*apkk1)*Fm1 + 
           amSig*(-2. + Elep*Fm1 - 1.*Enu*Fm1))
           + am*
           ((amSig*amSig)*Fm1 + 
           2.*(akpk - 1.*apkk1)*Fm1 + 
           4.*Ekaon*
           (1. + Elep*Fm1 - 1.*Enu*Fm1) + 
           2.*amSig*
           (-2. + Ekaon*Fm1 + Elep*Fm1 - 
           1.*Enu*Fm1)))*(t3*t3) + 
           C5*C6*Fm2*t4*
           (-1.*C5*C6*
           ((am*am*am)*Fm2 + 
           (am*am)*
           (-6. + 2.*amLam*Fm2 - 
           2.*Ekaon*Fm2 + Elep*Fm2 - 
           1.*Enu*Fm2) + 
           amLam*
           (2.*
           (-6.*Ekaon + akpk*Fm2 - 
           1.*apkk1*Fm2) + 
           amLam*(-6. + Elep*Fm2 - 1.*Enu*Fm2)
           ) + 
           am*
           ((amLam*amLam)*Fm2 + 
           2.*amLam*
           (-6. + Ekaon*Fm2 + Elep*Fm2 - 
           1.*Enu*Fm2) + 
           2.*
           ((akpk - 1.*apkk1)*Fm2 + 
           2.*Ekaon*
           (3. + Elep*Fm2 - 1.*Enu*Fm2))))*t4
           - 12.*C2*
           (t1 - 1.*C7*(am + amLam + Ekaon)*t5))
           - 6.*C1*C4*t3*
           (C5*C6*
           ((am*am*am)*Fm1*Fm2 - 
           1.*amSig*
           ((-1.*akpk + apkk1)*Fm1*Fm2 + 
           Ekaon*(3.*Fm1 + Fm2)) + 
           (am*am)*
           (-1.*Fm2 + 
           Fm1*
           (-3. + amLam*Fm2 + amSig*Fm2 - 
           2.*Ekaon*Fm2 + Elep*Fm2 - 
           1.*Enu*Fm2)) - 
           1.*amLam*
           ((-1.*akpk + apkk1)*Fm1*Fm2 + 
           Ekaon*(3.*Fm1 + Fm2) + 
           amSig*
           (Fm2 + 
           Fm1*(3. - 1.*Elep*Fm2 + Enu*Fm2)))
           + am*
           (amSig*
           (-1.*Fm2 + 
           Fm1*
           (-3. + Ekaon*Fm2 + Elep*Fm2 - 
           1.*Enu*Fm2)) + 
           amLam*
           (-1.*Fm2 + 
           Fm1*
           (-3. + amSig*Fm2 + Ekaon*Fm2 + 
           Elep*Fm2 - 1.*Enu*Fm2)) + 
           2.*
           ((akpk - 1.*apkk1)*Fm1*Fm2 + 
           Ekaon*
           (Fm2 + 
           Fm1*(3. + 2.*Elep*Fm2 - 2.*Enu*Fm2)
           ))))*t4 + 
           6.*C2*Fm1*
           (t1 - 1.*C7*(am + amSig + Ekaon)*t5)))
           + 2.*Ekaon*
           (-1.*(3.*amSig*C1*C4*Fm1*t3 + 
           amLam*C5*C6*Fm2*t4)*
           (6.*C2*t1 + 
           (akpk - 1.*apkk1)*
           (3.*amSig*C1*C4*Fm1*t3 + 
           amLam*C5*C6*Fm2*t4)) + 
           (am*am)*(3.*C1*C4*Fm1*t3 + C5*C6*Fm2*t4)*
           (akpk*(3.*C1*C4*Fm1*t3 + C5*C6*Fm2*t4) - 
           1.*apkk1*
           (3.*C1*C4*Fm1*t3 + C5*C6*Fm2*t4) + 
           2.*Ekaon*
           (3.*amSig*C1*C4*Fm1*t3 + 
           3.*C1*C4*(Elep - 1.*Enu)*Fm1*t3 + 
           C5*C6*(amLam + Elep - 1.*Enu)*Fm2*t4)
           ) - 
           6.*am*
           (6.*amSig*(C1*C1)*(C4*C4)*Ekaon*Fm1*(t3*t3) + 
           C5*C6*Fm2*t4*
           (-1.*C2*t1 + 2.*amLam*C5*C6*Ekaon*t4 + 
           2.*C2*C7*Ekaon*t5) + 
           C1*C4*t3*
           ((amLam + amSig)*C5*C6*Ekaon*
           (3.*Fm1 + Fm2)*t4 - 
           3.*C2*Fm1*(t1 - 2.*C7*Ekaon*t5))))))
           + 2.*(Enu*(-36.*apkk1*(t1*t1) - 72.*apkk1*C2*(t1*t1) - 
           36.*apkk1*(C2*C2)*(t1*t1) + 72.*am*Elep*(t1*t1) + 
           72.*am*(C2*C2)*Elep*(t1*t1) + 
           36.*(amk*amk)*apkk1*C1*C4*t1*t3 - 
           36.*(am*am)*apkk1*C1*C2*C4*t1*t3 + 
           36.*(amk*amk)*apkk1*C1*C2*C4*t1*t3 - 
           36.*am*amSig*apkk1*C1*C2*C4*t1*t3 + 
           36.*(am*am)*apkk1*C1*(C4*C4)*t1*t3 - 
           36.*(amk*amk)*apkk1*C1*(C4*C4)*t1*t3 + 
           36.*am*amSig*apkk1*C1*(C4*C4)*t1*t3 - 
           72.*(am*am)*apkk1*C1*C2*(C4*C4)*t1*t3 - 
           36.*(amk*amk)*apkk1*C1*C2*(C4*C4)*t1*t3 - 
           72.*am*amSig*apkk1*C1*C2*(C4*C4)*t1*t3 - 
           72.*am*apkk1*C1*C4*Ekaon*t1*t3 - 
           72.*am*apkk1*C1*C2*C4*Ekaon*t1*t3 + 
           72.*am*apkk1*C1*(C4*C4)*Ekaon*t1*t3 + 
           72.*am*apkk1*C1*C2*(C4*C4)*Ekaon*t1*t3 - 
           72.*am*(amk*amk)*C1*C2*C4*Elep*t1*t3 + 
           72.*am*(amk*amk)*C1*(C4*C4)*Elep*t1*t3 + 
           144.*(am*am)*C1*C2*C4*Ekaon*Elep*t1*t3 - 
           144.*(am*am)*C1*(C4*C4)*Ekaon*Elep*t1*t3 + 
           36.*am*(apkk1*apkk1)*C1*C4*Fm1*t1*t3 + 
           36.*amSig*(apkk1*apkk1)*C1*C4*Fm1*t1*t3 + 
           36.*am*(apkk1*apkk1)*C1*C2*C4*Fm1*t1*t3 + 
           36.*amSig*(apkk1*apkk1)*C1*C2*C4*Fm1*t1*t3 - 
           36.*(am*am)*apkk1*C1*C4*Elep*Fm1*t1*t3 - 
           36.*am*amSig*apkk1*C1*C4*Elep*Fm1*t1*t3 - 
           36.*(am*am)*apkk1*C1*C2*C4*Elep*Fm1*t1*t3 - 
           36.*am*amSig*apkk1*C1*C2*C4*Elep*Fm1*t1*t3 + 
           36.*(am*am)*apkk1*C1*C4*Enu*Fm1*t1*t3 + 
           36.*am*amSig*apkk1*C1*C4*Enu*Fm1*t1*t3 - 
           36.*(am*am)*apkk1*C1*C2*C4*Enu*Fm1*t1*t3 - 
           36.*am*amSig*apkk1*C1*C2*C4*Enu*Fm1*t1*t3 + 
           27.*(am*am)*(amk*amk)*apkk1*(C1*C1)*(C4*C4)*(t3*t3) - 
           9.*(amk*amk*amk*amk)*apkk1*(C1*C1)*(C4*C4)*(t3*t3) + 
           36.*am*(amk*amk)*amSig*apkk1*(C1*C1)*(C4*C4)*(t3*t3) + 
           9.*(amk*amk)*(amSig*amSig)*apkk1*(C1*C1)*(C4*C4)*(t3*t3) + 
           18.*(am*am)*(amk*amk)*apkk1*(C1*C1)*(C4*C4*C4)*(t3*t3) + 
           18.*(amk*amk*amk*amk)*apkk1*(C1*C1)*(C4*C4*C4)*(t3*t3) - 
           18.*(amk*amk)*(amSig*amSig)*apkk1*(C1*C1)*(C4*C4*C4)*(t3*t3) + 
           27.*(am*am)*(amk*amk)*apkk1*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) - 
           9.*(amk*amk*amk*amk)*apkk1*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) + 
           36.*am*(amk*amk)*amSig*apkk1*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) + 
           9.*(amk*amk)*(amSig*amSig)*apkk1*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) - 
           18.*(am*am*am)*apkk1*(C1*C1)*(C4*C4)*Ekaon*(t3*t3) + 
           36.*am*(amk*amk)*apkk1*(C1*C1)*(C4*C4)*Ekaon*(t3*t3) + 
           18.*am*(amSig*amSig)*apkk1*(C1*C1)*(C4*C4)*Ekaon*
           (t3*t3) - 
           36.*(am*am*am)*apkk1*(C1*C1)*(C4*C4*C4)*Ekaon*(t3*t3) - 
           72.*am*(amk*amk)*apkk1*(C1*C1)*(C4*C4*C4)*Ekaon*(t3*t3) + 
           36.*am*(amSig*amSig)*apkk1*(C1*C1)*(C4*C4*C4)*Ekaon*
           (t3*t3) - 
           18.*(am*am*am)*apkk1*(C1*C1)*(C4*C4*C4*C4)*Ekaon*(t3*t3) + 
           36.*am*(amk*amk)*apkk1*(C1*C1)*(C4*C4*C4*C4)*Ekaon*(t3*t3) + 
           18.*am*(amSig*amSig)*apkk1*(C1*C1)*(C4*C4*C4*C4)*Ekaon*
           (t3*t3) - 
           36.*(am*am)*apkk1*(C1*C1)*(C4*C4)*(Ekaon*Ekaon)*(t3*t3) + 
           72.*(am*am)*apkk1*(C1*C1)*(C4*C4*C4)*(Ekaon*Ekaon)*(t3*t3) - 
           36.*(am*am)*apkk1*(C1*C1)*(C4*C4*C4*C4)*(Ekaon*Ekaon)*(t3*t3) - 
           18.*(am*am*am)*(amk*amk)*(C1*C1)*(C4*C4)*Elep*(t3*t3) + 
           18.*am*(amk*amk*amk*amk)*(C1*C1)*(C4*C4)*Elep*(t3*t3) - 
           36.*(am*am)*(amk*amk)*amSig*(C1*C1)*(C4*C4)*Elep*
           (t3*t3) - 
           18.*am*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4)*Elep*
           (t3*t3) - 
           18.*(am*am*am)*(amk*amk)*(C1*C1)*(C4*C4*C4*C4)*Elep*(t3*t3) + 
           18.*am*(amk*amk*amk*amk)*(C1*C1)*(C4*C4*C4*C4)*Elep*(t3*t3) - 
           36.*(am*am)*(amk*amk)*amSig*(C1*C1)*(C4*C4*C4*C4)*Elep*
           (t3*t3) - 
           18.*am*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4*C4*C4)*Elep*
           (t3*t3) - 
           72.*(am*am)*(amk*amk)*(C1*C1)*(C4*C4)*Ekaon*Elep*
           (t3*t3) - 
           72.*(am*am)*(amk*amk)*(C1*C1)*(C4*C4*C4*C4)*Ekaon*Elep*
           (t3*t3) + 
           72.*(am*am*am)*(C1*C1)*(C4*C4)*(Ekaon*Ekaon)*Elep*(t3*t3) + 
           72.*(am*am*am)*(C1*C1)*(C4*C4*C4*C4)*(Ekaon*Ekaon)*Elep*(t3*t3) + 
           36.*(amk*amk)*apkk1*C5*C6*t1*t4 - 
           36.*(am*am)*apkk1*C2*C5*C6*t1*t4 + 
           36.*(amk*amk)*apkk1*C2*C5*C6*t1*t4 - 
           36.*am*amLam*apkk1*C2*C5*C6*t1*t4 - 
           12.*(am*am)*apkk1*(C5*C5)*C6*t1*t4 + 
           12.*(amk*amk)*apkk1*(C5*C5)*C6*t1*t4 - 
           12.*am*amLam*apkk1*(C5*C5)*C6*t1*t4 + 
           24.*(am*am)*apkk1*C2*(C5*C5)*C6*t1*t4 + 
           12.*(amk*amk)*apkk1*C2*(C5*C5)*C6*t1*t4 + 
           24.*am*amLam*apkk1*C2*(C5*C5)*C6*t1*t4 - 
           72.*am*apkk1*C5*C6*Ekaon*t1*t4 - 
           72.*am*apkk1*C2*C5*C6*Ekaon*t1*t4 - 
           24.*am*apkk1*(C5*C5)*C6*Ekaon*t1*t4 - 
           24.*am*apkk1*C2*(C5*C5)*C6*Ekaon*t1*t4 - 
           72.*am*(amk*amk)*C2*C5*C6*Elep*t1*t4 - 
           24.*am*(amk*amk)*(C5*C5)*C6*Elep*t1*t4 + 
           144.*(am*am)*C2*C5*C6*Ekaon*Elep*t1*t4 + 
           48.*(am*am)*(C5*C5)*C6*Ekaon*Elep*t1*t4 + 
           12.*am*(apkk1*apkk1)*C5*C6*Fm2*t1*t4 + 
           12.*amLam*(apkk1*apkk1)*C5*C6*Fm2*t1*t4 + 
           12.*am*(apkk1*apkk1)*C2*C5*C6*Fm2*t1*t4 + 
           12.*amLam*(apkk1*apkk1)*C2*C5*C6*Fm2*t1*t4 - 
           12.*(am*am)*apkk1*C5*C6*Elep*Fm2*t1*t4 - 
           12.*am*amLam*apkk1*C5*C6*Elep*Fm2*t1*t4 - 
           12.*(am*am)*apkk1*C2*C5*C6*Elep*Fm2*t1*t4 - 
           12.*am*amLam*apkk1*C2*C5*C6*Elep*Fm2*t1*t4 + 
           12.*(am*am)*apkk1*C5*C6*Enu*Fm2*t1*t4 + 
           12.*am*amLam*apkk1*C5*C6*Enu*Fm2*t1*t4 - 
           12.*(am*am)*apkk1*C2*C5*C6*Enu*Fm2*t1*t4 - 
           12.*am*amLam*apkk1*C2*C5*C6*Enu*Fm2*t1*t4 + 
           54.*(am*am)*(amk*amk)*apkk1*C1*C4*C5*C6*t3*t4 - 
           18.*(amk*amk*amk*amk)*apkk1*C1*C4*C5*C6*t3*t4 + 
           36.*am*(amk*amk)*amLam*apkk1*C1*C4*C5*C6*t3*t4 + 
           36.*am*(amk*amk)*amSig*apkk1*C1*C4*C5*C6*t3*t4 + 
           18.*(amk*amk)*amLam*amSig*apkk1*C1*C4*C5*C6*t3*
           t4 + 18.*(am*am)*(amk*amk)*apkk1*C1*(C4*C4)*C5*C6*
           t3*t4 + 
           18.*(amk*amk*amk*amk)*apkk1*C1*(C4*C4)*C5*C6*t3*t4 - 
           18.*am*(amk*amk)*amLam*apkk1*C1*(C4*C4)*C5*C6*t3*
           t4 + 18.*am*(amk*amk)*amSig*apkk1*C1*(C4*C4)*C5*
           C6*t3*t4 - 
           18.*(amk*amk)*amLam*amSig*apkk1*C1*(C4*C4)*C5*C6*
           t3*t4 - 
           6.*(am*am)*(amk*amk)*apkk1*C1*C4*(C5*C5)*C6*t3*t4 - 
           6.*(amk*amk*amk*amk)*apkk1*C1*C4*(C5*C5)*C6*t3*t4 - 
           6.*am*(amk*amk)*amLam*apkk1*C1*C4*(C5*C5)*C6*t3*
           t4 + 6.*am*(amk*amk)*amSig*apkk1*C1*C4*(C5*C5)*C6*
           t3*t4 + 
           6.*(amk*amk)*amLam*amSig*apkk1*C1*C4*(C5*C5)*C6*t3*
           t4 - 18.*(am*am)*(amk*amk)*apkk1*C1*(C4*C4)*(C5*C5)*
           C6*t3*t4 + 
           6.*(amk*amk*amk*amk)*apkk1*C1*(C4*C4)*(C5*C5)*C6*t3*t4 - 
           12.*am*(amk*amk)*amLam*apkk1*C1*(C4*C4)*(C5*C5)*C6*
           t3*t4 - 
           12.*am*(amk*amk)*amSig*apkk1*C1*(C4*C4)*(C5*C5)*C6*
           t3*t4 - 
           6.*(amk*amk)*amLam*amSig*apkk1*C1*(C4*C4)*(C5*C5)*C6*
           t3*t4 - 
           36.*(am*am*am)*apkk1*C1*C4*C5*C6*Ekaon*t3*t4 + 
           72.*am*(amk*amk)*apkk1*C1*C4*C5*C6*Ekaon*t3*t4 + 
           36.*am*amLam*amSig*apkk1*C1*C4*C5*C6*Ekaon*t3*
           t4 - 36.*(am*am*am)*apkk1*C1*(C4*C4)*C5*C6*Ekaon*t3*
           t4 - 72.*am*(amk*amk)*apkk1*C1*(C4*C4)*C5*C6*
           Ekaon*t3*t4 + 
           36.*(am*am)*amLam*apkk1*C1*(C4*C4)*C5*C6*Ekaon*t3*
           t4 - 36.*(am*am)*amSig*apkk1*C1*(C4*C4)*C5*C6*
           Ekaon*t3*t4 + 
           36.*am*amLam*amSig*apkk1*C1*(C4*C4)*C5*C6*Ekaon*
           t3*t4 + 
           12.*(am*am*am)*apkk1*C1*C4*(C5*C5)*C6*Ekaon*t3*t4 + 
           24.*am*(amk*amk)*apkk1*C1*C4*(C5*C5)*C6*Ekaon*t3*
           t4 + 12.*(am*am)*amLam*apkk1*C1*C4*(C5*C5)*C6*
           Ekaon*t3*t4 - 
           12.*(am*am)*amSig*apkk1*C1*C4*(C5*C5)*C6*Ekaon*t3*
           t4 - 12.*am*amLam*amSig*apkk1*C1*C4*(C5*C5)*C6*
           Ekaon*t3*t4 + 
           12.*(am*am*am)*apkk1*C1*(C4*C4)*(C5*C5)*C6*Ekaon*t3*
           t4 - 24.*am*(amk*amk)*apkk1*C1*(C4*C4)*(C5*C5)*C6*
           Ekaon*t3*t4 - 
           12.*am*amLam*amSig*apkk1*C1*(C4*C4)*(C5*C5)*C6*
           Ekaon*t3*t4 - 
           72.*(am*am)*apkk1*C1*C4*C5*C6*(Ekaon*Ekaon)*t3*t4 + 
           72.*(am*am)*apkk1*C1*(C4*C4)*C5*C6*(Ekaon*Ekaon)*t3*
           t4 - 24.*(am*am)*apkk1*C1*C4*(C5*C5)*C6*(Ekaon*Ekaon)*
           t3*t4 + 
           24.*(am*am)*apkk1*C1*(C4*C4)*(C5*C5)*C6*(Ekaon*Ekaon)*t3*
           t4 - 36.*(am*am*am)*(amk*amk)*C1*C4*C5*C6*Elep*t3*
           t4 + 36.*am*(amk*amk*amk*amk)*C1*C4*C5*C6*Elep*t3*t4 - 
           36.*(am*am)*(amk*amk)*amLam*C1*C4*C5*C6*Elep*t3*
           t4 - 36.*(am*am)*(amk*amk)*amSig*C1*C4*C5*C6*Elep*
           t3*t4 - 
           36.*am*(amk*amk)*amLam*amSig*C1*C4*C5*C6*Elep*t3*
           t4 + 12.*(am*am*am)*(amk*amk)*C1*(C4*C4)*(C5*C5)*C6*Elep*
           t3*t4 - 
           12.*am*(amk*amk*amk*amk)*C1*(C4*C4)*(C5*C5)*C6*Elep*t3*t4 + 
           12.*(am*am)*(amk*amk)*amLam*C1*(C4*C4)*(C5*C5)*C6*Elep*
           t3*t4 + 
           12.*(am*am)*(amk*amk)*amSig*C1*(C4*C4)*(C5*C5)*C6*Elep*
           t3*t4 + 
           12.*am*(amk*amk)*amLam*amSig*C1*(C4*C4)*(C5*C5)*C6*
           Elep*t3*t4 - 
           144.*(am*am)*(amk*amk)*C1*C4*C5*C6*Ekaon*Elep*t3*
           t4 + 48.*(am*am)*(amk*amk)*C1*(C4*C4)*(C5*C5)*C6*
           Ekaon*Elep*t3*t4 + 
           144.*(am*am*am)*C1*C4*C5*C6*(Ekaon*Ekaon)*Elep*t3*t4 - 
           48.*(am*am*am)*C1*(C4*C4)*(C5*C5)*C6*(Ekaon*Ekaon)*Elep*t3*
           t4 + 18.*(amk*amk)*amLam*(apkk1*apkk1)*C1*C4*C5*C6*
           Fm1*t3*t4 - 
           18.*(amk*amk)*amSig*(apkk1*apkk1)*C1*C4*C5*C6*Fm1*t3*
           t4 + 6.*(amk*amk)*amLam*(apkk1*apkk1)*C1*C4*(C5*C5)*C6*
           Fm1*t3*t4 - 
           6.*(amk*amk)*amSig*(apkk1*apkk1)*C1*C4*(C5*C5)*C6*Fm1*
           t3*t4 - 
           36.*am*amLam*(apkk1*apkk1)*C1*C4*C5*C6*Ekaon*Fm1*
           t3*t4 + 
           36.*am*amSig*(apkk1*apkk1)*C1*C4*C5*C6*Ekaon*Fm1*
           t3*t4 - 
           12.*am*amLam*(apkk1*apkk1)*C1*C4*(C5*C5)*C6*Ekaon*
           Fm1*t3*t4 + 
           12.*am*amSig*(apkk1*apkk1)*C1*C4*(C5*C5)*C6*Ekaon*
           Fm1*t3*t4 - 
           18.*am*(amk*amk)*amLam*apkk1*C1*C4*C5*C6*Elep*
           Fm1*t3*t4 + 
           18.*am*(amk*amk)*amSig*apkk1*C1*C4*C5*C6*Elep*
           Fm1*t3*t4 - 
           6.*am*(amk*amk)*amLam*apkk1*C1*C4*(C5*C5)*C6*Elep*
           Fm1*t3*t4 + 
           6.*am*(amk*amk)*amSig*apkk1*C1*C4*(C5*C5)*C6*Elep*
           Fm1*t3*t4 + 
           36.*(am*am)*amLam*apkk1*C1*C4*C5*C6*Ekaon*Elep*
           Fm1*t3*t4 - 
           36.*(am*am)*amSig*apkk1*C1*C4*C5*C6*Ekaon*Elep*
           Fm1*t3*t4 + 
           12.*(am*am)*amLam*apkk1*C1*C4*(C5*C5)*C6*Ekaon*
           Elep*Fm1*t3*t4 - 
           12.*(am*am)*amSig*apkk1*C1*C4*(C5*C5)*C6*Ekaon*
           Elep*Fm1*t3*t4 - 
           18.*am*(amk*amk)*amLam*apkk1*C1*C4*C5*C6*Enu*Fm1*
           t3*t4 + 
           18.*am*(amk*amk)*amSig*apkk1*C1*C4*C5*C6*Enu*Fm1*
           t3*t4 + 
           6.*am*(amk*amk)*amLam*apkk1*C1*C4*(C5*C5)*C6*Enu*
           Fm1*t3*t4 - 
           6.*am*(amk*amk)*amSig*apkk1*C1*C4*(C5*C5)*C6*Enu*
           Fm1*t3*t4 + 
           36.*(am*am)*amLam*apkk1*C1*C4*C5*C6*Ekaon*Enu*
           Fm1*t3*t4 - 
           36.*(am*am)*amSig*apkk1*C1*C4*C5*C6*Ekaon*Enu*
           Fm1*t3*t4 - 
           12.*(am*am)*amLam*apkk1*C1*C4*(C5*C5)*C6*Ekaon*
           Enu*Fm1*t3*t4 + 
           12.*(am*am)*amSig*apkk1*C1*C4*(C5*C5)*C6*Ekaon*
           Enu*Fm1*t3*t4 - 
           6.*(amk*amk)*amLam*(apkk1*apkk1)*C1*C4*C5*C6*Fm2*t3*
           t4 + 6.*(amk*amk)*amSig*(apkk1*apkk1)*C1*C4*C5*C6*
           Fm2*t3*t4 + 
           6.*(amk*amk)*amLam*(apkk1*apkk1)*C1*(C4*C4)*C5*C6*Fm2*
           t3*t4 - 
           6.*(amk*amk)*amSig*(apkk1*apkk1)*C1*(C4*C4)*C5*C6*Fm2*
           t3*t4 + 
           12.*am*amLam*(apkk1*apkk1)*C1*C4*C5*C6*Ekaon*Fm2*
           t3*t4 - 
           12.*am*amSig*(apkk1*apkk1)*C1*C4*C5*C6*Ekaon*Fm2*
           t3*t4 - 
           12.*am*amLam*(apkk1*apkk1)*C1*(C4*C4)*C5*C6*Ekaon*
           Fm2*t3*t4 + 
           12.*am*amSig*(apkk1*apkk1)*C1*(C4*C4)*C5*C6*Ekaon*
           Fm2*t3*t4 + 
           6.*am*(amk*amk)*amLam*apkk1*C1*C4*C5*C6*Elep*Fm2*
           t3*t4 - 
           6.*am*(amk*amk)*amSig*apkk1*C1*C4*C5*C6*Elep*Fm2*
           t3*t4 - 
           6.*am*(amk*amk)*amLam*apkk1*C1*(C4*C4)*C5*C6*Elep*
           Fm2*t3*t4 + 
           6.*am*(amk*amk)*amSig*apkk1*C1*(C4*C4)*C5*C6*Elep*
           Fm2*t3*t4 - 
           12.*(am*am)*amLam*apkk1*C1*C4*C5*C6*Ekaon*Elep*
           Fm2*t3*t4 + 
           12.*(am*am)*amSig*apkk1*C1*C4*C5*C6*Ekaon*Elep*
           Fm2*t3*t4 + 
           12.*(am*am)*amLam*apkk1*C1*(C4*C4)*C5*C6*Ekaon*
           Elep*Fm2*t3*t4 - 
           12.*(am*am)*amSig*apkk1*C1*(C4*C4)*C5*C6*Ekaon*
           Elep*Fm2*t3*t4 + 
           6.*am*(amk*amk)*amLam*apkk1*C1*C4*C5*C6*Enu*Fm2*
           t3*t4 - 
           6.*am*(amk*amk)*amSig*apkk1*C1*C4*C5*C6*Enu*Fm2*
           t3*t4 + 
           6.*am*(amk*amk)*amLam*apkk1*C1*(C4*C4)*C5*C6*Enu*
           Fm2*t3*t4 - 
           6.*am*(amk*amk)*amSig*apkk1*C1*(C4*C4)*C5*C6*Enu*
           Fm2*t3*t4 - 
           12.*(am*am)*amLam*apkk1*C1*C4*C5*C6*Ekaon*Enu*
           Fm2*t3*t4 + 
           12.*(am*am)*amSig*apkk1*C1*C4*C5*C6*Ekaon*Enu*
           Fm2*t3*t4 - 
           12.*(am*am)*amLam*apkk1*C1*(C4*C4)*C5*C6*Ekaon*
           Enu*Fm2*t3*t4 + 
           12.*(am*am)*amSig*apkk1*C1*(C4*C4)*C5*C6*Ekaon*
           Enu*Fm2*t3*t4 + 
           27.*(am*am)*(amk*amk)*apkk1*(C5*C5)*(C6*C6)*(t4*t4) - 
           9.*(amk*amk*amk*amk)*apkk1*(C5*C5)*(C6*C6)*(t4*t4) + 
           36.*am*(amk*amk)*amLam*apkk1*(C5*C5)*(C6*C6)*(t4*t4) + 
           9.*(amk*amk)*(amLam*amLam)*apkk1*(C5*C5)*(C6*C6)*(t4*t4) - 
           6.*(am*am)*(amk*amk)*apkk1*(C5*C5*C5)*(C6*C6)*(t4*t4) - 
           6.*(amk*amk*amk*amk)*apkk1*(C5*C5*C5)*(C6*C6)*(t4*t4) + 
           6.*(amk*amk)*(amLam*amLam)*apkk1*(C5*C5*C5)*(C6*C6)*(t4*t4) + 
           3.*(am*am)*(amk*amk)*apkk1*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) - 
           1.*(amk*amk*amk*amk)*apkk1*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) + 
           4.*am*(amk*amk)*amLam*apkk1*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) + 
           (amk*amk)*(amLam*amLam)*apkk1*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) - 
           18.*(am*am*am)*apkk1*(C5*C5)*(C6*C6)*Ekaon*(t4*t4) + 
           36.*am*(amk*amk)*apkk1*(C5*C5)*(C6*C6)*Ekaon*(t4*t4) + 
           18.*am*(amLam*amLam)*apkk1*(C5*C5)*(C6*C6)*Ekaon*
           (t4*t4) + 
           12.*(am*am*am)*apkk1*(C5*C5*C5)*(C6*C6)*Ekaon*(t4*t4) + 
           24.*am*(amk*amk)*apkk1*(C5*C5*C5)*(C6*C6)*Ekaon*(t4*t4) - 
           12.*am*(amLam*amLam)*apkk1*(C5*C5*C5)*(C6*C6)*Ekaon*
           (t4*t4) - 
           2.*(am*am*am)*apkk1*(C5*C5*C5*C5)*(C6*C6)*Ekaon*(t4*t4) + 
           4.*am*(amk*amk)*apkk1*(C5*C5*C5*C5)*(C6*C6)*Ekaon*(t4*t4) + 
           2.*am*(amLam*amLam)*apkk1*(C5*C5*C5*C5)*(C6*C6)*Ekaon*
           (t4*t4) - 
           36.*(am*am)*apkk1*(C5*C5)*(C6*C6)*(Ekaon*Ekaon)*(t4*t4) - 
           24.*(am*am)*apkk1*(C5*C5*C5)*(C6*C6)*(Ekaon*Ekaon)*(t4*t4) - 
           4.*(am*am)*apkk1*(C5*C5*C5*C5)*(C6*C6)*(Ekaon*Ekaon)*(t4*t4) - 
           18.*(am*am*am)*(amk*amk)*(C5*C5)*(C6*C6)*Elep*(t4*t4) + 
           18.*am*(amk*amk*amk*amk)*(C5*C5)*(C6*C6)*Elep*(t4*t4) - 
           36.*(am*am)*(amk*amk)*amLam*(C5*C5)*(C6*C6)*Elep*
           (t4*t4) - 
           18.*am*(amk*amk)*(amLam*amLam)*(C5*C5)*(C6*C6)*Elep*
           (t4*t4) - 
           2.*(am*am*am)*(amk*amk)*(C5*C5*C5*C5)*(C6*C6)*Elep*(t4*t4) + 
           2.*am*(amk*amk*amk*amk)*(C5*C5*C5*C5)*(C6*C6)*Elep*(t4*t4) - 
           4.*(am*am)*(amk*amk)*amLam*(C5*C5*C5*C5)*(C6*C6)*Elep*
           (t4*t4) - 
           2.*am*(amk*amk)*(amLam*amLam)*(C5*C5*C5*C5)*(C6*C6)*Elep*
           (t4*t4) - 
           72.*(am*am)*(amk*amk)*(C5*C5)*(C6*C6)*Ekaon*Elep*
           (t4*t4) - 
           8.*(am*am)*(amk*amk)*(C5*C5*C5*C5)*(C6*C6)*Ekaon*Elep*
           (t4*t4) + 
           72.*(am*am*am)*(C5*C5)*(C6*C6)*(Ekaon*Ekaon)*Elep*(t4*t4) + 
           8.*(am*am*am)*(C5*C5*C5*C5)*(C6*C6)*(Ekaon*Ekaon)*Elep*(t4*t4) + 
           36.*am*(amk*amk)*apkk1*C1*C2*C4*C7*t3*t5 + 
           36.*(amk*amk)*amSig*apkk1*C1*C2*C4*C7*t3*t5 + 
           36.*am*(apkk1*apkk1)*C1*C2*C4*C7*t3*t5 + 
           36.*amSig*(apkk1*apkk1)*C1*C2*C4*C7*t3*t5 + 
           36.*am*(apkk1*apkk1)*C1*C2*(C4*C4)*C7*t3*t5 + 
           36.*amSig*(apkk1*apkk1)*C1*C2*(C4*C4)*C7*t3*t5 + 
           36.*(amk*amk)*(apkk1*apkk1)*C1*C2*C4*C7*Fm1*t3*t5 - 
           72.*am*(apkk1*apkk1)*C1*C2*C4*C7*Ekaon*Fm1*t3*t5 + 
           36.*am*(amk*amk)*apkk1*C2*C5*C6*C7*t4*t5 + 
           36.*(amk*amk)*amLam*apkk1*C2*C5*C6*C7*t4*t5 + 
           36.*am*(apkk1*apkk1)*C2*C5*C6*C7*t4*t5 + 
           36.*amLam*(apkk1*apkk1)*C2*C5*C6*C7*t4*t5 - 
           12.*am*(apkk1*apkk1)*C2*(C5*C5)*C6*C7*t4*t5 - 
           12.*amLam*(apkk1*apkk1)*C2*(C5*C5)*C6*C7*t4*t5 + 
           12.*(amk*amk)*(apkk1*apkk1)*C2*C5*C6*C7*Fm2*t4*t5 - 
           24.*am*(apkk1*apkk1)*C2*C5*C6*C7*Ekaon*Fm2*t4*t5 + 
           (aml*aml*aml*aml)*(-6.*(amk*amk)*(C5*C5*C5)*C6*t1*t4 + 
           12.*am*(C5*C5*C5)*C6*Ekaon*t1*t4 + 
           3.*(am*am)*(amk*amk)*C1*(C4*C4)*(C5*C5*C5)*C6*t3*t4 - 
           3.*(amk*amk*amk*amk)*C1*(C4*C4)*(C5*C5*C5)*C6*t3*t4 + 
           3.*am*(amk*amk)*amLam*C1*(C4*C4)*(C5*C5*C5)*C6*t3*
           t4 + 
           3.*am*(amk*amk)*amSig*C1*(C4*C4)*(C5*C5*C5)*C6*t3*
           t4 + 
           3.*(amk*amk)*amLam*amSig*C1*(C4*C4)*(C5*C5*C5)*C6*t3*
           t4 + 
           12.*am*(amk*amk)*C1*(C4*C4)*(C5*C5*C5)*C6*Ekaon*t3*
           t4 - 
           12.*(am*am)*C1*(C4*C4)*(C5*C5*C5)*C6*(Ekaon*Ekaon)*t3*
           t4 - 1.*(am*am)*(amk*amk)*(C5*C5*C5*C5*C5)*(C6*C6)*(t4*t4) + 
           (amk*amk*amk*amk)*(C5*C5*C5*C5*C5)*(C6*C6)*(t4*t4) - 
           2.*am*(amk*amk)*amLam*(C5*C5*C5*C5*C5)*(C6*C6)*(t4*t4) - 
           1.*(amk*amk)*(amLam*amLam)*(C5*C5*C5*C5*C5)*(C6*C6)*(t4*t4) - 
           4.*am*(amk*amk)*(C5*C5*C5*C5*C5)*(C6*C6)*Ekaon*(t4*t4) + 
           4.*(am*am)*(C5*C5*C5*C5*C5)*(C6*C6)*(Ekaon*Ekaon)*(t4*t4) + 
           3.*C3*
           (6.*t1*
           (t2 + 
           C1*(C4*C4)*((amk*amk) - 2.*am*Ekaon)*t3)
           + (amk*amk*amk*amk)*C1*(C4*C4)*t3*
           (3.*C1*(C4*C4)*t3 - 1.*(C5*C5)*C6*t4) - 
           2.*am*Ekaon*
           (t2 - 2.*am*C1*(C4*C4)*Ekaon*t3)*
           (3.*C1*(C4*C4)*t3 - 1.*(C5*C5)*C6*t4) + 
           (amk*amk)*
           (-3.*(C1*C1)*(C4*C4*C4*C4)*
           ((am*am) + (amSig*amSig) + 
           2.*am*(amSig + 2.*Ekaon))*(t3*t3) - 
           1.*(C5*C5)*C6*t2*t4 + 
           C1*(C4*C4)*t3*
           (3.*t2 + 
           (C5*C5)*C6*
           ((am*am) + amLam*amSig + 
           am*(amLam + amSig + 4.*Ekaon))*t4))
           ) + 18.*(amk*amk)*C1*C2*C4*C7*Fm1*t3*t5 - 
           36.*am*C1*C2*C4*C7*Ekaon*Fm1*t3*t5 + 
           6.*(amk*amk)*C2*C5*C6*C7*Fm2*t4*t5 - 
           12.*am*C2*C5*C6*C7*Ekaon*Fm2*t4*t5) + 
           (aml*aml)*(-36.*((1. + C2)*(1. + C2))*(t1*t1) + 
           (amk*amk*amk*amk)*
           (9.*(C1*C1)*(C4*C4)*
           (-1. + 2.*C4 + 
           (C4*C4)*
           (-1. + apkk1*C3 + 
           2.*am*C3*(-1.*Elep + Enu)) - 
           2.*am*Fm1*(-1. + Enu*Fm1))*(t3*t3) + 
           (C5*C5)*C6*t4*
           (3.*C3*t2 + 
           C6*
           (-9. - 6.*C5 - 1.*(C5*C5) + 
           (C5*C5*C5)*
           (apkk1 + 2.*am*(-1.*Elep + Enu)) - 
           2.*am*Fm2*(-3. + Enu*Fm2))*t4) - 
           3.*C1*C4*t3*
           (C5*C6*
           (6. + 2.*C5 + 
           C4*
           (-6. - 2.*C5 + 
           (C5*C5)*
           (apkk1 - 2.*am*Elep + 2.*am*Enu))
           - 6.*am*Fm1 + 3.*amLam*Fm1 - 
           3.*amSig*Fm1 - 2.*am*Fm2 - 
           1.*amLam*Fm2 + amSig*Fm2 + 
           4.*am*Enu*Fm1*Fm2)*t4 + 
           C3*C4*
           (3.*t2 + 
           (C5*C5)*C6*
           (apkk1 - 2.*am*Elep + 2.*am*Enu)*t4
           ))) - 
           6.*t1*
           (-1.*apkk1*(2. + C2)*
           (3.*amSig*C1*C4*Fm1*t3 + 
           amLam*C5*C6*Fm2*t4) + 
           (amk*amk)*
           (3.*C1*C4*
           (2.*(-1. + C4) + 
           C2*
           (-2. + 2.*C4 + 3.*am*Fm1 + 
           amSig*Fm1))*t3 + 
           3.*C3*
           (t2 - 
           1.*C1*(C4*C4)*
           (apkk1 + 2.*am*(-1.*Elep + Enu))*t3
           ) + 
           C5*C6*
           (-6. - 2.*C5 + 
           (C5*C5)*
           (apkk1 - 2.*am*Elep + 2.*am*Enu) + 
           C2*
           (-6. - 2.*C5 + 3.*am*Fm2 + 
           amLam*Fm2))*t4) - 
           1.*am*
           (3.*C1*C4*
           (apkk1*(2. + C2)*Fm1 + 
           2.*Ekaon*
           (-2. + 2.*C2*(-1. + C4) + 2.*C4 + 
           amSig*Fm1))*t3 - 
           6.*C3*
           (-1.*Ekaon*t2 + Elep*t2 - 
           1.*Enu*t2 + 
           amSig*apkk1*C1*(C4*C4)*t3 + 
           apkk1*C1*(C4*C4)*Ekaon*t3) + 
           C5*C6*
           (-2.*
           (6. + 2.*C5 - 1.*apkk1*(C5*C5) + 
           2.*C2*(3. + C5))*Ekaon + 
           apkk1*(2. + C2)*Fm2 + 
           2.*amLam*(apkk1*(C5*C5) + Ekaon*Fm2))
           *t4) - 
           2.*(am*am)*
           (apkk1*
           (-3.*C1*C3*(C4*C4)*t3 + (C5*C5*C5)*C6*t4)
           + Ekaon*
           (3.*C1*C4*
           (2.*C3*C4*(Elep - 1.*Enu) + Fm1 + 
           2.*C2*Fm1)*t3 + 
           C5*C6*
           (-2.*(C5*C5)*(Elep - 1.*Enu) + Fm2 + 
           2.*C2*Fm2)*t4))) - 
           2.*((am*am*am)*Ekaon*
           (-1.*apkk1*
           (3.*C1*(C4*C4)*t3 - 1.*(C5*C5)*C6*t4)*
           (3.*C1*C3*(C4*C4)*t3 - 1.*(C5*C5*C5)*C6*t4)
           + 2.*Ekaon*
           (9.*(C1*C1)*(C4*C4)*
           (2.*C3*(C4*C4)*(Elep - 1.*Enu) + 
           Fm1*(-1. + 2.*C4 + 2.*Enu*Fm1))*
           (t3*t3) - 
           3.*C1*C4*C5*C6*
           (2.*C3*C4*C5*(Elep - 1.*Enu) + 
           3.*Fm1 + 2.*C5*Fm1 + 
           2.*C4*
           ((C5*C5)*(Elep - 1.*Enu) - 1.*Fm2) + 
           Fm2 - 4.*Enu*Fm1*Fm2)*t3*t4 + 
           (C5*C5)*(C6*C6)*
           (2.*(C5*C5*C5)*(Elep - 1.*Enu) - 
           2.*C5*Fm2 + Fm2*(-3. + 2.*Enu*Fm2))
           *(t4*t4))) + 
           (am*am)*Ekaon*
           (-18.*(C1*C1)*(C4*C4)*Ekaon*
           (-1. + (-1. + apkk1*C3)*(C4*C4) - 
           1.*amSig*Fm1 + 
           C4*(2. - 2.*amSig*Fm1))*(t3*t3) + 
           2.*(C5*C5)*C6*t4*
           (-3.*C3*(Ekaon - 1.*Elep + Enu)*
           t2 + 
           C6*Ekaon*
           (9. + (C5*C5) - 1.*apkk1*(C5*C5*C5) + 
           3.*amLam*Fm2 + 
           C5*(6. - 2.*amLam*Fm2))*t4) + 
           3.*C1*C4*t3*
           (C5*C6*
           (2.*
           (2.*(3. + C5) + 
           C4*(-6. - 2.*C5 + apkk1*(C5*C5)))*
           Ekaon + 
           amLam*
           (apkk1*C4*(C5*C5) + 
           2.*Ekaon*
           (3.*Fm1 - 1.*C5*Fm1 + C4*Fm2)) + 
           amSig*
           (-1.*apkk1*C4*(C5*C5) + 
           2.*Ekaon*
           (-1.*C5*Fm1 + Fm2 + C4*Fm2)))*t4 + 
           C3*C4*
           (6.*Ekaon*t2 - 6.*Elep*t2 + 
           6.*Enu*t2 - 
           1.*amLam*apkk1*(C5*C5)*C6*t4 + 
           amSig*apkk1*(C5*C5)*C6*t4 + 
           2.*apkk1*(C5*C5)*C6*Ekaon*t4))) - 
           9.*apkk1*C2*C7*
           (amSig*C1*C4*t3 + amLam*C5*C6*t4)*t5 + 
           am*apkk1*
           (9.*(amSig*amSig)*(C1*C1)*C3*(C4*C4*C4*C4)*Ekaon*
           (t3*t3) - 
           3.*amSig*C1*C4*C5*C6*Ekaon*
           (amLam*C4*C5*(C3 + C5) + 3.*Fm1 + 
           2.*C5*Fm1 - 1.*Fm2 + 2.*C4*Fm2)*t3*
           t4 + 
           3.*amLam*C1*C4*C5*C6*Ekaon*
           ((3. + 2.*C5)*Fm1 + 
           (-1. + 2.*C4)*Fm2)*t3*t4 + 
           (amLam*amLam)*(C5*C5*C5*C5*C5)*(C6*C6)*Ekaon*(t4*t4) + 
           9.*C2*C7*
           (C1*C4*(-1. + 6.*Ekaon*Fm1)*t3 + 
           C5*C6*(-1. + 2.*Ekaon*Fm2)*t4)*t5))
           + (amk*amk)*
           (9.*(amSig*amSig)*(C1*C1)*(C4*C4)*(t3*t3) - 
           18.*(amSig*amSig)*(C1*C1)*(C4*C4*C4)*(t3*t3) + 
           9.*(amSig*amSig)*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) - 
           9.*(amSig*amSig)*apkk1*(C1*C1)*C3*(C4*C4*C4*C4)*
           (t3*t3) + 
           18.*amLam*amSig*C1*C4*C5*C6*t3*t4 - 
           18.*amLam*amSig*C1*(C4*C4)*C5*C6*t3*t4 + 
           6.*amLam*amSig*C1*C4*(C5*C5)*C6*t3*t4 - 
           6.*amLam*amSig*C1*(C4*C4)*(C5*C5)*C6*t3*
           t4 + 
           3.*amLam*amSig*apkk1*C1*C3*(C4*C4)*(C5*C5)*
           C6*t3*t4 + 
           3.*amLam*amSig*apkk1*C1*(C4*C4)*(C5*C5*C5)*C6*
           t3*t4 + 
           9.*amLam*apkk1*C1*C4*C5*C6*Fm1*t3*t4 - 
           9.*amSig*apkk1*C1*C4*C5*C6*Fm1*t3*t4 + 
           6.*amLam*apkk1*C1*C4*(C5*C5)*C6*Fm1*t3*
           t4 - 
           6.*amSig*apkk1*C1*C4*(C5*C5)*C6*Fm1*t3*
           t4 - 
           3.*amLam*apkk1*C1*C4*C5*C6*Fm2*t3*t4 + 
           3.*amSig*apkk1*C1*C4*C5*C6*Fm2*t3*t4 + 
           6.*amLam*apkk1*C1*(C4*C4)*C5*C6*Fm2*t3*
           t4 - 
           6.*amSig*apkk1*C1*(C4*C4)*C5*C6*Fm2*t3*
           t4 + 9.*(amLam*amLam)*(C5*C5)*(C6*C6)*(t4*t4) + 
           6.*(amLam*amLam)*(C5*C5*C5)*(C6*C6)*(t4*t4) + 
           (amLam*amLam)*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) - 
           1.*(amLam*amLam)*apkk1*(C5*C5*C5*C5*C5)*(C6*C6)*(t4*t4) + 
           2.*(am*am*am)*
           (9.*(C1*C1)*(C4*C4)*
           (C3*(C4*C4)*(Elep - 1.*Enu) + 
           Fm1*(2.*C4 + Enu*Fm1))*(t3*t3) + 
           3.*C1*C4*C5*C6*
           (-1.*C4*(C5*C5)*Elep + C4*(C5*C5)*Enu + 
           C3*C4*C5*(-1.*Elep + Enu) - 
           2.*C5*Fm1 + 2.*C4*Fm2 + 
           2.*Enu*Fm1*Fm2)*t3*t4 + 
           (C5*C5)*(C6*C6)*
           ((C5*C5*C5)*(Elep - 1.*Enu) - 
           2.*C5*Fm2 + Enu*(Fm2*Fm2))*(t4*t4)) + 
           (am*am)*
           (-9.*(C1*C1)*(C4*C4)*
           (-1. + 
           (C4*C4)*
           (-1. + 3.*apkk1*C3 - 
           8.*C3*Ekaon*Elep + 
           8.*C3*Ekaon*Enu + 
           4.*amSig*C3*(-1.*Elep + Enu)) - 
           4.*amSig*Enu*(Fm1*Fm1) + 
           C4*
           (2. - 8.*amSig*Fm1 - 4.*Ekaon*Fm1)
           + 2.*Ekaon*Fm1*(3. - 4.*Enu*Fm1))*
           (t3*t3) + 
           3.*C1*C4*C5*C6*
           (C4*
           (-6. + 
           (C5*C5)*
           (3.*apkk1 - 
           2.*(amLam + amSig + 4.*Ekaon)*
           (Elep - 1.*Enu)) + 
           C5*
           (-2. + 3.*apkk1*C3 - 
           2.*amLam*C3*Elep - 
           2.*amSig*C3*Elep - 
           8.*C3*Ekaon*Elep + 
           2.*amLam*C3*Enu + 
           2.*amSig*C3*Enu + 8.*C3*Ekaon*Enu)
           + 4.*amLam*Fm2 + 4.*amSig*Fm2 + 
           4.*Ekaon*Fm2) + 
           2.*
           (3. - 9.*Ekaon*Fm1 - 
           1.*C5*
           (-1. + 2.*amLam*Fm1 + 
           2.*amSig*Fm1 + 2.*Ekaon*Fm1) - 
           3.*Ekaon*Fm2 + 
           2.*amLam*Enu*Fm1*Fm2 + 
           2.*amSig*Enu*Fm1*Fm2 + 
           8.*Ekaon*Enu*Fm1*Fm2))*t3*t4 + 
           (C5*C5)*(C6*C6)*
           (9. + (C5*C5) + 
           (C5*C5*C5)*
           (-3.*apkk1 + 
           4.*(amLam + 2.*Ekaon)*
           (Elep - 1.*Enu)) + 
           4.*amLam*Enu*(Fm2*Fm2) + 
           C5*
           (6. - 8.*amLam*Fm2 - 4.*Ekaon*Fm2)
           + 2.*Ekaon*Fm2*(-9. + 4.*Enu*Fm2))*
           (t4*t4)) + 
           18.*amSig*C1*C2*C4*C7*t3*t5 - 
           36.*amSig*C1*C2*(C4*C4)*C7*t3*t5 + 
           54.*apkk1*C1*C2*C4*C7*Fm1*t3*t5 + 
           18.*amLam*C2*C5*C6*C7*t4*t5 + 
           12.*amLam*C2*(C5*C5)*C6*C7*t4*t5 + 
           18.*apkk1*C2*C5*C6*C7*Fm2*t4*t5 + 
           am*
           (18.*(C1*C1)*(C4*C4)*
           (2.*
           (1. - 2.*C4 + 
           (1. - 1.*apkk1*C3)*(C4*C4))*Ekaon + 
           amSig*
           (1. + (1. - 2.*apkk1*C3)*(C4*C4) + 
           Ekaon*Fm1 + 2.*C4*(-1. + Ekaon*Fm1)
           ) + 
           (amSig*amSig)*
           (C3*(C4*C4)*(Elep - 1.*Enu) + 
           Fm1*(2.*C4 + Enu*Fm1)))*(t3*t3) + 
           3.*C1*C4*t3*
           (6.*amSig*C5*C6*t4 - 
           6.*amSig*C4*C5*C6*t4 + 
           2.*amSig*(C5*C5)*C6*t4 - 
           2.*amSig*C4*(C5*C5)*C6*t4 + 
           amSig*apkk1*C4*(C5*C5*C5)*C6*t4 + 
           24.*C5*C6*Ekaon*t4 - 
           24.*C4*C5*C6*Ekaon*t4 + 
           8.*(C5*C5)*C6*Ekaon*t4 - 
           8.*C4*(C5*C5)*C6*Ekaon*t4 + 
           4.*apkk1*C4*(C5*C5*C5)*C6*Ekaon*t4 - 
           6.*amSig*C5*C6*Ekaon*Fm1*t4 - 
           2.*amSig*(C5*C5)*C6*Ekaon*Fm1*t4 + 
           4.*amSig*C5*C6*Ekaon*Fm2*t4 + 
           2.*amSig*C4*C5*C6*Ekaon*Fm2*t4 + 
           amLam*C5*C6*
           (C4*
           (-6. - 2.*C5 + 
           (C5*C5)*
           (3.*apkk1 - 2.*amSig*Elep + 
           2.*amSig*Enu) + 4.*amSig*Fm2 + 
           2.*Ekaon*Fm2) + 
           2.*
           (3. + C5 - 2.*amSig*C5*Fm1 + 
           6.*Ekaon*Fm1 - 1.*C5*Ekaon*Fm1 - 
           1.*Ekaon*Fm2 + 2.*amSig*Enu*Fm1*Fm2
           ))*t4 + 
           C3*C4*
           (-6.*Elep*t2 + 6.*Enu*t2 + 
           amLam*apkk1*(C5*C5)*C6*t4 + 
           3.*amSig*apkk1*(C5*C5)*C6*t4 - 
           2.*amLam*amSig*(C5*C5)*C6*Elep*t4 + 
           2.*amLam*amSig*(C5*C5)*C6*Enu*t4 + 
           4.*Ekaon*
           (3.*t2 + apkk1*(C5*C5)*C6*t4)) + 
           6.*C2*C7*t5 - 12.*C2*C4*C7*t5) + 
           2.*C5*C6*t4*
           (-3.*C3*C5*
           (2.*Ekaon - 1.*Elep + Enu)*t2 + 
           18.*C5*C6*Ekaon*t4 + 
           12.*(C5*C5)*C6*Ekaon*t4 + 
           2.*(C5*C5*C5)*C6*Ekaon*t4 - 
           2.*apkk1*(C5*C5*C5*C5)*C6*Ekaon*t4 + 
           (amLam*amLam)*C5*C6*
           ((C5*C5*C5)*(Elep - 1.*Enu) - 
           2.*C5*Fm2 + Enu*(Fm2*Fm2))*t4 + 
           amLam*C5*C6*
           (9. + (C5*C5) - 2.*apkk1*(C5*C5*C5) + 
           3.*Ekaon*Fm2 + 
           C5*(6. - 2.*Ekaon*Fm2))*t4 + 
           9.*C2*C7*t5 + 6.*C2*C5*C7*t5)))))
           + 2.*(akpk*akpk)*
           (2.*(am*am)*(aml*aml)*Ekaon*
           (9.*(C1*C1)*(C4*C4)*(C3*(C4*C4) - 1.*(Fm1*Fm1))*
           (t3*t3) - 
           3.*C1*C4*C5*C6*
           (C3*C4*C5 + C4*(C5*C5) + 2.*Fm1*Fm2)*t3*t4
           + (C5*C5)*(C6*C6)*((C5*C5*C5) - 1.*(Fm2*Fm2))*(t4*t4))
           + 3.*Elep*
           (-2.*(amk*amk)*C2*C7*
           (3.*C1*C4*Fm1*t3 + C5*C6*Fm2*t4)*t5 + 
           amSig*C1*C4*t3*
           (-1.*(amk*amk)*(1. + C4)*C5*C6*Fm2*t4 + 
           Fm1*
           (-6.*(-1. + C2)*t1 - 
           1.*(amk*amk)*(-3. + C5)*C5*C6*t4) + 
           6.*C2*(-1. + C4)*C7*t5) + 
           amLam*C5*C6*t4*
           ((amk*amk)*C1*C4*(-3. + C5)*Fm1*t3 + 
           Fm2*
           (-2.*(-1. + C2)*t1 + 
           (amk*amk)*C1*C4*(1. + C4)*t3) - 
           2.*C2*(3. + C5)*C7*t5)) + 
           (aml*aml)*(-18.*(amSig*amSig)*(C1*C1)*(C4*C4)*Ekaon*
           (C3*(C4*C4) - 1.*(Fm1*Fm1))*(t3*t3) + 
           3.*amSig*C1*C4*t3*
           (2.*amLam*C4*(C5*C5*C5)*C6*Ekaon*t4 + 
           4.*amLam*C5*C6*Ekaon*Fm1*Fm2*t4 + 
           C3*C4*
           (-6.*t1 - 3.*t2 - 
           6.*(amk*amk)*C1*(C4*C4)*t3 + 
           (amk*amk)*(C5*C5)*C6*t4 + 
           2.*amLam*(C5*C5)*C6*Ekaon*t4) + 
           (amk*amk)*
           (6.*C1*C4*(Fm1*Fm1)*t3 + C4*(C5*C5*C5)*C6*t4 + 
           2.*C5*C6*Fm1*Fm2*t4) + 
           6.*C2*C7*Fm1*t5) + 
           amLam*C5*C6*t4*
           (3.*(C5*C5)*(2.*t1 + (amk*amk)*C1*(C4*C4)*t3) - 
           2.*(C5*C5*C5*C5)*C6*((amk*amk) + amLam*Ekaon)*t4 + 
           C5*
           (3.*C3*(t2 + (amk*amk)*C1*(C4*C4)*t3) + 
           2.*C6*((amk*amk) + amLam*Ekaon)*(Fm2*Fm2)*
           t4) + 
           6.*Fm2*((amk*amk)*C1*C4*Fm1*t3 + C2*C7*t5))
           ) + am*
           ((aml*aml)*
           (-18.*(amk*amk)*(C1*C1)*(C4*C4)*
           (C3*(C4*C4) - 1.*(Fm1*Fm1))*(t3*t3) + 
           3.*C1*C4*t3*
           (2.*(amk*amk)*C5*C6*
           (C4*(C5*C5) + 2.*Fm1*Fm2)*t4 + 
           C3*C4*
           (-6.*t1 - 3.*t2 + 
           2.*(amk*amk)*(C5*C5)*C6*t4) + 
           6.*C2*C7*Fm1*t5) + 
           C5*C6*t4*
           (6.*(C5*C5)*t1 + 3.*C3*C5*t2 - 
           2.*(amk*amk)*(C5*C5*C5*C5)*C6*t4 + 
           2.*(amk*amk)*C5*C6*(Fm2*Fm2)*t4 + 
           6.*C2*C7*Fm2*t5)) + 
           6.*Elep*
           (-1.*C5*C6*t4*
           (C2*(3. + C5)*C7*t5 + 
           Fm2*
           ((-1. + C2)*t1 - 2.*C2*C7*Ekaon*t5))
           + C1*C4*t3*
           (amSig*C5*C6*Ekaon*Fm2*t4 + 
           amSig*C4*C5*C6*Ekaon*Fm2*t4 - 
           1.*amLam*(1. + C4)*C5*C6*Ekaon*Fm2*
           t4 - 3.*C2*C7*t5 + 3.*C2*C4*C7*t5 + 
           Fm1*
           (-3.*(-1. + C2)*t1 + 
           Ekaon*
           (-1.*amLam*(-3. + C5)*C5*C6*t4 + 
           amSig*(-3. + C5)*C5*C6*t4 + 
           6.*C2*C7*t5)))))) + 
           akpk*(-36.*(aml*aml)*C3*Enu*t1*t2 + 
           36.*am*(aml*aml)*C1*C4*t1*t3 + 
           36.*(aml*aml)*amSig*C1*C4*t1*t3 + 
           36.*am*(aml*aml)*C1*C2*C4*t1*t3 + 
           36.*(aml*aml)*amSig*C1*C2*C4*t1*t3 + 
           72.*am*apkk1*C1*C2*C4*t1*t3 + 
           72.*amSig*apkk1*C1*C2*C4*t1*t3 - 
           36.*am*(aml*aml)*C1*(C4*C4)*t1*t3 - 
           36.*(aml*aml)*amSig*C1*(C4*C4)*t1*t3 - 
           72.*am*apkk1*C1*(C4*C4)*t1*t3 - 
           72.*amSig*apkk1*C1*(C4*C4)*t1*t3 - 
           36.*am*(aml*aml)*C1*C2*(C4*C4)*t1*t3 - 
           36.*(aml*aml)*amSig*C1*C2*(C4*C4)*t1*t3 + 
           18.*am*(aml*aml*aml*aml)*C1*C3*(C4*C4)*t1*t3 + 
           18.*(aml*aml*aml*aml)*amSig*C1*C3*(C4*C4)*t1*t3 + 
           36.*am*(aml*aml)*apkk1*C1*C3*(C4*C4)*t1*t3 + 
           36.*(aml*aml)*amSig*apkk1*C1*C3*(C4*C4)*t1*t3 + 
           36.*(am*am)*(aml*aml)*C1*C3*(C4*C4)*Enu*t1*t3 - 
           36.*(amk*amk)*(aml*aml)*C1*C3*(C4*C4)*Enu*t1*t3 + 
           36.*am*(aml*aml)*amSig*C1*C3*(C4*C4)*Enu*t1*t3 + 
           72.*am*(aml*aml)*C1*C3*(C4*C4)*Ekaon*Enu*t1*t3 - 
           36.*(am*am)*(aml*aml)*C1*C4*Fm1*t1*t3 - 
           36.*(amk*amk)*(aml*aml)*C1*C4*Fm1*t1*t3 - 
           36.*am*(aml*aml)*amSig*C1*C4*Fm1*t1*t3 - 
           36.*(am*am)*(aml*aml)*C1*C2*C4*Fm1*t1*t3 + 
           18.*(amk*amk)*(aml*aml)*C1*C2*C4*Fm1*t1*t3 - 
           36.*am*(aml*aml)*amSig*C1*C2*C4*Fm1*t1*t3 + 
           36.*am*(aml*aml)*C1*C4*Ekaon*Fm1*t1*t3 - 
           36.*(aml*aml)*amSig*C1*C4*Ekaon*Fm1*t1*t3 + 
           36.*(aml*aml)*amSig*C1*C2*C4*Ekaon*Fm1*t1*t3 - 
           36.*am*apkk1*C1*C4*Enu*Fm1*t1*t3 - 
           36.*amSig*apkk1*C1*C4*Enu*Fm1*t1*t3 + 
           36.*am*apkk1*C1*C2*C4*Enu*Fm1*t1*t3 + 
           36.*amSig*apkk1*C1*C2*C4*Enu*Fm1*t1*t3 - 
           9.*am*(amk*amk)*(aml*aml)*C1*C3*(C4*C4)*t2*t3 + 
           9.*am*(aml*aml*aml*aml)*C1*C3*(C4*C4)*t2*t3 - 
           9.*(amk*amk)*(aml*aml)*amSig*C1*C3*(C4*C4)*t2*t3 + 
           9.*(aml*aml*aml*aml)*amSig*C1*C3*(C4*C4)*t2*t3 + 
           18.*am*(aml*aml)*apkk1*C1*C3*(C4*C4)*t2*t3 + 
           18.*(aml*aml)*amSig*apkk1*C1*C3*(C4*C4)*t2*t3 + 
           18.*(am*am)*(aml*aml)*C1*C3*(C4*C4)*Ekaon*t2*t3 + 
           18.*am*(aml*aml)*amSig*C1*C3*(C4*C4)*Ekaon*t2*t3 + 
           18.*(am*am)*(aml*aml)*C1*C3*(C4*C4)*Enu*t2*t3 - 
           18.*(amk*amk)*(aml*aml)*C1*C3*(C4*C4)*Enu*t2*t3 + 
           18.*am*(aml*aml)*amSig*C1*C3*(C4*C4)*Enu*t2*t3 + 
           36.*am*(aml*aml)*C1*C3*(C4*C4)*Ekaon*Enu*t2*t3 - 
           18.*am*(amk*amk)*(aml*aml)*(C1*C1)*(C4*C4)*(t3*t3) - 
           18.*(amk*amk)*(aml*aml)*amSig*(C1*C1)*(C4*C4)*(t3*t3) - 
           36.*am*(amk*amk)*apkk1*(C1*C1)*(C4*C4)*(t3*t3) - 
           36.*(amk*amk)*amSig*apkk1*(C1*C1)*(C4*C4)*(t3*t3) + 
           36.*am*(amk*amk)*(aml*aml)*(C1*C1)*(C4*C4*C4)*(t3*t3) + 
           36.*(amk*amk)*(aml*aml)*amSig*(C1*C1)*(C4*C4*C4)*(t3*t3) - 
           18.*am*(amk*amk)*(aml*aml)*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) - 
           18.*(amk*amk)*(aml*aml)*amSig*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) - 
           36.*am*(amk*amk)*apkk1*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) - 
           36.*(amk*amk)*amSig*apkk1*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) + 
           18.*am*(amk*amk)*(aml*aml*aml*aml)*(C1*C1)*C3*(C4*C4*C4*C4)*(t3*t3) + 
           18.*(amk*amk)*(aml*aml*aml*aml)*amSig*(C1*C1)*C3*(C4*C4*C4*C4)*
           (t3*t3) + 
           36.*am*(amk*amk)*(aml*aml)*apkk1*(C1*C1)*C3*(C4*C4*C4*C4)*
           (t3*t3) + 
           36.*(amk*amk)*(aml*aml)*amSig*apkk1*(C1*C1)*C3*(C4*C4*C4*C4)*
           (t3*t3) + 
           18.*(am*am)*(aml*aml)*(C1*C1)*(C4*C4)*Ekaon*(t3*t3) - 
           18.*(aml*aml)*(amSig*amSig)*(C1*C1)*(C4*C4)*Ekaon*(t3*t3) + 
           36.*(am*am)*apkk1*(C1*C1)*(C4*C4)*Ekaon*(t3*t3) - 
           36.*(amSig*amSig)*apkk1*(C1*C1)*(C4*C4)*Ekaon*(t3*t3) - 
           36.*(am*am)*(aml*aml)*(C1*C1)*(C4*C4*C4)*Ekaon*(t3*t3) + 
           36.*(aml*aml)*(amSig*amSig)*(C1*C1)*(C4*C4*C4)*Ekaon*(t3*t3) + 
           18.*(am*am)*(aml*aml)*(C1*C1)*(C4*C4*C4*C4)*Ekaon*(t3*t3) - 
           18.*(aml*aml)*(amSig*amSig)*(C1*C1)*(C4*C4*C4*C4)*Ekaon*(t3*t3) + 
           36.*(am*am)*apkk1*(C1*C1)*(C4*C4*C4*C4)*Ekaon*(t3*t3) - 
           36.*(amSig*amSig)*apkk1*(C1*C1)*(C4*C4*C4*C4)*Ekaon*(t3*t3) - 
           18.*(am*am)*(aml*aml*aml*aml)*(C1*C1)*C3*(C4*C4*C4*C4)*Ekaon*(t3*t3) + 
           18.*(aml*aml*aml*aml)*(amSig*amSig)*(C1*C1)*C3*(C4*C4*C4*C4)*Ekaon*
           (t3*t3) - 
           36.*(am*am)*(aml*aml)*apkk1*(C1*C1)*C3*(C4*C4*C4*C4)*Ekaon*
           (t3*t3) + 
           36.*(aml*aml)*(amSig*amSig)*apkk1*(C1*C1)*C3*(C4*C4*C4*C4)*
           Ekaon*(t3*t3) + 
           54.*(am*am)*(amk*amk)*(aml*aml)*(C1*C1)*C3*(C4*C4*C4*C4)*Enu*
           (t3*t3) - 
           18.*(amk*amk*amk*amk)*(aml*aml)*(C1*C1)*C3*(C4*C4*C4*C4)*Enu*(t3*t3) + 
           72.*am*(amk*amk)*(aml*aml)*amSig*(C1*C1)*C3*(C4*C4*C4*C4)*Enu*
           (t3*t3) + 
           18.*(amk*amk)*(aml*aml)*(amSig*amSig)*(C1*C1)*C3*(C4*C4*C4*C4)*Enu*
           (t3*t3) - 
           36.*(am*am*am)*(aml*aml)*(C1*C1)*C3*(C4*C4*C4*C4)*Ekaon*Enu*
           (t3*t3) + 
           72.*am*(amk*amk)*(aml*aml)*(C1*C1)*C3*(C4*C4*C4*C4)*Ekaon*Enu*
           (t3*t3) + 
           36.*am*(aml*aml)*(amSig*amSig)*(C1*C1)*C3*(C4*C4*C4*C4)*Ekaon*
           Enu*(t3*t3) - 
           72.*(am*am)*(aml*aml)*(C1*C1)*C3*(C4*C4*C4*C4)*(Ekaon*Ekaon)*Enu*
           (t3*t3) + 
           9.*(am*am)*(amk*amk)*(aml*aml)*(C1*C1)*(C4*C4)*Fm1*
           (t3*t3) - 
           9.*(amk*amk*amk*amk)*(aml*aml)*(C1*C1)*(C4*C4)*Fm1*(t3*t3) - 
           9.*(amk*amk)*(aml*aml)*(amSig*amSig)*(C1*C1)*(C4*C4)*Fm1*
           (t3*t3) - 
           54.*(am*am)*(amk*amk)*(aml*aml)*(C1*C1)*(C4*C4*C4)*Fm1*
           (t3*t3) - 
           18.*(amk*amk*amk*amk)*(aml*aml)*(C1*C1)*(C4*C4*C4)*Fm1*(t3*t3) - 
           72.*am*(amk*amk)*(aml*aml)*amSig*(C1*C1)*(C4*C4*C4)*Fm1*
           (t3*t3) - 
           18.*(amk*amk)*(aml*aml)*(amSig*amSig)*(C1*C1)*(C4*C4*C4)*Fm1*
           (t3*t3) - 
           18.*(am*am*am)*(aml*aml)*(C1*C1)*(C4*C4)*Ekaon*Fm1*
           (t3*t3) + 
           18.*am*(amk*amk)*(aml*aml)*(C1*C1)*(C4*C4)*Ekaon*Fm1*
           (t3*t3) - 
           18.*(amk*amk)*(aml*aml)*amSig*(C1*C1)*(C4*C4)*Ekaon*Fm1*
           (t3*t3) + 
           18.*am*(aml*aml)*(amSig*amSig)*(C1*C1)*(C4*C4)*Ekaon*Fm1*
           (t3*t3) + 
           36.*(am*am*am)*(aml*aml)*(C1*C1)*(C4*C4*C4)*Ekaon*Fm1*
           (t3*t3) + 
           36.*am*(amk*amk)*(aml*aml)*(C1*C1)*(C4*C4*C4)*Ekaon*Fm1*
           (t3*t3) - 
           36.*(amk*amk)*(aml*aml)*amSig*(C1*C1)*(C4*C4*C4)*Ekaon*Fm1*
           (t3*t3) - 
           36.*am*(aml*aml)*(amSig*amSig)*(C1*C1)*(C4*C4*C4)*Ekaon*Fm1*
           (t3*t3) + 
           36.*am*(aml*aml)*amSig*(C1*C1)*(C4*C4)*(Ekaon*Ekaon)*Fm1*
           (t3*t3) + 
           72.*am*(aml*aml)*amSig*(C1*C1)*(C4*C4*C4)*(Ekaon*Ekaon)*Fm1*
           (t3*t3) - 
           54.*(am*am)*(amk*amk)*(aml*aml)*(C1*C1)*(C4*C4)*Enu*
           (Fm1*Fm1)*(t3*t3) + 
           18.*(amk*amk*amk*amk)*(aml*aml)*(C1*C1)*(C4*C4)*Enu*(Fm1*Fm1)*
           (t3*t3) - 
           72.*am*(amk*amk)*(aml*aml)*amSig*(C1*C1)*(C4*C4)*Enu*
           (Fm1*Fm1)*(t3*t3) - 
           18.*(amk*amk)*(aml*aml)*(amSig*amSig)*(C1*C1)*(C4*C4)*Enu*
           (Fm1*Fm1)*(t3*t3) + 
           36.*(am*am*am)*(aml*aml)*(C1*C1)*(C4*C4)*Ekaon*Enu*(Fm1*Fm1)*
           (t3*t3) - 
           72.*am*(amk*amk)*(aml*aml)*(C1*C1)*(C4*C4)*Ekaon*Enu*
           (Fm1*Fm1)*(t3*t3) - 
           36.*am*(aml*aml)*(amSig*amSig)*(C1*C1)*(C4*C4)*Ekaon*Enu*
           (Fm1*Fm1)*(t3*t3) + 
           72.*(am*am)*(aml*aml)*(C1*C1)*(C4*C4)*(Ekaon*Ekaon)*Enu*
           (Fm1*Fm1)*(t3*t3) + 36.*am*(aml*aml)*C5*C6*t1*t4 + 
           36.*(aml*aml)*amLam*C5*C6*t1*t4 + 
           36.*am*(aml*aml)*C2*C5*C6*t1*t4 + 
           36.*(aml*aml)*amLam*C2*C5*C6*t1*t4 + 
           72.*am*apkk1*C2*C5*C6*t1*t4 + 
           72.*amLam*apkk1*C2*C5*C6*t1*t4 + 
           12.*am*(aml*aml)*(C5*C5)*C6*t1*t4 + 
           12.*(aml*aml)*amLam*(C5*C5)*C6*t1*t4 + 
           24.*am*apkk1*(C5*C5)*C6*t1*t4 + 
           24.*amLam*apkk1*(C5*C5)*C6*t1*t4 + 
           12.*am*(aml*aml)*C2*(C5*C5)*C6*t1*t4 + 
           12.*(aml*aml)*amLam*C2*(C5*C5)*C6*t1*t4 - 
           6.*am*(aml*aml*aml*aml)*(C5*C5*C5)*C6*t1*t4 - 
           6.*(aml*aml*aml*aml)*amLam*(C5*C5*C5)*C6*t1*t4 - 
           12.*am*(aml*aml)*apkk1*(C5*C5*C5)*C6*t1*t4 - 
           12.*(aml*aml)*amLam*apkk1*(C5*C5*C5)*C6*t1*t4 - 
           12.*(am*am)*(aml*aml)*(C5*C5*C5)*C6*Enu*t1*t4 + 
           12.*(amk*amk)*(aml*aml)*(C5*C5*C5)*C6*Enu*t1*t4 - 
           12.*am*(aml*aml)*amLam*(C5*C5*C5)*C6*Enu*t1*t4 - 
           24.*am*(aml*aml)*(C5*C5*C5)*C6*Ekaon*Enu*t1*t4 - 
           12.*(am*am)*(aml*aml)*C5*C6*Fm2*t1*t4 - 
           12.*(amk*amk)*(aml*aml)*C5*C6*Fm2*t1*t4 - 
           12.*am*(aml*aml)*amLam*C5*C6*Fm2*t1*t4 - 
           12.*(am*am)*(aml*aml)*C2*C5*C6*Fm2*t1*t4 + 
           6.*(amk*amk)*(aml*aml)*C2*C5*C6*Fm2*t1*t4 - 
           12.*am*(aml*aml)*amLam*C2*C5*C6*Fm2*t1*t4 + 
           12.*am*(aml*aml)*C5*C6*Ekaon*Fm2*t1*t4 - 
           12.*(aml*aml)*amLam*C5*C6*Ekaon*Fm2*t1*t4 + 
           12.*(aml*aml)*amLam*C2*C5*C6*Ekaon*Fm2*t1*t4 - 
           12.*am*apkk1*C5*C6*Enu*Fm2*t1*t4 - 
           12.*amLam*apkk1*C5*C6*Enu*Fm2*t1*t4 + 
           12.*am*apkk1*C2*C5*C6*Enu*Fm2*t1*t4 + 
           12.*amLam*apkk1*C2*C5*C6*Enu*Fm2*t1*t4 + 
           3.*am*(amk*amk)*(aml*aml)*C3*(C5*C5)*C6*t2*t4 - 
           3.*am*(aml*aml*aml*aml)*C3*(C5*C5)*C6*t2*t4 + 
           3.*(amk*amk)*(aml*aml)*amLam*C3*(C5*C5)*C6*t2*t4 - 
           3.*(aml*aml*aml*aml)*amLam*C3*(C5*C5)*C6*t2*t4 - 
           6.*am*(aml*aml)*apkk1*C3*(C5*C5)*C6*t2*t4 - 
           6.*(aml*aml)*amLam*apkk1*C3*(C5*C5)*C6*t2*t4 - 
           6.*(am*am)*(aml*aml)*C3*(C5*C5)*C6*Ekaon*t2*t4 - 
           6.*am*(aml*aml)*amLam*C3*(C5*C5)*C6*Ekaon*t2*t4 - 
           6.*(am*am)*(aml*aml)*C3*(C5*C5)*C6*Enu*t2*t4 + 
           6.*(amk*amk)*(aml*aml)*C3*(C5*C5)*C6*Enu*t2*t4 - 
           6.*am*(aml*aml)*amLam*C3*(C5*C5)*C6*Enu*t2*t4 - 
           12.*am*(aml*aml)*C3*(C5*C5)*C6*Ekaon*Enu*t2*t4 - 
           36.*am*(amk*amk)*(aml*aml)*C1*C4*C5*C6*t3*t4 - 
           18.*(amk*amk)*(aml*aml)*amLam*C1*C4*C5*C6*t3*t4 - 
           18.*(amk*amk)*(aml*aml)*amSig*C1*C4*C5*C6*t3*t4 - 
           72.*am*(amk*amk)*apkk1*C1*C4*C5*C6*t3*t4 - 
           36.*(amk*amk)*amLam*apkk1*C1*C4*C5*C6*t3*t4 - 
           36.*(amk*amk)*amSig*apkk1*C1*C4*C5*C6*t3*t4 + 
           36.*am*(amk*amk)*(aml*aml)*C1*(C4*C4)*C5*C6*t3*t4 + 
           18.*(amk*amk)*(aml*aml)*amLam*C1*(C4*C4)*C5*C6*t3*
           t4 + 18.*(amk*amk)*(aml*aml)*amSig*C1*(C4*C4)*C5*C6*
           t3*t4 - 
           12.*am*(amk*amk)*(aml*aml)*C1*C4*(C5*C5)*C6*t3*t4 - 
           6.*(amk*amk)*(aml*aml)*amLam*C1*C4*(C5*C5)*C6*t3*t4 - 
           6.*(amk*amk)*(aml*aml)*amSig*C1*C4*(C5*C5)*C6*t3*t4 + 
           12.*am*(amk*amk)*(aml*aml)*C1*(C4*C4)*(C5*C5)*C6*t3*
           t4 + 6.*(amk*amk)*(aml*aml)*amLam*C1*(C4*C4)*(C5*C5)*
           C6*t3*t4 + 
           6.*(amk*amk)*(aml*aml)*amSig*C1*(C4*C4)*(C5*C5)*C6*t3*
           t4 + 24.*am*(amk*amk)*apkk1*C1*(C4*C4)*(C5*C5)*C6*
           t3*t4 + 
           12.*(amk*amk)*amLam*apkk1*C1*(C4*C4)*(C5*C5)*C6*t3*
           t4 + 12.*(amk*amk)*amSig*apkk1*C1*(C4*C4)*(C5*C5)*
           C6*t3*t4 - 
           6.*am*(amk*amk)*(aml*aml*aml*aml)*C1*C3*(C4*C4)*(C5*C5)*C6*t3*
           t4 - 3.*(amk*amk)*(aml*aml*aml*aml)*amLam*C1*C3*(C4*C4)*
           (C5*C5)*C6*t3*t4 - 
           3.*(amk*amk)*(aml*aml*aml*aml)*amSig*C1*C3*(C4*C4)*(C5*C5)*C6*
           t3*t4 - 
           12.*am*(amk*amk)*(aml*aml)*apkk1*C1*C3*(C4*C4)*(C5*C5)*
           C6*t3*t4 - 
           6.*(amk*amk)*(aml*aml)*amLam*apkk1*C1*C3*(C4*C4)*
           (C5*C5)*C6*t3*t4 - 
           6.*(amk*amk)*(aml*aml)*amSig*apkk1*C1*C3*(C4*C4)*
           (C5*C5)*C6*t3*t4 - 
           6.*am*(amk*amk)*(aml*aml*aml*aml)*C1*(C4*C4)*(C5*C5*C5)*C6*t3*t4 - 
           3.*(amk*amk)*(aml*aml*aml*aml)*amLam*C1*(C4*C4)*(C5*C5*C5)*C6*t3*
           t4 - 3.*(amk*amk)*(aml*aml*aml*aml)*amSig*C1*(C4*C4)*(C5*C5*C5)*
           C6*t3*t4 - 
           12.*am*(amk*amk)*(aml*aml)*apkk1*C1*(C4*C4)*(C5*C5*C5)*C6*
           t3*t4 - 
           6.*(amk*amk)*(aml*aml)*amLam*apkk1*C1*(C4*C4)*(C5*C5*C5)*
           C6*t3*t4 - 
           6.*(amk*amk)*(aml*aml)*amSig*apkk1*C1*(C4*C4)*(C5*C5*C5)*
           C6*t3*t4 + 
           36.*(am*am)*(aml*aml)*C1*C4*C5*C6*Ekaon*t3*t4 - 
           36.*(aml*aml)*amLam*amSig*C1*C4*C5*C6*Ekaon*t3*
           t4 + 72.*(am*am)*apkk1*C1*C4*C5*C6*Ekaon*t3*
           t4 - 72.*amLam*amSig*apkk1*C1*C4*C5*C6*Ekaon*
           t3*t4 - 
           36.*(am*am)*(aml*aml)*C1*(C4*C4)*C5*C6*Ekaon*t3*t4 + 
           36.*(aml*aml)*amLam*amSig*C1*(C4*C4)*C5*C6*Ekaon*
           t3*t4 + 
           12.*(am*am)*(aml*aml)*C1*C4*(C5*C5)*C6*Ekaon*t3*t4 - 
           12.*(aml*aml)*amLam*amSig*C1*C4*(C5*C5)*C6*Ekaon*
           t3*t4 - 
           12.*(am*am)*(aml*aml)*C1*(C4*C4)*(C5*C5)*C6*Ekaon*t3*
           t4 + 12.*(aml*aml)*amLam*amSig*C1*(C4*C4)*(C5*C5)*
           C6*Ekaon*t3*t4 - 
           24.*(am*am)*apkk1*C1*(C4*C4)*(C5*C5)*C6*Ekaon*t3*
           t4 + 24.*amLam*amSig*apkk1*C1*(C4*C4)*(C5*C5)*C6*
           Ekaon*t3*t4 + 
           6.*(am*am)*(aml*aml*aml*aml)*C1*C3*(C4*C4)*(C5*C5)*C6*Ekaon*t3*
           t4 - 6.*(aml*aml*aml*aml)*amLam*amSig*C1*C3*(C4*C4)*(C5*C5)*
           C6*Ekaon*t3*t4 + 
           12.*(am*am)*(aml*aml)*apkk1*C1*C3*(C4*C4)*(C5*C5)*C6*
           Ekaon*t3*t4 - 
           12.*(aml*aml)*amLam*amSig*apkk1*C1*C3*(C4*C4)*
           (C5*C5)*C6*Ekaon*t3*t4 + 
           6.*(am*am)*(aml*aml*aml*aml)*C1*(C4*C4)*(C5*C5*C5)*C6*Ekaon*t3*
           t4 - 6.*(aml*aml*aml*aml)*amLam*amSig*C1*(C4*C4)*(C5*C5*C5)*C6*
           Ekaon*t3*t4 + 
           12.*(am*am)*(aml*aml)*apkk1*C1*(C4*C4)*(C5*C5*C5)*C6*
           Ekaon*t3*t4 - 
           12.*(aml*aml)*amLam*amSig*apkk1*C1*(C4*C4)*(C5*C5*C5)*
           C6*Ekaon*t3*t4 - 
           18.*(am*am)*(amk*amk)*(aml*aml)*C1*C3*(C4*C4)*(C5*C5)*C6*
           Enu*t3*t4 + 
           6.*(amk*amk*amk*amk)*(aml*aml)*C1*C3*(C4*C4)*(C5*C5)*C6*Enu*t3*
           t4 - 12.*am*(amk*amk)*(aml*aml)*amLam*C1*C3*(C4*C4)*
           (C5*C5)*C6*Enu*t3*t4 - 
           12.*am*(amk*amk)*(aml*aml)*amSig*C1*C3*(C4*C4)*(C5*C5)*
           C6*Enu*t3*t4 - 
           6.*(amk*amk)*(aml*aml)*amLam*amSig*C1*C3*(C4*C4)*
           (C5*C5)*C6*Enu*t3*t4 - 
           18.*(am*am)*(amk*amk)*(aml*aml)*C1*(C4*C4)*(C5*C5*C5)*C6*Enu*
           t3*t4 + 
           6.*(amk*amk*amk*amk)*(aml*aml)*C1*(C4*C4)*(C5*C5*C5)*C6*Enu*t3*
           t4 - 12.*am*(amk*amk)*(aml*aml)*amLam*C1*(C4*C4)*
           (C5*C5*C5)*C6*Enu*t3*t4 - 
           12.*am*(amk*amk)*(aml*aml)*amSig*C1*(C4*C4)*(C5*C5*C5)*C6*
           Enu*t3*t4 - 
           6.*(amk*amk)*(aml*aml)*amLam*amSig*C1*(C4*C4)*(C5*C5*C5)*
           C6*Enu*t3*t4 + 
           12.*(am*am*am)*(aml*aml)*C1*C3*(C4*C4)*(C5*C5)*C6*Ekaon*
           Enu*t3*t4 - 
           24.*am*(amk*amk)*(aml*aml)*C1*C3*(C4*C4)*(C5*C5)*C6*
           Ekaon*Enu*t3*t4 - 
           12.*am*(aml*aml)*amLam*amSig*C1*C3*(C4*C4)*(C5*C5)*
           C6*Ekaon*Enu*t3*t4 + 
           12.*(am*am*am)*(aml*aml)*C1*(C4*C4)*(C5*C5*C5)*C6*Ekaon*Enu*
           t3*t4 - 
           24.*am*(amk*amk)*(aml*aml)*C1*(C4*C4)*(C5*C5*C5)*C6*Ekaon*
           Enu*t3*t4 - 
           12.*am*(aml*aml)*amLam*amSig*C1*(C4*C4)*(C5*C5*C5)*C6*
           Ekaon*Enu*t3*t4 + 
           24.*(am*am)*(aml*aml)*C1*C3*(C4*C4)*(C5*C5)*C6*
           (Ekaon*Ekaon)*Enu*t3*t4 + 
           24.*(am*am)*(aml*aml)*C1*(C4*C4)*(C5*C5*C5)*C6*(Ekaon*Ekaon)*
           Enu*t3*t4 + 
           9.*(am*am)*(amk*amk)*(aml*aml)*C1*C4*C5*C6*Fm1*t3*
           t4 - 9.*(amk*amk*amk*amk)*(aml*aml)*C1*C4*C5*C6*Fm1*t3*
           t4 - 9.*am*(amk*amk)*(aml*aml)*amLam*C1*C4*C5*C6*
           Fm1*t3*t4 + 
           9.*am*(amk*amk)*(aml*aml)*amSig*C1*C4*C5*C6*Fm1*t3*
           t4 - 9.*(amk*amk)*(aml*aml)*amLam*amSig*C1*C4*C5*
           C6*Fm1*t3*t4 + 
           18.*(am*am)*(amk*amk)*(aml*aml)*C1*C4*(C5*C5)*C6*Fm1*t3*
           t4 + 6.*(amk*amk*amk*amk)*(aml*aml)*C1*C4*(C5*C5)*C6*Fm1*t3*
           t4 + 12.*am*(amk*amk)*(aml*aml)*amLam*C1*C4*(C5*C5)*
           C6*Fm1*t3*t4 + 
           12.*am*(amk*amk)*(aml*aml)*amSig*C1*C4*(C5*C5)*C6*Fm1*
           t3*t4 + 
           6.*(amk*amk)*(aml*aml)*amLam*amSig*C1*C4*(C5*C5)*C6*
           Fm1*t3*t4 - 
           18.*(am*am*am)*(aml*aml)*C1*C4*C5*C6*Ekaon*Fm1*t3*
           t4 + 18.*am*(amk*amk)*(aml*aml)*C1*C4*C5*C6*Ekaon*
           Fm1*t3*t4 + 
           18.*(am*am)*(aml*aml)*amLam*C1*C4*C5*C6*Ekaon*Fm1*
           t3*t4 - 
           18.*(am*am)*(aml*aml)*amSig*C1*C4*C5*C6*Ekaon*Fm1*
           t3*t4 - 
           18.*(amk*amk)*(aml*aml)*amSig*C1*C4*C5*C6*Ekaon*Fm1*
           t3*t4 + 
           18.*am*(aml*aml)*amLam*amSig*C1*C4*C5*C6*Ekaon*
           Fm1*t3*t4 - 
           12.*(am*am*am)*(aml*aml)*C1*C4*(C5*C5)*C6*Ekaon*Fm1*t3*
           t4 - 12.*am*(amk*amk)*(aml*aml)*C1*C4*(C5*C5)*C6*
           Ekaon*Fm1*t3*t4 + 
           6.*(amk*amk)*(aml*aml)*amLam*C1*C4*(C5*C5)*C6*Ekaon*
           Fm1*t3*t4 + 
           6.*(amk*amk)*(aml*aml)*amSig*C1*C4*(C5*C5)*C6*Ekaon*
           Fm1*t3*t4 + 
           12.*am*(aml*aml)*amLam*amSig*C1*C4*(C5*C5)*C6*
           Ekaon*Fm1*t3*t4 + 
           36.*am*(aml*aml)*amSig*C1*C4*C5*C6*(Ekaon*Ekaon)*Fm1*
           t3*t4 - 
           12.*am*(aml*aml)*amLam*C1*C4*(C5*C5)*C6*(Ekaon*Ekaon)*
           Fm1*t3*t4 - 
           12.*am*(aml*aml)*amSig*C1*C4*(C5*C5)*C6*(Ekaon*Ekaon)*
           Fm1*t3*t4 + 
           18.*(amk*amk)*amLam*apkk1*C1*C4*C5*C6*Enu*Fm1*t3*
           t4 - 18.*(amk*amk)*amSig*apkk1*C1*C4*C5*C6*Enu*
           Fm1*t3*t4 - 
           6.*(amk*amk)*amLam*apkk1*C1*C4*(C5*C5)*C6*Enu*Fm1*
           t3*t4 + 
           6.*(amk*amk)*amSig*apkk1*C1*C4*(C5*C5)*C6*Enu*Fm1*
           t3*t4 - 
           36.*am*amLam*apkk1*C1*C4*C5*C6*Ekaon*Enu*Fm1*
           t3*t4 + 
           36.*am*amSig*apkk1*C1*C4*C5*C6*Ekaon*Enu*Fm1*
           t3*t4 + 
           12.*am*amLam*apkk1*C1*C4*(C5*C5)*C6*Ekaon*Enu*
           Fm1*t3*t4 - 
           12.*am*amSig*apkk1*C1*C4*(C5*C5)*C6*Ekaon*Enu*
           Fm1*t3*t4 + 
           3.*(am*am)*(amk*amk)*(aml*aml)*C1*C4*C5*C6*Fm2*t3*
           t4 - 3.*(amk*amk*amk*amk)*(aml*aml)*C1*C4*C5*C6*Fm2*t3*
           t4 + 3.*am*(amk*amk)*(aml*aml)*amLam*C1*C4*C5*C6*
           Fm2*t3*t4 - 
           3.*am*(amk*amk)*(aml*aml)*amSig*C1*C4*C5*C6*Fm2*t3*
           t4 - 3.*(amk*amk)*(aml*aml)*amLam*amSig*C1*C4*C5*
           C6*Fm2*t3*t4 - 
           18.*(am*am)*(amk*amk)*(aml*aml)*C1*(C4*C4)*C5*C6*Fm2*t3*
           t4 - 6.*(amk*amk*amk*amk)*(aml*aml)*C1*(C4*C4)*C5*C6*Fm2*t3*
           t4 - 12.*am*(amk*amk)*(aml*aml)*amLam*C1*(C4*C4)*C5*
           C6*Fm2*t3*t4 - 
           12.*am*(amk*amk)*(aml*aml)*amSig*C1*(C4*C4)*C5*C6*Fm2*
           t3*t4 - 
           6.*(amk*amk)*(aml*aml)*amLam*amSig*C1*(C4*C4)*C5*C6*
           Fm2*t3*t4 - 
           6.*(am*am*am)*(aml*aml)*C1*C4*C5*C6*Ekaon*Fm2*t3*t4 + 
           6.*am*(amk*amk)*(aml*aml)*C1*C4*C5*C6*Ekaon*Fm2*t3*
           t4 - 6.*(am*am)*(aml*aml)*amLam*C1*C4*C5*C6*Ekaon*
           Fm2*t3*t4 - 
           6.*(amk*amk)*(aml*aml)*amLam*C1*C4*C5*C6*Ekaon*Fm2*
           t3*t4 + 
           6.*(am*am)*(aml*aml)*amSig*C1*C4*C5*C6*Ekaon*Fm2*
           t3*t4 + 
           6.*am*(aml*aml)*amLam*amSig*C1*C4*C5*C6*Ekaon*
           Fm2*t3*t4 + 
           12.*(am*am*am)*(aml*aml)*C1*(C4*C4)*C5*C6*Ekaon*Fm2*t3*
           t4 + 12.*am*(amk*amk)*(aml*aml)*C1*(C4*C4)*C5*C6*
           Ekaon*Fm2*t3*t4 - 
           6.*(amk*amk)*(aml*aml)*amLam*C1*(C4*C4)*C5*C6*Ekaon*
           Fm2*t3*t4 - 
           6.*(amk*amk)*(aml*aml)*amSig*C1*(C4*C4)*C5*C6*Ekaon*
           Fm2*t3*t4 - 
           12.*am*(aml*aml)*amLam*amSig*C1*(C4*C4)*C5*C6*
           Ekaon*Fm2*t3*t4 + 
           12.*am*(aml*aml)*amLam*C1*C4*C5*C6*(Ekaon*Ekaon)*Fm2*
           t3*t4 + 
           12.*am*(aml*aml)*amLam*C1*(C4*C4)*C5*C6*(Ekaon*Ekaon)*
           Fm2*t3*t4 + 
           12.*am*(aml*aml)*amSig*C1*(C4*C4)*C5*C6*(Ekaon*Ekaon)*
           Fm2*t3*t4 - 
           6.*(amk*amk)*amLam*apkk1*C1*C4*C5*C6*Enu*Fm2*t3*
           t4 + 6.*(amk*amk)*amSig*apkk1*C1*C4*C5*C6*Enu*
           Fm2*t3*t4 - 
           6.*(amk*amk)*amLam*apkk1*C1*(C4*C4)*C5*C6*Enu*Fm2*
           t3*t4 + 
           6.*(amk*amk)*amSig*apkk1*C1*(C4*C4)*C5*C6*Enu*Fm2*
           t3*t4 + 
           12.*am*amLam*apkk1*C1*C4*C5*C6*Ekaon*Enu*Fm2*
           t3*t4 - 
           12.*am*amSig*apkk1*C1*C4*C5*C6*Ekaon*Enu*Fm2*
           t3*t4 + 
           12.*am*amLam*apkk1*C1*(C4*C4)*C5*C6*Ekaon*Enu*
           Fm2*t3*t4 - 
           12.*am*amSig*apkk1*C1*(C4*C4)*C5*C6*Ekaon*Enu*
           Fm2*t3*t4 - 
           36.*(am*am)*(amk*amk)*(aml*aml)*C1*C4*C5*C6*Enu*Fm1*
           Fm2*t3*t4 + 
           12.*(amk*amk*amk*amk)*(aml*aml)*C1*C4*C5*C6*Enu*Fm1*Fm2*t3*
           t4 - 24.*am*(amk*amk)*(aml*aml)*amLam*C1*C4*C5*C6*
           Enu*Fm1*Fm2*t3*t4 - 
           24.*am*(amk*amk)*(aml*aml)*amSig*C1*C4*C5*C6*Enu*
           Fm1*Fm2*t3*t4 - 
           12.*(amk*amk)*(aml*aml)*amLam*amSig*C1*C4*C5*C6*Enu*
           Fm1*Fm2*t3*t4 + 
           24.*(am*am*am)*(aml*aml)*C1*C4*C5*C6*Ekaon*Enu*Fm1*
           Fm2*t3*t4 - 
           48.*am*(amk*amk)*(aml*aml)*C1*C4*C5*C6*Ekaon*Enu*
           Fm1*Fm2*t3*t4 - 
           24.*am*(aml*aml)*amLam*amSig*C1*C4*C5*C6*Ekaon*
           Enu*Fm1*Fm2*t3*t4 + 
           48.*(am*am)*(aml*aml)*C1*C4*C5*C6*(Ekaon*Ekaon)*Enu*Fm1*
           Fm2*t3*t4 - 
           18.*am*(amk*amk)*(aml*aml)*(C5*C5)*(C6*C6)*(t4*t4) - 
           18.*(amk*amk)*(aml*aml)*amLam*(C5*C5)*(C6*C6)*(t4*t4) - 
           36.*am*(amk*amk)*apkk1*(C5*C5)*(C6*C6)*(t4*t4) - 
           36.*(amk*amk)*amLam*apkk1*(C5*C5)*(C6*C6)*(t4*t4) - 
           12.*am*(amk*amk)*(aml*aml)*(C5*C5*C5)*(C6*C6)*(t4*t4) - 
           12.*(amk*amk)*(aml*aml)*amLam*(C5*C5*C5)*(C6*C6)*(t4*t4) - 
           2.*am*(amk*amk)*(aml*aml)*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) - 
           2.*(amk*amk)*(aml*aml)*amLam*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) - 
           4.*am*(amk*amk)*apkk1*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) - 
           4.*(amk*amk)*amLam*apkk1*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) + 
           2.*am*(amk*amk)*(aml*aml*aml*aml)*(C5*C5*C5*C5*C5)*(C6*C6)*(t4*t4) + 
           2.*(amk*amk)*(aml*aml*aml*aml)*amLam*(C5*C5*C5*C5*C5)*(C6*C6)*(t4*t4) + 
           4.*am*(amk*amk)*(aml*aml)*apkk1*(C5*C5*C5*C5*C5)*(C6*C6)*(t4*t4) + 
           4.*(amk*amk)*(aml*aml)*amLam*apkk1*(C5*C5*C5*C5*C5)*(C6*C6)*
           (t4*t4) + 
           18.*(am*am)*(aml*aml)*(C5*C5)*(C6*C6)*Ekaon*(t4*t4) - 
           18.*(aml*aml)*(amLam*amLam)*(C5*C5)*(C6*C6)*Ekaon*(t4*t4) + 
           36.*(am*am)*apkk1*(C5*C5)*(C6*C6)*Ekaon*(t4*t4) - 
           36.*(amLam*amLam)*apkk1*(C5*C5)*(C6*C6)*Ekaon*(t4*t4) + 
           12.*(am*am)*(aml*aml)*(C5*C5*C5)*(C6*C6)*Ekaon*(t4*t4) - 
           12.*(aml*aml)*(amLam*amLam)*(C5*C5*C5)*(C6*C6)*Ekaon*(t4*t4) + 
           2.*(am*am)*(aml*aml)*(C5*C5*C5*C5)*(C6*C6)*Ekaon*(t4*t4) - 
           2.*(aml*aml)*(amLam*amLam)*(C5*C5*C5*C5)*(C6*C6)*Ekaon*(t4*t4) + 
           4.*(am*am)*apkk1*(C5*C5*C5*C5)*(C6*C6)*Ekaon*(t4*t4) - 
           4.*(amLam*amLam)*apkk1*(C5*C5*C5*C5)*(C6*C6)*Ekaon*(t4*t4) - 
           2.*(am*am)*(aml*aml*aml*aml)*(C5*C5*C5*C5*C5)*(C6*C6)*Ekaon*(t4*t4) + 
           2.*(aml*aml*aml*aml)*(amLam*amLam)*(C5*C5*C5*C5*C5)*(C6*C6)*Ekaon*(t4*t4) - 
           4.*(am*am)*(aml*aml)*apkk1*(C5*C5*C5*C5*C5)*(C6*C6)*Ekaon*
           (t4*t4) + 
           4.*(aml*aml)*(amLam*amLam)*apkk1*(C5*C5*C5*C5*C5)*(C6*C6)*Ekaon*
           (t4*t4) + 
           6.*(am*am)*(amk*amk)*(aml*aml)*(C5*C5*C5*C5*C5)*(C6*C6)*Enu*
           (t4*t4) - 
           2.*(amk*amk*amk*amk)*(aml*aml)*(C5*C5*C5*C5*C5)*(C6*C6)*Enu*(t4*t4) + 
           8.*am*(amk*amk)*(aml*aml)*amLam*(C5*C5*C5*C5*C5)*(C6*C6)*Enu*
           (t4*t4) + 
           2.*(amk*amk)*(aml*aml)*(amLam*amLam)*(C5*C5*C5*C5*C5)*(C6*C6)*Enu*
           (t4*t4) - 
           4.*(am*am*am)*(aml*aml)*(C5*C5*C5*C5*C5)*(C6*C6)*Ekaon*Enu*(t4*t4) + 
           8.*am*(amk*amk)*(aml*aml)*(C5*C5*C5*C5*C5)*(C6*C6)*Ekaon*Enu*
           (t4*t4) + 
           4.*am*(aml*aml)*(amLam*amLam)*(C5*C5*C5*C5*C5)*(C6*C6)*Ekaon*Enu*
           (t4*t4) - 
           8.*(am*am)*(aml*aml)*(C5*C5*C5*C5*C5)*(C6*C6)*(Ekaon*Ekaon)*Enu*
           (t4*t4) + 
           3.*(am*am)*(amk*amk)*(aml*aml)*(C5*C5)*(C6*C6)*Fm2*
           (t4*t4) - 
           3.*(amk*amk*amk*amk)*(aml*aml)*(C5*C5)*(C6*C6)*Fm2*(t4*t4) - 
           3.*(amk*amk)*(aml*aml)*(amLam*amLam)*(C5*C5)*(C6*C6)*Fm2*
           (t4*t4) + 
           6.*(am*am)*(amk*amk)*(aml*aml)*(C5*C5*C5)*(C6*C6)*Fm2*
           (t4*t4) + 
           2.*(amk*amk*amk*amk)*(aml*aml)*(C5*C5*C5)*(C6*C6)*Fm2*(t4*t4) + 
           8.*am*(amk*amk)*(aml*aml)*amLam*(C5*C5*C5)*(C6*C6)*Fm2*
           (t4*t4) + 
           2.*(amk*amk)*(aml*aml)*(amLam*amLam)*(C5*C5*C5)*(C6*C6)*Fm2*
           (t4*t4) - 
           6.*(am*am*am)*(aml*aml)*(C5*C5)*(C6*C6)*Ekaon*Fm2*(t4*t4) + 
           6.*am*(amk*amk)*(aml*aml)*(C5*C5)*(C6*C6)*Ekaon*Fm2*
           (t4*t4) - 
           6.*(amk*amk)*(aml*aml)*amLam*(C5*C5)*(C6*C6)*Ekaon*Fm2*
           (t4*t4) + 
           6.*am*(aml*aml)*(amLam*amLam)*(C5*C5)*(C6*C6)*Ekaon*Fm2*
           (t4*t4) - 
           4.*(am*am*am)*(aml*aml)*(C5*C5*C5)*(C6*C6)*Ekaon*Fm2*(t4*t4) - 
           4.*am*(amk*amk)*(aml*aml)*(C5*C5*C5)*(C6*C6)*Ekaon*Fm2*
           (t4*t4) + 
           4.*(amk*amk)*(aml*aml)*amLam*(C5*C5*C5)*(C6*C6)*Ekaon*Fm2*
           (t4*t4) + 
           4.*am*(aml*aml)*(amLam*amLam)*(C5*C5*C5)*(C6*C6)*Ekaon*Fm2*
           (t4*t4) + 
           12.*am*(aml*aml)*amLam*(C5*C5)*(C6*C6)*(Ekaon*Ekaon)*Fm2*
           (t4*t4) - 
           8.*am*(aml*aml)*amLam*(C5*C5*C5)*(C6*C6)*(Ekaon*Ekaon)*Fm2*
           (t4*t4) - 
           6.*(am*am)*(amk*amk)*(aml*aml)*(C5*C5)*(C6*C6)*Enu*(Fm2*Fm2)*
           (t4*t4) + 
           2.*(amk*amk*amk*amk)*(aml*aml)*(C5*C5)*(C6*C6)*Enu*(Fm2*Fm2)*
           (t4*t4) - 
           8.*am*(amk*amk)*(aml*aml)*amLam*(C5*C5)*(C6*C6)*Enu*
           (Fm2*Fm2)*(t4*t4) - 
           2.*(amk*amk)*(aml*aml)*(amLam*amLam)*(C5*C5)*(C6*C6)*Enu*
           (Fm2*Fm2)*(t4*t4) + 
           4.*(am*am*am)*(aml*aml)*(C5*C5)*(C6*C6)*Ekaon*Enu*(Fm2*Fm2)*
           (t4*t4) - 
           8.*am*(amk*amk)*(aml*aml)*(C5*C5)*(C6*C6)*Ekaon*Enu*
           (Fm2*Fm2)*(t4*t4) - 
           4.*am*(aml*aml)*(amLam*amLam)*(C5*C5)*(C6*C6)*Ekaon*Enu*
           (Fm2*Fm2)*(t4*t4) + 
           8.*(am*am)*(aml*aml)*(C5*C5)*(C6*C6)*(Ekaon*Ekaon)*Enu*
           (Fm2*Fm2)*(t4*t4) + 
           6.*am*(Elep*Elep)*
           (amLam*C5*C6*
           ((amk*amk)*C1*C4*(3. + C5)*Fm1*t3 + 
           Fm2*
           (2.*(1. + C2)*t1 + 
           (amk*amk)*C1*(-1. + C4)*C4*t3))*t4 + 
           amSig*C1*C4*t3*
           (-1.*(amk*amk)*(-1. + C4)*C5*C6*Fm2*t4 + 
           Fm1*
           (6.*(1. + C2)*t1 - 
           1.*(amk*amk)*C5*(3. + C5)*C6*t4)) + 
           2.*am*
           ((1. + C2)*C5*C6*Fm2*t1*t4 + 
           C1*C4*t3*
           (-1.*(amLam - 1.*amSig)*(-1. + C4)*C5*
           C6*Ekaon*Fm2*t4 + 
           Fm1*
           (3.*(1. + C2)*t1 - 
           1.*(amLam - 1.*amSig)*C5*(3. + C5)*
           C6*Ekaon*t4)))) + 
           108.*(aml*aml)*(C2*C2)*C7*t1*t5 + 
           144.*apkk1*(C2*C2)*C7*t1*t5 - 
           54.*(amk*amk)*(aml*aml)*C1*C2*C4*C7*t3*t5 - 
           72.*(amk*amk)*apkk1*C1*C2*C4*C7*t3*t5 + 
           36.*am*(aml*aml)*C1*C2*C4*C7*Ekaon*t3*t5 - 
           72.*(aml*aml)*amSig*C1*C2*C4*C7*Ekaon*t3*t5 - 
           144.*amSig*apkk1*C1*C2*C4*C7*Ekaon*t3*t5 + 
           36.*am*(aml*aml)*C1*C2*(C4*C4)*C7*Ekaon*t3*t5 + 
           36.*(aml*aml)*amSig*C1*C2*(C4*C4)*C7*Ekaon*t3*t5 + 
           36.*am*apkk1*C1*C2*C4*C7*Enu*t3*t5 + 
           36.*amSig*apkk1*C1*C2*C4*C7*Enu*t3*t5 - 
           36.*am*apkk1*C1*C2*(C4*C4)*C7*Enu*t3*t5 - 
           36.*amSig*apkk1*C1*C2*(C4*C4)*C7*Enu*t3*t5 + 
           18.*am*(aml*aml*aml*aml)*C1*C2*C4*C7*Fm1*t3*t5 + 
           18.*(aml*aml*aml*aml)*amSig*C1*C2*C4*C7*Fm1*t3*t5 + 
           36.*am*(aml*aml)*apkk1*C1*C2*C4*C7*Fm1*t3*t5 + 
           36.*(aml*aml)*amSig*apkk1*C1*C2*C4*C7*Fm1*t3*
           t5 + 36.*(amk*amk)*(aml*aml)*C1*C2*C4*C7*Enu*Fm1*
           t3*t5 + 
           36.*(amk*amk)*apkk1*C1*C2*C4*C7*Enu*Fm1*t3*t5 - 
           72.*am*(aml*aml)*C1*C2*C4*C7*Ekaon*Enu*Fm1*t3*
           t5 - 72.*am*apkk1*C1*C2*C4*C7*Ekaon*Enu*Fm1*
           t3*t5 - 
           54.*(amk*amk)*(aml*aml)*C2*C5*C6*C7*t4*t5 - 
           72.*(amk*amk)*apkk1*C2*C5*C6*C7*t4*t5 + 
           36.*am*(aml*aml)*C2*C5*C6*C7*Ekaon*t4*t5 - 
           72.*(aml*aml)*amLam*C2*C5*C6*C7*Ekaon*t4*t5 - 
           144.*amLam*apkk1*C2*C5*C6*C7*Ekaon*t4*t5 - 
           12.*am*(aml*aml)*C2*(C5*C5)*C6*C7*Ekaon*t4*t5 - 
           12.*(aml*aml)*amLam*C2*(C5*C5)*C6*C7*Ekaon*t4*t5 + 
           36.*am*apkk1*C2*C5*C6*C7*Enu*t4*t5 + 
           36.*amLam*apkk1*C2*C5*C6*C7*Enu*t4*t5 + 
           12.*am*apkk1*C2*(C5*C5)*C6*C7*Enu*t4*t5 + 
           12.*amLam*apkk1*C2*(C5*C5)*C6*C7*Enu*t4*t5 + 
           6.*am*(aml*aml*aml*aml)*C2*C5*C6*C7*Fm2*t4*t5 + 
           6.*(aml*aml*aml*aml)*amLam*C2*C5*C6*C7*Fm2*t4*t5 + 
           12.*am*(aml*aml)*apkk1*C2*C5*C6*C7*Fm2*t4*t5 + 
           12.*(aml*aml)*amLam*apkk1*C2*C5*C6*C7*Fm2*t4*
           t5 + 12.*(amk*amk)*(aml*aml)*C2*C5*C6*C7*Enu*Fm2*
           t4*t5 + 
           12.*(amk*amk)*apkk1*C2*C5*C6*C7*Enu*Fm2*t4*t5 - 
           24.*am*(aml*aml)*C2*C5*C6*C7*Ekaon*Enu*Fm2*t4*
           t5 - 24.*am*apkk1*C2*C5*C6*C7*Ekaon*Enu*Fm2*
           t4*t5 - 72.*(aml*aml)*(C2*C2)*(C7*C7)*Ekaon*(t5*t5) - 
           144.*apkk1*(C2*C2)*(C7*C7)*Ekaon*(t5*t5) + 
           72.*(aml*aml)*(C2*C2)*(C7*C7)*Enu*(t5*t5) + 
           144.*apkk1*(C2*C2)*(C7*C7)*Enu*(t5*t5) - 
           1.*Elep*(36.*((-1. + C2)*(-1. + C2))*(t1*t1) - 
           18.*(amk*amk)*(aml*aml)*C1*C3*(C4*C4)*t2*t3 + 
           9.*(amk*amk*amk*amk)*(C1*C1)*(C4*C4)*(t3*t3) - 
           9.*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4)*(t3*t3) + 
           18.*(amk*amk*amk*amk)*(C1*C1)*(C4*C4*C4)*(t3*t3) - 
           18.*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4*C4)*(t3*t3) + 
           9.*(amk*amk*amk*amk)*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) - 
           9.*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) - 
           9.*(amk*amk*amk*amk)*(aml*aml)*(C1*C1)*C3*(C4*C4*C4*C4)*(t3*t3) + 
           9.*(amk*amk)*(aml*aml)*(amSig*amSig)*(C1*C1)*C3*(C4*C4*C4*C4)*
           (t3*t3) + 
           6.*(amk*amk)*(aml*aml)*C3*(C5*C5)*C6*t2*t4 + 
           18.*(amk*amk*amk*amk)*C1*C4*C5*C6*t3*t4 - 
           18.*(amk*amk)*amLam*amSig*C1*C4*C5*C6*t3*t4 + 
           18.*(amk*amk*amk*amk)*C1*(C4*C4)*C5*C6*t3*t4 - 
           18.*(amk*amk)*amLam*amSig*C1*(C4*C4)*C5*C6*t3*
           t4 - 6.*(amk*amk*amk*amk)*C1*C4*(C5*C5)*C6*t3*t4 + 
           6.*(amk*amk)*amLam*amSig*C1*C4*(C5*C5)*C6*t3*
           t4 - 6.*(amk*amk*amk*amk)*C1*(C4*C4)*(C5*C5)*C6*t3*t4 + 
           6.*(amk*amk)*amLam*amSig*C1*(C4*C4)*(C5*C5)*C6*t3*
           t4 + 
           3.*(amk*amk*amk*amk)*(aml*aml)*C1*C3*(C4*C4)*(C5*C5)*C6*t3*
           t4 - 
           3.*(amk*amk)*(aml*aml)*amLam*amSig*C1*C3*(C4*C4)*
           (C5*C5)*C6*t3*t4 + 
           3.*(amk*amk*amk*amk)*(aml*aml)*C1*(C4*C4)*(C5*C5*C5)*C6*t3*t4 - 
           3.*(amk*amk)*(aml*aml)*amLam*amSig*C1*(C4*C4)*
           (C5*C5*C5)*C6*t3*t4 + 
           9.*(amk*amk)*(aml*aml)*amLam*C1*C4*C5*C6*Fm1*t3*
           t4 - 
           9.*(amk*amk)*(aml*aml)*amSig*C1*C4*C5*C6*Fm1*t3*
           t4 + 
           18.*(amk*amk)*amLam*apkk1*C1*C4*C5*C6*Fm1*t3*
           t4 - 
           18.*(amk*amk)*amSig*apkk1*C1*C4*C5*C6*Fm1*t3*
           t4 + 
           6.*(amk*amk)*(aml*aml)*amLam*C1*C4*(C5*C5)*C6*Fm1*
           t3*t4 - 
           6.*(amk*amk)*(aml*aml)*amSig*C1*C4*(C5*C5)*C6*Fm1*
           t3*t4 + 
           6.*(amk*amk)*amLam*apkk1*C1*C4*(C5*C5)*C6*Fm1*
           t3*t4 - 
           6.*(amk*amk)*amSig*apkk1*C1*C4*(C5*C5)*C6*Fm1*
           t3*t4 - 
           3.*(amk*amk)*(aml*aml)*amLam*C1*C4*C5*C6*Fm2*t3*
           t4 + 
           3.*(amk*amk)*(aml*aml)*amSig*C1*C4*C5*C6*Fm2*t3*
           t4 - 
           6.*(amk*amk)*amLam*apkk1*C1*C4*C5*C6*Fm2*t3*
           t4 + 
           6.*(amk*amk)*amSig*apkk1*C1*C4*C5*C6*Fm2*t3*
           t4 + 
           6.*(amk*amk)*(aml*aml)*amLam*C1*(C4*C4)*C5*C6*Fm2*
           t3*t4 - 
           6.*(amk*amk)*(aml*aml)*amSig*C1*(C4*C4)*C5*C6*Fm2*
           t3*t4 + 
           6.*(amk*amk)*amLam*apkk1*C1*(C4*C4)*C5*C6*Fm2*
           t3*t4 - 
           6.*(amk*amk)*amSig*apkk1*C1*(C4*C4)*C5*C6*Fm2*
           t3*t4 + 9.*(amk*amk*amk*amk)*(C5*C5)*(C6*C6)*(t4*t4) - 
           9.*(amk*amk)*(amLam*amLam)*(C5*C5)*(C6*C6)*(t4*t4) - 
           6.*(amk*amk*amk*amk)*(C5*C5*C5)*(C6*C6)*(t4*t4) + 
           6.*(amk*amk)*(amLam*amLam)*(C5*C5*C5)*(C6*C6)*(t4*t4) + 
           (amk*amk*amk*amk)*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) - 
           1.*(amk*amk)*(amLam*amLam)*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) - 
           1.*(amk*amk*amk*amk)*(aml*aml)*(C5*C5*C5*C5*C5)*(C6*C6)*(t4*t4) + 
           (amk*amk)*(aml*aml)*(amLam*amLam)*(C5*C5*C5*C5*C5)*(C6*C6)*(t4*t4) - 
           2.*(am*am*am)*Ekaon*
           (9.*(C1*C1)*(C4*C4)*
           (-1. + 2.*C4 + 
           (-1. + (aml*aml)*C3)*(C4*C4))*(t3*t3) - 
           3.*C1*C4*C5*
           (2.*(3. + C5) + 
           C4*
           (-6. - 2.*C5 + (aml*aml)*C3*C5 + 
           (aml*aml)*(C5*C5)))*C6*t3*t4 + 
           (C5*C5)*
           (-9. - 6.*C5 - 1.*(C5*C5) + 
           (aml*aml)*(C5*C5*C5))*(C6*C6)*(t4*t4)) - 
           6.*t1*
           ((aml*aml)*
           (-6.*amSig*C1*C4*Fm1*t3 - 
           3.*amSig*C1*C2*C4*Fm1*t3 + 
           3.*C3*
           (2.*t2 + 
           C1*(C4*C4)*((amk*amk) - 2.*am*Ekaon)*t3)
           - 1.*(amk*amk)*(C5*C5*C5)*C6*t4 - 
           2.*amLam*C5*C6*Fm2*t4 - 
           1.*amLam*C2*C5*C6*Fm2*t4 - 
           1.*am*
           (3.*C1*(2. + C2)*C4*Fm1*t3 + 
           C5*C6*
           (-2.*(C5*C5)*Ekaon + (2. + C2)*Fm2)*
           t4)) + 
           2.*
           ((amk*amk)*(-1. + C2)*
           (3.*C1*C4*(1. + C4)*t3 - 
           1.*(-3. + C5)*C5*C6*t4) - 
           1.*apkk1*(1. + C2)*
           (3.*amSig*C1*C4*Fm1*t3 + 
           amLam*C5*C6*Fm2*t4) + 
           (am*am)*
           (3.*C1*C4*
           (C4 - 1.*Enu*Fm1 + 
           C2*(-1. + 2.*C4 + Enu*Fm1))*t3 - 
           1.*C5*C6*
           (C5 + Enu*Fm2 + 
           C2*(3. + 2.*C5 - 1.*Enu*Fm2))*t4)
           + am*
           (-3.*C1*C4*
           (2.*(-1. + C2)*(1. + C4)*Ekaon + 
           apkk1*(1. + C2)*Fm1)*t3 + 
           3.*amSig*C1*C4*
           (C4 - 1.*Enu*Fm1 + 
           C2*(-1. + 2.*C4 + Enu*Fm1))*t3 - 
           1.*C5*C6*
           (-2.*(-1. + C2)*(-3. + C5)*Ekaon + 
           apkk1*(1. + C2)*Fm2 + 
           amLam*
           (C5 + Enu*Fm2 + 
           C2*(3. + 2.*C5 - 1.*Enu*Fm2)))*t4))
           ) - 
           1.*(am*am)*
           (3.*(amk*amk)*
           (3.*(C1*C1)*(C4*C4)*
           (3. - 2.*C4 + 3.*(C4*C4))*(t3*t3) + 
           2.*C1*C4*C5*
           (9. + C5 - 3.*C4*(1. + C5))*C6*t3*t4
           + (C5*C5)*(9. + 2.*C5 + (C5*C5))*(C6*C6)*
           (t4*t4)) - 
           4.*Ekaon*
           (9.*(C1*C1)*(C4*C4)*((1.+C4)*(1.+C4))*Ekaon*
           (t3*t3) - 
           3.*C1*C4*C5*C6*
           (2.*(1. + C4)*(-3. + C5)*Ekaon + 
           amSig*
           (3.*C4 + C5 + 3.*Enu*Fm1 - 
           1.*C5*Enu*Fm1 - 1.*Enu*Fm2 - 
           1.*C4*Enu*Fm2) + 
           amLam*
           (-1.*C5 - 3.*Enu*Fm1 + 
           C5*Enu*Fm1 + Enu*Fm2 + 
           C4*(-3. + Enu*Fm2)))*t3*t4 + 
           ((-3. + C5)*(-3. + C5))*(C5*C5)*(C6*C6)*Ekaon*(t4*t4)
           ) + 
           (aml*aml)*
           (-9.*(C1*C1)*C3*(C4*C4*C4*C4)*
           (3.*(amk*amk) - 4.*(Ekaon*Ekaon))*(t3*t3) + 
           (C5*C5)*C6*t4*
           (6.*C3*t2 + 
           (C5*C5*C5)*C6*
           (-3.*(amk*amk) + 4.*(Ekaon*Ekaon))*t4) - 
           3.*C1*(C4*C4)*t3*
           ((C5*C5*C5)*C6*
           (-3.*(amk*amk) + 
           2.*Ekaon*
           (-1.*amLam + amSig + 2.*Ekaon))*t4
           + C3*
           (6.*t2 + 
           (C5*C5)*C6*
           (-3.*(amk*amk) + 
           2.*Ekaon*
           (amLam - 1.*amSig + 2.*Ekaon))*t4))
           )) - 
           36.*(amk*amk)*amSig*C1*C2*C4*C7*t3*t5 + 
           18.*(aml*aml)*amSig*C1*C2*C4*C7*t3*t5 + 
           36.*amSig*apkk1*C1*C2*C4*C7*t3*t5 + 
           36.*amSig*apkk1*C1*C2*(C4*C4)*C7*t3*t5 + 
           18.*(amk*amk)*(aml*aml)*C1*C2*C4*C7*Fm1*t3*t5 + 
           36.*(amk*amk)*apkk1*C1*C2*C4*C7*Fm1*t3*t5 - 
           36.*(amk*amk)*amLam*C2*C5*C6*C7*t4*t5 + 
           18.*(aml*aml)*amLam*C2*C5*C6*C7*t4*t5 + 
           36.*amLam*apkk1*C2*C5*C6*C7*t4*t5 - 
           12.*amLam*apkk1*C2*(C5*C5)*C6*C7*t4*t5 + 
           6.*(amk*amk)*(aml*aml)*C2*C5*C6*C7*Fm2*t4*t5 + 
           12.*(amk*amk)*apkk1*C2*C5*C6*C7*Fm2*t4*t5 + 
           72.*(aml*aml)*(C2*C2)*(C7*C7)*(t5*t5) + 
           144.*apkk1*(C2*C2)*(C7*C7)*(t5*t5) - 
           1.*am*
           ((aml*aml)*
           (-36.*(amk*amk)*(C1*C1)*C3*(C4*C4*C4*C4)*Ekaon*
           (t3*t3) - 
           18.*(amSig*amSig)*(C1*C1)*C3*(C4*C4*C4*C4)*Ekaon*
           (t3*t3) - 
           3.*amSig*C1*C4*t3*
           (C5*C6*
           (-3.*(amk*amk)*C4*(C5*C5) - 
           2.*Ekaon*
           (amLam*C4*(C5*C5) - 3.*Fm1 - 
           2.*C5*Fm1 + Fm2 - 2.*C4*Fm2))*t4 + 
           C3*C4*
           (6.*t2 - 
           2.*amLam*(C5*C5)*C6*Ekaon*t4 + 
           (amk*amk)*
           (12.*C1*(C4*C4)*t3 - 1.*(C5*C5)*C6*t4))
           ) + 
           3.*C1*C4*t3*
           ((amk*amk)*C4*(C5*C5*C5)*C6*
           (amLam + 4.*Ekaon)*t4 + 
           6.*amLam*C5*C6*Ekaon*Fm1*t4 + 
           4.*amLam*(C5*C5)*C6*Ekaon*Fm1*t4 - 
           2.*amLam*C5*C6*Ekaon*Fm2*t4 + 
           4.*amLam*C4*C5*C6*Ekaon*Fm2*t4 + 
           C3*C4*
           (3.*(amk*amk)*amLam*(C5*C5)*C6*t4 + 
           4.*Ekaon*
           (-3.*t2 + (amk*amk)*(C5*C5)*C6*t4)) - 
           6.*C2*C7*t5 + 
           12.*C2*C7*Ekaon*Fm1*t5) - 
           2.*C5*C6*t4*
           (-3.*amLam*C3*C5*t2 - 
           6.*C3*C5*Ekaon*t2 + 
           2.*(amk*amk)*amLam*(C5*C5*C5*C5)*C6*t4 + 
           2.*(amk*amk)*(C5*C5*C5*C5)*C6*Ekaon*t4 + 
           (amLam*amLam)*(C5*C5*C5*C5)*C6*Ekaon*t4 + 
           9.*C2*C7*t5 - 6.*C2*C7*Ekaon*Fm2*t5
           )) + 
           2.*
           (9.*(amSig*amSig)*(C1*C1)*(C4*C4)*Ekaon*(t3*t3) - 
           18.*(amSig*amSig)*(C1*C1)*(C4*C4*C4)*Ekaon*
           (t3*t3) + 
           9.*(amSig*amSig)*(C1*C1)*(C4*C4*C4*C4)*Ekaon*
           (t3*t3) + 
           18.*amLam*amSig*C1*C4*C5*C6*Ekaon*t3*
           t4 - 
           18.*amLam*amSig*C1*(C4*C4)*C5*C6*Ekaon*
           t3*t4 + 
           6.*amLam*amSig*C1*C4*(C5*C5)*C6*Ekaon*
           t3*t4 - 
           6.*amLam*amSig*C1*(C4*C4)*(C5*C5)*C6*
           Ekaon*t3*t4 + 
           18.*amLam*apkk1*C1*C4*C5*C6*Ekaon*
           Fm1*t3*t4 - 
           18.*amSig*apkk1*C1*C4*C5*C6*Ekaon*
           Fm1*t3*t4 + 
           6.*amLam*apkk1*C1*C4*(C5*C5)*C6*Ekaon*
           Fm1*t3*t4 - 
           6.*amSig*apkk1*C1*C4*(C5*C5)*C6*Ekaon*
           Fm1*t3*t4 - 
           6.*amLam*apkk1*C1*C4*C5*C6*Ekaon*Fm2*
           t3*t4 + 
           6.*amSig*apkk1*C1*C4*C5*C6*Ekaon*Fm2*
           t3*t4 + 
           6.*amLam*apkk1*C1*(C4*C4)*C5*C6*Ekaon*
           Fm2*t3*t4 - 
           6.*amSig*apkk1*C1*(C4*C4)*C5*C6*Ekaon*
           Fm2*t3*t4 + 
           9.*(amLam*amLam)*(C5*C5)*(C6*C6)*Ekaon*
           (t4*t4) + 
           6.*(amLam*amLam)*(C5*C5*C5)*(C6*C6)*Ekaon*
           (t4*t4) + 
           (amLam*amLam)*(C5*C5*C5*C5)*(C6*C6)*Ekaon*(t4*t4) - 
           18.*apkk1*C1*C2*C4*C7*t3*t5 - 
           18.*apkk1*C1*C2*(C4*C4)*C7*t3*t5 + 
           36.*apkk1*C1*C2*C4*C7*Ekaon*Fm1*t3*
           t5 - 18.*apkk1*C2*C5*C6*C7*t4*t5 + 
           6.*apkk1*C2*(C5*C5)*C6*C7*t4*t5 + 
           12.*apkk1*C2*C5*C6*C7*Ekaon*Fm2*t4*
           t5 + 
           (amk*amk)*
           (18.*(C1*C1)*(C4*C4)*((1.+C4)*(1.+C4))*Ekaon*
           (t3*t3) + 
           3.*amSig*C1*C4*t3*
           (6.*C1*C4*(1. + (C4*C4))*t3 + 
           C5*C6*
           (6. - 1.*C5 - 3.*Enu*Fm1 + 
           C5*Enu*Fm1 + Enu*Fm2 + 
           C4*(-3. - 2.*C5 + Enu*Fm2))*t4) + 
           2.*C5*C6*t4*
           (amLam*C5*(9. + (C5*C5))*C6*t4 + 
           9.*C5*C6*Ekaon*t4 - 
           6.*(C5*C5)*C6*Ekaon*t4 + 
           (C5*C5*C5)*C6*Ekaon*t4 + 9.*C2*C7*t5) - 
           3.*C1*C4*t3*
           (amLam*C5*C6*
           (-6. - 1.*C5 - 3.*Enu*Fm1 + 
           C5*Enu*Fm1 + Enu*Fm2 + 
           C4*(-3. + 2.*C5 + Enu*Fm2))*t4 + 
           2.*
           (-6.*(1. + C4)*C5*C6*Ekaon*t4 + 
           2.*(1. + C4)*(C5*C5)*C6*Ekaon*t4 - 
           3.*C2*C7*t5))))))))));
  
  if (sol <= 0.0) {
    LOG( "SKXSec", pERROR )
      << "Matrix element is 0 or negative! This is bad and really should throw an exception!";
  }  
  
  return sol;
}


//____________________________________________________________________________
double AlamSimoAtharVacasSKPXSec2014::Amatrix_PP(double theta, double phikq) const
{
  double sol = 0.;
  
  double akk1=0., zdotq=0., qdotpk=0., akcrosk1=0., qcrospk=0.;
  double zdotpk=0., azpk=0., aqkaon=0., akpk=0., apkk1=0.;
  double C1=0., C2=0., C3=0., C4=0., C5=0., C6=0., C7=0., C8=0., C9=0.;
  double aq2=0., gform=0., con=0., t1=0., t2=0., t3=0., t4=0., t5=0., t6=0.;
  
  akk1=Enu*Elep-Enu*alepvec*cos(theta);
  zdotq=(Enu*Enu-alepvec*alepvec)/2.0;
  qdotpk=aqvec*pkvec*angkq;
  akcrosk1 = Enu*alepvec*sin(theta);
  qcrospk=aqvec*pkvec*sqrt(1.0-angkq*angkq);
  zdotpk=(akcrosk1*qcrospk*cos(phikq)+zdotq*qdotpk)/(aqvec*aqvec);
  azpk=Ekaon*(Enu+Elep)/2.0-zdotpk;
  aqkaon=aq0*Ekaon-aqvec*pkvec*angkq;
  akpk=azpk + aqkaon/2.0;
  apkk1=azpk - aqkaon/2.0;                                          
  C1=1.0/(am*am+amk*amk-2.0*am*Ekaon-amSig*amSig);
  C2=d+f;
  C3=1./(aml*aml-2.0*akk1-amk*amk);
  C4=d-f;
  C5=d+3.*f;
  C6=1.0/(am*am+amk*amk-2.0*am*Ekaon-amLam*amLam);       
  C7=2.0*am/(aml*aml-2.0*akk1+amk*amk-2.*aqkaon-ampi*ampi);
  C8=2.0*am/(aml*aml-2.0*akk1+amk*amk-2.*aqkaon-amEta*amEta);
  C9=d - 3.*f;
  aq2=aml*aml-2.0*akk1;              
  gform=1.0/pow(1.0-aq2/(1.1*1.1),4);
      
  con=g*g*Vus*Vus/(4.0*fpi*fpi);      
      
  t1=1.0;
  t2=1.0;               // !Full Term
  t3=1.0;
  t4=1.0;
  t5=1.0;
  t6=1.0;

  sol  =   con*gform*(0.4444444444444444*am*
           (akk1*(-288.*Elep*(t1*t1) + 288.*Enu*(t1*t1) + 
           576.*Elep*f*(t1*t1) + 576.*Enu*f*(t1*t1) - 
           288.*Elep*(f*f)*(t1*t1) + 288.*Enu*(f*f)*(t1*t1) - 
           288.*(aml*aml)*C3*Ekaon*t1*t2 - 
           288.*(aml*aml)*C3*Enu*t1*t2 + 
           72.*akpk*(aml*aml)*(C3*C3)*Ekaon*(t2*t2) - 
           36.*(amk*amk)*(aml*aml)*(C3*C3)*Ekaon*(t2*t2) + 
           108.*(aml*aml*aml*aml)*(C3*C3)*Ekaon*(t2*t2) - 
           72.*(aml*aml)*apkk1*(C3*C3)*Ekaon*(t2*t2) + 
           72.*akpk*(aml*aml)*(C3*C3)*Elep*(t2*t2) + 
           108.*(amk*amk)*(aml*aml)*(C3*C3)*Elep*(t2*t2) - 
           36.*(aml*aml*aml*aml)*(C3*C3)*Elep*(t2*t2) - 
           72.*(aml*aml)*apkk1*(C3*C3)*Elep*(t2*t2) - 
           72.*akpk*(aml*aml)*(C3*C3)*Enu*(t2*t2) - 
           108.*(amk*amk)*(aml*aml)*(C3*C3)*Enu*(t2*t2) + 
           36.*(aml*aml*aml*aml)*(C3*C3)*Enu*(t2*t2) + 
           72.*(aml*aml)*apkk1*(C3*C3)*Enu*(t2*t2) + 
           144.*akpk*amSig*C1*C4*t1*t3 + 
           144.*amSig*apkk1*C1*C4*t1*t3 - 
           144.*akpk*amSig*C1*(C4*C4)*t1*t3 + 
           144.*amSig*apkk1*C1*(C4*C4)*t1*t3 + 
           144.*akpk*(aml*aml)*amSig*C1*C3*(C4*C4)*t1*t3 + 
           72.*(amk*amk)*(aml*aml)*amSig*C1*C3*(C4*C4)*t1*t3 + 
           72.*(amk*amk)*(aml*aml)*C1*C3*(C4*C4)*Ekaon*t1*t3 + 
           144.*(amk*amk)*C1*C4*Elep*t1*t3 + 
           144.*(amk*amk)*C1*(C4*C4)*Elep*t1*t3 + 
           144.*(amk*amk)*C1*C4*Enu*t1*t3 - 
           144.*(amk*amk)*C1*(C4*C4)*Enu*t1*t3 + 
           144.*(amk*amk)*(aml*aml)*C1*C3*(C4*C4)*Enu*t1*t3 + 
           144.*akpk*amSig*C1*C4*f*t1*t3 - 
           144.*amSig*apkk1*C1*C4*f*t1*t3 - 
           144.*akpk*amSig*C1*(C4*C4)*f*t1*t3 - 
           144.*amSig*apkk1*C1*(C4*C4)*f*t1*t3 - 
           144.*(amk*amk)*C1*C4*Elep*f*t1*t3 - 
           144.*(amk*amk)*C1*(C4*C4)*Elep*f*t1*t3 + 
           144.*(amk*amk)*C1*C4*Enu*f*t1*t3 - 
           144.*(amk*amk)*C1*(C4*C4)*Enu*f*t1*t3 - 
           144.*akpk*(amk*amk)*C1*C4*Fm1*t1*t3 - 
           144.*(amk*amk)*apkk1*C1*C4*Fm1*t1*t3 - 
           144.*akpk*amSig*C1*C4*Ekaon*Fm1*t1*t3 - 
           144.*amSig*apkk1*C1*C4*Ekaon*Fm1*t1*t3 - 
           288.*akpk*amSig*C1*C4*Elep*Fm1*t1*t3 + 
           288.*amSig*apkk1*C1*C4*Enu*Fm1*t1*t3 + 
           144.*akpk*(amk*amk)*C1*C4*f*Fm1*t1*t3 - 
           216.*(amk*amk)*(aml*aml)*C1*C4*f*Fm1*t1*t3 - 
           144.*(amk*amk)*apkk1*C1*C4*f*Fm1*t1*t3 + 
           144.*akpk*amSig*C1*C4*Ekaon*f*Fm1*t1*t3 - 
           216.*(aml*aml)*amSig*C1*C4*Ekaon*f*Fm1*t1*t3 - 
           144.*amSig*apkk1*C1*C4*Ekaon*f*Fm1*t1*t3 + 
           72.*akpk*(aml*aml)*amSig*C1*C3*(C4*C4)*t2*t3 + 
           72.*(amk*amk)*(aml*aml)*amSig*C1*C3*(C4*C4)*t2*t3 + 
           72.*(akpk*akpk)*(aml*aml)*amSig*C1*(C3*C3)*(C4*C4)*t2*t3 + 
           36.*akpk*(amk*amk)*(aml*aml)*amSig*C1*(C3*C3)*(C4*C4)*t2*
           t3 - 36.*akpk*(aml*aml*aml*aml)*amSig*C1*(C3*C3)*(C4*C4)*t2*
           t3 - 72.*(amk*amk)*(aml*aml*aml*aml)*amSig*C1*(C3*C3)*(C4*C4)*t2*
           t3 - 144.*akpk*(aml*aml)*amSig*apkk1*C1*(C3*C3)*
           (C4*C4)*t2*t3 - 
           36.*(amk*amk)*(aml*aml)*amSig*apkk1*C1*(C3*C3)*(C4*C4)*t2*
           t3 + 36.*(aml*aml*aml*aml)*amSig*apkk1*C1*(C3*C3)*(C4*C4)*t2*
           t3 + 72.*(aml*aml)*amSig*(apkk1*apkk1)*C1*(C3*C3)*(C4*C4)*
           t2*t3 + 72.*(amk*amk)*(aml*aml)*C1*C3*(C4*C4)*Ekaon*t2*
           t3 - 72.*(amk*amk)*(aml*aml*aml*aml)*C1*(C3*C3)*(C4*C4)*Ekaon*t2*
           t3 - 72.*akpk*(amk*amk)*(aml*aml)*C1*(C3*C3)*(C4*C4)*Elep*
           t2*t3 - 36.*(amk*amk*amk*amk)*(aml*aml)*C1*(C3*C3)*(C4*C4)*Elep*
           t2*t3 + 36.*(amk*amk)*(aml*aml*aml*aml)*C1*(C3*C3)*(C4*C4)*Elep*
           t2*t3 + 72.*(amk*amk)*(aml*aml)*apkk1*C1*(C3*C3)*(C4*C4)*
           Elep*t2*t3 + 
           72.*(amk*amk)*(aml*aml)*C1*C3*(C4*C4)*Enu*t2*t3 + 
           72.*akpk*(amk*amk)*(aml*aml)*C1*(C3*C3)*(C4*C4)*Enu*t2*
           t3 + 36.*(amk*amk*amk*amk)*(aml*aml)*C1*(C3*C3)*(C4*C4)*Enu*t2*
           t3 - 36.*(amk*amk)*(aml*aml*aml*aml)*C1*(C3*C3)*(C4*C4)*Enu*t2*
           t3 - 72.*(amk*amk)*(aml*aml)*apkk1*C1*(C3*C3)*(C4*C4)*Enu*
           t2*t3 + 36.*akpk*(amk*amk)*amSig*(C1*C1)*(C4*C4)*
           (t3*t3) - 36.*(amk*amk)*amSig*apkk1*(C1*C1)*(C4*C4)*
           (t3*t3) - 72.*akpk*(amk*amk)*amSig*(C1*C1)*(C4*C4*C4)*
           (t3*t3) - 72.*(amk*amk)*amSig*apkk1*(C1*C1)*(C4*C4*C4)*
           (t3*t3) + 36.*akpk*(amk*amk)*amSig*(C1*C1)*(C4*C4*C4*C4)*
           (t3*t3) - 36.*(amk*amk)*amSig*apkk1*(C1*C1)*(C4*C4*C4*C4)*
           (t3*t3) - 72.*akpk*(amk*amk)*(aml*aml)*amSig*(C1*C1)*C3*
           (C4*C4*C4*C4)*(t3*t3) - 
           36.*(amk*amk*amk*amk)*(aml*aml)*amSig*(C1*C1)*C3*(C4*C4*C4*C4)*(t3*t3) - 
           36.*(akpk*akpk)*(amk*amk)*(aml*aml)*amSig*(C1*C1)*(C3*C3)*
           (C4*C4*C4*C4)*(t3*t3) + 
           18.*akpk*(amk*amk)*(aml*aml*aml*aml)*amSig*(C1*C1)*(C3*C3)*(C4*C4*C4*C4)*
           (t3*t3) + 18.*(amk*amk*amk*amk)*(aml*aml*aml*aml)*amSig*(C1*C1)*(C3*C3)*
           (C4*C4*C4*C4)*(t3*t3) + 
           72.*akpk*(amk*amk)*(aml*aml)*amSig*apkk1*(C1*C1)*(C3*C3)*
           (C4*C4*C4*C4)*(t3*t3) - 
           18.*(amk*amk)*(aml*aml*aml*aml)*amSig*apkk1*(C1*C1)*(C3*C3)*(C4*C4*C4*C4)*
           (t3*t3) - 36.*(amk*amk)*(aml*aml)*amSig*(apkk1*apkk1)*(C1*C1)*
           (C3*C3)*(C4*C4*C4*C4)*(t3*t3) + 
           36.*akpk*(amSig*amSig)*(C1*C1)*(C4*C4)*Ekaon*(t3*t3) - 
           36.*(amSig*amSig)*apkk1*(C1*C1)*(C4*C4)*Ekaon*(t3*t3) - 
           72.*akpk*(amSig*amSig)*(C1*C1)*(C4*C4*C4)*Ekaon*(t3*t3) - 
           72.*(amSig*amSig)*apkk1*(C1*C1)*(C4*C4*C4)*Ekaon*(t3*t3) + 
           36.*akpk*(amSig*amSig)*(C1*C1)*(C4*C4*C4*C4)*Ekaon*(t3*t3) - 
           36.*(amSig*amSig)*apkk1*(C1*C1)*(C4*C4*C4*C4)*Ekaon*(t3*t3) - 
           18.*(amk*amk*amk*amk)*(aml*aml)*(C1*C1)*C3*(C4*C4*C4*C4)*Ekaon*(t3*t3) - 
           72.*akpk*(aml*aml)*(amSig*amSig)*(C1*C1)*C3*(C4*C4*C4*C4)*Ekaon*
           (t3*t3) - 18.*(amk*amk)*(aml*aml)*(amSig*amSig)*(C1*C1)*C3*
           (C4*C4*C4*C4)*Ekaon*(t3*t3) + 
           9.*(amk*amk*amk*amk)*(aml*aml*aml*aml)*(C1*C1)*(C3*C3)*(C4*C4*C4*C4)*Ekaon*(t3*t3) - 
           36.*(akpk*akpk)*(aml*aml)*(amSig*amSig)*(C1*C1)*(C3*C3)*(C4*C4*C4*C4)*
           Ekaon*(t3*t3) + 
           18.*akpk*(aml*aml*aml*aml)*(amSig*amSig)*(C1*C1)*(C3*C3)*(C4*C4*C4*C4)*Ekaon*
           (t3*t3) + 9.*(amk*amk)*(aml*aml*aml*aml)*(amSig*amSig)*(C1*C1)*(C3*C3)*
           (C4*C4*C4*C4)*Ekaon*(t3*t3) + 
           72.*akpk*(aml*aml)*(amSig*amSig)*apkk1*(C1*C1)*(C3*C3)*(C4*C4*C4*C4)*
           Ekaon*(t3*t3) - 
           18.*(aml*aml*aml*aml)*(amSig*amSig)*apkk1*(C1*C1)*(C3*C3)*(C4*C4*C4*C4)*
           Ekaon*(t3*t3) - 
           36.*(aml*aml)*(amSig*amSig)*(apkk1*apkk1)*(C1*C1)*(C3*C3)*(C4*C4*C4*C4)*
           Ekaon*(t3*t3) - 
           18.*(amk*amk*amk*amk)*(C1*C1)*(C4*C4)*Elep*(t3*t3) + 
           18.*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4)*Elep*(t3*t3) - 
           36.*(amk*amk*amk*amk)*(C1*C1)*(C4*C4*C4)*Elep*(t3*t3) + 
           36.*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4*C4)*Elep*(t3*t3) - 
           18.*(amk*amk*amk*amk)*(C1*C1)*(C4*C4*C4*C4)*Elep*(t3*t3) + 
           18.*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4*C4*C4)*Elep*(t3*t3) + 
           18.*akpk*(amk*amk*amk*amk)*(aml*aml)*(C1*C1)*(C3*C3)*(C4*C4*C4*C4)*Elep*
           (t3*t3) - 9.*(amk*amk*amk*amk)*(aml*aml*aml*aml)*(C1*C1)*(C3*C3)*(C4*C4*C4*C4)*Elep*
           (t3*t3) - 18.*akpk*(amk*amk)*(aml*aml)*(amSig*amSig)*(C1*C1)*
           (C3*C3)*(C4*C4*C4*C4)*Elep*(t3*t3) + 
           9.*(amk*amk)*(aml*aml*aml*aml)*(amSig*amSig)*(C1*C1)*(C3*C3)*(C4*C4*C4*C4)*Elep*
           (t3*t3) - 18.*(amk*amk*amk*amk)*(aml*aml)*apkk1*(C1*C1)*(C3*C3)*
           (C4*C4*C4*C4)*Elep*(t3*t3) + 
           18.*(amk*amk)*(aml*aml)*(amSig*amSig)*apkk1*(C1*C1)*(C3*C3)*
           (C4*C4*C4*C4)*Elep*(t3*t3) + 
           18.*(amk*amk*amk*amk)*(C1*C1)*(C4*C4)*Enu*(t3*t3) - 
           18.*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4)*Enu*(t3*t3) - 
           36.*(amk*amk*amk*amk)*(C1*C1)*(C4*C4*C4)*Enu*(t3*t3) + 
           36.*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4*C4)*Enu*(t3*t3) + 
           18.*(amk*amk*amk*amk)*(C1*C1)*(C4*C4*C4*C4)*Enu*(t3*t3) - 
           18.*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4*C4*C4)*Enu*(t3*t3) - 
           36.*(amk*amk*amk*amk)*(aml*aml)*(C1*C1)*C3*(C4*C4*C4*C4)*Enu*(t3*t3) + 
           36.*(amk*amk)*(aml*aml)*(amSig*amSig)*(C1*C1)*C3*(C4*C4*C4*C4)*Enu*
           (t3*t3) - 18.*akpk*(amk*amk*amk*amk)*(aml*aml)*(C1*C1)*(C3*C3)*
           (C4*C4*C4*C4)*Enu*(t3*t3) + 
           9.*(amk*amk*amk*amk)*(aml*aml*aml*aml)*(C1*C1)*(C3*C3)*(C4*C4*C4*C4)*Enu*(t3*t3) + 
           18.*akpk*(amk*amk)*(aml*aml)*(amSig*amSig)*(C1*C1)*(C3*C3)*
           (C4*C4*C4*C4)*Enu*(t3*t3) - 
           9.*(amk*amk)*(aml*aml*aml*aml)*(amSig*amSig)*(C1*C1)*(C3*C3)*(C4*C4*C4*C4)*Enu*
           (t3*t3) + 18.*(amk*amk*amk*amk)*(aml*aml)*apkk1*(C1*C1)*(C3*C3)*
           (C4*C4*C4*C4)*Enu*(t3*t3) - 
           18.*(amk*amk)*(aml*aml)*(amSig*amSig)*apkk1*(C1*C1)*(C3*C3)*
           (C4*C4*C4*C4)*Enu*(t3*t3) + 
           36.*akpk*(amk*amk*amk*amk)*(C1*C1)*(C4*C4)*Fm1*(t3*t3) - 
           54.*(amk*amk*amk*amk)*(aml*aml)*(C1*C1)*(C4*C4)*Fm1*(t3*t3) + 
           36.*akpk*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4)*Fm1*(t3*t3) - 
           54.*(amk*amk)*(aml*aml)*(amSig*amSig)*(C1*C1)*(C4*C4)*Fm1*
           (t3*t3) - 36.*(amk*amk*amk*amk)*apkk1*(C1*C1)*(C4*C4)*Fm1*
           (t3*t3) - 36.*(amk*amk)*(amSig*amSig)*apkk1*(C1*C1)*(C4*C4)*
           Fm1*(t3*t3) + 
           36.*akpk*(amk*amk*amk*amk)*(C1*C1)*(C4*C4*C4)*Fm1*(t3*t3) + 
           36.*akpk*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4*C4)*Fm1*(t3*t3) + 
           36.*(amk*amk*amk*amk)*apkk1*(C1*C1)*(C4*C4*C4)*Fm1*(t3*t3) + 
           36.*(amk*amk)*(amSig*amSig)*apkk1*(C1*C1)*(C4*C4*C4)*Fm1*
           (t3*t3) + 72.*akpk*(amk*amk)*amSig*(C1*C1)*(C4*C4)*Ekaon*
           Fm1*(t3*t3) - 
           108.*(amk*amk)*(aml*aml)*amSig*(C1*C1)*(C4*C4)*Ekaon*Fm1*
           (t3*t3) - 72.*(amk*amk)*amSig*apkk1*(C1*C1)*(C4*C4)*
           Ekaon*Fm1*(t3*t3) + 
           72.*akpk*(amk*amk)*amSig*(C1*C1)*(C4*C4*C4)*Ekaon*Fm1*
           (t3*t3) + 72.*(amk*amk)*amSig*apkk1*(C1*C1)*(C4*C4*C4)*
           Ekaon*Fm1*(t3*t3) - 
           72.*(akpk*akpk)*(amk*amk)*amSig*(C1*C1)*(C4*C4)*(Fm1*Fm1)*
           (t3*t3) + 54.*akpk*(amk*amk)*(aml*aml)*amSig*(C1*C1)*
           (C4*C4)*(Fm1*Fm1)*(t3*t3) + 
           18.*(amk*amk*amk*amk)*(aml*aml)*amSig*(C1*C1)*(C4*C4)*(Fm1*Fm1)*
           (t3*t3) - 54.*(amk*amk)*(aml*aml)*amSig*apkk1*(C1*C1)*
           (C4*C4)*(Fm1*Fm1)*(t3*t3) - 
           72.*(amk*amk)*amSig*(apkk1*apkk1)*(C1*C1)*(C4*C4)*(Fm1*Fm1)*
           (t3*t3) + 9.*(amk*amk*amk*amk)*(aml*aml)*(C1*C1)*(C4*C4)*Ekaon*
           (Fm1*Fm1)*(t3*t3) - 
           72.*(akpk*akpk)*(amSig*amSig)*(C1*C1)*(C4*C4)*Ekaon*(Fm1*Fm1)*
           (t3*t3) + 54.*akpk*(aml*aml)*(amSig*amSig)*(C1*C1)*(C4*C4)*
           Ekaon*(Fm1*Fm1)*(t3*t3) + 
           9.*(amk*amk)*(aml*aml)*(amSig*amSig)*(C1*C1)*(C4*C4)*Ekaon*
           (Fm1*Fm1)*(t3*t3) - 
           54.*(aml*aml)*(amSig*amSig)*apkk1*(C1*C1)*(C4*C4)*Ekaon*
           (Fm1*Fm1)*(t3*t3) - 
           72.*(amSig*amSig)*(apkk1*apkk1)*(C1*C1)*(C4*C4)*Ekaon*(Fm1*Fm1)*
           (t3*t3) - 27.*(amk*amk*amk*amk)*(aml*aml)*(C1*C1)*(C4*C4)*Elep*
           (Fm1*Fm1)*(t3*t3) + 
           27.*(amk*amk)*(aml*aml)*(amSig*amSig)*(C1*C1)*(C4*C4)*Elep*
           (Fm1*Fm1)*(t3*t3) - 
           36.*(amk*amk*amk*amk)*apkk1*(C1*C1)*(C4*C4)*Elep*(Fm1*Fm1)*(t3*t3) + 
           36.*(amk*amk)*(amSig*amSig)*apkk1*(C1*C1)*(C4*C4)*Elep*
           (Fm1*Fm1)*(t3*t3) - 
           36.*akpk*(amk*amk*amk*amk)*(C1*C1)*(C4*C4)*Enu*(Fm1*Fm1)*(t3*t3) + 
           27.*(amk*amk*amk*amk)*(aml*aml)*(C1*C1)*(C4*C4)*Enu*(Fm1*Fm1)*(t3*t3) + 
           36.*akpk*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4)*Enu*(Fm1*Fm1)*
           (t3*t3) - 27.*(amk*amk)*(aml*aml)*(amSig*amSig)*(C1*C1)*(C4*C4)*
           Enu*(Fm1*Fm1)*(t3*t3) - 
           144.*akpk*amLam*C5*C6*t1*t4 - 
           144.*amLam*apkk1*C5*C6*t1*t4 - 
           48.*akpk*amLam*(C5*C5)*C6*t1*t4 + 
           48.*amLam*apkk1*(C5*C5)*C6*t1*t4 + 
           48.*akpk*(aml*aml)*amLam*C3*(C5*C5)*C6*t1*t4 + 
           24.*(amk*amk)*(aml*aml)*amLam*C3*(C5*C5)*C6*t1*t4 + 
           24.*(amk*amk)*(aml*aml)*C3*(C5*C5)*C6*Ekaon*t1*t4 - 
           144.*(amk*amk)*C5*C6*Elep*t1*t4 + 
           48.*(amk*amk)*(C5*C5)*C6*Elep*t1*t4 - 
           144.*(amk*amk)*C5*C6*Enu*t1*t4 - 
           48.*(amk*amk)*(C5*C5)*C6*Enu*t1*t4 + 
           48.*(amk*amk)*(aml*aml)*C3*(C5*C5)*C6*Enu*t1*t4 - 
           144.*akpk*amLam*C5*C6*f*t1*t4 + 
           144.*amLam*apkk1*C5*C6*f*t1*t4 - 
           48.*akpk*amLam*(C5*C5)*C6*f*t1*t4 - 
           48.*amLam*apkk1*(C5*C5)*C6*f*t1*t4 + 
           144.*(amk*amk)*C5*C6*Elep*f*t1*t4 - 
           48.*(amk*amk)*(C5*C5)*C6*Elep*f*t1*t4 - 
           144.*(amk*amk)*C5*C6*Enu*f*t1*t4 - 
           48.*(amk*amk)*(C5*C5)*C6*Enu*f*t1*t4 + 
           48.*akpk*(amk*amk)*C5*C6*Fm2*t1*t4 + 
           48.*(amk*amk)*apkk1*C5*C6*Fm2*t1*t4 + 
           48.*akpk*amLam*C5*C6*Ekaon*Fm2*t1*t4 + 
           48.*amLam*apkk1*C5*C6*Ekaon*Fm2*t1*t4 + 
           96.*akpk*amLam*C5*C6*Elep*Fm2*t1*t4 - 
           96.*amLam*apkk1*C5*C6*Enu*Fm2*t1*t4 - 
           48.*akpk*(amk*amk)*C5*C6*f*Fm2*t1*t4 + 
           72.*(amk*amk)*(aml*aml)*C5*C6*f*Fm2*t1*t4 + 
           48.*(amk*amk)*apkk1*C5*C6*f*Fm2*t1*t4 - 
           48.*akpk*amLam*C5*C6*Ekaon*f*Fm2*t1*t4 + 
           72.*(aml*aml)*amLam*C5*C6*Ekaon*f*Fm2*t1*t4 + 
           48.*amLam*apkk1*C5*C6*Ekaon*f*Fm2*t1*t4 + 
           24.*akpk*(aml*aml)*amLam*C3*(C5*C5)*C6*t2*t4 + 
           24.*(amk*amk)*(aml*aml)*amLam*C3*(C5*C5)*C6*t2*t4 + 
           24.*(akpk*akpk)*(aml*aml)*amLam*(C3*C3)*(C5*C5)*C6*t2*t4 + 
           12.*akpk*(amk*amk)*(aml*aml)*amLam*(C3*C3)*(C5*C5)*C6*t2*
           t4 - 12.*akpk*(aml*aml*aml*aml)*amLam*(C3*C3)*(C5*C5)*C6*t2*
           t4 - 24.*(amk*amk)*(aml*aml*aml*aml)*amLam*(C3*C3)*(C5*C5)*C6*t2*
           t4 - 48.*akpk*(aml*aml)*amLam*apkk1*(C3*C3)*(C5*C5)*C6*
           t2*t4 - 12.*(amk*amk)*(aml*aml)*amLam*apkk1*(C3*C3)*
           (C5*C5)*C6*t2*t4 + 
           12.*(aml*aml*aml*aml)*amLam*apkk1*(C3*C3)*(C5*C5)*C6*t2*t4 + 
           24.*(aml*aml)*amLam*(apkk1*apkk1)*(C3*C3)*(C5*C5)*C6*t2*t4 + 
           24.*(amk*amk)*(aml*aml)*C3*(C5*C5)*C6*Ekaon*t2*t4 - 
           24.*(amk*amk)*(aml*aml*aml*aml)*(C3*C3)*(C5*C5)*C6*Ekaon*t2*t4 - 
           24.*akpk*(amk*amk)*(aml*aml)*(C3*C3)*(C5*C5)*C6*Elep*t2*
           t4 - 12.*(amk*amk*amk*amk)*(aml*aml)*(C3*C3)*(C5*C5)*C6*Elep*t2*
           t4 + 12.*(amk*amk)*(aml*aml*aml*aml)*(C3*C3)*(C5*C5)*C6*Elep*t2*
           t4 + 24.*(amk*amk)*(aml*aml)*apkk1*(C3*C3)*(C5*C5)*C6*
           Elep*t2*t4 + 
           24.*(amk*amk)*(aml*aml)*C3*(C5*C5)*C6*Enu*t2*t4 + 
           24.*akpk*(amk*amk)*(aml*aml)*(C3*C3)*(C5*C5)*C6*Enu*t2*
           t4 + 12.*(amk*amk*amk*amk)*(aml*aml)*(C3*C3)*(C5*C5)*C6*Enu*t2*
           t4 - 12.*(amk*amk)*(aml*aml*aml*aml)*(C3*C3)*(C5*C5)*C6*Enu*t2*
           t4 - 24.*(amk*amk)*(aml*aml)*apkk1*(C3*C3)*(C5*C5)*C6*Enu*
           t2*t4 - 36.*akpk*(amk*amk)*amLam*C1*C4*C5*C6*t3*
           t4 - 36.*akpk*(amk*amk)*amSig*C1*C4*C5*C6*t3*t4 + 
           36.*(amk*amk)*amLam*apkk1*C1*C4*C5*C6*t3*t4 + 
           36.*(amk*amk)*amSig*apkk1*C1*C4*C5*C6*t3*t4 + 
           36.*akpk*(amk*amk)*amLam*C1*(C4*C4)*C5*C6*t3*t4 + 
           36.*akpk*(amk*amk)*amSig*C1*(C4*C4)*C5*C6*t3*t4 + 
           36.*(amk*amk)*amLam*apkk1*C1*(C4*C4)*C5*C6*t3*t4 + 
           36.*(amk*amk)*amSig*apkk1*C1*(C4*C4)*C5*C6*t3*t4 - 
           12.*akpk*(amk*amk)*amLam*C1*C4*(C5*C5)*C6*t3*t4 - 
           12.*akpk*(amk*amk)*amSig*C1*C4*(C5*C5)*C6*t3*t4 - 
           12.*(amk*amk)*amLam*apkk1*C1*C4*(C5*C5)*C6*t3*t4 - 
           12.*(amk*amk)*amSig*apkk1*C1*C4*(C5*C5)*C6*t3*t4 + 
           12.*akpk*(amk*amk)*amLam*C1*(C4*C4)*(C5*C5)*C6*t3*t4 + 
           12.*akpk*(amk*amk)*amSig*C1*(C4*C4)*(C5*C5)*C6*t3*t4 - 
           12.*(amk*amk)*amLam*apkk1*C1*(C4*C4)*(C5*C5)*C6*t3*t4 - 
           12.*(amk*amk)*amSig*apkk1*C1*(C4*C4)*(C5*C5)*C6*t3*t4 - 
           24.*akpk*(amk*amk)*(aml*aml)*amLam*C1*C3*(C4*C4)*(C5*C5)*
           C6*t3*t4 - 
           12.*(amk*amk*amk*amk)*(aml*aml)*amLam*C1*C3*(C4*C4)*(C5*C5)*C6*t3*
           t4 - 24.*akpk*(amk*amk)*(aml*aml)*amSig*C1*C3*(C4*C4)*
           (C5*C5)*C6*t3*t4 - 
           12.*(amk*amk*amk*amk)*(aml*aml)*amSig*C1*C3*(C4*C4)*(C5*C5)*C6*t3*
           t4 - 12.*(akpk*akpk)*(amk*amk)*(aml*aml)*amLam*C1*(C3*C3)*
           (C4*C4)*(C5*C5)*C6*t3*t4 + 
           6.*akpk*(amk*amk)*(aml*aml*aml*aml)*amLam*C1*(C3*C3)*(C4*C4)*(C5*C5)*
           C6*t3*t4 + 
           6.*(amk*amk*amk*amk)*(aml*aml*aml*aml)*amLam*C1*(C3*C3)*(C4*C4)*(C5*C5)*C6*
           t3*t4 - 12.*(akpk*akpk)*(amk*amk)*(aml*aml)*amSig*C1*
           (C3*C3)*(C4*C4)*(C5*C5)*C6*t3*t4 + 
           6.*akpk*(amk*amk)*(aml*aml*aml*aml)*amSig*C1*(C3*C3)*(C4*C4)*(C5*C5)*
           C6*t3*t4 + 
           6.*(amk*amk*amk*amk)*(aml*aml*aml*aml)*amSig*C1*(C3*C3)*(C4*C4)*(C5*C5)*C6*
           t3*t4 + 24.*akpk*(amk*amk)*(aml*aml)*amLam*apkk1*C1*
           (C3*C3)*(C4*C4)*(C5*C5)*C6*t3*t4 - 
           6.*(amk*amk)*(aml*aml*aml*aml)*amLam*apkk1*C1*(C3*C3)*(C4*C4)*
           (C5*C5)*C6*t3*t4 + 
           24.*akpk*(amk*amk)*(aml*aml)*amSig*apkk1*C1*(C3*C3)*
           (C4*C4)*(C5*C5)*C6*t3*t4 - 
           6.*(amk*amk)*(aml*aml*aml*aml)*amSig*apkk1*C1*(C3*C3)*(C4*C4)*
           (C5*C5)*C6*t3*t4 - 
           12.*(amk*amk)*(aml*aml)*amLam*(apkk1*apkk1)*C1*(C3*C3)*(C4*C4)*
           (C5*C5)*C6*t3*t4 - 
           12.*(amk*amk)*(aml*aml)*amSig*(apkk1*apkk1)*C1*(C3*C3)*(C4*C4)*
           (C5*C5)*C6*t3*t4 - 
           72.*akpk*amLam*amSig*C1*C4*C5*C6*Ekaon*t3*t4 + 
           72.*amLam*amSig*apkk1*C1*C4*C5*C6*Ekaon*t3*t4 + 
           72.*akpk*amLam*amSig*C1*(C4*C4)*C5*C6*Ekaon*t3*
           t4 + 72.*amLam*amSig*apkk1*C1*(C4*C4)*C5*C6*Ekaon*
           t3*t4 - 24.*akpk*amLam*amSig*C1*C4*(C5*C5)*C6*
           Ekaon*t3*t4 - 
           24.*amLam*amSig*apkk1*C1*C4*(C5*C5)*C6*Ekaon*t3*
           t4 + 24.*akpk*amLam*amSig*C1*(C4*C4)*(C5*C5)*C6*
           Ekaon*t3*t4 - 
           24.*amLam*amSig*apkk1*C1*(C4*C4)*(C5*C5)*C6*Ekaon*t3*
           t4 - 12.*(amk*amk*amk*amk)*(aml*aml)*C1*C3*(C4*C4)*(C5*C5)*C6*
           Ekaon*t3*t4 - 
           48.*akpk*(aml*aml)*amLam*amSig*C1*C3*(C4*C4)*(C5*C5)*C6*
           Ekaon*t3*t4 - 
           12.*(amk*amk)*(aml*aml)*amLam*amSig*C1*C3*(C4*C4)*(C5*C5)*
           C6*Ekaon*t3*t4 + 
           6.*(amk*amk*amk*amk)*(aml*aml*aml*aml)*C1*(C3*C3)*(C4*C4)*(C5*C5)*C6*Ekaon*
           t3*t4 - 24.*(akpk*akpk)*(aml*aml)*amLam*amSig*C1*(C3*C3)*
           (C4*C4)*(C5*C5)*C6*Ekaon*t3*t4 + 
           12.*akpk*(aml*aml*aml*aml)*amLam*amSig*C1*(C3*C3)*(C4*C4)*(C5*C5)*
           C6*Ekaon*t3*t4 + 
           6.*(amk*amk)*(aml*aml*aml*aml)*amLam*amSig*C1*(C3*C3)*(C4*C4)*
           (C5*C5)*C6*Ekaon*t3*t4 + 
           48.*akpk*(aml*aml)*amLam*amSig*apkk1*C1*(C3*C3)*(C4*C4)*
           (C5*C5)*C6*Ekaon*t3*t4 - 
           12.*(aml*aml*aml*aml)*amLam*amSig*apkk1*C1*(C3*C3)*(C4*C4)*
           (C5*C5)*C6*Ekaon*t3*t4 - 
           24.*(aml*aml)*amLam*amSig*(apkk1*apkk1)*C1*(C3*C3)*(C4*C4)*
           (C5*C5)*C6*Ekaon*t3*t4 + 
           36.*(amk*amk*amk*amk)*C1*C4*C5*C6*Elep*t3*t4 - 
           36.*(amk*amk)*amLam*amSig*C1*C4*C5*C6*Elep*t3*t4 + 
           36.*(amk*amk*amk*amk)*C1*(C4*C4)*C5*C6*Elep*t3*t4 - 
           36.*(amk*amk)*amLam*amSig*C1*(C4*C4)*C5*C6*Elep*t3*
           t4 - 12.*(amk*amk*amk*amk)*C1*C4*(C5*C5)*C6*Elep*t3*t4 + 
           12.*(amk*amk)*amLam*amSig*C1*C4*(C5*C5)*C6*Elep*t3*
           t4 - 12.*(amk*amk*amk*amk)*C1*(C4*C4)*(C5*C5)*C6*Elep*t3*t4 + 
           12.*(amk*amk)*amLam*amSig*C1*(C4*C4)*(C5*C5)*C6*Elep*t3*
           t4 + 12.*akpk*(amk*amk*amk*amk)*(aml*aml)*C1*(C3*C3)*(C4*C4)*
           (C5*C5)*C6*Elep*t3*t4 - 
           6.*(amk*amk*amk*amk)*(aml*aml*aml*aml)*C1*(C3*C3)*(C4*C4)*(C5*C5)*C6*Elep*t3*
           t4 - 12.*akpk*(amk*amk)*(aml*aml)*amLam*amSig*C1*
           (C3*C3)*(C4*C4)*(C5*C5)*C6*Elep*t3*t4 + 
           6.*(amk*amk)*(aml*aml*aml*aml)*amLam*amSig*C1*(C3*C3)*(C4*C4)*
           (C5*C5)*C6*Elep*t3*t4 - 
           12.*(amk*amk*amk*amk)*(aml*aml)*apkk1*C1*(C3*C3)*(C4*C4)*(C5*C5)*C6*
           Elep*t3*t4 + 
           12.*(amk*amk)*(aml*aml)*amLam*amSig*apkk1*C1*(C3*C3)*
           (C4*C4)*(C5*C5)*C6*Elep*t3*t4 - 
           36.*(amk*amk*amk*amk)*C1*C4*C5*C6*Enu*t3*t4 + 
           36.*(amk*amk)*amLam*amSig*C1*C4*C5*C6*Enu*t3*t4 + 
           36.*(amk*amk*amk*amk)*C1*(C4*C4)*C5*C6*Enu*t3*t4 - 
           36.*(amk*amk)*amLam*amSig*C1*(C4*C4)*C5*C6*Enu*t3*
           t4 - 12.*(amk*amk*amk*amk)*C1*C4*(C5*C5)*C6*Enu*t3*t4 + 
           12.*(amk*amk)*amLam*amSig*C1*C4*(C5*C5)*C6*Enu*t3*
           t4 + 12.*(amk*amk*amk*amk)*C1*(C4*C4)*(C5*C5)*C6*Enu*t3*t4 - 
           12.*(amk*amk)*amLam*amSig*C1*(C4*C4)*(C5*C5)*C6*Enu*t3*
           t4 - 24.*(amk*amk*amk*amk)*(aml*aml)*C1*C3*(C4*C4)*(C5*C5)*C6*Enu*
           t3*t4 + 24.*(amk*amk)*(aml*aml)*amLam*amSig*C1*C3*
           (C4*C4)*(C5*C5)*C6*Enu*t3*t4 - 
           12.*akpk*(amk*amk*amk*amk)*(aml*aml)*C1*(C3*C3)*(C4*C4)*(C5*C5)*C6*
           Enu*t3*t4 + 
           6.*(amk*amk*amk*amk)*(aml*aml*aml*aml)*C1*(C3*C3)*(C4*C4)*(C5*C5)*C6*Enu*t3*
           t4 + 12.*akpk*(amk*amk)*(aml*aml)*amLam*amSig*C1*
           (C3*C3)*(C4*C4)*(C5*C5)*C6*Enu*t3*t4 - 
           6.*(amk*amk)*(aml*aml*aml*aml)*amLam*amSig*C1*(C3*C3)*(C4*C4)*
           (C5*C5)*C6*Enu*t3*t4 + 
           12.*(amk*amk*amk*amk)*(aml*aml)*apkk1*C1*(C3*C3)*(C4*C4)*(C5*C5)*C6*
           Enu*t3*t4 - 
           12.*(amk*amk)*(aml*aml)*amLam*amSig*apkk1*C1*(C3*C3)*
           (C4*C4)*(C5*C5)*C6*Enu*t3*t4 - 
           36.*akpk*(amk*amk*amk*amk)*C1*C4*C5*C6*Fm1*t3*t4 + 
           54.*(amk*amk*amk*amk)*(aml*aml)*C1*C4*C5*C6*Fm1*t3*t4 - 
           36.*akpk*(amk*amk)*amLam*amSig*C1*C4*C5*C6*Fm1*t3*
           t4 + 54.*(amk*amk)*(aml*aml)*amLam*amSig*C1*C4*C5*C6*
           Fm1*t3*t4 + 
           36.*(amk*amk*amk*amk)*apkk1*C1*C4*C5*C6*Fm1*t3*t4 + 
           36.*(amk*amk)*amLam*amSig*apkk1*C1*C4*C5*C6*Fm1*t3*
           t4 + 12.*akpk*(amk*amk*amk*amk)*C1*C4*(C5*C5)*C6*Fm1*t3*t4 + 
           12.*akpk*(amk*amk)*amLam*amSig*C1*C4*(C5*C5)*C6*Fm1*
           t3*t4 + 12.*(amk*amk*amk*amk)*apkk1*C1*C4*(C5*C5)*C6*Fm1*t3*
           t4 + 12.*(amk*amk)*amLam*amSig*apkk1*C1*C4*(C5*C5)*
           C6*Fm1*t3*t4 - 
           36.*akpk*(amk*amk)*amLam*C1*C4*C5*C6*Ekaon*Fm1*t3*
           t4 + 54.*(amk*amk)*(aml*aml)*amLam*C1*C4*C5*C6*Ekaon*
           Fm1*t3*t4 - 
           36.*akpk*(amk*amk)*amSig*C1*C4*C5*C6*Ekaon*Fm1*t3*
           t4 + 54.*(amk*amk)*(aml*aml)*amSig*C1*C4*C5*C6*Ekaon*
           Fm1*t3*t4 + 
           36.*(amk*amk)*amLam*apkk1*C1*C4*C5*C6*Ekaon*Fm1*t3*
           t4 + 36.*(amk*amk)*amSig*apkk1*C1*C4*C5*C6*Ekaon*
           Fm1*t3*t4 + 
           12.*akpk*(amk*amk)*amLam*C1*C4*(C5*C5)*C6*Ekaon*Fm1*
           t3*t4 + 12.*akpk*(amk*amk)*amSig*C1*C4*(C5*C5)*C6*
           Ekaon*Fm1*t3*t4 + 
           12.*(amk*amk)*amLam*apkk1*C1*C4*(C5*C5)*C6*Ekaon*Fm1*
           t3*t4 + 12.*(amk*amk)*amSig*apkk1*C1*C4*(C5*C5)*C6*
           Ekaon*Fm1*t3*t4 - 
           24.*akpk*(amk*amk)*amLam*C1*C4*(C5*C5)*C6*Elep*Fm1*t3*
           t4 + 24.*akpk*(amk*amk)*amSig*C1*C4*(C5*C5)*C6*Elep*
           Fm1*t3*t4 + 
           24.*(amk*amk)*amLam*apkk1*C1*C4*(C5*C5)*C6*Enu*Fm1*t3*
           t4 - 24.*(amk*amk)*amSig*apkk1*C1*C4*(C5*C5)*C6*Enu*
           Fm1*t3*t4 - 
           12.*akpk*(amk*amk*amk*amk)*C1*C4*C5*C6*Fm2*t3*t4 + 
           18.*(amk*amk*amk*amk)*(aml*aml)*C1*C4*C5*C6*Fm2*t3*t4 - 
           12.*akpk*(amk*amk)*amLam*amSig*C1*C4*C5*C6*Fm2*t3*
           t4 + 18.*(amk*amk)*(aml*aml)*amLam*amSig*C1*C4*C5*C6*
           Fm2*t3*t4 + 
           12.*(amk*amk*amk*amk)*apkk1*C1*C4*C5*C6*Fm2*t3*t4 + 
           12.*(amk*amk)*amLam*amSig*apkk1*C1*C4*C5*C6*Fm2*t3*
           t4 - 12.*akpk*(amk*amk*amk*amk)*C1*(C4*C4)*C5*C6*Fm2*t3*t4 - 
           12.*akpk*(amk*amk)*amLam*amSig*C1*(C4*C4)*C5*C6*Fm2*
           t3*t4 - 12.*(amk*amk*amk*amk)*apkk1*C1*(C4*C4)*C5*C6*Fm2*t3*
           t4 - 12.*(amk*amk)*amLam*amSig*apkk1*C1*(C4*C4)*C5*
           C6*Fm2*t3*t4 - 
           12.*akpk*(amk*amk)*amLam*C1*C4*C5*C6*Ekaon*Fm2*t3*
           t4 + 18.*(amk*amk)*(aml*aml)*amLam*C1*C4*C5*C6*Ekaon*
           Fm2*t3*t4 - 
           12.*akpk*(amk*amk)*amSig*C1*C4*C5*C6*Ekaon*Fm2*t3*
           t4 + 18.*(amk*amk)*(aml*aml)*amSig*C1*C4*C5*C6*Ekaon*
           Fm2*t3*t4 + 
           12.*(amk*amk)*amLam*apkk1*C1*C4*C5*C6*Ekaon*Fm2*t3*
           t4 + 12.*(amk*amk)*amSig*apkk1*C1*C4*C5*C6*Ekaon*
           Fm2*t3*t4 - 
           12.*akpk*(amk*amk)*amLam*C1*(C4*C4)*C5*C6*Ekaon*Fm2*
           t3*t4 - 12.*akpk*(amk*amk)*amSig*C1*(C4*C4)*C5*C6*
           Ekaon*Fm2*t3*t4 - 
           12.*(amk*amk)*amLam*apkk1*C1*(C4*C4)*C5*C6*Ekaon*Fm2*
           t3*t4 - 12.*(amk*amk)*amSig*apkk1*C1*(C4*C4)*C5*C6*
           Ekaon*Fm2*t3*t4 - 
           24.*akpk*(amk*amk)*amLam*C1*(C4*C4)*C5*C6*Elep*Fm2*t3*
           t4 + 24.*akpk*(amk*amk)*amSig*C1*(C4*C4)*C5*C6*Elep*
           Fm2*t3*t4 + 
           24.*(amk*amk)*amLam*apkk1*C1*(C4*C4)*C5*C6*Enu*Fm2*t3*
           t4 - 24.*(amk*amk)*amSig*apkk1*C1*(C4*C4)*C5*C6*Enu*
           Fm2*t3*t4 + 
           24.*(akpk*akpk)*(amk*amk)*amLam*C1*C4*C5*C6*Fm1*Fm2*t3*
           t4 - 18.*akpk*(amk*amk)*(aml*aml)*amLam*C1*C4*C5*C6*
           Fm1*Fm2*t3*t4 - 
           6.*(amk*amk*amk*amk)*(aml*aml)*amLam*C1*C4*C5*C6*Fm1*Fm2*t3*
           t4 + 24.*(akpk*akpk)*(amk*amk)*amSig*C1*C4*C5*C6*Fm1*
           Fm2*t3*t4 - 
           18.*akpk*(amk*amk)*(aml*aml)*amSig*C1*C4*C5*C6*Fm1*Fm2*
           t3*t4 - 6.*(amk*amk*amk*amk)*(aml*aml)*amSig*C1*C4*C5*C6*Fm1*
           Fm2*t3*t4 + 
           18.*(amk*amk)*(aml*aml)*amLam*apkk1*C1*C4*C5*C6*Fm1*
           Fm2*t3*t4 + 
           18.*(amk*amk)*(aml*aml)*amSig*apkk1*C1*C4*C5*C6*Fm1*
           Fm2*t3*t4 + 
           24.*(amk*amk)*amLam*(apkk1*apkk1)*C1*C4*C5*C6*Fm1*Fm2*t3*
           t4 + 24.*(amk*amk)*amSig*(apkk1*apkk1)*C1*C4*C5*C6*Fm1*
           Fm2*t3*t4 - 
           6.*(amk*amk*amk*amk)*(aml*aml)*C1*C4*C5*C6*Ekaon*Fm1*Fm2*t3*
           t4 + 48.*(akpk*akpk)*amLam*amSig*C1*C4*C5*C6*Ekaon*
           Fm1*Fm2*t3*t4 - 
           36.*akpk*(aml*aml)*amLam*amSig*C1*C4*C5*C6*Ekaon*
           Fm1*Fm2*t3*t4 - 
           6.*(amk*amk)*(aml*aml)*amLam*amSig*C1*C4*C5*C6*Ekaon*
           Fm1*Fm2*t3*t4 + 
           36.*(aml*aml)*amLam*amSig*apkk1*C1*C4*C5*C6*Ekaon*
           Fm1*Fm2*t3*t4 + 
           48.*amLam*amSig*(apkk1*apkk1)*C1*C4*C5*C6*Ekaon*Fm1*
           Fm2*t3*t4 + 
           18.*(amk*amk*amk*amk)*(aml*aml)*C1*C4*C5*C6*Elep*Fm1*Fm2*t3*
           t4 - 18.*(amk*amk)*(aml*aml)*amLam*amSig*C1*C4*C5*C6*
           Elep*Fm1*Fm2*t3*t4 + 
           24.*(amk*amk*amk*amk)*apkk1*C1*C4*C5*C6*Elep*Fm1*Fm2*t3*
           t4 - 24.*(amk*amk)*amLam*amSig*apkk1*C1*C4*C5*C6*
           Elep*Fm1*Fm2*t3*t4 + 
           24.*akpk*(amk*amk*amk*amk)*C1*C4*C5*C6*Enu*Fm1*Fm2*t3*t4 - 
           18.*(amk*amk*amk*amk)*(aml*aml)*C1*C4*C5*C6*Enu*Fm1*Fm2*t3*
           t4 - 24.*akpk*(amk*amk)*amLam*amSig*C1*C4*C5*C6*
           Enu*Fm1*Fm2*t3*t4 + 
           18.*(amk*amk)*(aml*aml)*amLam*amSig*C1*C4*C5*C6*Enu*
           Fm1*Fm2*t3*t4 + 
           36.*akpk*(amk*amk)*amLam*(C5*C5)*(C6*C6)*(t4*t4) - 
           36.*(amk*amk)*amLam*apkk1*(C5*C5)*(C6*C6)*(t4*t4) + 
           24.*akpk*(amk*amk)*amLam*(C5*C5*C5)*(C6*C6)*(t4*t4) + 
           24.*(amk*amk)*amLam*apkk1*(C5*C5*C5)*(C6*C6)*(t4*t4) + 
           4.*akpk*(amk*amk)*amLam*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) - 
           4.*(amk*amk)*amLam*apkk1*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) - 
           8.*akpk*(amk*amk)*(aml*aml)*amLam*C3*(C5*C5*C5*C5)*(C6*C6)*
           (t4*t4) - 4.*(amk*amk*amk*amk)*(aml*aml)*amLam*C3*(C5*C5*C5*C5)*(C6*C6)*
           (t4*t4) - 4.*(akpk*akpk)*(amk*amk)*(aml*aml)*amLam*(C3*C3)*
           (C5*C5*C5*C5)*(C6*C6)*(t4*t4) + 
           2.*akpk*(amk*amk)*(aml*aml*aml*aml)*amLam*(C3*C3)*(C5*C5*C5*C5)*(C6*C6)*
           (t4*t4) + 2.*(amk*amk*amk*amk)*(aml*aml*aml*aml)*amLam*(C3*C3)*(C5*C5*C5*C5)*
           (C6*C6)*(t4*t4) + 
           8.*akpk*(amk*amk)*(aml*aml)*amLam*apkk1*(C3*C3)*(C5*C5*C5*C5)*
           (C6*C6)*(t4*t4) - 
           2.*(amk*amk)*(aml*aml*aml*aml)*amLam*apkk1*(C3*C3)*(C5*C5*C5*C5)*(C6*C6)*
           (t4*t4) - 4.*(amk*amk)*(aml*aml)*amLam*(apkk1*apkk1)*(C3*C3)*
           (C5*C5*C5*C5)*(C6*C6)*(t4*t4) + 
           36.*akpk*(amLam*amLam)*(C5*C5)*(C6*C6)*Ekaon*(t4*t4) - 
           36.*(amLam*amLam)*apkk1*(C5*C5)*(C6*C6)*Ekaon*(t4*t4) + 
           24.*akpk*(amLam*amLam)*(C5*C5*C5)*(C6*C6)*Ekaon*(t4*t4) + 
           24.*(amLam*amLam)*apkk1*(C5*C5*C5)*(C6*C6)*Ekaon*(t4*t4) + 
           4.*akpk*(amLam*amLam)*(C5*C5*C5*C5)*(C6*C6)*Ekaon*(t4*t4) - 
           4.*(amLam*amLam)*apkk1*(C5*C5*C5*C5)*(C6*C6)*Ekaon*(t4*t4) - 
           2.*(amk*amk*amk*amk)*(aml*aml)*C3*(C5*C5*C5*C5)*(C6*C6)*Ekaon*(t4*t4) - 
           8.*akpk*(aml*aml)*(amLam*amLam)*C3*(C5*C5*C5*C5)*(C6*C6)*Ekaon*
           (t4*t4) - 2.*(amk*amk)*(aml*aml)*(amLam*amLam)*C3*(C5*C5*C5*C5)*
           (C6*C6)*Ekaon*(t4*t4) + 
           (amk*amk*amk*amk)*(aml*aml*aml*aml)*(C3*C3)*(C5*C5*C5*C5)*(C6*C6)*Ekaon*(t4*t4) - 
           4.*(akpk*akpk)*(aml*aml)*(amLam*amLam)*(C3*C3)*(C5*C5*C5*C5)*(C6*C6)*
           Ekaon*(t4*t4) + 
           2.*akpk*(aml*aml*aml*aml)*(amLam*amLam)*(C3*C3)*(C5*C5*C5*C5)*(C6*C6)*Ekaon*
           (t4*t4) + (amk*amk)*(aml*aml*aml*aml)*(amLam*amLam)*(C3*C3)*(C5*C5*C5*C5)*
           (C6*C6)*Ekaon*(t4*t4) + 
           8.*akpk*(aml*aml)*(amLam*amLam)*apkk1*(C3*C3)*(C5*C5*C5*C5)*(C6*C6)*
           Ekaon*(t4*t4) - 
           2.*(aml*aml*aml*aml)*(amLam*amLam)*apkk1*(C3*C3)*(C5*C5*C5*C5)*(C6*C6)*Ekaon*
           (t4*t4) - 4.*(aml*aml)*(amLam*amLam)*(apkk1*apkk1)*(C3*C3)*(C5*C5*C5*C5)*
           (C6*C6)*Ekaon*(t4*t4) - 
           18.*(amk*amk*amk*amk)*(C5*C5)*(C6*C6)*Elep*(t4*t4) + 
           18.*(amk*amk)*(amLam*amLam)*(C5*C5)*(C6*C6)*Elep*(t4*t4) + 
           12.*(amk*amk*amk*amk)*(C5*C5*C5)*(C6*C6)*Elep*(t4*t4) - 
           12.*(amk*amk)*(amLam*amLam)*(C5*C5*C5)*(C6*C6)*Elep*(t4*t4) - 
           2.*(amk*amk*amk*amk)*(C5*C5*C5*C5)*(C6*C6)*Elep*(t4*t4) + 
           2.*(amk*amk)*(amLam*amLam)*(C5*C5*C5*C5)*(C6*C6)*Elep*(t4*t4) + 
           2.*akpk*(amk*amk*amk*amk)*(aml*aml)*(C3*C3)*(C5*C5*C5*C5)*(C6*C6)*Elep*
           (t4*t4) - 1.*(amk*amk*amk*amk)*(aml*aml*aml*aml)*(C3*C3)*(C5*C5*C5*C5)*(C6*C6)*Elep*
           (t4*t4) - 2.*akpk*(amk*amk)*(aml*aml)*(amLam*amLam)*(C3*C3)*
           (C5*C5*C5*C5)*(C6*C6)*Elep*(t4*t4) + 
           (amk*amk)*(aml*aml*aml*aml)*(amLam*amLam)*(C3*C3)*(C5*C5*C5*C5)*(C6*C6)*Elep*
           (t4*t4) - 2.*(amk*amk*amk*amk)*(aml*aml)*apkk1*(C3*C3)*(C5*C5*C5*C5)*
           (C6*C6)*Elep*(t4*t4) + 
           2.*(amk*amk)*(aml*aml)*(amLam*amLam)*apkk1*(C3*C3)*(C5*C5*C5*C5)*
           (C6*C6)*Elep*(t4*t4) + 
           18.*(amk*amk*amk*amk)*(C5*C5)*(C6*C6)*Enu*(t4*t4) - 
           18.*(amk*amk)*(amLam*amLam)*(C5*C5)*(C6*C6)*Enu*(t4*t4) + 
           12.*(amk*amk*amk*amk)*(C5*C5*C5)*(C6*C6)*Enu*(t4*t4) - 
           12.*(amk*amk)*(amLam*amLam)*(C5*C5*C5)*(C6*C6)*Enu*(t4*t4) + 
           2.*(amk*amk*amk*amk)*(C5*C5*C5*C5)*(C6*C6)*Enu*(t4*t4) - 
           2.*(amk*amk)*(amLam*amLam)*(C5*C5*C5*C5)*(C6*C6)*Enu*(t4*t4) - 
           4.*(amk*amk*amk*amk)*(aml*aml)*C3*(C5*C5*C5*C5)*(C6*C6)*Enu*(t4*t4) + 
           4.*(amk*amk)*(aml*aml)*(amLam*amLam)*C3*(C5*C5*C5*C5)*(C6*C6)*Enu*
           (t4*t4) - 2.*akpk*(amk*amk*amk*amk)*(aml*aml)*(C3*C3)*(C5*C5*C5*C5)*(C6*C6)*
           Enu*(t4*t4) + 
           (amk*amk*amk*amk)*(aml*aml*aml*aml)*(C3*C3)*(C5*C5*C5*C5)*(C6*C6)*Enu*(t4*t4) + 
           2.*akpk*(amk*amk)*(aml*aml)*(amLam*amLam)*(C3*C3)*(C5*C5*C5*C5)*(C6*C6)*
           Enu*(t4*t4) - 
           1.*(amk*amk)*(aml*aml*aml*aml)*(amLam*amLam)*(C3*C3)*(C5*C5*C5*C5)*(C6*C6)*Enu*
           (t4*t4) + 2.*(amk*amk*amk*amk)*(aml*aml)*apkk1*(C3*C3)*(C5*C5*C5*C5)*
           (C6*C6)*Enu*(t4*t4) - 
           2.*(amk*amk)*(aml*aml)*(amLam*amLam)*apkk1*(C3*C3)*(C5*C5*C5*C5)*
           (C6*C6)*Enu*(t4*t4) + 
           12.*akpk*(amk*amk*amk*amk)*(C5*C5)*(C6*C6)*Fm2*(t4*t4) - 
           18.*(amk*amk*amk*amk)*(aml*aml)*(C5*C5)*(C6*C6)*Fm2*(t4*t4) + 
           12.*akpk*(amk*amk)*(amLam*amLam)*(C5*C5)*(C6*C6)*Fm2*(t4*t4) - 
           18.*(amk*amk)*(aml*aml)*(amLam*amLam)*(C5*C5)*(C6*C6)*Fm2*
           (t4*t4) - 12.*(amk*amk*amk*amk)*apkk1*(C5*C5)*(C6*C6)*Fm2*
           (t4*t4) - 12.*(amk*amk)*(amLam*amLam)*apkk1*(C5*C5)*(C6*C6)*
           Fm2*(t4*t4) - 
           4.*akpk*(amk*amk*amk*amk)*(C5*C5*C5)*(C6*C6)*Fm2*(t4*t4) - 
           4.*akpk*(amk*amk)*(amLam*amLam)*(C5*C5*C5)*(C6*C6)*Fm2*(t4*t4) - 
           4.*(amk*amk*amk*amk)*apkk1*(C5*C5*C5)*(C6*C6)*Fm2*(t4*t4) - 
           4.*(amk*amk)*(amLam*amLam)*apkk1*(C5*C5*C5)*(C6*C6)*Fm2*(t4*t4) + 
           24.*akpk*(amk*amk)*amLam*(C5*C5)*(C6*C6)*Ekaon*Fm2*
           (t4*t4) - 36.*(amk*amk)*(aml*aml)*amLam*(C5*C5)*(C6*C6)*
           Ekaon*Fm2*(t4*t4) - 
           24.*(amk*amk)*amLam*apkk1*(C5*C5)*(C6*C6)*Ekaon*Fm2*
           (t4*t4) - 8.*akpk*(amk*amk)*amLam*(C5*C5*C5)*(C6*C6)*Ekaon*
           Fm2*(t4*t4) - 
           8.*(amk*amk)*amLam*apkk1*(C5*C5*C5)*(C6*C6)*Ekaon*Fm2*
           (t4*t4) - 8.*(akpk*akpk)*(amk*amk)*amLam*(C5*C5)*(C6*C6)*
           (Fm2*Fm2)*(t4*t4) + 
           6.*akpk*(amk*amk)*(aml*aml)*amLam*(C5*C5)*(C6*C6)*(Fm2*Fm2)*
           (t4*t4) + 2.*(amk*amk*amk*amk)*(aml*aml)*amLam*(C5*C5)*(C6*C6)*
           (Fm2*Fm2)*(t4*t4) - 
           6.*(amk*amk)*(aml*aml)*amLam*apkk1*(C5*C5)*(C6*C6)*(Fm2*Fm2)*
           (t4*t4) - 8.*(amk*amk)*amLam*(apkk1*apkk1)*(C5*C5)*(C6*C6)*
           (Fm2*Fm2)*(t4*t4) + 
           (amk*amk*amk*amk)*(aml*aml)*(C5*C5)*(C6*C6)*Ekaon*(Fm2*Fm2)*(t4*t4) - 
           8.*(akpk*akpk)*(amLam*amLam)*(C5*C5)*(C6*C6)*Ekaon*(Fm2*Fm2)*
           (t4*t4) + 6.*akpk*(aml*aml)*(amLam*amLam)*(C5*C5)*(C6*C6)*
           Ekaon*(Fm2*Fm2)*(t4*t4) + 
           (amk*amk)*(aml*aml)*(amLam*amLam)*(C5*C5)*(C6*C6)*Ekaon*(Fm2*Fm2)*
           (t4*t4) - 6.*(aml*aml)*(amLam*amLam)*apkk1*(C5*C5)*(C6*C6)*
           Ekaon*(Fm2*Fm2)*(t4*t4) - 
           8.*(amLam*amLam)*(apkk1*apkk1)*(C5*C5)*(C6*C6)*Ekaon*(Fm2*Fm2)*
           (t4*t4) - 3.*(amk*amk*amk*amk)*(aml*aml)*(C5*C5)*(C6*C6)*Elep*
           (Fm2*Fm2)*(t4*t4) + 
           3.*(amk*amk)*(aml*aml)*(amLam*amLam)*(C5*C5)*(C6*C6)*Elep*
           (Fm2*Fm2)*(t4*t4) - 
           4.*(amk*amk*amk*amk)*apkk1*(C5*C5)*(C6*C6)*Elep*(Fm2*Fm2)*(t4*t4) + 
           4.*(amk*amk)*(amLam*amLam)*apkk1*(C5*C5)*(C6*C6)*Elep*(Fm2*Fm2)*
           (t4*t4) - 4.*akpk*(amk*amk*amk*amk)*(C5*C5)*(C6*C6)*Enu*(Fm2*Fm2)*
           (t4*t4) + 3.*(amk*amk*amk*amk)*(aml*aml)*(C5*C5)*(C6*C6)*Enu*(Fm2*Fm2)*
           (t4*t4) + 4.*akpk*(amk*amk)*(amLam*amLam)*(C5*C5)*(C6*C6)*Enu*
           (Fm2*Fm2)*(t4*t4) - 
           3.*(amk*amk)*(aml*aml)*(amLam*amLam)*(C5*C5)*(C6*C6)*Enu*(Fm2*Fm2)*
           (t4*t4) + 2.*(am*am*am)*
           ((amk*amk)*(9.*(C1*C1)*(C4*C4)*
           (1. + 
           (C4*C4)*
           (-1. + (aml*aml*aml*aml)*(C3*C3) - 
           1.*(aml*aml)*C3*
           (2. + C3*((Elep - 1.*Enu)*(Elep - 1.*Enu)))) - 
           4.*C4*(Elep + Enu)*Fm1 + 
           2.*(aml*aml)*(Fm1*Fm1) - 2.*(Elep*Elep)*(Fm1*Fm1) - 
           2.*(Enu*Enu)*(Fm1*Fm1))*(t3*t3) + 
           6.*C1*C4*C5*C6*
           (-3. - 2.*C5*(Elep + Enu)*Fm1 - 
           2.*(aml*aml)*Fm1*Fm2 + 
           2.*(Elep*Elep)*Fm1*Fm2 + 
           2.*(Enu*Enu)*Fm1*Fm2 + 
           C4*
           (C5*
           (-1. + (aml*aml*aml*aml)*(C3*C3) - 
           1.*(aml*aml)*C3*
           (2. + C3*((Elep - 1.*Enu)*(Elep - 1.*Enu)))) + 
           2.*(Elep + Enu)*Fm2))*t3*t4 + 
           (C5*C5)*(C6*C6)*
           (9. + 
           (C5*C5)*
           (-1. + (aml*aml*aml*aml)*(C3*C3) - 
           1.*(aml*aml)*C3*
           (2. + C3*((Elep - 1.*Enu)*(Elep - 1.*Enu)))) + 
           4.*C5*(Elep + Enu)*Fm2 + 
           2.*(aml*aml)*(Fm2*Fm2) - 2.*(Elep*Elep)*(Fm2*Fm2) - 
           2.*(Enu*Enu)*(Fm2*Fm2))*(t4*t4)) + 
           Ekaon*(-1.*(aml*aml*aml*aml)*(C3*C3)*Ekaon*
           ((3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4)*(3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4)) + 
           (aml*aml)*
           (9.*(C1*C1)*(C4*C4)*
           (2.*C3*(C4*C4)*Ekaon + 
           2.*(C3*C3)*(C4*C4)*(Elep - 1.*Enu)*
           (akpk - 1.*apkk1 + 2.*Ekaon*Elep - 
           2.*Ekaon*Enu) - 1.*Ekaon*(Fm1*Fm1))*
           (t3*t3) + 
           6.*C1*C4*C5*C6*
           (2.*C3*C4*C5*Ekaon + 
           2.*(C3*C3)*C4*C5*(Elep - 1.*Enu)*
           (akpk - 1.*apkk1 + 2.*Ekaon*Elep - 
           2.*Ekaon*Enu) + Ekaon*Fm1*Fm2)*t3*
           t4 + 
           (C5*C5)*(C6*C6)*
           (2.*C3*(C5*C5)*Ekaon + 
           2.*(C3*C3)*(C5*C5)*(Elep - 1.*Enu)*
           (akpk - 1.*apkk1 + 2.*Ekaon*Elep - 
           2.*Ekaon*Enu) - 1.*Ekaon*(Fm2*Fm2))*
           (t4*t4)) - 
           4.*(3.*C1*C4*Fm1*t3 - 1.*C5*C6*Fm2*t4)*
           (apkk1*
           (3.*C1*C4*(1. + C4 + Elep*Fm1)*t3 + 
           C5*C6*(-3. + C5 - 1.*Elep*Fm2)*t4) + 
           akpk*
           (3.*C1*C4*(-1. + C4 + Enu*Fm1)*t3 + 
           C5*C6*(3. + C5 - 1.*Enu*Fm2)*t4) - 
           2.*Ekaon*
           (3.*C1*C4*
           ((1. + C4)*Elep + (Elep*Elep)*Fm1 + 
           Enu*(-1. + C4 + Enu*Fm1))*t3 + 
           C5*C6*
           ((-3. + C5)*Elep - 1.*(Elep*Elep)*Fm2 + 
           Enu*(3. + C5 - 1.*Enu*Fm2))*t4))))
           + (am*am)*(9.*(C1*C1)*(C4*C4)*
           ((amk*amk)*
           (2.*Elep + 4.*C4*Elep + 2.*(C4*C4)*Elep - 
           6.*akpk*(aml*aml)*(C3*C3)*(C4*C4)*Elep + 
           (aml*aml*aml*aml)*(C3*C3)*(C4*C4)*Elep + 
           6.*(aml*aml)*apkk1*(C3*C3)*(C4*C4)*Elep - 
           2.*Enu + 4.*C4*Enu - 2.*(C4*C4)*Enu + 
           4.*(aml*aml)*C3*(C4*C4)*Enu + 
           6.*akpk*(aml*aml)*(C3*C3)*(C4*C4)*Enu - 
           1.*(aml*aml*aml*aml)*(C3*C3)*(C4*C4)*Enu - 
           6.*(aml*aml)*apkk1*(C3*C3)*(C4*C4)*Enu - 
           4.*akpk*Fm1 - 6.*(aml*aml)*Fm1 + 
           4.*apkk1*Fm1 + 12.*akpk*C4*Fm1 + 
           12.*apkk1*C4*Fm1 + 
           3.*(aml*aml)*Elep*(Fm1*Fm1) + 
           12.*apkk1*Elep*(Fm1*Fm1) + 
           12.*akpk*Enu*(Fm1*Fm1) - 
           3.*(aml*aml)*Enu*(Fm1*Fm1) + 
           4.*amSig*
           (1. + 
           (C4*C4)*
           (-1. + (aml*aml*aml*aml)*(C3*C3) - 
           1.*(aml*aml)*C3*
           (2. + C3*((Elep - 1.*Enu)*(Elep - 1.*Enu)))) - 
           4.*C4*(Elep + Enu)*Fm1 + 
           2.*(aml*aml)*(Fm1*Fm1) - 
           2.*(Elep*Elep)*(Fm1*Fm1) - 2.*(Enu*Enu)*(Fm1*Fm1))
           - 1.*Ekaon*
           (4. + 
           (C4*C4)*
           (-4. + 3.*(aml*aml*aml*aml)*(C3*C3) + 
           2.*(aml*aml)*C3*
           (-3. + 4.*C3*((Elep - 1.*Enu)*(Elep - 1.*Enu))))
           + 24.*Elep*Fm1 - 24.*Enu*Fm1 + 
           8.*C4*(Elep + Enu)*Fm1 + 
           7.*(aml*aml)*(Fm1*Fm1) + 
           16.*(Elep*Elep)*(Fm1*Fm1) + 
           16.*(Enu*Enu)*(Fm1*Fm1))) + 
           2.*Ekaon*
           (2.*(akpk*akpk)*
           ((aml*aml)*(C3*C3)*(C4*C4) + 2.*(Fm1*Fm1)) + 
           2.*(apkk1*apkk1)*
           ((aml*aml)*(C3*C3)*(C4*C4) + 2.*(Fm1*Fm1)) + 
           apkk1*
           (2. + 4.*C4 + 
           (C4*C4)*
           (2. + (aml*aml*aml*aml)*(C3*C3) + 
           4.*(aml*aml)*(C3*C3)*Ekaon*
           (-1.*Elep + Enu)) + 
           3.*(aml*aml)*(Fm1*Fm1) - 
           8.*Ekaon*Elep*(Fm1*Fm1)) - 
           1.*akpk*
           (2. - 4.*C4 + 
           (C4*C4)*
           (2. + (aml*aml*aml*aml)*(C3*C3) + 
           4.*(aml*aml)*C3*
           (-1. + apkk1*C3 + 
           C3*Ekaon*(-1.*Elep + Enu))) + 
           3.*(aml*aml)*(Fm1*Fm1) + 
           8.*Ekaon*Enu*(Fm1*Fm1)) + 
           2.*Ekaon*
           ((aml*aml*aml*aml)*(C3*C3)*(C4*C4)*
           (Ekaon - 1.*Elep + Enu) - 
           2.*
           (((1.+C4)*(1.+C4))*Elep - 
           1.*((-1.+C4)*(-1.+C4))*Enu) + 
           amSig*
           (-2. + 
           (2. + 2.*(aml*aml)*C3 - 
           1.*(aml*aml*aml*aml)*(C3*C3))*(C4*C4) - 
           4.*Elep*Fm1 + 4.*Enu*Fm1 + 
           4.*C4*(Elep + Enu)*Fm1 - 
           3.*(aml*aml)*(Fm1*Fm1)) + 
           (aml*aml)*
           (-2.*C3*(C4*C4)*(Ekaon + 2.*Enu) + 
           (Ekaon - 3.*Elep + 3.*Enu)*(Fm1*Fm1)))
           ))*(t3*t3) + 
           6.*C1*C4*t3*
           (-24.*akpk*Fm1*t1 - 24.*apkk1*Fm1*t1 - 
           24.*akpk*f*Fm1*t1 + 24.*apkk1*f*Fm1*t1 - 
           6.*(amk*amk)*amLam*C5*C6*t4 - 
           6.*(amk*amk)*amSig*C5*C6*t4 - 
           6.*(amk*amk)*C5*C6*Elep*t4 + 
           2.*(amk*amk)*(C5*C5)*C6*Elep*t4 + 
           6.*(amk*amk)*C5*C6*Enu*t4 + 
           2.*(amk*amk)*(C5*C5)*C6*Enu*t4 + 
           6.*akpk*(amk*amk)*C5*C6*Fm1*t4 + 
           9.*(amk*amk)*(aml*aml)*C5*C6*Fm1*t4 - 
           6.*(amk*amk)*apkk1*C5*C6*Fm1*t4 + 
           6.*akpk*(amk*amk)*(C5*C5)*C6*Fm1*t4 + 
           6.*(amk*amk)*apkk1*(C5*C5)*C6*Fm1*t4 - 
           4.*(amk*amk)*amLam*(C5*C5)*C6*Elep*Fm1*t4 - 
           4.*(amk*amk)*amSig*(C5*C5)*C6*Elep*Fm1*t4 - 
           4.*(amk*amk)*amLam*(C5*C5)*C6*Enu*Fm1*t4 - 
           4.*(amk*amk)*amSig*(C5*C5)*C6*Enu*Fm1*t4 + 
           2.*akpk*(amk*amk)*C5*C6*Fm2*t4 + 
           3.*(amk*amk)*(aml*aml)*C5*C6*Fm2*t4 - 
           2.*(amk*amk)*apkk1*C5*C6*Fm2*t4 - 
           4.*(amk*amk)*(aml*aml)*amLam*C5*C6*Fm1*Fm2*t4 - 
           4.*(amk*amk)*(aml*aml)*amSig*C5*C6*Fm1*Fm2*t4 - 
           4.*(aml*aml)*C5*C6*(Ekaon*Ekaon*Ekaon)*Fm1*Fm2*t4 - 
           3.*(amk*amk)*(aml*aml)*C5*C6*Elep*Fm1*Fm2*t4 - 
           12.*(amk*amk)*apkk1*C5*C6*Elep*Fm1*Fm2*t4 + 
           4.*(amk*amk)*amLam*C5*C6*(Elep*Elep)*Fm1*Fm2*t4 + 
           4.*(amk*amk)*amSig*C5*C6*(Elep*Elep)*Fm1*Fm2*t4 - 
           12.*akpk*(amk*amk)*C5*C6*Enu*Fm1*Fm2*t4 + 
           3.*(amk*amk)*(aml*aml)*C5*C6*Enu*Fm1*Fm2*t4 + 
           4.*(amk*amk)*amLam*C5*C6*(Enu*Enu)*Fm1*Fm2*t4 + 
           4.*(amk*amk)*amSig*C5*C6*(Enu*Enu)*Fm1*Fm2*t4 + 
           2.*C5*C6*(Ekaon*Ekaon)*
           (amLam*
           (6. + 
           2.*Enu*((-3. + C5)*Fm1 - 1.*Fm2) + 
           3.*(aml*aml)*Fm1*Fm2 + 
           2.*Elep*((3. + C5)*Fm1 + Fm2)) + 
           amSig*
           (6. + 
           2.*Enu*((-3. + C5)*Fm1 - 1.*Fm2) + 
           3.*(aml*aml)*Fm1*Fm2 + 
           2.*Elep*((3. + C5)*Fm1 + Fm2)) + 
           2.*
           (-1.*Enu*
           (6. + 2.*C5 - 4.*akpk*Fm1*Fm2 + 
           3.*(aml*aml)*Fm1*Fm2) + 
           Elep*
           (6. - 2.*C5 + 3.*(aml*aml)*Fm1*Fm2 + 
           4.*apkk1*Fm1*Fm2)))*t4 + 
           Ekaon*
           (24.*Enu*Fm1*t1 - 
           24.*f*(1. + 3.*Elep*Fm1 - 3.*Enu*Fm1)*
           t1 + 12.*akpk*C5*C6*t4 + 
           12.*(amk*amk)*C5*C6*t4 - 
           12.*apkk1*C5*C6*t4 + 
           4.*akpk*(C5*C5)*C6*t4 + 
           4.*apkk1*(C5*C5)*C6*t4 - 
           36.*(amk*amk)*C5*C6*Enu*Fm1*t4 - 
           4.*(amk*amk)*(C5*C5)*C6*Enu*Fm1*t4 - 
           12.*(amk*amk)*C5*C6*Enu*Fm2*t4 - 
           8.*(akpk*akpk)*C5*C6*Fm1*Fm2*t4 + 
           6.*akpk*(aml*aml)*C5*C6*Fm1*Fm2*t4 + 
           7.*(amk*amk)*(aml*aml)*C5*C6*Fm1*Fm2*t4 - 
           6.*(aml*aml)*apkk1*C5*C6*Fm1*Fm2*t4 - 
           8.*(apkk1*apkk1)*C5*C6*Fm1*Fm2*t4 + 
           16.*(amk*amk)*C5*C6*(Elep*Elep)*Fm1*Fm2*t4 + 
           16.*(amk*amk)*C5*C6*(Enu*Enu)*Fm1*Fm2*t4 + 
           4.*Elep*
           (3.*(amk*amk)*C5*C6*Fm2*t4 + 
           Fm1*
           (6.*t1 - 
           1.*(amk*amk)*(-9. + C5)*C5*C6*t4))) + 
           C4*(-12.*(aml*aml)*apkk1*(C3*C3)*Elep*t2 + 
           12.*(aml*aml)*apkk1*(C3*C3)*Enu*t2 - 
           2.*(amk*amk)*amLam*(C5*C5)*C6*t4 - 
           2.*(amk*amk)*amSig*(C5*C5)*C6*t4 - 
           4.*(amk*amk)*(aml*aml)*amLam*C3*(C5*C5)*C6*t4 - 
           4.*(amk*amk)*(aml*aml)*amSig*C3*(C5*C5)*C6*t4 + 
           2.*(amk*amk)*(aml*aml*aml*aml)*amLam*(C3*C3)*(C5*C5)*C6*
           t4 + 
           2.*(amk*amk)*(aml*aml*aml*aml)*amSig*(C3*C3)*(C5*C5)*C6*
           t4 + 
           4.*(aml*aml)*C3*(-2. + (aml*aml)*C3)*(C5*C5)*C6*
           (Ekaon*Ekaon*Ekaon)*t4 - 
           6.*(amk*amk)*C5*C6*Elep*t4 + 
           2.*(amk*amk)*(C5*C5)*C6*Elep*t4 + 
           (amk*amk)*(aml*aml*aml*aml)*(C3*C3)*(C5*C5)*C6*Elep*t4 + 
           6.*(amk*amk)*(aml*aml)*apkk1*(C3*C3)*(C5*C5)*C6*
           Elep*t4 - 
           2.*(amk*amk)*(aml*aml)*amLam*(C3*C3)*(C5*C5)*C6*
           (Elep*Elep)*t4 - 
           2.*(amk*amk)*(aml*aml)*amSig*(C3*C3)*(C5*C5)*C6*
           (Elep*Elep)*t4 - 6.*(amk*amk)*C5*C6*Enu*t4 - 
           2.*(amk*amk)*(C5*C5)*C6*Enu*t4 + 
           4.*(amk*amk)*(aml*aml)*C3*(C5*C5)*C6*Enu*t4 - 
           1.*(amk*amk)*(aml*aml*aml*aml)*(C3*C3)*(C5*C5)*C6*Enu*
           t4 - 
           6.*(amk*amk)*(aml*aml)*apkk1*(C3*C3)*(C5*C5)*C6*
           Enu*t4 + 
           4.*(amk*amk)*(aml*aml)*amLam*(C3*C3)*(C5*C5)*C6*
           Elep*Enu*t4 + 
           4.*(amk*amk)*(aml*aml)*amSig*(C3*C3)*(C5*C5)*C6*
           Elep*Enu*t4 - 
           2.*(amk*amk)*(aml*aml)*amLam*(C3*C3)*(C5*C5)*C6*
           (Enu*Enu)*t4 - 
           2.*(amk*amk)*(aml*aml)*amSig*(C3*C3)*(C5*C5)*C6*
           (Enu*Enu)*t4 - 
           6.*(amk*amk)*apkk1*C5*C6*Fm2*t4 + 
           4.*(amk*amk)*amLam*C5*C6*Elep*Fm2*t4 + 
           4.*(amk*amk)*amSig*C5*C6*Elep*Fm2*t4 + 
           4.*(amk*amk)*amLam*C5*C6*Enu*Fm2*t4 + 
           4.*(amk*amk)*amSig*C5*C6*Enu*Fm2*t4 + 
           6.*akpk*
           (-1.*(amk*amk)*C5*C6*Fm2*t4 + 
           (aml*aml)*(C3*C3)*(Elep - 1.*Enu)*
           (2.*t2 - 1.*(amk*amk)*(C5*C5)*C6*t4)) - 
           2.*(Ekaon*Ekaon)*
           ((aml*aml*aml*aml)*(C3*C3)*(C5*C5)*C6*
           (amLam + amSig + 2.*Elep - 2.*Enu)*
           t4 - 
           2.*C5*C6*
           (6.*Elep - 2.*C5*Elep + 6.*Enu + 
           2.*C5*Enu + 
           amLam*(C5 - 1.*(Elep + Enu)*Fm2) + 
           amSig*(C5 - 1.*(Elep + Enu)*Fm2))*
           t4 + 
           2.*(aml*aml)*C3*
           (-1.*(C5*C5)*C6*
           (amLam + amSig - 4.*Enu)*t4 + 
           2.*C3*(Elep - 1.*Enu)*
           (3.*t2 + 
           (-1.*akpk + apkk1)*(C5*C5)*C6*t4)))
           + Ekaon*
           (-24.*t1 + 12.*(aml*aml)*apkk1*(C3*C3)*t2 + 
           24.*(aml*aml)*(C3*C3)*(Elep*Elep)*t2 - 
           48.*(aml*aml)*(C3*C3)*Elep*Enu*t2 + 
           24.*(aml*aml)*(C3*C3)*(Enu*Enu)*t2 - 
           12.*apkk1*C5*C6*t4 + 
           4.*(amk*amk)*(C5*C5)*C6*t4 + 
           4.*apkk1*(C5*C5)*C6*t4 + 
           6.*(amk*amk)*(aml*aml)*C3*(C5*C5)*C6*t4 + 
           4.*(akpk*akpk)*(aml*aml)*(C3*C3)*(C5*C5)*C6*
           t4 - 
           3.*(amk*amk)*(aml*aml*aml*aml)*(C3*C3)*(C5*C5)*C6*t4 + 
           2.*(aml*aml*aml*aml)*apkk1*(C3*C3)*(C5*C5)*C6*t4 + 
           4.*(aml*aml)*(apkk1*apkk1)*(C3*C3)*(C5*C5)*C6*
           t4 - 
           8.*(amk*amk)*(aml*aml)*(C3*C3)*(C5*C5)*C6*
           (Elep*Elep)*t4 + 
           16.*(amk*amk)*(aml*aml)*(C3*C3)*(C5*C5)*C6*
           Elep*Enu*t4 - 
           8.*(amk*amk)*(aml*aml)*(C3*C3)*(C5*C5)*C6*
           (Enu*Enu)*t4 + 
           4.*(amk*amk)*C5*C6*Elep*Fm2*t4 + 
           4.*(amk*amk)*C5*C6*Enu*Fm2*t4 - 
           2.*akpk*
           ((aml*aml*aml*aml)*(C3*C3)*(C5*C5)*C6*t4 + 
           2.*C5*(3. + C5)*C6*t4 + 
           2.*(aml*aml)*C3*
           (3.*C3*t2 - 2.*(C5*C5)*C6*t4 + 
           2.*apkk1*C3*(C5*C5)*C6*t4))))) + 
           C5*C6*t4*
           (48.*((akpk + apkk1 + akpk*f - 1.*apkk1*f)*
           Fm2 + 
           Ekaon*
           (-1.*(Elep + Enu)*Fm2 + 
           3.*f*(1. + Elep*Fm2 - 1.*Enu*Fm2)))*
           t1 + 
           (C5*C5*C5)*C6*
           ((amk*amk)*
           (2.*Elep - 6.*akpk*(aml*aml)*(C3*C3)*Elep + 
           (aml*aml*aml*aml)*(C3*C3)*Elep + 
           6.*(aml*aml)*apkk1*(C3*C3)*Elep + 
           Ekaon*
           (4. - 3.*(aml*aml*aml*aml)*(C3*C3) + 
           2.*(aml*aml)*C3*
           (3. - 4.*C3*((Elep - 1.*Enu)*(Elep - 1.*Enu)))) + 
           4.*amLam*
           (-1. + (aml*aml*aml*aml)*(C3*C3) - 
           1.*(aml*aml)*C3*
           (2. + C3*((Elep - 1.*Enu)*(Elep - 1.*Enu)))) - 
           2.*Enu + 4.*(aml*aml)*C3*Enu + 
           6.*akpk*(aml*aml)*(C3*C3)*Enu - 
           1.*(aml*aml*aml*aml)*(C3*C3)*Enu - 
           6.*(aml*aml)*apkk1*(C3*C3)*Enu) + 
           2.*Ekaon*
           (2.*(akpk*akpk)*(aml*aml)*(C3*C3) + 
           2.*(aml*aml)*(apkk1*apkk1)*(C3*C3) + 
           apkk1*
           (2. + (aml*aml*aml*aml)*(C3*C3) + 
           4.*(aml*aml)*(C3*C3)*Ekaon*
           (-1.*Elep + Enu)) + 
           2.*Ekaon*
           (amLam*
           (2. + 2.*(aml*aml)*C3 - 
           1.*(aml*aml*aml*aml)*(C3*C3)) - 2.*Elep + 
           2.*Enu + 
           (aml*aml*aml*aml)*(C3*C3)*
           (Ekaon - 1.*Elep + Enu) - 
           2.*(aml*aml)*C3*(Ekaon + 2.*Enu)) - 
           1.*akpk*
           (2. + (aml*aml*aml*aml)*(C3*C3) + 
           4.*(aml*aml)*C3*
           (-1. + apkk1*C3 + 
           C3*Ekaon*(-1.*Elep + Enu)))))*t4 - 
           4.*(C5*C5)*C6*
           (3.*akpk*(2.*Ekaon + (amk*amk)*Fm2) + 
           3.*apkk1*(2.*Ekaon + (amk*amk)*Fm2) - 
           1.*(Elep + Enu)*
           (-4.*(Ekaon*Ekaon)*(-3. + amLam*Fm2) + 
           (amk*amk)*
           (-3. + 4.*amLam*Fm2 + 2.*Ekaon*Fm2))
           )*t4 + 
           C5*(-24.*(aml*aml)*apkk1*(C3*C3)*Elep*t2 + 
           24.*(aml*aml)*apkk1*(C3*C3)*Enu*t2 + 
           36.*(amk*amk)*amLam*C6*t4 + 
           18.*(amk*amk)*C6*Elep*t4 - 
           18.*(amk*amk)*C6*Enu*t4 - 
           18.*(amk*amk)*(aml*aml)*C6*Fm2*t4 + 
           12.*(amk*amk)*apkk1*C6*Fm2*t4 + 
           8.*(amk*amk)*(aml*aml)*amLam*C6*(Fm2*Fm2)*t4 + 
           4.*(aml*aml)*C6*(Ekaon*Ekaon*Ekaon)*(Fm2*Fm2)*t4 + 
           3.*(amk*amk)*(aml*aml)*C6*Elep*(Fm2*Fm2)*t4 + 
           12.*(amk*amk)*apkk1*C6*Elep*(Fm2*Fm2)*t4 - 
           8.*(amk*amk)*amLam*C6*(Elep*Elep)*(Fm2*Fm2)*t4 - 
           3.*(amk*amk)*(aml*aml)*C6*Enu*(Fm2*Fm2)*t4 - 
           8.*(amk*amk)*amLam*C6*(Enu*Enu)*(Fm2*Fm2)*t4 + 
           12.*akpk*
           (2.*(aml*aml)*(C3*C3)*(Elep - 1.*Enu)*t2 + 
           (amk*amk)*C6*Fm2*(-1. + Enu*Fm2)*t4) - 
           4.*(Ekaon*Ekaon)*
           (2.*C6*
           (9.*Elep - 9.*Enu + 
           2.*apkk1*Elep*(Fm2*Fm2) + 
           2.*akpk*Enu*(Fm2*Fm2) + 
           amLam*
           (9. + 6.*Elep*Fm2 - 6.*Enu*Fm2))*t4
           + 3.*(aml*aml)*
           (4.*(C3*C3)*(Elep - 1.*Enu)*t2 + 
           C6*(amLam + Elep - 1.*Enu)*(Fm2*Fm2)*
           t4)) - 
           1.*Ekaon*
           (48.*t1 - 24.*(aml*aml)*apkk1*(C3*C3)*t2 - 
           48.*(aml*aml)*(C3*C3)*(Elep*Elep)*t2 + 
           96.*(aml*aml)*(C3*C3)*Elep*Enu*t2 - 
           48.*(aml*aml)*(C3*C3)*(Enu*Enu)*t2 + 
           36.*(amk*amk)*C6*t4 - 36.*apkk1*C6*t4 + 
           72.*(amk*amk)*C6*Elep*Fm2*t4 - 
           72.*(amk*amk)*C6*Enu*Fm2*t4 - 
           8.*(akpk*akpk)*C6*(Fm2*Fm2)*t4 + 
           7.*(amk*amk)*(aml*aml)*C6*(Fm2*Fm2)*t4 - 
           6.*(aml*aml)*apkk1*C6*(Fm2*Fm2)*t4 - 
           8.*(apkk1*apkk1)*C6*(Fm2*Fm2)*t4 + 
           16.*(amk*amk)*C6*(Elep*Elep)*(Fm2*Fm2)*t4 + 
           16.*(amk*amk)*C6*(Enu*Enu)*(Fm2*Fm2)*t4 + 
           6.*akpk*
           (6.*C6*t4 + 
           (aml*aml)*(4.*(C3*C3)*t2 + C6*(Fm2*Fm2)*t4)
           ))))) - 
           144.*(amk*amk)*C2*C7*f*t1*t5 + 
           72.*(aml*aml)*C2*C7*f*t1*t5 - 
           36.*(amk*amk*amk*amk)*C1*C2*C4*C7*t3*t5 + 
           18.*(amk*amk)*(aml*aml)*C1*C2*C4*C7*t3*t5 - 
           36.*akpk*amSig*C1*C2*C4*C7*Ekaon*t3*t5 - 
           36.*(amk*amk)*amSig*C1*C2*C4*C7*Ekaon*t3*t5 + 
           18.*(aml*aml)*amSig*C1*C2*C4*C7*Ekaon*t3*t5 + 
           36.*amSig*apkk1*C1*C2*C4*C7*Ekaon*t3*t5 + 
           36.*akpk*amSig*C1*C2*(C4*C4)*C7*Ekaon*t3*t5 + 
           36.*amSig*apkk1*C1*C2*(C4*C4)*C7*Ekaon*t3*t5 - 
           36.*(amk*amk)*amSig*C1*C2*C4*C7*Elep*t3*t5 - 
           36.*(amk*amk)*amSig*C1*C2*(C4*C4)*C7*Elep*t3*t5 + 
           36.*(amk*amk)*amSig*C1*C2*C4*C7*Enu*t3*t5 - 
           36.*(amk*amk)*amSig*C1*C2*(C4*C4)*C7*Enu*t3*t5 + 
           18.*akpk*(aml*aml)*amSig*C1*C2*C4*C7*Fm1*t3*t5 + 
           36.*(amk*amk)*(aml*aml)*amSig*C1*C2*C4*C7*Fm1*t3*t5 + 
           144.*akpk*amSig*apkk1*C1*C2*C4*C7*Fm1*t3*t5 + 
           18.*(aml*aml)*amSig*apkk1*C1*C2*C4*C7*Fm1*t3*t5 - 
           36.*akpk*(amk*amk)*C1*C2*C4*C7*Ekaon*Fm1*t3*t5 + 
           36.*(amk*amk)*(aml*aml)*C1*C2*C4*C7*Ekaon*Fm1*t3*t5 + 
           36.*(amk*amk)*apkk1*C1*C2*C4*C7*Ekaon*Fm1*t3*t5 + 
           72.*akpk*(amk*amk)*C1*C2*C4*C7*Elep*Fm1*t3*t5 - 
           36.*(amk*amk*amk*amk)*C1*C2*C4*C7*Elep*Fm1*t3*t5 + 
           18.*(amk*amk)*(aml*aml)*C1*C2*C4*C7*Elep*Fm1*t3*t5 + 
           36.*(amk*amk*amk*amk)*C1*C2*C4*C7*Enu*Fm1*t3*t5 + 
           18.*(amk*amk)*(aml*aml)*C1*C2*C4*C7*Enu*Fm1*t3*t5 + 
           72.*(amk*amk)*apkk1*C1*C2*C4*C7*Enu*Fm1*t3*t5 + 
           36.*(amk*amk*amk*amk)*C2*C5*C6*C7*t4*t5 - 
           18.*(amk*amk)*(aml*aml)*C2*C5*C6*C7*t4*t5 + 
           36.*akpk*amLam*C2*C5*C6*C7*Ekaon*t4*t5 + 
           36.*(amk*amk)*amLam*C2*C5*C6*C7*Ekaon*t4*t5 - 
           18.*(aml*aml)*amLam*C2*C5*C6*C7*Ekaon*t4*t5 - 
           36.*amLam*apkk1*C2*C5*C6*C7*Ekaon*t4*t5 + 
           12.*akpk*amLam*C2*(C5*C5)*C6*C7*Ekaon*t4*t5 + 
           12.*amLam*apkk1*C2*(C5*C5)*C6*C7*Ekaon*t4*t5 + 
           36.*(amk*amk)*amLam*C2*C5*C6*C7*Elep*t4*t5 - 
           12.*(amk*amk)*amLam*C2*(C5*C5)*C6*C7*Elep*t4*t5 - 
           36.*(amk*amk)*amLam*C2*C5*C6*C7*Enu*t4*t5 - 
           12.*(amk*amk)*amLam*C2*(C5*C5)*C6*C7*Enu*t4*t5 - 
           6.*akpk*(aml*aml)*amLam*C2*C5*C6*C7*Fm2*t4*t5 - 
           12.*(amk*amk)*(aml*aml)*amLam*C2*C5*C6*C7*Fm2*t4*t5 - 
           48.*akpk*amLam*apkk1*C2*C5*C6*C7*Fm2*t4*t5 - 
           6.*(aml*aml)*amLam*apkk1*C2*C5*C6*C7*Fm2*t4*t5 + 
           12.*akpk*(amk*amk)*C2*C5*C6*C7*Ekaon*Fm2*t4*t5 - 
           12.*(amk*amk)*(aml*aml)*C2*C5*C6*C7*Ekaon*Fm2*t4*t5 - 
           12.*(amk*amk)*apkk1*C2*C5*C6*C7*Ekaon*Fm2*t4*t5 - 
           24.*akpk*(amk*amk)*C2*C5*C6*C7*Elep*Fm2*t4*t5 + 
           12.*(amk*amk*amk*amk)*C2*C5*C6*C7*Elep*Fm2*t4*t5 - 
           6.*(amk*amk)*(aml*aml)*C2*C5*C6*C7*Elep*Fm2*t4*t5 - 
           12.*(amk*amk*amk*amk)*C2*C5*C6*C7*Enu*Fm2*t4*t5 - 
           6.*(amk*amk)*(aml*aml)*C2*C5*C6*C7*Enu*Fm2*t4*t5 - 
           24.*(amk*amk)*apkk1*C2*C5*C6*C7*Enu*Fm2*t4*t5 + 
           36.*(amk*amk)*(C2*C2)*(C7*C7)*Ekaon*(t5*t5) - 
           9.*(aml*aml)*(C2*C2)*(C7*C7)*Ekaon*(t5*t5) + 
           36.*(amk*amk)*(C2*C2)*(C7*C7)*Elep*(t5*t5) - 
           9.*(aml*aml)*(C2*C2)*(C7*C7)*Elep*(t5*t5) - 
           36.*(amk*amk)*(C2*C2)*(C7*C7)*Enu*(t5*t5) + 
           9.*(aml*aml)*(C2*C2)*(C7*C7)*Enu*(t5*t5) + 
           144.*(amk*amk)*C8*C9*f*t1*t6 - 
           72.*(aml*aml)*C8*C9*f*t1*t6 + 
           36.*(amk*amk*amk*amk)*C1*C4*C8*C9*t3*t6 - 
           18.*(amk*amk)*(aml*aml)*C1*C4*C8*C9*t3*t6 + 
           36.*akpk*amSig*C1*C4*C8*C9*Ekaon*t3*t6 + 
           36.*(amk*amk)*amSig*C1*C4*C8*C9*Ekaon*t3*t6 - 
           18.*(aml*aml)*amSig*C1*C4*C8*C9*Ekaon*t3*t6 - 
           36.*amSig*apkk1*C1*C4*C8*C9*Ekaon*t3*t6 - 
           36.*akpk*amSig*C1*(C4*C4)*C8*C9*Ekaon*t3*t6 - 
           36.*amSig*apkk1*C1*(C4*C4)*C8*C9*Ekaon*t3*t6 + 
           36.*(amk*amk)*amSig*C1*C4*C8*C9*Elep*t3*t6 + 
           36.*(amk*amk)*amSig*C1*(C4*C4)*C8*C9*Elep*t3*t6 - 
           36.*(amk*amk)*amSig*C1*C4*C8*C9*Enu*t3*t6 + 
           36.*(amk*amk)*amSig*C1*(C4*C4)*C8*C9*Enu*t3*t6 - 
           18.*akpk*(aml*aml)*amSig*C1*C4*C8*C9*Fm1*t3*t6 - 
           36.*(amk*amk)*(aml*aml)*amSig*C1*C4*C8*C9*Fm1*t3*t6 - 
           144.*akpk*amSig*apkk1*C1*C4*C8*C9*Fm1*t3*t6 - 
           18.*(aml*aml)*amSig*apkk1*C1*C4*C8*C9*Fm1*t3*t6 + 
           36.*akpk*(amk*amk)*C1*C4*C8*C9*Ekaon*Fm1*t3*t6 - 
           36.*(amk*amk)*(aml*aml)*C1*C4*C8*C9*Ekaon*Fm1*t3*t6 - 
           36.*(amk*amk)*apkk1*C1*C4*C8*C9*Ekaon*Fm1*t3*t6 - 
           72.*akpk*(amk*amk)*C1*C4*C8*C9*Elep*Fm1*t3*t6 + 
           36.*(amk*amk*amk*amk)*C1*C4*C8*C9*Elep*Fm1*t3*t6 - 
           18.*(amk*amk)*(aml*aml)*C1*C4*C8*C9*Elep*Fm1*t3*t6 - 
           36.*(amk*amk*amk*amk)*C1*C4*C8*C9*Enu*Fm1*t3*t6 - 
           18.*(amk*amk)*(aml*aml)*C1*C4*C8*C9*Enu*Fm1*t3*t6 - 
           72.*(amk*amk)*apkk1*C1*C4*C8*C9*Enu*Fm1*t3*t6 - 
           36.*(amk*amk*amk*amk)*C5*C6*C8*C9*t4*t6 + 
           18.*(amk*amk)*(aml*aml)*C5*C6*C8*C9*t4*t6 - 
           36.*akpk*amLam*C5*C6*C8*C9*Ekaon*t4*t6 - 
           36.*(amk*amk)*amLam*C5*C6*C8*C9*Ekaon*t4*t6 + 
           18.*(aml*aml)*amLam*C5*C6*C8*C9*Ekaon*t4*t6 + 
           36.*amLam*apkk1*C5*C6*C8*C9*Ekaon*t4*t6 - 
           12.*akpk*amLam*(C5*C5)*C6*C8*C9*Ekaon*t4*t6 - 
           12.*amLam*apkk1*(C5*C5)*C6*C8*C9*Ekaon*t4*t6 - 
           36.*(amk*amk)*amLam*C5*C6*C8*C9*Elep*t4*t6 + 
           12.*(amk*amk)*amLam*(C5*C5)*C6*C8*C9*Elep*t4*t6 + 
           36.*(amk*amk)*amLam*C5*C6*C8*C9*Enu*t4*t6 + 
           12.*(amk*amk)*amLam*(C5*C5)*C6*C8*C9*Enu*t4*t6 + 
           6.*akpk*(aml*aml)*amLam*C5*C6*C8*C9*Fm2*t4*t6 + 
           12.*(amk*amk)*(aml*aml)*amLam*C5*C6*C8*C9*Fm2*t4*t6 + 
           48.*akpk*amLam*apkk1*C5*C6*C8*C9*Fm2*t4*t6 + 
           6.*(aml*aml)*amLam*apkk1*C5*C6*C8*C9*Fm2*t4*t6 - 
           12.*akpk*(amk*amk)*C5*C6*C8*C9*Ekaon*Fm2*t4*t6 + 
           12.*(amk*amk)*(aml*aml)*C5*C6*C8*C9*Ekaon*Fm2*t4*t6 + 
           12.*(amk*amk)*apkk1*C5*C6*C8*C9*Ekaon*Fm2*t4*t6 + 
           24.*akpk*(amk*amk)*C5*C6*C8*C9*Elep*Fm2*t4*t6 - 
           12.*(amk*amk*amk*amk)*C5*C6*C8*C9*Elep*Fm2*t4*t6 + 
           6.*(amk*amk)*(aml*aml)*C5*C6*C8*C9*Elep*Fm2*t4*t6 + 
           12.*(amk*amk*amk*amk)*C5*C6*C8*C9*Enu*Fm2*t4*t6 + 
           6.*(amk*amk)*(aml*aml)*C5*C6*C8*C9*Enu*Fm2*t4*t6 + 
           24.*(amk*amk)*apkk1*C5*C6*C8*C9*Enu*Fm2*t4*t6 - 
           72.*(amk*amk)*C2*C7*C8*C9*Ekaon*t5*t6 + 
           18.*(aml*aml)*C2*C7*C8*C9*Ekaon*t5*t6 - 
           72.*(amk*amk)*C2*C7*C8*C9*Elep*t5*t6 + 
           18.*(aml*aml)*C2*C7*C8*C9*Elep*t5*t6 + 
           72.*(amk*amk)*C2*C7*C8*C9*Enu*t5*t6 - 
           18.*(aml*aml)*C2*C7*C8*C9*Enu*t5*t6 + 
           36.*(amk*amk)*(C8*C8)*(C9*C9)*Ekaon*(t6*t6) - 
           9.*(aml*aml)*(C8*C8)*(C9*C9)*Ekaon*(t6*t6) + 
           36.*(amk*amk)*(C8*C8)*(C9*C9)*Elep*(t6*t6) - 
           9.*(aml*aml)*(C8*C8)*(C9*C9)*Elep*(t6*t6) - 
           36.*(amk*amk)*(C8*C8)*(C9*C9)*Enu*(t6*t6) + 
           9.*(aml*aml)*(C8*C8)*(C9*C9)*Enu*(t6*t6) + 
           2.*am*(144.*(-1. + (f*f))*(t1*t1) + 
           18.*akpk*(amk*amk)*(C1*C1)*(C4*C4)*(t3*t3) + 
           9.*(amk*amk*amk*amk)*(C1*C1)*(C4*C4)*(t3*t3) + 
           9.*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4)*(t3*t3) - 
           18.*(amk*amk)*apkk1*(C1*C1)*(C4*C4)*(t3*t3) - 
           36.*akpk*(amk*amk)*(C1*C1)*(C4*C4*C4)*(t3*t3) - 
           36.*(amk*amk)*apkk1*(C1*C1)*(C4*C4*C4)*(t3*t3) + 
           18.*akpk*(amk*amk)*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) - 
           9.*(amk*amk*amk*amk)*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) - 
           9.*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) - 
           18.*(amk*amk)*apkk1*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) + 
           18.*(amk*amk)*amSig*(C1*C1)*(C4*C4)*Ekaon*(t3*t3) - 
           18.*(amk*amk)*amSig*(C1*C1)*(C4*C4*C4*C4)*Ekaon*(t3*t3) + 
           18.*(amk*amk)*amSig*(C1*C1)*(C4*C4)*Elep*(t3*t3) + 
           36.*(amk*amk)*amSig*(C1*C1)*(C4*C4*C4)*Elep*(t3*t3) + 
           18.*(amk*amk)*amSig*(C1*C1)*(C4*C4*C4*C4)*Elep*(t3*t3) + 
           36.*(amk*amk)*(C1*C1)*(C4*C4)*Ekaon*Elep*(t3*t3) + 
           72.*(amk*amk)*(C1*C1)*(C4*C4*C4)*Ekaon*Elep*(t3*t3) + 
           36.*(amk*amk)*(C1*C1)*(C4*C4*C4*C4)*Ekaon*Elep*(t3*t3) - 
           18.*(amk*amk)*amSig*(C1*C1)*(C4*C4)*Enu*(t3*t3) + 
           36.*(amk*amk)*amSig*(C1*C1)*(C4*C4*C4)*Enu*(t3*t3) - 
           18.*(amk*amk)*amSig*(C1*C1)*(C4*C4*C4*C4)*Enu*(t3*t3) - 
           36.*(amk*amk)*(C1*C1)*(C4*C4)*Ekaon*Enu*(t3*t3) + 
           72.*(amk*amk)*(C1*C1)*(C4*C4*C4)*Ekaon*Enu*(t3*t3) - 
           36.*(amk*amk)*(C1*C1)*(C4*C4*C4*C4)*Ekaon*Enu*(t3*t3) + 
           72.*akpk*(amk*amk)*amSig*(C1*C1)*(C4*C4*C4)*Fm1*(t3*t3) + 
           72.*(amk*amk)*amSig*apkk1*(C1*C1)*(C4*C4*C4)*Fm1*
           (t3*t3) - 
           36.*akpk*(amk*amk)*(C1*C1)*(C4*C4)*Ekaon*Fm1*(t3*t3) - 
           36.*akpk*(amSig*amSig)*(C1*C1)*(C4*C4)*Ekaon*Fm1*
           (t3*t3) + 
           36.*(amk*amk)*apkk1*(C1*C1)*(C4*C4)*Ekaon*Fm1*
           (t3*t3) + 
           36.*(amSig*amSig)*apkk1*(C1*C1)*(C4*C4)*Ekaon*Fm1*
           (t3*t3) - 
           36.*akpk*(amk*amk)*(C1*C1)*(C4*C4*C4)*Ekaon*Fm1*(t3*t3) + 
           36.*akpk*(amSig*amSig)*(C1*C1)*(C4*C4*C4)*Ekaon*Fm1*
           (t3*t3) - 
           36.*(amk*amk)*apkk1*(C1*C1)*(C4*C4*C4)*Ekaon*Fm1*
           (t3*t3) + 
           36.*(amSig*amSig)*apkk1*(C1*C1)*(C4*C4*C4)*Ekaon*Fm1*
           (t3*t3) - 
           72.*akpk*amSig*(C1*C1)*(C4*C4)*(Ekaon*Ekaon)*Fm1*
           (t3*t3) + 
           72.*amSig*apkk1*(C1*C1)*(C4*C4)*(Ekaon*Ekaon)*Fm1*
           (t3*t3) - 
           72.*akpk*amSig*(C1*C1)*(C4*C4*C4)*(Ekaon*Ekaon)*Fm1*
           (t3*t3) - 
           72.*amSig*apkk1*(C1*C1)*(C4*C4*C4)*(Ekaon*Ekaon)*Fm1*
           (t3*t3) + 
           36.*(amk*amk*amk*amk)*(C1*C1)*(C4*C4)*Elep*Fm1*(t3*t3) - 
           36.*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4*C4)*Elep*Fm1*
           (t3*t3) + 
           36.*(amk*amk)*amSig*(C1*C1)*(C4*C4)*Ekaon*Elep*Fm1*
           (t3*t3) - 
           36.*(amk*amk)*amSig*(C1*C1)*(C4*C4*C4)*Ekaon*Elep*Fm1*
           (t3*t3) - 
           36.*(amk*amk*amk*amk)*(C1*C1)*(C4*C4)*Enu*Fm1*(t3*t3) - 
           36.*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4*C4)*Enu*Fm1*
           (t3*t3) - 
           36.*(amk*amk)*amSig*(C1*C1)*(C4*C4)*Ekaon*Enu*Fm1*
           (t3*t3) - 
           36.*(amk*amk)*amSig*(C1*C1)*(C4*C4*C4)*Ekaon*Enu*Fm1*
           (t3*t3) - 
           36.*(akpk*akpk)*(amk*amk)*(C1*C1)*(C4*C4)*(Fm1*Fm1)*(t3*t3) - 
           36.*(amk*amk)*(apkk1*apkk1)*(C1*C1)*(C4*C4)*(Fm1*Fm1)*
           (t3*t3) + 
           72.*(amk*amk)*amSig*apkk1*(C1*C1)*(C4*C4)*Elep*
           (Fm1*Fm1)*(t3*t3) + 
           72.*(amk*amk)*apkk1*(C1*C1)*(C4*C4)*Ekaon*Elep*
           (Fm1*Fm1)*(t3*t3) + 
           36.*(amSig*amSig)*apkk1*(C1*C1)*(C4*C4)*Ekaon*Elep*
           (Fm1*Fm1)*(t3*t3) + 
           18.*(amk*amk*amk*amk)*(C1*C1)*(C4*C4)*(Elep*Elep)*(Fm1*Fm1)*(t3*t3) - 
           18.*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4)*(Elep*Elep)*
           (Fm1*Fm1)*(t3*t3) + 
           72.*akpk*(amk*amk)*amSig*(C1*C1)*(C4*C4)*Enu*(Fm1*Fm1)*
           (t3*t3) + 
           72.*akpk*(amk*amk)*(C1*C1)*(C4*C4)*Ekaon*Enu*(Fm1*Fm1)*
           (t3*t3) + 
           36.*akpk*(amSig*amSig)*(C1*C1)*(C4*C4)*Ekaon*Enu*
           (Fm1*Fm1)*(t3*t3) + 
           18.*(amk*amk*amk*amk)*(C1*C1)*(C4*C4)*(Enu*Enu)*(Fm1*Fm1)*(t3*t3) - 
           18.*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4)*(Enu*Enu)*(Fm1*Fm1)*
           (t3*t3) - 36.*akpk*(amk*amk)*C1*C4*C5*C6*t3*t4 - 
           18.*(amk*amk*amk*amk)*C1*C4*C5*C6*t3*t4 - 
           18.*(amk*amk)*amLam*amSig*C1*C4*C5*C6*t3*t4 + 
           36.*(amk*amk)*apkk1*C1*C4*C5*C6*t3*t4 + 
           36.*akpk*(amk*amk)*C1*(C4*C4)*C5*C6*t3*t4 + 
           36.*(amk*amk)*apkk1*C1*(C4*C4)*C5*C6*t3*t4 - 
           12.*akpk*(amk*amk)*C1*C4*(C5*C5)*C6*t3*t4 - 
           12.*(amk*amk)*apkk1*C1*C4*(C5*C5)*C6*t3*t4 + 
           12.*akpk*(amk*amk)*C1*(C4*C4)*(C5*C5)*C6*t3*t4 - 
           6.*(amk*amk*amk*amk)*C1*(C4*C4)*(C5*C5)*C6*t3*t4 - 
           6.*(amk*amk)*amLam*amSig*C1*(C4*C4)*(C5*C5)*C6*t3*
           t4 - 12.*(amk*amk)*apkk1*C1*(C4*C4)*(C5*C5)*C6*t3*
           t4 - 18.*(amk*amk)*amLam*C1*C4*C5*C6*Ekaon*t3*
           t4 - 18.*(amk*amk)*amSig*C1*C4*C5*C6*Ekaon*t3*
           t4 - 6.*(amk*amk)*amLam*C1*(C4*C4)*(C5*C5)*C6*Ekaon*
           t3*t4 - 
           6.*(amk*amk)*amSig*C1*(C4*C4)*(C5*C5)*C6*Ekaon*t3*
           t4 - 18.*(amk*amk)*amLam*C1*C4*C5*C6*Elep*t3*
           t4 - 18.*(amk*amk)*amSig*C1*C4*C5*C6*Elep*t3*
           t4 - 18.*(amk*amk)*amLam*C1*(C4*C4)*C5*C6*Elep*t3*
           t4 - 18.*(amk*amk)*amSig*C1*(C4*C4)*C5*C6*Elep*t3*
           t4 + 6.*(amk*amk)*amLam*C1*C4*(C5*C5)*C6*Elep*t3*
           t4 + 6.*(amk*amk)*amSig*C1*C4*(C5*C5)*C6*Elep*t3*
           t4 + 6.*(amk*amk)*amLam*C1*(C4*C4)*(C5*C5)*C6*Elep*
           t3*t4 + 
           6.*(amk*amk)*amSig*C1*(C4*C4)*(C5*C5)*C6*Elep*t3*
           t4 - 72.*(amk*amk)*C1*C4*C5*C6*Ekaon*Elep*t3*
           t4 - 72.*(amk*amk)*C1*(C4*C4)*C5*C6*Ekaon*Elep*t3*
           t4 + 24.*(amk*amk)*C1*C4*(C5*C5)*C6*Ekaon*Elep*t3*
           t4 + 24.*(amk*amk)*C1*(C4*C4)*(C5*C5)*C6*Ekaon*Elep*
           t3*t4 + 
           18.*(amk*amk)*amLam*C1*C4*C5*C6*Enu*t3*t4 + 
           18.*(amk*amk)*amSig*C1*C4*C5*C6*Enu*t3*t4 - 
           18.*(amk*amk)*amLam*C1*(C4*C4)*C5*C6*Enu*t3*t4 - 
           18.*(amk*amk)*amSig*C1*(C4*C4)*C5*C6*Enu*t3*t4 + 
           6.*(amk*amk)*amLam*C1*C4*(C5*C5)*C6*Enu*t3*t4 + 
           6.*(amk*amk)*amSig*C1*C4*(C5*C5)*C6*Enu*t3*t4 - 
           6.*(amk*amk)*amLam*C1*(C4*C4)*(C5*C5)*C6*Enu*t3*t4 - 
           6.*(amk*amk)*amSig*C1*(C4*C4)*(C5*C5)*C6*Enu*t3*t4 + 
           72.*(amk*amk)*C1*C4*C5*C6*Ekaon*Enu*t3*t4 - 
           72.*(amk*amk)*C1*(C4*C4)*C5*C6*Ekaon*Enu*t3*t4 + 
           24.*(amk*amk)*C1*C4*(C5*C5)*C6*Ekaon*Enu*t3*t4 - 
           24.*(amk*amk)*C1*(C4*C4)*(C5*C5)*C6*Ekaon*Enu*t3*
           t4 + 12.*akpk*(amk*amk)*amLam*C1*C4*(C5*C5)*C6*
           Fm1*t3*t4 + 
           12.*akpk*(amk*amk)*amSig*C1*C4*(C5*C5)*C6*Fm1*t3*
           t4 + 12.*(amk*amk)*amLam*apkk1*C1*C4*(C5*C5)*C6*
           Fm1*t3*t4 + 
           12.*(amk*amk)*amSig*apkk1*C1*C4*(C5*C5)*C6*Fm1*t3*
           t4 + 36.*akpk*(amk*amk)*C1*C4*C5*C6*Ekaon*Fm1*
           t3*t4 + 
           36.*akpk*amLam*amSig*C1*C4*C5*C6*Ekaon*Fm1*t3*
           t4 - 36.*(amk*amk)*apkk1*C1*C4*C5*C6*Ekaon*Fm1*
           t3*t4 - 
           36.*amLam*amSig*apkk1*C1*C4*C5*C6*Ekaon*Fm1*
           t3*t4 - 
           12.*akpk*(amk*amk)*C1*C4*(C5*C5)*C6*Ekaon*Fm1*t3*
           t4 + 12.*akpk*amLam*amSig*C1*C4*(C5*C5)*C6*
           Ekaon*Fm1*t3*t4 - 
           12.*(amk*amk)*apkk1*C1*C4*(C5*C5)*C6*Ekaon*Fm1*t3*
           t4 + 12.*amLam*amSig*apkk1*C1*C4*(C5*C5)*C6*
           Ekaon*Fm1*t3*t4 + 
           36.*akpk*amLam*C1*C4*C5*C6*(Ekaon*Ekaon)*Fm1*t3*
           t4 + 36.*akpk*amSig*C1*C4*C5*C6*(Ekaon*Ekaon)*Fm1*
           t3*t4 - 
           36.*amLam*apkk1*C1*C4*C5*C6*(Ekaon*Ekaon)*Fm1*t3*
           t4 - 36.*amSig*apkk1*C1*C4*C5*C6*(Ekaon*Ekaon)*
           Fm1*t3*t4 - 
           12.*akpk*amLam*C1*C4*(C5*C5)*C6*(Ekaon*Ekaon)*Fm1*t3*
           t4 - 12.*akpk*amSig*C1*C4*(C5*C5)*C6*(Ekaon*Ekaon)*
           Fm1*t3*t4 - 
           12.*amLam*apkk1*C1*C4*(C5*C5)*C6*(Ekaon*Ekaon)*Fm1*
           t3*t4 - 
           12.*amSig*apkk1*C1*C4*(C5*C5)*C6*(Ekaon*Ekaon)*Fm1*
           t3*t4 - 
           36.*(amk*amk*amk*amk)*C1*C4*C5*C6*Elep*Fm1*t3*t4 - 
           12.*(amk*amk)*amLam*amSig*C1*C4*(C5*C5)*C6*Elep*
           Fm1*t3*t4 - 
           18.*(amk*amk)*amLam*C1*C4*C5*C6*Ekaon*Elep*Fm1*
           t3*t4 - 
           18.*(amk*amk)*amSig*C1*C4*C5*C6*Ekaon*Elep*Fm1*
           t3*t4 + 
           24.*akpk*amLam*C1*C4*(C5*C5)*C6*Ekaon*Elep*Fm1*
           t3*t4 - 
           6.*(amk*amk)*amLam*C1*C4*(C5*C5)*C6*Ekaon*Elep*Fm1*
           t3*t4 - 
           24.*akpk*amSig*C1*C4*(C5*C5)*C6*Ekaon*Elep*Fm1*
           t3*t4 - 
           6.*(amk*amk)*amSig*C1*C4*(C5*C5)*C6*Ekaon*Elep*Fm1*
           t3*t4 + 
           36.*(amk*amk*amk*amk)*C1*C4*C5*C6*Enu*Fm1*t3*t4 - 
           12.*(amk*amk)*amLam*amSig*C1*C4*(C5*C5)*C6*Enu*Fm1*
           t3*t4 + 
           18.*(amk*amk)*amLam*C1*C4*C5*C6*Ekaon*Enu*Fm1*t3*
           t4 + 18.*(amk*amk)*amSig*C1*C4*C5*C6*Ekaon*Enu*
           Fm1*t3*t4 - 
           6.*(amk*amk)*amLam*C1*C4*(C5*C5)*C6*Ekaon*Enu*Fm1*
           t3*t4 - 
           6.*(amk*amk)*amSig*C1*C4*(C5*C5)*C6*Ekaon*Enu*Fm1*
           t3*t4 - 
           24.*amLam*apkk1*C1*C4*(C5*C5)*C6*Ekaon*Enu*Fm1*
           t3*t4 + 
           24.*amSig*apkk1*C1*C4*(C5*C5)*C6*Ekaon*Enu*Fm1*
           t3*t4 - 
           12.*akpk*(amk*amk)*amLam*C1*(C4*C4)*C5*C6*Fm2*t3*
           t4 - 12.*akpk*(amk*amk)*amSig*C1*(C4*C4)*C5*C6*
           Fm2*t3*t4 - 
           12.*(amk*amk)*amLam*apkk1*C1*(C4*C4)*C5*C6*Fm2*t3*
           t4 - 12.*(amk*amk)*amSig*apkk1*C1*(C4*C4)*C5*C6*
           Fm2*t3*t4 + 
           12.*akpk*(amk*amk)*C1*C4*C5*C6*Ekaon*Fm2*t3*t4 + 
           12.*akpk*amLam*amSig*C1*C4*C5*C6*Ekaon*Fm2*t3*
           t4 - 12.*(amk*amk)*apkk1*C1*C4*C5*C6*Ekaon*Fm2*
           t3*t4 - 
           12.*amLam*amSig*apkk1*C1*C4*C5*C6*Ekaon*Fm2*
           t3*t4 + 
           12.*akpk*(amk*amk)*C1*(C4*C4)*C5*C6*Ekaon*Fm2*t3*
           t4 - 12.*akpk*amLam*amSig*C1*(C4*C4)*C5*C6*
           Ekaon*Fm2*t3*t4 + 
           12.*(amk*amk)*apkk1*C1*(C4*C4)*C5*C6*Ekaon*Fm2*t3*
           t4 - 12.*amLam*amSig*apkk1*C1*(C4*C4)*C5*C6*
           Ekaon*Fm2*t3*t4 + 
           12.*akpk*amLam*C1*C4*C5*C6*(Ekaon*Ekaon)*Fm2*t3*
           t4 + 12.*akpk*amSig*C1*C4*C5*C6*(Ekaon*Ekaon)*Fm2*
           t3*t4 - 
           12.*amLam*apkk1*C1*C4*C5*C6*(Ekaon*Ekaon)*Fm2*t3*
           t4 - 12.*amSig*apkk1*C1*C4*C5*C6*(Ekaon*Ekaon)*
           Fm2*t3*t4 + 
           12.*akpk*amLam*C1*(C4*C4)*C5*C6*(Ekaon*Ekaon)*Fm2*t3*
           t4 + 12.*akpk*amSig*C1*(C4*C4)*C5*C6*(Ekaon*Ekaon)*
           Fm2*t3*t4 + 
           12.*amLam*apkk1*C1*(C4*C4)*C5*C6*(Ekaon*Ekaon)*Fm2*
           t3*t4 + 
           12.*amSig*apkk1*C1*(C4*C4)*C5*C6*(Ekaon*Ekaon)*Fm2*
           t3*t4 - 
           12.*(amk*amk*amk*amk)*C1*C4*C5*C6*Elep*Fm2*t3*t4 + 
           12.*(amk*amk)*amLam*amSig*C1*(C4*C4)*C5*C6*Elep*
           Fm2*t3*t4 - 
           6.*(amk*amk)*amLam*C1*C4*C5*C6*Ekaon*Elep*Fm2*t3*
           t4 - 6.*(amk*amk)*amSig*C1*C4*C5*C6*Ekaon*Elep*
           Fm2*t3*t4 + 
           24.*akpk*amLam*C1*(C4*C4)*C5*C6*Ekaon*Elep*Fm2*
           t3*t4 + 
           6.*(amk*amk)*amLam*C1*(C4*C4)*C5*C6*Ekaon*Elep*Fm2*
           t3*t4 - 
           24.*akpk*amSig*C1*(C4*C4)*C5*C6*Ekaon*Elep*Fm2*
           t3*t4 + 
           6.*(amk*amk)*amSig*C1*(C4*C4)*C5*C6*Ekaon*Elep*Fm2*
           t3*t4 + 
           12.*(amk*amk*amk*amk)*C1*C4*C5*C6*Enu*Fm2*t3*t4 + 
           12.*(amk*amk)*amLam*amSig*C1*(C4*C4)*C5*C6*Enu*Fm2*
           t3*t4 + 
           6.*(amk*amk)*amLam*C1*C4*C5*C6*Ekaon*Enu*Fm2*t3*
           t4 + 6.*(amk*amk)*amSig*C1*C4*C5*C6*Ekaon*Enu*
           Fm2*t3*t4 + 
           6.*(amk*amk)*amLam*C1*(C4*C4)*C5*C6*Ekaon*Enu*Fm2*
           t3*t4 + 
           6.*(amk*amk)*amSig*C1*(C4*C4)*C5*C6*Ekaon*Enu*Fm2*
           t3*t4 - 
           24.*amLam*apkk1*C1*(C4*C4)*C5*C6*Ekaon*Enu*Fm2*
           t3*t4 + 
           24.*amSig*apkk1*C1*(C4*C4)*C5*C6*Ekaon*Enu*Fm2*
           t3*t4 + 
           24.*(akpk*akpk)*(amk*amk)*C1*C4*C5*C6*Fm1*Fm2*t3*
           t4 + 24.*(amk*amk)*(apkk1*apkk1)*C1*C4*C5*C6*Fm1*Fm2*
           t3*t4 - 
           24.*(amk*amk)*amLam*apkk1*C1*C4*C5*C6*Elep*Fm1*
           Fm2*t3*t4 - 
           24.*(amk*amk)*amSig*apkk1*C1*C4*C5*C6*Elep*Fm1*
           Fm2*t3*t4 - 
           48.*(amk*amk)*apkk1*C1*C4*C5*C6*Ekaon*Elep*Fm1*
           Fm2*t3*t4 - 
           24.*amLam*amSig*apkk1*C1*C4*C5*C6*Ekaon*Elep*
           Fm1*Fm2*t3*t4 - 
           12.*(amk*amk*amk*amk)*C1*C4*C5*C6*(Elep*Elep)*Fm1*Fm2*t3*
           t4 + 12.*(amk*amk)*amLam*amSig*C1*C4*C5*C6*
           (Elep*Elep)*Fm1*Fm2*t3*t4 - 
           24.*akpk*(amk*amk)*amLam*C1*C4*C5*C6*Enu*Fm1*Fm2*
           t3*t4 - 
           24.*akpk*(amk*amk)*amSig*C1*C4*C5*C6*Enu*Fm1*Fm2*
           t3*t4 - 
           48.*akpk*(amk*amk)*C1*C4*C5*C6*Ekaon*Enu*Fm1*Fm2*
           t3*t4 - 
           24.*akpk*amLam*amSig*C1*C4*C5*C6*Ekaon*Enu*
           Fm1*Fm2*t3*t4 - 
           12.*(amk*amk*amk*amk)*C1*C4*C5*C6*(Enu*Enu)*Fm1*Fm2*t3*t4 + 
           12.*(amk*amk)*amLam*amSig*C1*C4*C5*C6*(Enu*Enu)*Fm1*
           Fm2*t3*t4 + 
           18.*akpk*(amk*amk)*(C5*C5)*(C6*C6)*(t4*t4) + 
           9.*(amk*amk*amk*amk)*(C5*C5)*(C6*C6)*(t4*t4) + 
           9.*(amk*amk)*(amLam*amLam)*(C5*C5)*(C6*C6)*(t4*t4) - 
           18.*(amk*amk)*apkk1*(C5*C5)*(C6*C6)*(t4*t4) + 
           12.*akpk*(amk*amk)*(C5*C5*C5)*(C6*C6)*(t4*t4) + 
           12.*(amk*amk)*apkk1*(C5*C5*C5)*(C6*C6)*(t4*t4) + 
           2.*akpk*(amk*amk)*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) - 
           1.*(amk*amk*amk*amk)*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) - 
           1.*(amk*amk)*(amLam*amLam)*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) - 
           2.*(amk*amk)*apkk1*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) + 
           18.*(amk*amk)*amLam*(C5*C5)*(C6*C6)*Ekaon*(t4*t4) - 
           2.*(amk*amk)*amLam*(C5*C5*C5*C5)*(C6*C6)*Ekaon*(t4*t4) + 
           18.*(amk*amk)*amLam*(C5*C5)*(C6*C6)*Elep*(t4*t4) - 
           12.*(amk*amk)*amLam*(C5*C5*C5)*(C6*C6)*Elep*(t4*t4) + 
           2.*(amk*amk)*amLam*(C5*C5*C5*C5)*(C6*C6)*Elep*(t4*t4) + 
           36.*(amk*amk)*(C5*C5)*(C6*C6)*Ekaon*Elep*(t4*t4) - 
           24.*(amk*amk)*(C5*C5*C5)*(C6*C6)*Ekaon*Elep*(t4*t4) + 
           4.*(amk*amk)*(C5*C5*C5*C5)*(C6*C6)*Ekaon*Elep*(t4*t4) - 
           18.*(amk*amk)*amLam*(C5*C5)*(C6*C6)*Enu*(t4*t4) - 
           12.*(amk*amk)*amLam*(C5*C5*C5)*(C6*C6)*Enu*(t4*t4) - 
           2.*(amk*amk)*amLam*(C5*C5*C5*C5)*(C6*C6)*Enu*(t4*t4) - 
           36.*(amk*amk)*(C5*C5)*(C6*C6)*Ekaon*Enu*(t4*t4) - 
           24.*(amk*amk)*(C5*C5*C5)*(C6*C6)*Ekaon*Enu*(t4*t4) - 
           4.*(amk*amk)*(C5*C5*C5*C5)*(C6*C6)*Ekaon*Enu*(t4*t4) - 
           8.*akpk*(amk*amk)*amLam*(C5*C5*C5)*(C6*C6)*Fm2*(t4*t4) - 
           8.*(amk*amk)*amLam*apkk1*(C5*C5*C5)*(C6*C6)*Fm2*(t4*t4) - 
           12.*akpk*(amk*amk)*(C5*C5)*(C6*C6)*Ekaon*Fm2*(t4*t4) - 
           12.*akpk*(amLam*amLam)*(C5*C5)*(C6*C6)*Ekaon*Fm2*
           (t4*t4) + 
           12.*(amk*amk)*apkk1*(C5*C5)*(C6*C6)*Ekaon*Fm2*
           (t4*t4) + 
           12.*(amLam*amLam)*apkk1*(C5*C5)*(C6*C6)*Ekaon*Fm2*
           (t4*t4) + 
           4.*akpk*(amk*amk)*(C5*C5*C5)*(C6*C6)*Ekaon*Fm2*(t4*t4) - 
           4.*akpk*(amLam*amLam)*(C5*C5*C5)*(C6*C6)*Ekaon*Fm2*
           (t4*t4) + 
           4.*(amk*amk)*apkk1*(C5*C5*C5)*(C6*C6)*Ekaon*Fm2*(t4*t4) - 
           4.*(amLam*amLam)*apkk1*(C5*C5*C5)*(C6*C6)*Ekaon*Fm2*
           (t4*t4) - 
           24.*akpk*amLam*(C5*C5)*(C6*C6)*(Ekaon*Ekaon)*Fm2*
           (t4*t4) + 
           24.*amLam*apkk1*(C5*C5)*(C6*C6)*(Ekaon*Ekaon)*Fm2*
           (t4*t4) + 
           8.*akpk*amLam*(C5*C5*C5)*(C6*C6)*(Ekaon*Ekaon)*Fm2*
           (t4*t4) + 
           8.*amLam*apkk1*(C5*C5*C5)*(C6*C6)*(Ekaon*Ekaon)*Fm2*
           (t4*t4) + 
           12.*(amk*amk*amk*amk)*(C5*C5)*(C6*C6)*Elep*Fm2*(t4*t4) + 
           4.*(amk*amk)*(amLam*amLam)*(C5*C5*C5)*(C6*C6)*Elep*Fm2*
           (t4*t4) + 
           12.*(amk*amk)*amLam*(C5*C5)*(C6*C6)*Ekaon*Elep*Fm2*
           (t4*t4) + 
           4.*(amk*amk)*amLam*(C5*C5*C5)*(C6*C6)*Ekaon*Elep*Fm2*
           (t4*t4) - 
           12.*(amk*amk*amk*amk)*(C5*C5)*(C6*C6)*Enu*Fm2*(t4*t4) + 
           4.*(amk*amk)*(amLam*amLam)*(C5*C5*C5)*(C6*C6)*Enu*Fm2*
           (t4*t4) - 
           12.*(amk*amk)*amLam*(C5*C5)*(C6*C6)*Ekaon*Enu*Fm2*
           (t4*t4) + 
           4.*(amk*amk)*amLam*(C5*C5*C5)*(C6*C6)*Ekaon*Enu*Fm2*
           (t4*t4) - 
           4.*(akpk*akpk)*(amk*amk)*(C5*C5)*(C6*C6)*(Fm2*Fm2)*(t4*t4) - 
           4.*(amk*amk)*(apkk1*apkk1)*(C5*C5)*(C6*C6)*(Fm2*Fm2)*(t4*t4) + 
           8.*(amk*amk)*amLam*apkk1*(C5*C5)*(C6*C6)*Elep*(Fm2*Fm2)*
           (t4*t4) + 
           8.*(amk*amk)*apkk1*(C5*C5)*(C6*C6)*Ekaon*Elep*(Fm2*Fm2)*
           (t4*t4) + 
           4.*(amLam*amLam)*apkk1*(C5*C5)*(C6*C6)*Ekaon*Elep*
           (Fm2*Fm2)*(t4*t4) + 
           2.*(amk*amk*amk*amk)*(C5*C5)*(C6*C6)*(Elep*Elep)*(Fm2*Fm2)*(t4*t4) - 
           2.*(amk*amk)*(amLam*amLam)*(C5*C5)*(C6*C6)*(Elep*Elep)*(Fm2*Fm2)*
           (t4*t4) + 
           8.*akpk*(amk*amk)*amLam*(C5*C5)*(C6*C6)*Enu*(Fm2*Fm2)*
           (t4*t4) + 
           8.*akpk*(amk*amk)*(C5*C5)*(C6*C6)*Ekaon*Enu*(Fm2*Fm2)*
           (t4*t4) + 
           4.*akpk*(amLam*amLam)*(C5*C5)*(C6*C6)*Ekaon*Enu*(Fm2*Fm2)*
           (t4*t4) + 
           2.*(amk*amk*amk*amk)*(C5*C5)*(C6*C6)*(Enu*Enu)*(Fm2*Fm2)*(t4*t4) - 
           2.*(amk*amk)*(amLam*amLam)*(C5*C5)*(C6*C6)*(Enu*Enu)*(Fm2*Fm2)*
           (t4*t4) + 
           12.*t1*(6.*(amk*amk)*C1*(C4*C4)*t3 + 
           3.*(amk*amk)*(aml*aml)*C1*C3*(C4*C4)*t3 + 
           6.*amSig*C1*(C4*C4)*Ekaon*t3 - 
           6.*(aml*aml)*C1*C3*(C4*C4)*(Ekaon*Ekaon)*t3 - 
           12.*C1*C4*Ekaon*Elep*t3 - 
           12.*C1*(C4*C4)*Ekaon*Elep*t3 - 
           12.*C1*C4*Ekaon*Enu*t3 + 
           12.*C1*(C4*C4)*Ekaon*Enu*t3 - 
           12.*(aml*aml)*C1*C3*(C4*C4)*Ekaon*Enu*t3 + 
           6.*(amk*amk)*C1*C4*f*t3 + 
           6.*amSig*C1*C4*Ekaon*f*t3 + 
           12.*C1*C4*Ekaon*Elep*f*t3 + 
           12.*C1*(C4*C4)*Ekaon*Elep*f*t3 - 
           12.*C1*C4*Ekaon*Enu*f*t3 + 
           12.*C1*(C4*C4)*Ekaon*Enu*f*t3 + 
           6.*amSig*C1*C4*Ekaon*Elep*Fm1*t3 + 
           6.*amSig*C1*C4*Ekaon*Enu*Fm1*t3 + 
           9.*(aml*aml)*C1*C4*Ekaon*f*Fm1*t3 + 
           12.*(amk*amk)*C1*C4*Elep*f*Fm1*t3 + 
           6.*amSig*C1*C4*Ekaon*Elep*f*Fm1*t3 - 
           12.*(amk*amk)*C1*C4*Enu*f*Fm1*t3 - 
           6.*amSig*C1*C4*Ekaon*Enu*f*Fm1*t3 + 
           2.*(amk*amk)*(C5*C5)*C6*t4 + 
           (amk*amk)*(aml*aml)*C3*(C5*C5)*C6*t4 + 
           2.*amLam*(C5*C5)*C6*Ekaon*t4 - 
           2.*(aml*aml)*C3*(C5*C5)*C6*(Ekaon*Ekaon)*t4 + 
           12.*C5*C6*Ekaon*Elep*t4 - 
           4.*(C5*C5)*C6*Ekaon*Elep*t4 + 
           12.*C5*C6*Ekaon*Enu*t4 + 
           4.*(C5*C5)*C6*Ekaon*Enu*t4 - 
           4.*(aml*aml)*C3*(C5*C5)*C6*Ekaon*Enu*t4 - 
           6.*(amk*amk)*C5*C6*f*t4 - 
           6.*amLam*C5*C6*Ekaon*f*t4 - 
           12.*C5*C6*Ekaon*Elep*f*t4 + 
           4.*(C5*C5)*C6*Ekaon*Elep*f*t4 + 
           12.*C5*C6*Ekaon*Enu*f*t4 + 
           4.*(C5*C5)*C6*Ekaon*Enu*f*t4 - 
           2.*amLam*C5*C6*Ekaon*Elep*Fm2*t4 - 
           2.*amLam*C5*C6*Ekaon*Enu*Fm2*t4 - 
           3.*(aml*aml)*C5*C6*Ekaon*f*Fm2*t4 - 
           4.*(amk*amk)*C5*C6*Elep*f*Fm2*t4 - 
           2.*amLam*C5*C6*Ekaon*Elep*f*Fm2*t4 + 
           4.*(amk*amk)*C5*C6*Enu*f*Fm2*t4 + 
           2.*amLam*C5*C6*Ekaon*Enu*f*Fm2*t4 + 
           2.*akpk*
           (3.*C1*C4*
           (1. + C4*(-1. + (aml*aml)*C3 - 1.*f) - 
           1.*amSig*Fm1 + Ekaon*Fm1 - 
           2.*Elep*Fm1 - 
           1.*f*(-1. + amSig*Fm1 + Ekaon*Fm1))*
           t3 + 
           C5*C6*
           (-3. + C5*(-1. + (aml*aml)*C3 - 1.*f) + 
           amLam*Fm2 - 1.*Ekaon*Fm2 + 
           2.*Elep*Fm2 + 
           f*(-3. + amLam*Fm2 + Ekaon*Fm2))*t4)
           - 2.*apkk1*
           (3.*C1*C4*
           (-1. + C4*(-1. + f) + amSig*Fm1 - 
           1.*Ekaon*Fm1 - 2.*Enu*Fm1 - 
           1.*f*(-1. + amSig*Fm1 + Ekaon*Fm1))*
           t3 + 
           C5*C6*
           (3. + C5*(-1. + f) - 1.*amLam*Fm2 + 
           Ekaon*Fm2 + 2.*Enu*Fm2 + 
           f*(-3. + amLam*Fm2 + Ekaon*Fm2))*t4))
           + (aml*aml*aml*aml)*(C3*C3)*
           (18.*apkk1*C1*(C4*C4)*t2*t3 + 
           72.*C1*(C4*C4)*(Ekaon*Ekaon)*t2*t3 - 
           36.*C1*(C4*C4)*Ekaon*Elep*t2*t3 + 
           36.*C1*(C4*C4)*Ekaon*Enu*t2*t3 - 
           9.*(amSig*amSig)*(C1*C1)*(C4*C4*C4*C4)*(Ekaon*Ekaon)*(t3*t3) + 
           6.*apkk1*(C5*C5)*C6*t2*t4 + 
           24.*(C5*C5)*C6*(Ekaon*Ekaon)*t2*t4 - 
           12.*(C5*C5)*C6*Ekaon*Elep*t2*t4 + 
           12.*(C5*C5)*C6*Ekaon*Enu*t2*t4 - 
           6.*amLam*amSig*C1*(C4*C4)*(C5*C5)*C6*(Ekaon*Ekaon)*
           t3*t4 - 
           1.*(amLam*amLam)*(C5*C5*C5*C5)*(C6*C6)*(Ekaon*Ekaon)*(t4*t4) + 
           (amk*amk*amk*amk)*((3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4)*(3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4)) + 
           akpk*(3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4)*
           (-6.*t2 + 
           (amk*amk)*(3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4))
           + (amk*amk)*
           (9.*(C1*C1)*(C4*C4*C4*C4)*
           ((amSig*amSig) - 1.*apkk1 - 
           1.*amSig*(Ekaon - 1.*Elep + Enu) - 
           2.*Ekaon*(Ekaon - 1.*Elep + Enu))*
           (t3*t3) + 
           (C5*C5)*C6*t4*
           (-12.*t2 - 
           1.*(C5*C5)*C6*
           (-1.*(amLam*amLam) + apkk1 + 
           amLam*(Ekaon - 1.*Elep + Enu) + 
           2.*Ekaon*(Ekaon - 1.*Elep + Enu))*
           t4) - 
           3.*C1*(C4*C4)*t3*
           (12.*t2 + 
           (C5*C5)*C6*
           (2.*apkk1 + 
           (amSig + 4.*Ekaon)*
           (Ekaon - 1.*Elep + Enu) + 
           amLam*
           (-2.*amSig + Ekaon - 1.*Elep + Enu)
           )*t4))) - 
           18.*akpk*C1*C2*C4*C7*Ekaon*t3*t5 + 
           18.*(amk*amk)*C1*C2*C4*C7*Ekaon*t3*t5 + 
           18.*apkk1*C1*C2*C4*C7*Ekaon*t3*t5 + 
           18.*akpk*C1*C2*(C4*C4)*C7*Ekaon*t3*t5 + 
           18.*apkk1*C1*C2*(C4*C4)*C7*Ekaon*t3*t5 - 
           18.*(amk*amk)*C1*C2*C4*C7*Elep*t3*t5 - 
           18.*(amk*amk)*C1*C2*(C4*C4)*C7*Elep*t3*t5 + 
           18.*(amk*amk)*C1*C2*C4*C7*Enu*t3*t5 - 
           18.*(amk*amk)*C1*C2*(C4*C4)*C7*Enu*t3*t5 + 
           72.*akpk*apkk1*C1*C2*C4*C7*Fm1*t3*t5 + 
           36.*akpk*C1*C2*C4*C7*(Ekaon*Ekaon)*Fm1*t3*t5 - 
           36.*apkk1*C1*C2*C4*C7*(Ekaon*Ekaon)*Fm1*t3*t5 - 
           72.*akpk*C1*C2*C4*C7*Ekaon*Elep*Fm1*t3*t5 + 
           36.*(amk*amk)*C1*C2*C4*C7*Ekaon*Elep*Fm1*t3*t5 - 
           36.*(amk*amk)*C1*C2*C4*C7*Ekaon*Enu*Fm1*t3*t5 - 
           72.*apkk1*C1*C2*C4*C7*Ekaon*Enu*Fm1*t3*t5 + 
           18.*akpk*C2*C5*C6*C7*Ekaon*t4*t5 - 
           18.*(amk*amk)*C2*C5*C6*C7*Ekaon*t4*t5 - 
           18.*apkk1*C2*C5*C6*C7*Ekaon*t4*t5 + 
           6.*akpk*C2*(C5*C5)*C6*C7*Ekaon*t4*t5 + 
           6.*apkk1*C2*(C5*C5)*C6*C7*Ekaon*t4*t5 + 
           18.*(amk*amk)*C2*C5*C6*C7*Elep*t4*t5 - 
           6.*(amk*amk)*C2*(C5*C5)*C6*C7*Elep*t4*t5 - 
           18.*(amk*amk)*C2*C5*C6*C7*Enu*t4*t5 - 
           6.*(amk*amk)*C2*(C5*C5)*C6*C7*Enu*t4*t5 - 
           24.*akpk*apkk1*C2*C5*C6*C7*Fm2*t4*t5 - 
           12.*akpk*C2*C5*C6*C7*(Ekaon*Ekaon)*Fm2*t4*t5 + 
           12.*apkk1*C2*C5*C6*C7*(Ekaon*Ekaon)*Fm2*t4*t5 + 
           24.*akpk*C2*C5*C6*C7*Ekaon*Elep*Fm2*t4*t5 - 
           12.*(amk*amk)*C2*C5*C6*C7*Ekaon*Elep*Fm2*t4*t5 + 
           12.*(amk*amk)*C2*C5*C6*C7*Ekaon*Enu*Fm2*t4*t5 + 
           24.*apkk1*C2*C5*C6*C7*Ekaon*Enu*Fm2*t4*t5 + 
           18.*akpk*C1*C4*C8*C9*Ekaon*t3*t6 - 
           18.*(amk*amk)*C1*C4*C8*C9*Ekaon*t3*t6 - 
           18.*apkk1*C1*C4*C8*C9*Ekaon*t3*t6 - 
           18.*akpk*C1*(C4*C4)*C8*C9*Ekaon*t3*t6 - 
           18.*apkk1*C1*(C4*C4)*C8*C9*Ekaon*t3*t6 + 
           18.*(amk*amk)*C1*C4*C8*C9*Elep*t3*t6 + 
           18.*(amk*amk)*C1*(C4*C4)*C8*C9*Elep*t3*t6 - 
           18.*(amk*amk)*C1*C4*C8*C9*Enu*t3*t6 + 
           18.*(amk*amk)*C1*(C4*C4)*C8*C9*Enu*t3*t6 - 
           72.*akpk*apkk1*C1*C4*C8*C9*Fm1*t3*t6 - 
           36.*akpk*C1*C4*C8*C9*(Ekaon*Ekaon)*Fm1*t3*t6 + 
           36.*apkk1*C1*C4*C8*C9*(Ekaon*Ekaon)*Fm1*t3*t6 + 
           72.*akpk*C1*C4*C8*C9*Ekaon*Elep*Fm1*t3*t6 - 
           36.*(amk*amk)*C1*C4*C8*C9*Ekaon*Elep*Fm1*t3*t6 + 
           36.*(amk*amk)*C1*C4*C8*C9*Ekaon*Enu*Fm1*t3*t6 + 
           72.*apkk1*C1*C4*C8*C9*Ekaon*Enu*Fm1*t3*t6 - 
           18.*akpk*C5*C6*C8*C9*Ekaon*t4*t6 + 
           18.*(amk*amk)*C5*C6*C8*C9*Ekaon*t4*t6 + 
           18.*apkk1*C5*C6*C8*C9*Ekaon*t4*t6 - 
           6.*akpk*(C5*C5)*C6*C8*C9*Ekaon*t4*t6 - 
           6.*apkk1*(C5*C5)*C6*C8*C9*Ekaon*t4*t6 - 
           18.*(amk*amk)*C5*C6*C8*C9*Elep*t4*t6 + 
           6.*(amk*amk)*(C5*C5)*C6*C8*C9*Elep*t4*t6 + 
           18.*(amk*amk)*C5*C6*C8*C9*Enu*t4*t6 + 
           6.*(amk*amk)*(C5*C5)*C6*C8*C9*Enu*t4*t6 + 
           24.*akpk*apkk1*C5*C6*C8*C9*Fm2*t4*t6 + 
           12.*akpk*C5*C6*C8*C9*(Ekaon*Ekaon)*Fm2*t4*t6 - 
           12.*apkk1*C5*C6*C8*C9*(Ekaon*Ekaon)*Fm2*t4*t6 - 
           24.*akpk*C5*C6*C8*C9*Ekaon*Elep*Fm2*t4*t6 + 
           12.*(amk*amk)*C5*C6*C8*C9*Ekaon*Elep*Fm2*t4*t6 - 
           12.*(amk*amk)*C5*C6*C8*C9*Ekaon*Enu*Fm2*t4*t6 - 
           24.*apkk1*C5*C6*C8*C9*Ekaon*Enu*Fm2*t4*t6 + 
           (aml*aml)*(-54.*(amk*amk)*amSig*(C1*C1)*(C4*C4)*Fm1*
           (t3*t3) + 
           54.*(amk*amk)*(C1*C1)*(C4*C4)*Ekaon*Fm1*(t3*t3) + 
           108.*amSig*(C1*C1)*(C4*C4)*(Ekaon*Ekaon)*Fm1*
           (t3*t3) + 
           27.*akpk*(amk*amk)*(C1*C1)*(C4*C4)*(Fm1*Fm1)*(t3*t3) + 
           18.*(amk*amk*amk*amk)*(C1*C1)*(C4*C4)*(Fm1*Fm1)*(t3*t3) + 
           18.*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4)*(Fm1*Fm1)*
           (t3*t3) - 
           27.*(amk*amk)*apkk1*(C1*C1)*(C4*C4)*(Fm1*Fm1)*
           (t3*t3) + 
           9.*(amk*amk)*amSig*(C1*C1)*(C4*C4)*Ekaon*(Fm1*Fm1)*
           (t3*t3) - 
           18.*(amk*amk)*(C1*C1)*(C4*C4)*(Ekaon*Ekaon)*(Fm1*Fm1)*
           (t3*t3) - 
           9.*(amSig*amSig)*(C1*C1)*(C4*C4)*(Ekaon*Ekaon)*(Fm1*Fm1)*
           (t3*t3) + 
           27.*(amk*amk)*amSig*(C1*C1)*(C4*C4)*Elep*(Fm1*Fm1)*
           (t3*t3) + 
           54.*(amk*amk)*(C1*C1)*(C4*C4)*Ekaon*Elep*(Fm1*Fm1)*
           (t3*t3) - 
           27.*(amk*amk)*amSig*(C1*C1)*(C4*C4)*Enu*(Fm1*Fm1)*
           (t3*t3) - 
           54.*(amk*amk)*(C1*C1)*(C4*C4)*Ekaon*Enu*(Fm1*Fm1)*
           (t3*t3) + 
           27.*(amk*amk)*amLam*C1*C4*C5*C6*Fm1*t3*t4 + 
           27.*(amk*amk)*amSig*C1*C4*C5*C6*Fm1*t3*t4 - 
           54.*(amk*amk)*C1*C4*C5*C6*Ekaon*Fm1*t3*t4 - 
           54.*amLam*C1*C4*C5*C6*(Ekaon*Ekaon)*Fm1*t3*t4 - 
           54.*amSig*C1*C4*C5*C6*(Ekaon*Ekaon)*Fm1*t3*t4 + 
           9.*(amk*amk)*amLam*C1*C4*C5*C6*Fm2*t3*t4 + 
           9.*(amk*amk)*amSig*C1*C4*C5*C6*Fm2*t3*t4 - 
           18.*(amk*amk)*C1*C4*C5*C6*Ekaon*Fm2*t3*t4 - 
           18.*amLam*C1*C4*C5*C6*(Ekaon*Ekaon)*Fm2*t3*t4 - 
           18.*amSig*C1*C4*C5*C6*(Ekaon*Ekaon)*Fm2*t3*t4 - 
           18.*akpk*(amk*amk)*C1*C4*C5*C6*Fm1*Fm2*t3*
           t4 - 
           12.*(amk*amk*amk*amk)*C1*C4*C5*C6*Fm1*Fm2*t3*t4 - 
           12.*(amk*amk)*amLam*amSig*C1*C4*C5*C6*Fm1*Fm2*
           t3*t4 + 
           18.*(amk*amk)*apkk1*C1*C4*C5*C6*Fm1*Fm2*t3*
           t4 - 
           3.*(amk*amk)*amLam*C1*C4*C5*C6*Ekaon*Fm1*Fm2*
           t3*t4 - 
           3.*(amk*amk)*amSig*C1*C4*C5*C6*Ekaon*Fm1*Fm2*
           t3*t4 + 
           12.*(amk*amk)*C1*C4*C5*C6*(Ekaon*Ekaon)*Fm1*Fm2*t3*
           t4 + 
           6.*amLam*amSig*C1*C4*C5*C6*(Ekaon*Ekaon)*Fm1*
           Fm2*t3*t4 - 
           9.*(amk*amk)*amLam*C1*C4*C5*C6*Elep*Fm1*Fm2*
           t3*t4 - 
           9.*(amk*amk)*amSig*C1*C4*C5*C6*Elep*Fm1*Fm2*
           t3*t4 - 
           36.*(amk*amk)*C1*C4*C5*C6*Ekaon*Elep*Fm1*Fm2*
           t3*t4 + 
           9.*(amk*amk)*amLam*C1*C4*C5*C6*Enu*Fm1*Fm2*t3*
           t4 + 
           9.*(amk*amk)*amSig*C1*C4*C5*C6*Enu*Fm1*Fm2*t3*
           t4 + 
           36.*(amk*amk)*C1*C4*C5*C6*Ekaon*Enu*Fm1*Fm2*
           t3*t4 - 
           18.*(amk*amk)*amLam*(C5*C5)*(C6*C6)*Fm2*(t4*t4) + 
           18.*(amk*amk)*(C5*C5)*(C6*C6)*Ekaon*Fm2*(t4*t4) + 
           36.*amLam*(C5*C5)*(C6*C6)*(Ekaon*Ekaon)*Fm2*(t4*t4) + 
           3.*akpk*(amk*amk)*(C5*C5)*(C6*C6)*(Fm2*Fm2)*(t4*t4) + 
           2.*(amk*amk*amk*amk)*(C5*C5)*(C6*C6)*(Fm2*Fm2)*(t4*t4) + 
           2.*(amk*amk)*(amLam*amLam)*(C5*C5)*(C6*C6)*(Fm2*Fm2)*
           (t4*t4) - 
           3.*(amk*amk)*apkk1*(C5*C5)*(C6*C6)*(Fm2*Fm2)*(t4*t4) + 
           (amk*amk)*amLam*(C5*C5)*(C6*C6)*Ekaon*(Fm2*Fm2)*
           (t4*t4) - 
           2.*(amk*amk)*(C5*C5)*(C6*C6)*(Ekaon*Ekaon)*(Fm2*Fm2)*
           (t4*t4) - 
           1.*(amLam*amLam)*(C5*C5)*(C6*C6)*(Ekaon*Ekaon)*(Fm2*Fm2)*
           (t4*t4) + 
           3.*(amk*amk)*amLam*(C5*C5)*(C6*C6)*Elep*(Fm2*Fm2)*
           (t4*t4) + 
           6.*(amk*amk)*(C5*C5)*(C6*C6)*Ekaon*Elep*(Fm2*Fm2)*
           (t4*t4) - 
           3.*(amk*amk)*amLam*(C5*C5)*(C6*C6)*Enu*(Fm2*Fm2)*
           (t4*t4) - 
           6.*(amk*amk)*(C5*C5)*(C6*C6)*Ekaon*Enu*(Fm2*Fm2)*
           (t4*t4) - 
           2.*C3*
           ((amk*amk*amk*amk)*
           ((3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4)*(3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4)) + 
           2.*akpk*(3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4)*
           (-3.*t2 + 
           (amk*amk)*(3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4)
           ) + 
           Ekaon*
           (-9.*(amSig*amSig)*(C1*C1)*(C4*C4*C4*C4)*Ekaon*
           (t3*t3) + 
           (C5*C5)*C6*t4*
           (12.*Ekaon*t2 + 12.*Enu*t2 - 
           1.*(amLam*amLam)*(C5*C5)*C6*Ekaon*t4) + 
           6.*C1*(C4*C4)*t3*
           (6.*Ekaon*t2 + 6.*Enu*t2 - 
           1.*amLam*amSig*(C5*C5)*C6*Ekaon*t4))
           + (amk*amk)*
           (9.*(C1*C1)*(C4*C4*C4*C4)*
           ((amSig*amSig) - 
           1.*amSig*(Ekaon + 2.*Enu) - 
           2.*Ekaon*(Ekaon + 2.*Enu))*(t3*t3) + 
           (C5*C5)*C6*t4*
           (-6.*t2 + 
           (C5*C5)*C6*
           ((amLam*amLam) - 
           1.*amLam*(Ekaon + 2.*Enu) - 
           2.*Ekaon*(Ekaon + 2.*Enu))*t4) - 
           3.*C1*(C4*C4)*t3*
           (6.*t2 + 
           (C5*C5)*C6*
           ((amSig + 4.*Ekaon)*
           (Ekaon + 2.*Enu) + 
           amLam*(-2.*amSig + Ekaon + 2.*Enu))
           *t4))) + 
           (C3*C3)*
           (36.*(Ekaon*Ekaon)*(t2*t2) + 36.*(Enu*Enu)*(t2*t2) + 
           36.*(akpk*akpk)*C1*(C4*C4)*t2*t3 + 
           18.*akpk*(amk*amk)*C1*(C4*C4)*t2*t3 - 
           72.*akpk*apkk1*C1*(C4*C4)*t2*t3 - 
           18.*(amk*amk)*apkk1*C1*(C4*C4)*t2*t3 + 
           36.*(apkk1*apkk1)*C1*(C4*C4)*t2*t3 - 
           36.*akpk*amSig*C1*(C4*C4)*Enu*t2*t3 + 
           36.*amSig*apkk1*C1*(C4*C4)*Enu*t2*t3 - 
           36.*(amk*amk)*C1*(C4*C4)*(Enu*Enu)*t2*t3 - 
           18.*(akpk*akpk)*(amk*amk)*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) + 
           36.*akpk*(amk*amk)*apkk1*(C1*C1)*(C4*C4*C4*C4)*
           (t3*t3) - 
           18.*(amk*amk)*(apkk1*apkk1)*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) + 
           36.*akpk*(amk*amk)*amSig*(C1*C1)*(C4*C4*C4*C4)*Enu*
           (t3*t3) - 
           36.*(amk*amk)*amSig*apkk1*(C1*C1)*(C4*C4*C4*C4)*Enu*
           (t3*t3) + 
           9.*(amk*amk*amk*amk)*(C1*C1)*(C4*C4*C4*C4)*(Enu*Enu)*(t3*t3) - 
           9.*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4*C4*C4)*(Enu*Enu)*
           (t3*t3) + 12.*(akpk*akpk)*(C5*C5)*C6*t2*t4 + 
           6.*akpk*(amk*amk)*(C5*C5)*C6*t2*t4 - 
           24.*akpk*apkk1*(C5*C5)*C6*t2*t4 - 
           6.*(amk*amk)*apkk1*(C5*C5)*C6*t2*t4 + 
           12.*(apkk1*apkk1)*(C5*C5)*C6*t2*t4 - 
           12.*akpk*amLam*(C5*C5)*C6*Enu*t2*t4 + 
           12.*amLam*apkk1*(C5*C5)*C6*Enu*t2*t4 - 
           12.*(amk*amk)*(C5*C5)*C6*(Enu*Enu)*t2*t4 - 
           12.*(akpk*akpk)*(amk*amk)*C1*(C4*C4)*(C5*C5)*C6*t3*
           t4 + 
           24.*akpk*(amk*amk)*apkk1*C1*(C4*C4)*(C5*C5)*C6*
           t3*t4 - 
           12.*(amk*amk)*(apkk1*apkk1)*C1*(C4*C4)*(C5*C5)*C6*
           t3*t4 + 
           12.*akpk*(amk*amk)*amLam*C1*(C4*C4)*(C5*C5)*C6*
           Enu*t3*t4 + 
           12.*akpk*(amk*amk)*amSig*C1*(C4*C4)*(C5*C5)*C6*
           Enu*t3*t4 - 
           12.*(amk*amk)*amLam*apkk1*C1*(C4*C4)*(C5*C5)*
           C6*Enu*t3*t4 - 
           12.*(amk*amk)*amSig*apkk1*C1*(C4*C4)*(C5*C5)*
           C6*Enu*t3*t4 + 
           6.*(amk*amk*amk*amk)*C1*(C4*C4)*(C5*C5)*C6*(Enu*Enu)*t3*
           t4 - 
           6.*(amk*amk)*amLam*amSig*C1*(C4*C4)*(C5*C5)*C6*
           (Enu*Enu)*t3*t4 - 
           2.*(akpk*akpk)*(amk*amk)*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) + 
           4.*akpk*(amk*amk)*apkk1*(C5*C5*C5*C5)*(C6*C6)*
           (t4*t4) - 
           2.*(amk*amk)*(apkk1*apkk1)*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) + 
           4.*akpk*(amk*amk)*amLam*(C5*C5*C5*C5)*(C6*C6)*Enu*
           (t4*t4) - 
           4.*(amk*amk)*amLam*apkk1*(C5*C5*C5*C5)*(C6*C6)*Enu*
           (t4*t4) + 
           (amk*amk*amk*amk)*(C5*C5*C5*C5)*(C6*C6)*(Enu*Enu)*(t4*t4) - 
           1.*(amk*amk)*(amLam*amLam)*(C5*C5*C5*C5)*(C6*C6)*(Enu*Enu)*
           (t4*t4) + 
           (Elep*Elep)*
           (36.*(t2*t2) - 
           12.*(amk*amk)*t2*
           (3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4) + 
           (amk*amk)*
           ((amk*amk)*
           ((3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4)*(3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4))
           - 1.*
           ((3.*amSig*C1*(C4*C4)*t3 + amLam*(C5*C5)*C6*t4)*(3.*amSig*C1*(C4*C4)*t3 + amLam*(C5*C5)*C6*t4)))) - 
           2.*Elep*
           (2.*(akpk - 1.*apkk1)*
           (3.*amSig*C1*(C4*C4)*t3 + 
           amLam*(C5*C5)*C6*t4)*
           (-3.*t2 + 
           (amk*amk)*
           (3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4)) + 
           Enu*
           (36.*(t2*t2) - 
           12.*(amk*amk)*t2*
           (3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4) + 
           (amk*amk)*
           ((amk*amk)*
           ((3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4)*(3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4))
           - 1.*
           ((3.*amSig*C1*(C4*C4)*t3 + amLam*(C5*C5)*C6*t4) * 
           (3.*amSig*C1*(C4*C4)*t3 + amLam*(C5*C5)*C6*t4))))) - 
           2.*Ekaon*
           (6.*(akpk - 1.*apkk1)*t2*
           (3.*amSig*C1*(C4*C4)*t3 + 
           amLam*(C5*C5)*C6*t4) + 
           Enu*
           (-36.*(t2*t2) + 
           12.*(akpk + (amk*amk) - 1.*apkk1)*t2*
           (3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4) - 
           1.*(akpk - 1.*apkk1)*
           (2.*(amk*amk)*
           ((3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4)*(3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4))
           + ((3.*amSig*C1*(C4*C4)*t3 + amLam*(C5*C5)*C6*t4)*(3.*amSig*C1*(C4*C4)*t3 + amLam*(C5*C5)*C6*t4)))) + 
           Elep*
           (36.*(t2*t2) - 
           12.*(akpk + (amk*amk) - 1.*apkk1)*t2*
           (3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4) + 
           (akpk - 1.*apkk1)*
           (2.*(amk*amk)*
           ((3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4)*(3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4))
           + ((3.*amSig*C1*(C4*C4)*t3 + amLam*(C5*C5)*C6*t4)*(3.*amSig*C1*(C4*C4)*t3 + amLam*(C5*C5)*C6*t4)))))) - 
           9.*C1*C2*C4*C7*Ekaon*t3*t5 + 
           9.*akpk*C1*C2*C4*C7*Fm1*t3*t5 + 
           18.*(amk*amk)*C1*C2*C4*C7*Fm1*t3*t5 + 
           9.*apkk1*C1*C2*C4*C7*Fm1*t3*t5 - 
           36.*C1*C2*C4*C7*(Ekaon*Ekaon)*Fm1*t3*t5 - 
           18.*C1*C2*C4*C7*Ekaon*Elep*Fm1*t3*t5 - 
           18.*C1*C2*C4*C7*Ekaon*Enu*Fm1*t3*t5 + 
           9.*C2*C5*C6*C7*Ekaon*t4*t5 - 
           3.*akpk*C2*C5*C6*C7*Fm2*t4*t5 - 
           6.*(amk*amk)*C2*C5*C6*C7*Fm2*t4*t5 - 
           3.*apkk1*C2*C5*C6*C7*Fm2*t4*t5 + 
           12.*C2*C5*C6*C7*(Ekaon*Ekaon)*Fm2*t4*t5 + 
           6.*C2*C5*C6*C7*Ekaon*Elep*Fm2*t4*t5 + 
           6.*C2*C5*C6*C7*Ekaon*Enu*Fm2*t4*t5 + 
           9.*C1*C4*C8*C9*Ekaon*t3*t6 - 
           9.*akpk*C1*C4*C8*C9*Fm1*t3*t6 - 
           18.*(amk*amk)*C1*C4*C8*C9*Fm1*t3*t6 - 
           9.*apkk1*C1*C4*C8*C9*Fm1*t3*t6 + 
           36.*C1*C4*C8*C9*(Ekaon*Ekaon)*Fm1*t3*t6 + 
           18.*C1*C4*C8*C9*Ekaon*Elep*Fm1*t3*t6 + 
           18.*C1*C4*C8*C9*Ekaon*Enu*Fm1*t3*t6 - 
           9.*C5*C6*C8*C9*Ekaon*t4*t6 + 
           3.*akpk*C5*C6*C8*C9*Fm2*t4*t6 + 
           6.*(amk*amk)*C5*C6*C8*C9*Fm2*t4*t6 + 
           3.*apkk1*C5*C6*C8*C9*Fm2*t4*t6 - 
           12.*C5*C6*C8*C9*(Ekaon*Ekaon)*Fm2*t4*t6 - 
           6.*C5*C6*C8*C9*Ekaon*Elep*Fm2*t4*t6 - 
           6.*C5*C6*C8*C9*Ekaon*Enu*Fm2*t4*t6))) - 
           2.*(akk1*akk1)*((aml*aml)*(C3*C3)*
           (36.*Enu*(t2*t2) - 36.*akpk*am*C1*(C4*C4)*t2*t3 - 
           72.*am*(amk*amk)*C1*(C4*C4)*t2*t3 - 
           36.*akpk*amSig*C1*(C4*C4)*t2*t3 - 
           72.*(amk*amk)*amSig*C1*(C4*C4)*t2*t3 + 
           36.*am*apkk1*C1*(C4*C4)*t2*t3 + 
           36.*amSig*apkk1*C1*(C4*C4)*t2*t3 - 
           36.*(amk*amk)*C1*(C4*C4)*Enu*t2*t3 + 
           18.*akpk*am*(amk*amk)*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) + 
           18.*(am*am*am)*(amk*amk)*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) + 
           18.*am*(amk*amk*amk*amk)*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) + 
           18.*akpk*(amk*amk)*amSig*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) + 
           36.*(am*am)*(amk*amk)*amSig*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) + 
           18.*(amk*amk*amk*amk)*amSig*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) + 
           18.*am*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) - 
           18.*am*(amk*amk)*apkk1*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) - 
           18.*(amk*amk)*amSig*apkk1*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) - 
           9.*(am*am)*(amk*amk)*(C1*C1)*(C4*C4*C4*C4)*Enu*(t3*t3) + 
           9.*(amk*amk*amk*amk)*(C1*C1)*(C4*C4*C4*C4)*Enu*(t3*t3) - 
           18.*am*(amk*amk)*amSig*(C1*C1)*(C4*C4*C4*C4)*Enu*(t3*t3) - 
           9.*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4*C4*C4)*Enu*(t3*t3) - 
           12.*akpk*am*(C5*C5)*C6*t2*t4 - 
           24.*am*(amk*amk)*(C5*C5)*C6*t2*t4 - 
           12.*akpk*amLam*(C5*C5)*C6*t2*t4 - 
           24.*(amk*amk)*amLam*(C5*C5)*C6*t2*t4 + 
           12.*am*apkk1*(C5*C5)*C6*t2*t4 + 
           12.*amLam*apkk1*(C5*C5)*C6*t2*t4 - 
           12.*(amk*amk)*(C5*C5)*C6*Enu*t2*t4 + 
           12.*akpk*am*(amk*amk)*C1*(C4*C4)*(C5*C5)*C6*t3*t4 + 
           12.*(am*am*am)*(amk*amk)*C1*(C4*C4)*(C5*C5)*C6*t3*t4 + 
           12.*am*(amk*amk*amk*amk)*C1*(C4*C4)*(C5*C5)*C6*t3*t4 + 
           6.*akpk*(amk*amk)*amLam*C1*(C4*C4)*(C5*C5)*C6*t3*
           t4 + 12.*(am*am)*(amk*amk)*amLam*C1*(C4*C4)*(C5*C5)*
           C6*t3*t4 + 
           6.*(amk*amk*amk*amk)*amLam*C1*(C4*C4)*(C5*C5)*C6*t3*t4 + 
           6.*akpk*(amk*amk)*amSig*C1*(C4*C4)*(C5*C5)*C6*t3*
           t4 + 12.*(am*am)*(amk*amk)*amSig*C1*(C4*C4)*(C5*C5)*
           C6*t3*t4 + 
           6.*(amk*amk*amk*amk)*amSig*C1*(C4*C4)*(C5*C5)*C6*t3*t4 + 
           12.*am*(amk*amk)*amLam*amSig*C1*(C4*C4)*(C5*C5)*C6*
           t3*t4 - 
           12.*am*(amk*amk)*apkk1*C1*(C4*C4)*(C5*C5)*C6*t3*t4 - 
           6.*(amk*amk)*amLam*apkk1*C1*(C4*C4)*(C5*C5)*C6*t3*
           t4 - 6.*(amk*amk)*amSig*apkk1*C1*(C4*C4)*(C5*C5)*C6*
           t3*t4 - 
           6.*(am*am)*(amk*amk)*C1*(C4*C4)*(C5*C5)*C6*Enu*t3*t4 + 
           6.*(amk*amk*amk*amk)*C1*(C4*C4)*(C5*C5)*C6*Enu*t3*t4 - 
           6.*am*(amk*amk)*amLam*C1*(C4*C4)*(C5*C5)*C6*Enu*t3*
           t4 - 6.*am*(amk*amk)*amSig*C1*(C4*C4)*(C5*C5)*C6*
           Enu*t3*t4 - 
           6.*(amk*amk)*amLam*amSig*C1*(C4*C4)*(C5*C5)*C6*Enu*
           t3*t4 + 
           2.*akpk*am*(amk*amk)*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) + 
           2.*(am*am*am)*(amk*amk)*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) + 
           2.*am*(amk*amk*amk*amk)*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) + 
           2.*akpk*(amk*amk)*amLam*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) + 
           4.*(am*am)*(amk*amk)*amLam*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) + 
           2.*(amk*amk*amk*amk)*amLam*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) + 
           2.*am*(amk*amk)*(amLam*amLam)*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) - 
           2.*am*(amk*amk)*apkk1*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) - 
           2.*(amk*amk)*amLam*apkk1*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) - 
           1.*(am*am)*(amk*amk)*(C5*C5*C5*C5)*(C6*C6)*Enu*(t4*t4) + 
           (amk*amk*amk*amk)*(C5*C5*C5*C5)*(C6*C6)*Enu*(t4*t4) - 
           2.*am*(amk*amk)*amLam*(C5*C5*C5*C5)*(C6*C6)*Enu*(t4*t4) - 
           1.*(amk*amk)*(amLam*amLam)*(C5*C5*C5*C5)*(C6*C6)*Enu*(t4*t4) + 
           4.*(am*am)*(Ekaon*Ekaon*Ekaon)*
           ((3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4)*(3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4)) - 
           2.*am*(Ekaon*Ekaon)*
           (9.*(C1*C1)*(C4*C4*C4*C4)*
           ((am*am) + 2.*(amk*amk) + (amSig*amSig) + 
           2.*am*(amSig + Elep - 1.*Enu))*(t3*t3) + 
           6.*C1*(C4*C4)*t3*
           (-12.*t2 + 
           (C5*C5)*C6*
           ((am*am) + 2.*(amk*amk) + amLam*amSig + 
           am*(amLam + amSig + 2.*Elep - 2.*Enu)
           )*t4) + 
           (C5*C5)*C6*t4*
           (-24.*t2 + 
           (C5*C5)*C6*
           ((am*am) + 2.*(amk*amk) + (amLam*amLam) + 
           2.*am*(amLam + Elep - 1.*Enu))*t4))
           + Ekaon*
           (108.*(t2*t2) - 
           27.*(am*am)*(amk*amk)*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) + 
           9.*(amk*amk*amk*amk)*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) - 
           18.*am*(amk*amk)*amSig*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) + 
           9.*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) + 
           18.*(am*am)*apkk1*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) - 
           18.*(amSig*amSig)*apkk1*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) + 
           36.*am*(amk*amk)*(C1*C1)*(C4*C4*C4*C4)*Elep*(t3*t3) - 
           36.*am*(amk*amk)*(C1*C1)*(C4*C4*C4*C4)*Enu*(t3*t3) - 
           18.*(am*am)*(amk*amk)*C1*(C4*C4)*(C5*C5)*C6*t3*t4 + 
           6.*(amk*amk*amk*amk)*C1*(C4*C4)*(C5*C5)*C6*t3*t4 - 
           6.*am*(amk*amk)*amLam*C1*(C4*C4)*(C5*C5)*C6*t3*
           t4 - 
           6.*am*(amk*amk)*amSig*C1*(C4*C4)*(C5*C5)*C6*t3*
           t4 + 
           6.*(amk*amk)*amLam*amSig*C1*(C4*C4)*(C5*C5)*C6*t3*
           t4 + 
           12.*(am*am)*apkk1*C1*(C4*C4)*(C5*C5)*C6*t3*t4 - 
           12.*amLam*amSig*apkk1*C1*(C4*C4)*(C5*C5)*C6*t3*
           t4 + 
           24.*am*(amk*amk)*C1*(C4*C4)*(C5*C5)*C6*Elep*t3*
           t4 - 
           24.*am*(amk*amk)*C1*(C4*C4)*(C5*C5)*C6*Enu*t3*
           t4 - 3.*(am*am)*(amk*amk)*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) + 
           (amk*amk*amk*amk)*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) - 
           2.*am*(amk*amk)*amLam*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) + 
           (amk*amk)*(amLam*amLam)*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) + 
           2.*(am*am)*apkk1*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) - 
           2.*(amLam*amLam)*apkk1*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) + 
           4.*am*(amk*amk)*(C5*C5*C5*C5)*(C6*C6)*Elep*(t4*t4) - 
           4.*am*(amk*amk)*(C5*C5*C5*C5)*(C6*C6)*Enu*(t4*t4) - 
           24.*((amk*amk) + am*(Elep - 1.*Enu))*t2*
           (3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4) - 
           2.*akpk*
           ((am*am)*
           ((3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4)*(3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4)) - 
           1.*
           ((3.*amSig*C1*(C4*C4)*t3 + amLam*(C5*C5)*C6*t4)*(3.*amSig*C1*(C4*C4)*t3 + amLam*(C5*C5)*C6*t4)))) + 
           Elep*(-36.*(t2*t2) + 
           12.*(amk*amk)*t2*
           (3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4) + 
           (amk*amk)*
           ((am*am)*
           ((3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4)*(3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4)) - 
           1.*(amk*amk)*
           ((3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4)*(3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4)) + 
           2.*am*(3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4)*
           (3.*amSig*C1*(C4*C4)*t3 + 
           amLam*(C5*C5)*C6*t4) + 
           ((3.*amSig*C1*(C4*C4)*t3 + amLam*(C5*C5)*C6*t4)*(3.*amSig*C1*(C4*C4)*t3 + amLam*(C5*C5)*C6*t4))))) + 
           2.*((amk*amk*amk*amk)*(3.*C1*C4*Fm1*t3 - 1.*C5*C6*Fm2*t4)*
           (3.*C1*C4*
           (-2. + am*Fm1 - 1.*Elep*Fm1 + Enu*Fm1)*t3
           + C5*C6*
           (6. - 1.*am*Fm2 + Elep*Fm2 - 1.*Enu*Fm2)*
           t4) + 
           (amk*amk)*(9.*(C1*C1)*(C4*C4)*Fm1*
           ((am*am*am)*Fm1 + 
           (am*am)*
           (-2. + 2.*amSig*Fm1 - 2.*Ekaon*Fm1 + 
           Elep*Fm1 - 1.*Enu*Fm1) + 
           amSig*
           (-4.*Ekaon + 
           2.*(akpk - 1.*apkk1)*Fm1 + 
           amSig*(-2. + Elep*Fm1 - 1.*Enu*Fm1))
           + am*
           ((amSig*amSig)*Fm1 + 
           2.*(akpk - 1.*apkk1)*Fm1 + 
           4.*Ekaon*
           (1. + Elep*Fm1 - 1.*Enu*Fm1) + 
           2.*amSig*
           (-2. + Ekaon*Fm1 + Elep*Fm1 - 
           1.*Enu*Fm1)))*(t3*t3) - 
           6.*C1*C4*t3*
           (12.*f*Fm1*t1 - 
           3.*amLam*amSig*C5*C6*Fm1*t4 - 
           3.*amLam*C5*C6*Ekaon*Fm1*t4 - 
           3.*amSig*C5*C6*Ekaon*Fm1*t4 - 
           1.*amLam*amSig*C5*C6*Fm2*t4 - 
           1.*amLam*C5*C6*Ekaon*Fm2*t4 - 
           1.*amSig*C5*C6*Ekaon*Fm2*t4 + 
           (am*am*am)*C5*C6*Fm1*Fm2*t4 + 
           akpk*amLam*C5*C6*Fm1*Fm2*t4 + 
           akpk*amSig*C5*C6*Fm1*Fm2*t4 - 
           1.*amLam*apkk1*C5*C6*Fm1*Fm2*t4 - 
           1.*amSig*apkk1*C5*C6*Fm1*Fm2*t4 + 
           amLam*amSig*C5*C6*Elep*Fm1*Fm2*t4 - 
           1.*amLam*amSig*C5*C6*Enu*Fm1*Fm2*t4 + 
           (am*am)*C5*C6*
           (-1.*Fm2 + 
           Fm1*
           (-3. + amLam*Fm2 + amSig*Fm2 - 
           2.*Ekaon*Fm2 + Elep*Fm2 - 
           1.*Enu*Fm2))*t4 - 
           3.*amSig*C2*C7*Fm1*t5 - 
           3.*C2*C7*Ekaon*Fm1*t5 + 
           3.*amSig*C8*C9*Fm1*t6 + 
           3.*C8*C9*Ekaon*Fm1*t6 + 
           am*
           (6.*C5*C6*Ekaon*Fm1*t4 + 
           2.*C5*C6*Ekaon*Fm2*t4 + 
           2.*akpk*C5*C6*Fm1*Fm2*t4 - 
           2.*apkk1*C5*C6*Fm1*Fm2*t4 + 
           4.*C5*C6*Ekaon*Elep*Fm1*Fm2*t4 - 
           4.*C5*C6*Ekaon*Enu*Fm1*Fm2*t4 + 
           amSig*C5*C6*
           (-1.*Fm2 + 
           Fm1*
           (-3. + Ekaon*Fm2 + Elep*Fm2 - 
           1.*Enu*Fm2))*t4 + 
           amLam*C5*C6*
           (-1.*Fm2 + 
           Fm1*
           (-3. + amSig*Fm2 + Ekaon*Fm2 + 
           Elep*Fm2 - 1.*Enu*Fm2))*t4 - 
           3.*C2*C7*Fm1*t5 + 3.*C8*C9*Fm1*t6))
           + C5*C6*Fm2*t4*
           (24.*f*t1 - 6.*(amLam*amLam)*C5*C6*t4 - 
           12.*amLam*C5*C6*Ekaon*t4 + 
           (am*am*am)*C5*C6*Fm2*t4 + 
           2.*akpk*amLam*C5*C6*Fm2*t4 - 
           2.*amLam*apkk1*C5*C6*Fm2*t4 + 
           (amLam*amLam)*C5*C6*Elep*Fm2*t4 - 
           1.*(amLam*amLam)*C5*C6*Enu*Fm2*t4 + 
           (am*am)*C5*C6*
           (-6. + 2.*amLam*Fm2 - 2.*Ekaon*Fm2 + 
           Elep*Fm2 - 1.*Enu*Fm2)*t4 - 
           6.*amLam*C2*C7*t5 - 6.*C2*C7*Ekaon*t5 + 
           6.*amLam*C8*C9*t6 + 6.*C8*C9*Ekaon*t6 + 
           am*
           ((amLam*amLam)*C5*C6*Fm2*t4 + 
           2.*amLam*C5*C6*
           (-6. + Ekaon*Fm2 + Elep*Fm2 - 
           1.*Enu*Fm2)*t4 + 
           2.*
           (C5*C6*
           ((akpk - 1.*apkk1)*Fm2 + 
           2.*Ekaon*
           (3. + Elep*Fm2 - 1.*Enu*Fm2))*t4 - 
           3.*C2*C7*t5 + 3.*C8*C9*t6)))) - 
           2.*Ekaon*
           (-1.*(3.*amSig*C1*C4*Fm1*t3 - 
           1.*amLam*C5*C6*Fm2*t4)*
           (-12.*f*t1 + 
           (akpk - 1.*apkk1)*
           (3.*amSig*C1*C4*Fm1*t3 - 
           1.*amLam*C5*C6*Fm2*t4)) + 
           (am*am)*(3.*C1*C4*Fm1*t3 - 1.*C5*C6*Fm2*t4)*
           (akpk*
           (3.*C1*C4*Fm1*t3 - 1.*C5*C6*Fm2*t4) + 
           apkk1*
           (-3.*C1*C4*Fm1*t3 + C5*C6*Fm2*t4) + 
           2.*Ekaon*
           (3.*amSig*C1*C4*Fm1*t3 + 
           3.*C1*C4*(Elep - 1.*Enu)*Fm1*t3 - 
           1.*C5*C6*(amLam + Elep - 1.*Enu)*Fm2*
           t4)) - 
           6.*am*
           (6.*amSig*(C1*C1)*(C4*C4)*Ekaon*Fm1*(t3*t3) + 
           C5*C6*Fm2*t4*
           (-2.*f*t1 + 
           Ekaon*
           (2.*amLam*C5*C6*t4 + C2*C7*t5 - 
           1.*C8*C9*t6)) + 
           C1*C4*t3*
           (6.*f*Fm1*t1 - 
           1.*Ekaon*
           (amLam*C5*C6*(3.*Fm1 + Fm2)*t4 + 
           amSig*C5*C6*(3.*Fm1 + Fm2)*t4 + 
           3.*Fm1*(C2*C7*t5 - 1.*C8*C9*t6)))))
           )) + 2.*(Enu*
           (-144.*apkk1*(t1*t1) + 288.*am*Elep*(t1*t1) - 
           288.*apkk1*f*(t1*t1) - 144.*apkk1*(f*f)*(t1*t1) + 
           288.*am*Elep*(f*f)*(t1*t1) - 
           72.*(amk*amk)*apkk1*C1*C4*t1*t3 - 
           72.*(am*am)*apkk1*C1*(C4*C4)*t1*t3 + 
           72.*(amk*amk)*apkk1*C1*(C4*C4)*t1*t3 - 
           72.*am*amSig*apkk1*C1*(C4*C4)*t1*t3 + 
           144.*am*apkk1*C1*C4*Ekaon*t1*t3 - 
           144.*am*apkk1*C1*(C4*C4)*Ekaon*t1*t3 - 
           144.*am*(amk*amk)*C1*(C4*C4)*Elep*t1*t3 + 
           288.*(am*am)*C1*(C4*C4)*Ekaon*Elep*t1*t3 + 
           72.*(am*am)*apkk1*C1*C4*f*t1*t3 - 
           72.*(amk*amk)*apkk1*C1*C4*f*t1*t3 + 
           72.*am*amSig*apkk1*C1*C4*f*t1*t3 + 
           144.*(am*am)*apkk1*C1*(C4*C4)*f*t1*t3 + 
           72.*(amk*amk)*apkk1*C1*(C4*C4)*f*t1*t3 + 
           144.*am*amSig*apkk1*C1*(C4*C4)*f*t1*t3 + 
           144.*am*apkk1*C1*C4*Ekaon*f*t1*t3 - 
           144.*am*apkk1*C1*(C4*C4)*Ekaon*f*t1*t3 + 
           144.*am*(amk*amk)*C1*C4*Elep*f*t1*t3 - 
           288.*(am*am)*C1*C4*Ekaon*Elep*f*t1*t3 - 
           72.*am*(apkk1*apkk1)*C1*C4*Fm1*t1*t3 - 
           72.*amSig*(apkk1*apkk1)*C1*C4*Fm1*t1*t3 + 
           72.*(am*am)*apkk1*C1*C4*Elep*Fm1*t1*t3 + 
           72.*am*amSig*apkk1*C1*C4*Elep*Fm1*t1*t3 - 
           72.*(am*am)*apkk1*C1*C4*Enu*Fm1*t1*t3 - 
           72.*am*amSig*apkk1*C1*C4*Enu*Fm1*t1*t3 - 
           72.*am*(apkk1*apkk1)*C1*C4*f*Fm1*t1*t3 - 
           72.*amSig*(apkk1*apkk1)*C1*C4*f*Fm1*t1*t3 + 
           72.*(am*am)*apkk1*C1*C4*Elep*f*Fm1*t1*t3 + 
           72.*am*amSig*apkk1*C1*C4*Elep*f*Fm1*t1*t3 + 
           72.*(am*am)*apkk1*C1*C4*Enu*f*Fm1*t1*t3 + 
           72.*am*amSig*apkk1*C1*C4*Enu*f*Fm1*t1*t3 + 
           27.*(am*am)*(amk*amk)*apkk1*(C1*C1)*(C4*C4)*(t3*t3) - 
           9.*(amk*amk*amk*amk)*apkk1*(C1*C1)*(C4*C4)*(t3*t3) + 
           36.*am*(amk*amk)*amSig*apkk1*(C1*C1)*(C4*C4)*(t3*t3) + 
           9.*(amk*amk)*(amSig*amSig)*apkk1*(C1*C1)*(C4*C4)*(t3*t3) + 
           18.*(am*am)*(amk*amk)*apkk1*(C1*C1)*(C4*C4*C4)*(t3*t3) + 
           18.*(amk*amk*amk*amk)*apkk1*(C1*C1)*(C4*C4*C4)*(t3*t3) - 
           18.*(amk*amk)*(amSig*amSig)*apkk1*(C1*C1)*(C4*C4*C4)*(t3*t3) + 
           27.*(am*am)*(amk*amk)*apkk1*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) - 
           9.*(amk*amk*amk*amk)*apkk1*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) + 
           36.*am*(amk*amk)*amSig*apkk1*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) + 
           9.*(amk*amk)*(amSig*amSig)*apkk1*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) - 
           18.*(am*am*am)*apkk1*(C1*C1)*(C4*C4)*Ekaon*(t3*t3) + 
           36.*am*(amk*amk)*apkk1*(C1*C1)*(C4*C4)*Ekaon*(t3*t3) + 
           18.*am*(amSig*amSig)*apkk1*(C1*C1)*(C4*C4)*Ekaon*
           (t3*t3) - 
           36.*(am*am*am)*apkk1*(C1*C1)*(C4*C4*C4)*Ekaon*(t3*t3) - 
           72.*am*(amk*amk)*apkk1*(C1*C1)*(C4*C4*C4)*Ekaon*(t3*t3) + 
           36.*am*(amSig*amSig)*apkk1*(C1*C1)*(C4*C4*C4)*Ekaon*
           (t3*t3) - 
           18.*(am*am*am)*apkk1*(C1*C1)*(C4*C4*C4*C4)*Ekaon*(t3*t3) + 
           36.*am*(amk*amk)*apkk1*(C1*C1)*(C4*C4*C4*C4)*Ekaon*(t3*t3) + 
           18.*am*(amSig*amSig)*apkk1*(C1*C1)*(C4*C4*C4*C4)*Ekaon*
           (t3*t3) - 
           36.*(am*am)*apkk1*(C1*C1)*(C4*C4)*(Ekaon*Ekaon)*(t3*t3) + 
           72.*(am*am)*apkk1*(C1*C1)*(C4*C4*C4)*(Ekaon*Ekaon)*(t3*t3) - 
           36.*(am*am)*apkk1*(C1*C1)*(C4*C4*C4*C4)*(Ekaon*Ekaon)*(t3*t3) - 
           18.*(am*am*am)*(amk*amk)*(C1*C1)*(C4*C4)*Elep*(t3*t3) + 
           18.*am*(amk*amk*amk*amk)*(C1*C1)*(C4*C4)*Elep*(t3*t3) - 
           36.*(am*am)*(amk*amk)*amSig*(C1*C1)*(C4*C4)*Elep*
           (t3*t3) - 
           18.*am*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4)*Elep*
           (t3*t3) - 
           18.*(am*am*am)*(amk*amk)*(C1*C1)*(C4*C4*C4*C4)*Elep*(t3*t3) + 
           18.*am*(amk*amk*amk*amk)*(C1*C1)*(C4*C4*C4*C4)*Elep*(t3*t3) - 
           36.*(am*am)*(amk*amk)*amSig*(C1*C1)*(C4*C4*C4*C4)*Elep*
           (t3*t3) - 
           18.*am*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4*C4*C4)*Elep*
           (t3*t3) - 
           72.*(am*am)*(amk*amk)*(C1*C1)*(C4*C4)*Ekaon*Elep*
           (t3*t3) - 
           72.*(am*am)*(amk*amk)*(C1*C1)*(C4*C4*C4*C4)*Ekaon*Elep*
           (t3*t3) + 
           72.*(am*am*am)*(C1*C1)*(C4*C4)*(Ekaon*Ekaon)*Elep*(t3*t3) + 
           72.*(am*am*am)*(C1*C1)*(C4*C4*C4*C4)*(Ekaon*Ekaon)*Elep*(t3*t3) + 
           72.*(amk*amk)*apkk1*C5*C6*t1*t4 - 
           24.*(am*am)*apkk1*(C5*C5)*C6*t1*t4 + 
           24.*(amk*amk)*apkk1*(C5*C5)*C6*t1*t4 - 
           24.*am*amLam*apkk1*(C5*C5)*C6*t1*t4 - 
           144.*am*apkk1*C5*C6*Ekaon*t1*t4 - 
           48.*am*apkk1*(C5*C5)*C6*Ekaon*t1*t4 - 
           48.*am*(amk*amk)*(C5*C5)*C6*Elep*t1*t4 + 
           96.*(am*am)*(C5*C5)*C6*Ekaon*Elep*t1*t4 - 
           72.*(am*am)*apkk1*C5*C6*f*t1*t4 + 
           72.*(amk*amk)*apkk1*C5*C6*f*t1*t4 - 
           72.*am*amLam*apkk1*C5*C6*f*t1*t4 + 
           48.*(am*am)*apkk1*(C5*C5)*C6*f*t1*t4 + 
           24.*(amk*amk)*apkk1*(C5*C5)*C6*f*t1*t4 + 
           48.*am*amLam*apkk1*(C5*C5)*C6*f*t1*t4 - 
           144.*am*apkk1*C5*C6*Ekaon*f*t1*t4 - 
           48.*am*apkk1*(C5*C5)*C6*Ekaon*f*t1*t4 - 
           144.*am*(amk*amk)*C5*C6*Elep*f*t1*t4 + 
           288.*(am*am)*C5*C6*Ekaon*Elep*f*t1*t4 + 
           24.*am*(apkk1*apkk1)*C5*C6*Fm2*t1*t4 + 
           24.*amLam*(apkk1*apkk1)*C5*C6*Fm2*t1*t4 - 
           24.*(am*am)*apkk1*C5*C6*Elep*Fm2*t1*t4 - 
           24.*am*amLam*apkk1*C5*C6*Elep*Fm2*t1*t4 + 
           24.*(am*am)*apkk1*C5*C6*Enu*Fm2*t1*t4 + 
           24.*am*amLam*apkk1*C5*C6*Enu*Fm2*t1*t4 + 
           24.*am*(apkk1*apkk1)*C5*C6*f*Fm2*t1*t4 + 
           24.*amLam*(apkk1*apkk1)*C5*C6*f*Fm2*t1*t4 - 
           24.*(am*am)*apkk1*C5*C6*Elep*f*Fm2*t1*t4 - 
           24.*am*amLam*apkk1*C5*C6*Elep*f*Fm2*t1*t4 - 
           24.*(am*am)*apkk1*C5*C6*Enu*f*Fm2*t1*t4 - 
           24.*am*amLam*apkk1*C5*C6*Enu*f*Fm2*t1*t4 - 
           54.*(am*am)*(amk*amk)*apkk1*C1*C4*C5*C6*t3*t4 + 
           18.*(amk*amk*amk*amk)*apkk1*C1*C4*C5*C6*t3*t4 - 
           36.*am*(amk*amk)*amLam*apkk1*C1*C4*C5*C6*t3*t4 - 
           36.*am*(amk*amk)*amSig*apkk1*C1*C4*C5*C6*t3*t4 - 
           18.*(amk*amk)*amLam*amSig*apkk1*C1*C4*C5*C6*t3*
           t4 - 18.*(am*am)*(amk*amk)*apkk1*C1*(C4*C4)*C5*C6*
           t3*t4 - 
           18.*(amk*amk*amk*amk)*apkk1*C1*(C4*C4)*C5*C6*t3*t4 + 
           18.*am*(amk*amk)*amLam*apkk1*C1*(C4*C4)*C5*C6*t3*
           t4 - 18.*am*(amk*amk)*amSig*apkk1*C1*(C4*C4)*C5*
           C6*t3*t4 + 
           18.*(amk*amk)*amLam*amSig*apkk1*C1*(C4*C4)*C5*C6*
           t3*t4 + 
           6.*(am*am)*(amk*amk)*apkk1*C1*C4*(C5*C5)*C6*t3*t4 + 
           6.*(amk*amk*amk*amk)*apkk1*C1*C4*(C5*C5)*C6*t3*t4 + 
           6.*am*(amk*amk)*amLam*apkk1*C1*C4*(C5*C5)*C6*t3*
           t4 - 6.*am*(amk*amk)*amSig*apkk1*C1*C4*(C5*C5)*C6*
           t3*t4 - 
           6.*(amk*amk)*amLam*amSig*apkk1*C1*C4*(C5*C5)*C6*t3*
           t4 + 18.*(am*am)*(amk*amk)*apkk1*C1*(C4*C4)*(C5*C5)*
           C6*t3*t4 - 
           6.*(amk*amk*amk*amk)*apkk1*C1*(C4*C4)*(C5*C5)*C6*t3*t4 + 
           12.*am*(amk*amk)*amLam*apkk1*C1*(C4*C4)*(C5*C5)*C6*
           t3*t4 + 
           12.*am*(amk*amk)*amSig*apkk1*C1*(C4*C4)*(C5*C5)*C6*
           t3*t4 + 
           6.*(amk*amk)*amLam*amSig*apkk1*C1*(C4*C4)*(C5*C5)*C6*
           t3*t4 + 
           36.*(am*am*am)*apkk1*C1*C4*C5*C6*Ekaon*t3*t4 - 
           72.*am*(amk*amk)*apkk1*C1*C4*C5*C6*Ekaon*t3*t4 - 
           36.*am*amLam*amSig*apkk1*C1*C4*C5*C6*Ekaon*t3*
           t4 + 36.*(am*am*am)*apkk1*C1*(C4*C4)*C5*C6*Ekaon*t3*
           t4 + 72.*am*(amk*amk)*apkk1*C1*(C4*C4)*C5*C6*
           Ekaon*t3*t4 - 
           36.*(am*am)*amLam*apkk1*C1*(C4*C4)*C5*C6*Ekaon*t3*
           t4 + 36.*(am*am)*amSig*apkk1*C1*(C4*C4)*C5*C6*
           Ekaon*t3*t4 - 
           36.*am*amLam*amSig*apkk1*C1*(C4*C4)*C5*C6*Ekaon*
           t3*t4 - 
           12.*(am*am*am)*apkk1*C1*C4*(C5*C5)*C6*Ekaon*t3*t4 - 
           24.*am*(amk*amk)*apkk1*C1*C4*(C5*C5)*C6*Ekaon*t3*
           t4 - 12.*(am*am)*amLam*apkk1*C1*C4*(C5*C5)*C6*
           Ekaon*t3*t4 + 
           12.*(am*am)*amSig*apkk1*C1*C4*(C5*C5)*C6*Ekaon*t3*
           t4 + 12.*am*amLam*amSig*apkk1*C1*C4*(C5*C5)*C6*
           Ekaon*t3*t4 - 
           12.*(am*am*am)*apkk1*C1*(C4*C4)*(C5*C5)*C6*Ekaon*t3*
           t4 + 24.*am*(amk*amk)*apkk1*C1*(C4*C4)*(C5*C5)*C6*
           Ekaon*t3*t4 + 
           12.*am*amLam*amSig*apkk1*C1*(C4*C4)*(C5*C5)*C6*
           Ekaon*t3*t4 + 
           72.*(am*am)*apkk1*C1*C4*C5*C6*(Ekaon*Ekaon)*t3*t4 - 
           72.*(am*am)*apkk1*C1*(C4*C4)*C5*C6*(Ekaon*Ekaon)*t3*
           t4 + 24.*(am*am)*apkk1*C1*C4*(C5*C5)*C6*(Ekaon*Ekaon)*
           t3*t4 - 
           24.*(am*am)*apkk1*C1*(C4*C4)*(C5*C5)*C6*(Ekaon*Ekaon)*t3*
           t4 + 36.*(am*am*am)*(amk*amk)*C1*C4*C5*C6*Elep*t3*
           t4 - 36.*am*(amk*amk*amk*amk)*C1*C4*C5*C6*Elep*t3*t4 + 
           36.*(am*am)*(amk*amk)*amLam*C1*C4*C5*C6*Elep*t3*
           t4 + 36.*(am*am)*(amk*amk)*amSig*C1*C4*C5*C6*Elep*
           t3*t4 + 
           36.*am*(amk*amk)*amLam*amSig*C1*C4*C5*C6*Elep*t3*
           t4 - 12.*(am*am*am)*(amk*amk)*C1*(C4*C4)*(C5*C5)*C6*Elep*
           t3*t4 + 
           12.*am*(amk*amk*amk*amk)*C1*(C4*C4)*(C5*C5)*C6*Elep*t3*t4 - 
           12.*(am*am)*(amk*amk)*amLam*C1*(C4*C4)*(C5*C5)*C6*Elep*
           t3*t4 - 
           12.*(am*am)*(amk*amk)*amSig*C1*(C4*C4)*(C5*C5)*C6*Elep*
           t3*t4 - 
           12.*am*(amk*amk)*amLam*amSig*C1*(C4*C4)*(C5*C5)*C6*
           Elep*t3*t4 + 
           144.*(am*am)*(amk*amk)*C1*C4*C5*C6*Ekaon*Elep*t3*
           t4 - 48.*(am*am)*(amk*amk)*C1*(C4*C4)*(C5*C5)*C6*
           Ekaon*Elep*t3*t4 - 
           144.*(am*am*am)*C1*C4*C5*C6*(Ekaon*Ekaon)*Elep*t3*t4 + 
           48.*(am*am*am)*C1*(C4*C4)*(C5*C5)*C6*(Ekaon*Ekaon)*Elep*t3*
           t4 - 18.*(amk*amk)*amLam*(apkk1*apkk1)*C1*C4*C5*C6*
           Fm1*t3*t4 + 
           18.*(amk*amk)*amSig*(apkk1*apkk1)*C1*C4*C5*C6*Fm1*t3*
           t4 - 6.*(amk*amk)*amLam*(apkk1*apkk1)*C1*C4*(C5*C5)*C6*
           Fm1*t3*t4 + 
           6.*(amk*amk)*amSig*(apkk1*apkk1)*C1*C4*(C5*C5)*C6*Fm1*
           t3*t4 + 
           36.*am*amLam*(apkk1*apkk1)*C1*C4*C5*C6*Ekaon*Fm1*
           t3*t4 - 
           36.*am*amSig*(apkk1*apkk1)*C1*C4*C5*C6*Ekaon*Fm1*
           t3*t4 + 
           12.*am*amLam*(apkk1*apkk1)*C1*C4*(C5*C5)*C6*Ekaon*
           Fm1*t3*t4 - 
           12.*am*amSig*(apkk1*apkk1)*C1*C4*(C5*C5)*C6*Ekaon*
           Fm1*t3*t4 + 
           18.*am*(amk*amk)*amLam*apkk1*C1*C4*C5*C6*Elep*
           Fm1*t3*t4 - 
           18.*am*(amk*amk)*amSig*apkk1*C1*C4*C5*C6*Elep*
           Fm1*t3*t4 + 
           6.*am*(amk*amk)*amLam*apkk1*C1*C4*(C5*C5)*C6*Elep*
           Fm1*t3*t4 - 
           6.*am*(amk*amk)*amSig*apkk1*C1*C4*(C5*C5)*C6*Elep*
           Fm1*t3*t4 - 
           36.*(am*am)*amLam*apkk1*C1*C4*C5*C6*Ekaon*Elep*
           Fm1*t3*t4 + 
           36.*(am*am)*amSig*apkk1*C1*C4*C5*C6*Ekaon*Elep*
           Fm1*t3*t4 - 
           12.*(am*am)*amLam*apkk1*C1*C4*(C5*C5)*C6*Ekaon*
           Elep*Fm1*t3*t4 + 
           12.*(am*am)*amSig*apkk1*C1*C4*(C5*C5)*C6*Ekaon*
           Elep*Fm1*t3*t4 + 
           18.*am*(amk*amk)*amLam*apkk1*C1*C4*C5*C6*Enu*Fm1*
           t3*t4 - 
           18.*am*(amk*amk)*amSig*apkk1*C1*C4*C5*C6*Enu*Fm1*
           t3*t4 - 
           6.*am*(amk*amk)*amLam*apkk1*C1*C4*(C5*C5)*C6*Enu*
           Fm1*t3*t4 + 
           6.*am*(amk*amk)*amSig*apkk1*C1*C4*(C5*C5)*C6*Enu*
           Fm1*t3*t4 - 
           36.*(am*am)*amLam*apkk1*C1*C4*C5*C6*Ekaon*Enu*
           Fm1*t3*t4 + 
           36.*(am*am)*amSig*apkk1*C1*C4*C5*C6*Ekaon*Enu*
           Fm1*t3*t4 + 
           12.*(am*am)*amLam*apkk1*C1*C4*(C5*C5)*C6*Ekaon*
           Enu*Fm1*t3*t4 - 
           12.*(am*am)*amSig*apkk1*C1*C4*(C5*C5)*C6*Ekaon*
           Enu*Fm1*t3*t4 + 
           6.*(amk*amk)*amLam*(apkk1*apkk1)*C1*C4*C5*C6*Fm2*t3*
           t4 - 6.*(amk*amk)*amSig*(apkk1*apkk1)*C1*C4*C5*C6*
           Fm2*t3*t4 - 
           6.*(amk*amk)*amLam*(apkk1*apkk1)*C1*(C4*C4)*C5*C6*Fm2*
           t3*t4 + 
           6.*(amk*amk)*amSig*(apkk1*apkk1)*C1*(C4*C4)*C5*C6*Fm2*
           t3*t4 - 
           12.*am*amLam*(apkk1*apkk1)*C1*C4*C5*C6*Ekaon*Fm2*
           t3*t4 + 
           12.*am*amSig*(apkk1*apkk1)*C1*C4*C5*C6*Ekaon*Fm2*
           t3*t4 + 
           12.*am*amLam*(apkk1*apkk1)*C1*(C4*C4)*C5*C6*Ekaon*
           Fm2*t3*t4 - 
           12.*am*amSig*(apkk1*apkk1)*C1*(C4*C4)*C5*C6*Ekaon*
           Fm2*t3*t4 - 
           6.*am*(amk*amk)*amLam*apkk1*C1*C4*C5*C6*Elep*Fm2*
           t3*t4 + 
           6.*am*(amk*amk)*amSig*apkk1*C1*C4*C5*C6*Elep*Fm2*
           t3*t4 + 
           6.*am*(amk*amk)*amLam*apkk1*C1*(C4*C4)*C5*C6*Elep*
           Fm2*t3*t4 - 
           6.*am*(amk*amk)*amSig*apkk1*C1*(C4*C4)*C5*C6*Elep*
           Fm2*t3*t4 + 
           12.*(am*am)*amLam*apkk1*C1*C4*C5*C6*Ekaon*Elep*
           Fm2*t3*t4 - 
           12.*(am*am)*amSig*apkk1*C1*C4*C5*C6*Ekaon*Elep*
           Fm2*t3*t4 - 
           12.*(am*am)*amLam*apkk1*C1*(C4*C4)*C5*C6*Ekaon*
           Elep*Fm2*t3*t4 + 
           12.*(am*am)*amSig*apkk1*C1*(C4*C4)*C5*C6*Ekaon*
           Elep*Fm2*t3*t4 - 
           6.*am*(amk*amk)*amLam*apkk1*C1*C4*C5*C6*Enu*Fm2*
           t3*t4 + 
           6.*am*(amk*amk)*amSig*apkk1*C1*C4*C5*C6*Enu*Fm2*
           t3*t4 - 
           6.*am*(amk*amk)*amLam*apkk1*C1*(C4*C4)*C5*C6*Enu*
           Fm2*t3*t4 + 
           6.*am*(amk*amk)*amSig*apkk1*C1*(C4*C4)*C5*C6*Enu*
           Fm2*t3*t4 + 
           12.*(am*am)*amLam*apkk1*C1*C4*C5*C6*Ekaon*Enu*
           Fm2*t3*t4 - 
           12.*(am*am)*amSig*apkk1*C1*C4*C5*C6*Ekaon*Enu*
           Fm2*t3*t4 + 
           12.*(am*am)*amLam*apkk1*C1*(C4*C4)*C5*C6*Ekaon*
           Enu*Fm2*t3*t4 - 
           12.*(am*am)*amSig*apkk1*C1*(C4*C4)*C5*C6*Ekaon*
           Enu*Fm2*t3*t4 + 
           27.*(am*am)*(amk*amk)*apkk1*(C5*C5)*(C6*C6)*(t4*t4) - 
           9.*(amk*amk*amk*amk)*apkk1*(C5*C5)*(C6*C6)*(t4*t4) + 
           36.*am*(amk*amk)*amLam*apkk1*(C5*C5)*(C6*C6)*(t4*t4) + 
           9.*(amk*amk)*(amLam*amLam)*apkk1*(C5*C5)*(C6*C6)*(t4*t4) - 
           6.*(am*am)*(amk*amk)*apkk1*(C5*C5*C5)*(C6*C6)*(t4*t4) - 
           6.*(amk*amk*amk*amk)*apkk1*(C5*C5*C5)*(C6*C6)*(t4*t4) + 
           6.*(amk*amk)*(amLam*amLam)*apkk1*(C5*C5*C5)*(C6*C6)*(t4*t4) + 
           3.*(am*am)*(amk*amk)*apkk1*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) - 
           1.*(amk*amk*amk*amk)*apkk1*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) + 
           4.*am*(amk*amk)*amLam*apkk1*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) + 
           (amk*amk)*(amLam*amLam)*apkk1*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) - 
           18.*(am*am*am)*apkk1*(C5*C5)*(C6*C6)*Ekaon*(t4*t4) + 
           36.*am*(amk*amk)*apkk1*(C5*C5)*(C6*C6)*Ekaon*(t4*t4) + 
           18.*am*(amLam*amLam)*apkk1*(C5*C5)*(C6*C6)*Ekaon*
           (t4*t4) + 
           12.*(am*am*am)*apkk1*(C5*C5*C5)*(C6*C6)*Ekaon*(t4*t4) + 
           24.*am*(amk*amk)*apkk1*(C5*C5*C5)*(C6*C6)*Ekaon*(t4*t4) - 
           12.*am*(amLam*amLam)*apkk1*(C5*C5*C5)*(C6*C6)*Ekaon*
           (t4*t4) - 
           2.*(am*am*am)*apkk1*(C5*C5*C5*C5)*(C6*C6)*Ekaon*(t4*t4) + 
           4.*am*(amk*amk)*apkk1*(C5*C5*C5*C5)*(C6*C6)*Ekaon*(t4*t4) + 
           2.*am*(amLam*amLam)*apkk1*(C5*C5*C5*C5)*(C6*C6)*Ekaon*
           (t4*t4) - 
           36.*(am*am)*apkk1*(C5*C5)*(C6*C6)*(Ekaon*Ekaon)*(t4*t4) - 
           24.*(am*am)*apkk1*(C5*C5*C5)*(C6*C6)*(Ekaon*Ekaon)*(t4*t4) - 
           4.*(am*am)*apkk1*(C5*C5*C5*C5)*(C6*C6)*(Ekaon*Ekaon)*(t4*t4) - 
           18.*(am*am*am)*(amk*amk)*(C5*C5)*(C6*C6)*Elep*(t4*t4) + 
           18.*am*(amk*amk*amk*amk)*(C5*C5)*(C6*C6)*Elep*(t4*t4) - 
           36.*(am*am)*(amk*amk)*amLam*(C5*C5)*(C6*C6)*Elep*
           (t4*t4) - 
           18.*am*(amk*amk)*(amLam*amLam)*(C5*C5)*(C6*C6)*Elep*
           (t4*t4) - 
           2.*(am*am*am)*(amk*amk)*(C5*C5*C5*C5)*(C6*C6)*Elep*(t4*t4) + 
           2.*am*(amk*amk*amk*amk)*(C5*C5*C5*C5)*(C6*C6)*Elep*(t4*t4) - 
           4.*(am*am)*(amk*amk)*amLam*(C5*C5*C5*C5)*(C6*C6)*Elep*
           (t4*t4) - 
           2.*am*(amk*amk)*(amLam*amLam)*(C5*C5*C5*C5)*(C6*C6)*Elep*
           (t4*t4) - 
           72.*(am*am)*(amk*amk)*(C5*C5)*(C6*C6)*Ekaon*Elep*
           (t4*t4) - 
           8.*(am*am)*(amk*amk)*(C5*C5*C5*C5)*(C6*C6)*Ekaon*Elep*
           (t4*t4) + 
           72.*(am*am*am)*(C5*C5)*(C6*C6)*(Ekaon*Ekaon)*Elep*(t4*t4) + 
           8.*(am*am*am)*(C5*C5*C5*C5)*(C6*C6)*(Ekaon*Ekaon)*Elep*(t4*t4) - 
           18.*am*(amk*amk)*apkk1*C1*C2*C4*C7*t3*t5 - 
           18.*(amk*amk)*amSig*apkk1*C1*C2*C4*C7*t3*t5 - 
           18.*am*(apkk1*apkk1)*C1*C2*C4*C7*t3*t5 - 
           18.*amSig*(apkk1*apkk1)*C1*C2*C4*C7*t3*t5 - 
           18.*am*(apkk1*apkk1)*C1*C2*(C4*C4)*C7*t3*t5 - 
           18.*amSig*(apkk1*apkk1)*C1*C2*(C4*C4)*C7*t3*t5 - 
           18.*(amk*amk)*(apkk1*apkk1)*C1*C2*C4*C7*Fm1*t3*t5 + 
           36.*am*(apkk1*apkk1)*C1*C2*C4*C7*Ekaon*Fm1*t3*t5 + 
           18.*am*(amk*amk)*apkk1*C2*C5*C6*C7*t4*t5 + 
           18.*(amk*amk)*amLam*apkk1*C2*C5*C6*C7*t4*t5 + 
           18.*am*(apkk1*apkk1)*C2*C5*C6*C7*t4*t5 + 
           18.*amLam*(apkk1*apkk1)*C2*C5*C6*C7*t4*t5 - 
           6.*am*(apkk1*apkk1)*C2*(C5*C5)*C6*C7*t4*t5 - 
           6.*amLam*(apkk1*apkk1)*C2*(C5*C5)*C6*C7*t4*t5 + 
           6.*(amk*amk)*(apkk1*apkk1)*C2*C5*C6*C7*Fm2*t4*t5 - 
           12.*am*(apkk1*apkk1)*C2*C5*C6*C7*Ekaon*Fm2*t4*t5 + 
           18.*am*(amk*amk)*apkk1*C1*C4*C8*C9*t3*t6 + 
           18.*(amk*amk)*amSig*apkk1*C1*C4*C8*C9*t3*t6 + 
           18.*am*(apkk1*apkk1)*C1*C4*C8*C9*t3*t6 + 
           18.*amSig*(apkk1*apkk1)*C1*C4*C8*C9*t3*t6 + 
           18.*am*(apkk1*apkk1)*C1*(C4*C4)*C8*C9*t3*t6 + 
           18.*amSig*(apkk1*apkk1)*C1*(C4*C4)*C8*C9*t3*t6 + 
           18.*(amk*amk)*(apkk1*apkk1)*C1*C4*C8*C9*Fm1*t3*t6 - 
           36.*am*(apkk1*apkk1)*C1*C4*C8*C9*Ekaon*Fm1*t3*t6 - 
           18.*am*(amk*amk)*apkk1*C5*C6*C8*C9*t4*t6 - 
           18.*(amk*amk)*amLam*apkk1*C5*C6*C8*C9*t4*t6 - 
           18.*am*(apkk1*apkk1)*C5*C6*C8*C9*t4*t6 - 
           18.*amLam*(apkk1*apkk1)*C5*C6*C8*C9*t4*t6 + 
           6.*am*(apkk1*apkk1)*(C5*C5)*C6*C8*C9*t4*t6 + 
           6.*amLam*(apkk1*apkk1)*(C5*C5)*C6*C8*C9*t4*t6 - 
           6.*(amk*amk)*(apkk1*apkk1)*C5*C6*C8*C9*Fm2*t4*t6 + 
           12.*am*(apkk1*apkk1)*C5*C6*C8*C9*Ekaon*Fm2*t4*t6 + 
           (aml*aml*aml*aml)*(C3*
           ((amk*amk*amk*amk)*
           ((3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4)*(3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4)) + 
           4.*am*Ekaon*
           (3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4)*
           (3.*t2 + 
           am*Ekaon*
           (3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4)) + 
           12.*t1*
           (6.*t2 - 
           1.*((amk*amk) - 2.*am*Ekaon)*
           (3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4)) - 
           1.*(amk*amk)*
           (9.*(C1*C1)*(C4*C4*C4*C4)*
           ((am*am) + (amSig*amSig) + 
           2.*am*(amSig + 2.*Ekaon))*(t3*t3) + 
           (C5*C5)*C6*t4*
           (6.*t2 + 
           (C5*C5)*C6*
           ((am*am) + (amLam*amLam) + 
           2.*am*(amLam + 2.*Ekaon))*t4) + 
           6.*C1*(C4*C4)*t3*
           (3.*t2 + 
           (C5*C5)*C6*
           ((am*am) + amLam*amSig + 
           am*(amLam + amSig + 4.*Ekaon))*t4))
           ) - 
           3.*((amk*amk) - 2.*am*Ekaon)*
           (3.*C1*C4*Fm1*t3 - 1.*C5*C6*Fm2*t4)*
           (C2*C7*t5 - 1.*C8*C9*t6)) + 
           (aml*aml)*(-144.*((1. + f)*(1. + f))*(t1*t1) + 
           72.*(am*am)*C1*C3*(C4*C4)*(Ekaon*Ekaon)*t2*t3 - 
           72.*(am*am)*C1*C3*(C4*C4)*Ekaon*Elep*t2*t3 + 
           72.*(am*am)*C1*C3*(C4*C4)*Ekaon*Enu*t2*t3 + 
           18.*(am*am*am)*apkk1*(C1*C1)*C3*(C4*C4*C4*C4)*Ekaon*
           (t3*t3) - 
           18.*am*(amSig*amSig)*apkk1*(C1*C1)*C3*(C4*C4*C4*C4)*Ekaon*
           (t3*t3) - 
           36.*(am*am)*(C1*C1)*(C4*C4)*(Ekaon*Ekaon)*(t3*t3) + 
           72.*(am*am)*(C1*C1)*(C4*C4*C4)*(Ekaon*Ekaon)*(t3*t3) - 
           36.*(am*am)*(C1*C1)*(C4*C4*C4*C4)*(Ekaon*Ekaon)*(t3*t3) + 
           36.*(am*am)*apkk1*(C1*C1)*C3*(C4*C4*C4*C4)*(Ekaon*Ekaon)*
           (t3*t3) - 
           72.*(am*am*am)*(C1*C1)*C3*(C4*C4*C4*C4)*(Ekaon*Ekaon)*Elep*
           (t3*t3) + 
           72.*(am*am*am)*(C1*C1)*C3*(C4*C4*C4*C4)*(Ekaon*Ekaon)*Enu*
           (t3*t3) + 
           36.*(am*am*am)*(C1*C1)*(C4*C4)*(Ekaon*Ekaon)*Fm1*(t3*t3) - 
           36.*(am*am)*amSig*(C1*C1)*(C4*C4)*(Ekaon*Ekaon)*Fm1*
           (t3*t3) - 
           72.*(am*am*am)*(C1*C1)*(C4*C4*C4)*(Ekaon*Ekaon)*Fm1*(t3*t3) - 
           72.*(am*am)*amSig*(C1*C1)*(C4*C4*C4)*(Ekaon*Ekaon)*Fm1*
           (t3*t3) - 
           72.*(am*am*am)*(C1*C1)*(C4*C4)*(Ekaon*Ekaon)*Enu*(Fm1*Fm1)*
           (t3*t3) + 
           24.*(am*am)*C3*(C5*C5)*C6*(Ekaon*Ekaon)*t2*t4 - 
           24.*(am*am)*C3*(C5*C5)*C6*Ekaon*Elep*t2*t4 + 
           24.*(am*am)*C3*(C5*C5)*C6*Ekaon*Enu*t2*t4 + 
           12.*(am*am*am)*apkk1*C1*C3*(C4*C4)*(C5*C5)*C6*Ekaon*
           t3*t4 - 
           12.*am*amLam*amSig*apkk1*C1*C3*(C4*C4)*(C5*C5)*
           C6*Ekaon*t3*t4 + 
           72.*(am*am)*C1*C4*C5*C6*(Ekaon*Ekaon)*t3*t4 - 
           72.*(am*am)*C1*(C4*C4)*C5*C6*(Ekaon*Ekaon)*t3*t4 + 
           24.*(am*am)*C1*C4*(C5*C5)*C6*(Ekaon*Ekaon)*t3*t4 - 
           24.*(am*am)*C1*(C4*C4)*(C5*C5)*C6*(Ekaon*Ekaon)*t3*
           t4 + 
           24.*(am*am)*apkk1*C1*C3*(C4*C4)*(C5*C5)*C6*
           (Ekaon*Ekaon)*t3*t4 - 
           48.*(am*am*am)*C1*C3*(C4*C4)*(C5*C5)*C6*(Ekaon*Ekaon)*
           Elep*t3*t4 + 
           48.*(am*am*am)*C1*C3*(C4*C4)*(C5*C5)*C6*(Ekaon*Ekaon)*
           Enu*t3*t4 + 
           18.*am*amLam*apkk1*C1*C4*C5*C6*Ekaon*Fm1*
           t3*t4 - 
           18.*am*amSig*apkk1*C1*C4*C5*C6*Ekaon*Fm1*
           t3*t4 + 
           12.*am*amLam*apkk1*C1*C4*(C5*C5)*C6*Ekaon*
           Fm1*t3*t4 - 
           12.*am*amSig*apkk1*C1*C4*(C5*C5)*C6*Ekaon*
           Fm1*t3*t4 - 
           36.*(am*am*am)*C1*C4*C5*C6*(Ekaon*Ekaon)*Fm1*t3*t4 + 
           36.*(am*am)*amLam*C1*C4*C5*C6*(Ekaon*Ekaon)*Fm1*
           t3*t4 - 
           24.*(am*am*am)*C1*C4*(C5*C5)*C6*(Ekaon*Ekaon)*Fm1*t3*
           t4 - 
           12.*(am*am)*amLam*C1*C4*(C5*C5)*C6*(Ekaon*Ekaon)*
           Fm1*t3*t4 - 
           12.*(am*am)*amSig*C1*C4*(C5*C5)*C6*(Ekaon*Ekaon)*
           Fm1*t3*t4 - 
           6.*am*amLam*apkk1*C1*C4*C5*C6*Ekaon*Fm2*t3*
           t4 + 
           6.*am*amSig*apkk1*C1*C4*C5*C6*Ekaon*Fm2*t3*
           t4 + 
           12.*am*amLam*apkk1*C1*(C4*C4)*C5*C6*Ekaon*
           Fm2*t3*t4 - 
           12.*am*amSig*apkk1*C1*(C4*C4)*C5*C6*Ekaon*
           Fm2*t3*t4 - 
           12.*(am*am*am)*C1*C4*C5*C6*(Ekaon*Ekaon)*Fm2*t3*t4 + 
           12.*(am*am)*amSig*C1*C4*C5*C6*(Ekaon*Ekaon)*Fm2*
           t3*t4 + 
           24.*(am*am*am)*C1*(C4*C4)*C5*C6*(Ekaon*Ekaon)*Fm2*t3*
           t4 + 
           12.*(am*am)*amLam*C1*(C4*C4)*C5*C6*(Ekaon*Ekaon)*
           Fm2*t3*t4 + 
           12.*(am*am)*amSig*C1*(C4*C4)*C5*C6*(Ekaon*Ekaon)*
           Fm2*t3*t4 + 
           48.*(am*am*am)*C1*C4*C5*C6*(Ekaon*Ekaon)*Enu*Fm1*Fm2*
           t3*t4 + 
           2.*(am*am*am)*apkk1*C3*(C5*C5*C5*C5)*(C6*C6)*Ekaon*
           (t4*t4) - 
           2.*am*(amLam*amLam)*apkk1*C3*(C5*C5*C5*C5)*(C6*C6)*Ekaon*
           (t4*t4) - 
           36.*(am*am)*(C5*C5)*(C6*C6)*(Ekaon*Ekaon)*(t4*t4) - 
           24.*(am*am)*(C5*C5*C5)*(C6*C6)*(Ekaon*Ekaon)*(t4*t4) - 
           4.*(am*am)*(C5*C5*C5*C5)*(C6*C6)*(Ekaon*Ekaon)*(t4*t4) + 
           4.*(am*am)*apkk1*C3*(C5*C5*C5*C5)*(C6*C6)*(Ekaon*Ekaon)*
           (t4*t4) - 
           8.*(am*am*am)*C3*(C5*C5*C5*C5)*(C6*C6)*(Ekaon*Ekaon)*Elep*
           (t4*t4) + 
           8.*(am*am*am)*C3*(C5*C5*C5*C5)*(C6*C6)*(Ekaon*Ekaon)*Enu*
           (t4*t4) + 
           12.*(am*am*am)*(C5*C5)*(C6*C6)*(Ekaon*Ekaon)*Fm2*(t4*t4) - 
           12.*(am*am)*amLam*(C5*C5)*(C6*C6)*(Ekaon*Ekaon)*Fm2*
           (t4*t4) + 
           8.*(am*am*am)*(C5*C5*C5)*(C6*C6)*(Ekaon*Ekaon)*Fm2*(t4*t4) + 
           8.*(am*am)*amLam*(C5*C5*C5)*(C6*C6)*(Ekaon*Ekaon)*Fm2*
           (t4*t4) - 
           8.*(am*am*am)*(C5*C5)*(C6*C6)*(Ekaon*Ekaon)*Enu*(Fm2*Fm2)*
           (t4*t4) + 
           (amk*amk*amk*amk)*
           (9.*(C1*C1)*(C4*C4)*
           (-1. + 2.*C4 + 
           (C4*C4)*
           (-1. + apkk1*C3 + 
           2.*am*C3*(-1.*Elep + Enu)) - 
           2.*am*Fm1*(-1. + Enu*Fm1))*(t3*t3) + 
           3.*C1*C4*t3*
           (C5*C6*
           (6. + 2.*C5 - 2.*C4*(3. + C5) - 
           6.*am*Fm1 + 3.*amLam*Fm1 - 
           3.*amSig*Fm1 - 2.*am*Fm2 - 
           1.*amLam*Fm2 + amSig*Fm2 + 
           4.*am*Enu*Fm1*Fm2)*t4 + 
           2.*C3*C4*
           (3.*t2 + 
           (C5*C5)*C6*
           (apkk1 - 2.*am*Elep + 2.*am*Enu)*t4
           )) + 
           (C5*C5)*C6*t4*
           (-1.*C6*
           (9. + 6.*C5 + (C5*C5) + 
           2.*am*Fm2*(-3. + Enu*Fm2))*t4 + 
           C3*
           (6.*t2 + 
           (C5*C5)*C6*
           (apkk1 + 2.*am*(-1.*Elep + Enu))*t4
           ))) - 
           12.*t1*
           (apkk1*(2. + f)*
           (3.*amSig*C1*C4*Fm1*t3 - 
           1.*amLam*C5*C6*Fm2*t4) - 
           2.*(am*am)*
           (apkk1*C3*
           (3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4) + 
           Ekaon*
           (-3.*C1*C4*
           (2.*C3*C4*(Elep - 1.*Enu) + Fm1 + 
           2.*f*Fm1)*t3 + 
           C5*C6*
           (2.*C3*C5*(-1.*Elep + Enu) + Fm2 + 
           2.*f*Fm2)*t4)) + 
           (amk*amk)*
           (-3.*C1*C4*
           (-2. + 2.*C4*(1. + f) + 
           f*(-2. + 3.*am*Fm1 + amSig*Fm1))*t3
           - 1.*C5*C6*
           (6. + 2.*C5*(1. + f) - 
           1.*f*(-6. + 3.*am*Fm2 + amLam*Fm2))
           *t4 + 
           C3*
           (6.*t2 + 
           (apkk1 + 2.*am*(-1.*Elep + Enu))*
           (3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4))) - 
           1.*am*
           (-3.*C1*C4*
           (apkk1*(2. + f)*Fm1 + 
           2.*Ekaon*
           (-2. - 2.*f + 2.*C4*(1. + f) + 
           amSig*Fm1))*t3 - 
           1.*C5*C6*
           (-1.*apkk1*(2. + f)*Fm2 + 
           2.*Ekaon*
           (6. + 6.*f + 2.*C5*(1. + f) - 
           1.*amLam*Fm2))*t4 + 
           2.*C3*
           (-6.*Elep*t2 + 6.*Enu*t2 + 
           3.*amSig*apkk1*C1*(C4*C4)*t3 + 
           amLam*apkk1*(C5*C5)*C6*t4 + 
           Ekaon*
           (6.*t2 + 3.*apkk1*C1*(C4*C4)*t3 + 
           apkk1*(C5*C5)*C6*t4)))) - 
           9.*am*apkk1*C1*C2*C4*C7*t3*t5 - 
           9.*amSig*apkk1*C1*C2*C4*C7*t3*t5 + 
           54.*am*apkk1*C1*C2*C4*C7*Ekaon*Fm1*t3*t5 + 
           9.*am*apkk1*C2*C5*C6*C7*t4*t5 + 
           9.*amLam*apkk1*C2*C5*C6*C7*t4*t5 - 
           18.*am*apkk1*C2*C5*C6*C7*Ekaon*Fm2*t4*t5 + 
           9.*am*apkk1*C1*C4*C8*C9*t3*t6 + 
           9.*amSig*apkk1*C1*C4*C8*C9*t3*t6 - 
           54.*am*apkk1*C1*C4*C8*C9*Ekaon*Fm1*t3*t6 - 
           9.*am*apkk1*C5*C6*C8*C9*t4*t6 - 
           9.*amLam*apkk1*C5*C6*C8*C9*t4*t6 + 
           18.*am*apkk1*C5*C6*C8*C9*Ekaon*Fm2*t4*t6 + 
           (amk*amk)*
           (-9.*(amSig*amSig)*(C1*C1)*(C4*C4)*
           (-1. + 2.*C4 + (-1. + apkk1*C3)*(C4*C4))*
           (t3*t3) - 
           9.*amLam*apkk1*C1*C4*C5*C6*Fm1*t3*t4 - 
           6.*amLam*apkk1*C1*C4*(C5*C5)*C6*Fm1*t3*
           t4 + 
           3.*amLam*apkk1*C1*C4*C5*C6*Fm2*t3*t4 - 
           6.*amLam*apkk1*C1*(C4*C4)*C5*C6*Fm2*t3*
           t4 + 9.*(amLam*amLam)*(C5*C5)*(C6*C6)*(t4*t4) + 
           6.*(amLam*amLam)*(C5*C5*C5)*(C6*C6)*(t4*t4) + 
           (amLam*amLam)*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) - 
           1.*(amLam*amLam)*apkk1*C3*(C5*C5*C5*C5)*(C6*C6)*
           (t4*t4) + 
           2.*(am*am*am)*
           (9.*(C1*C1)*(C4*C4)*
           (C3*(C4*C4)*(Elep - 1.*Enu) + 
           Fm1*(2.*C4 + Enu*Fm1))*(t3*t3) + 
           6.*C1*C4*C5*C6*
           (C3*C4*C5*(Elep - 1.*Enu) + 
           C5*Fm1 - 1.*(C4 + Enu*Fm1)*Fm2)*t3*
           t4 + 
           (C5*C5)*(C6*C6)*
           (C3*(C5*C5)*(Elep - 1.*Enu) + 
           Fm2*(-2.*C5 + Enu*Fm2))*(t4*t4)) + 
           (am*am)*
           (-9.*(C1*C1)*(C4*C4)*
           (-1. + 
           (C4*C4)*
           (-1. + 3.*apkk1*C3 - 
           8.*C3*Ekaon*Elep + 
           8.*C3*Ekaon*Enu + 
           4.*amSig*C3*(-1.*Elep + Enu)) - 
           4.*amSig*Enu*(Fm1*Fm1) + 
           C4*
           (2. - 8.*amSig*Fm1 - 4.*Ekaon*Fm1)
           + 2.*Ekaon*Fm1*(3. - 4.*Enu*Fm1))*
           (t3*t3) - 
           6.*C1*C4*C5*C6*
           (3. + C5 - 2.*amLam*C5*Fm1 - 
           2.*amSig*C5*Fm1 - 9.*Ekaon*Fm1 - 
           2.*C5*Ekaon*Fm1 - 3.*Ekaon*Fm2 + 
           2.*amLam*Enu*Fm1*Fm2 + 
           2.*amSig*Enu*Fm1*Fm2 + 
           8.*Ekaon*Enu*Fm1*Fm2 + 
           C4*
           (-3. + 
           C5*
           (-1. + 3.*apkk1*C3 - 
           2.*amSig*C3*Elep - 
           8.*C3*Ekaon*Elep + 
           2.*amSig*C3*Enu + 
           8.*C3*Ekaon*Enu + 
           2.*amLam*C3*(-1.*Elep + Enu)) + 
           2.*amLam*Fm2 + 2.*amSig*Fm2 + 
           2.*Ekaon*Fm2))*t3*t4 + 
           (C5*C5)*(C6*C6)*
           (9. + 
           (C5*C5)*
           (1. - 3.*apkk1*C3 + 
           8.*C3*Ekaon*Elep + 
           4.*amLam*C3*(Elep - 1.*Enu) - 
           8.*C3*Ekaon*Enu) + 
           4.*amLam*Enu*(Fm2*Fm2) + 
           C5*
           (6. - 8.*amLam*Fm2 - 4.*Ekaon*Fm2)
           + 2.*Ekaon*Fm2*(-9. + 4.*Enu*Fm2))*
           (t4*t4)) - 
           27.*apkk1*C1*C2*C4*C7*Fm1*t3*t5 + 
           9.*amLam*C2*C5*C6*C7*t4*t5 + 
           6.*amLam*C2*(C5*C5)*C6*C7*t4*t5 + 
           9.*apkk1*C2*C5*C6*C7*Fm2*t4*t5 + 
           27.*apkk1*C1*C4*C8*C9*Fm1*t3*t6 - 
           9.*amLam*C5*C6*C8*C9*t4*t6 - 
           6.*amLam*(C5*C5)*C6*C8*C9*t4*t6 - 
           9.*apkk1*C5*C6*C8*C9*Fm2*t4*t6 - 
           3.*amSig*C1*C4*t3*
           (2.*amLam*C5*
           (3. + C5 + 
           C4*(-3. - 1.*C5 + apkk1*C3*C5))*C6*
           t4 + 
           apkk1*C5*C6*
           (-3.*Fm1 - 2.*C5*Fm1 + Fm2 - 
           2.*C4*Fm2)*t4 - 
           3.*(-1. + 2.*C4)*
           (C2*C7*t5 - 1.*C8*C9*t6)) + 
           am*
           (18.*(C1*C1)*(C4*C4)*
           (2.*
           (1. - 2.*C4 + 
           (1. - 1.*apkk1*C3)*(C4*C4))*Ekaon + 
           amSig*
           (1. + (1. - 2.*apkk1*C3)*(C4*C4) + 
           Ekaon*Fm1 + 2.*C4*(-1. + Ekaon*Fm1)
           ) + 
           (amSig*amSig)*
           (C3*(C4*C4)*(Elep - 1.*Enu) + 
           Fm1*(2.*C4 + Enu*Fm1)))*(t3*t3) - 
           3.*C1*C4*t3*
           (6.*amSig*C5*C6*t4 - 
           6.*amSig*C4*C5*C6*t4 + 
           2.*amSig*(C5*C5)*C6*t4 - 
           2.*amSig*C4*(C5*C5)*C6*t4 + 
           24.*C5*C6*Ekaon*t4 - 
           24.*C4*C5*C6*Ekaon*t4 + 
           8.*(C5*C5)*C6*Ekaon*t4 - 
           8.*C4*(C5*C5)*C6*Ekaon*t4 - 
           6.*amSig*C5*C6*Ekaon*Fm1*t4 - 
           2.*amSig*(C5*C5)*C6*Ekaon*Fm1*t4 + 
           4.*amSig*C5*C6*Ekaon*Fm2*t4 + 
           2.*amSig*C4*C5*C6*Ekaon*Fm2*t4 - 
           2.*amLam*C5*C6*
           (-3. - 6.*Ekaon*Fm1 + 
           C5*
           (-1. + 2.*amSig*Fm1 + Ekaon*Fm1) + 
           Ekaon*Fm2 - 2.*amSig*Enu*Fm1*Fm2 + 
           C4*
           (3. + C5 - 2.*amSig*Fm2 - 
           1.*Ekaon*Fm2))*t4 + 
           4.*C3*C4*
           (6.*Ekaon*t2 - 3.*Elep*t2 + 
           3.*Enu*t2 + 
           amLam*apkk1*(C5*C5)*C6*t4 + 
           amSig*apkk1*(C5*C5)*C6*t4 + 
           2.*apkk1*(C5*C5)*C6*Ekaon*t4 - 
           1.*amLam*amSig*(C5*C5)*C6*Elep*t4 + 
           amLam*amSig*(C5*C5)*C6*Enu*t4) + 
           3.*C2*C7*t5 - 6.*C2*C4*C7*t5 - 
           3.*C8*C9*t6 + 6.*C4*C8*C9*t6) + 
           C5*C6*t4*
           (36.*C5*C6*Ekaon*t4 + 
           24.*(C5*C5)*C6*Ekaon*t4 + 
           4.*(C5*C5*C5)*C6*Ekaon*t4 + 
           2.*amLam*C5*C6*
           (9. + 6.*C5 + (C5*C5) + 
           3.*Ekaon*Fm2 - 2.*C5*Ekaon*Fm2)*t4
           + 2.*(amLam*amLam)*C5*C6*Fm2*
           (-2.*C5 + Enu*Fm2)*t4 - 
           2.*C3*C5*
           (-6.*Elep*t2 + 6.*Enu*t2 + 
           2.*amLam*apkk1*(C5*C5)*C6*t4 - 
           1.*(amLam*amLam)*(C5*C5)*C6*Elep*t4 + 
           (amLam*amLam)*(C5*C5)*C6*Enu*t4 + 
           2.*Ekaon*
           (6.*t2 + apkk1*(C5*C5)*C6*t4)) + 
           9.*C2*C7*t5 + 6.*C2*C5*C7*t5 - 
           9.*C8*C9*t6 - 6.*C5*C8*C9*t6)))))
           + 2.*(akpk*akpk)*
           (2.*(am*am)*(aml*aml)*Ekaon*
           (9.*(C1*C1)*(C4*C4)*(C3*(C4*C4) - 1.*(Fm1*Fm1))*
           (t3*t3) + 
           6.*C1*C4*C5*C6*(C3*C4*C5 + Fm1*Fm2)*t3*
           t4 + 
           (C5*C5)*(C6*C6)*(C3*(C5*C5) - 1.*(Fm2*Fm2))*(t4*t4))
           + 3.*Elep*
           ((amk*amk)*(3.*C1*C4*Fm1*t3 - 1.*C5*C6*Fm2*t4)*
           (C2*C7*t5 - 1.*C8*C9*t6) + 
           amSig*C1*C4*t3*
           ((amk*amk)*(1. + C4)*C5*C6*Fm2*t4 + 
           Fm1*
           (12.*(-1. + f)*t1 + 
           (amk*amk)*(-3. + C5)*C5*C6*t4) - 
           3.*(-1. + C4)*(C2*C7*t5 - 1.*C8*C9*t6))
           - 1.*amLam*C5*C6*t4*
           ((amk*amk)*C1*C4*(-3. + C5)*Fm1*t3 + 
           Fm2*
           (4.*(-1. + f)*t1 + 
           (amk*amk)*C1*C4*(1. + C4)*t3) + 
           (3. + C5)*(C2*C7*t5 - 1.*C8*C9*t6))) - 
           1.*(aml*aml)*
           (18.*(amSig*amSig)*(C1*C1)*(C4*C4)*Ekaon*
           (C3*(C4*C4) - 1.*(Fm1*Fm1))*(t3*t3) + 
           3.*amSig*C1*C4*t3*
           (2.*C3*C4*
           (-6.*t1 - 3.*t2 + 
           3.*(amk*amk)*C1*(C4*C4)*t3 + 
           (amk*amk)*(C5*C5)*C6*t4 + 
           2.*amLam*(C5*C5)*C6*Ekaon*t4) + 
           Fm1*
           (4.*amLam*C5*C6*Ekaon*Fm2*t4 + 
           (amk*amk)*
           (-6.*C1*C4*Fm1*t3 + 2.*C5*C6*Fm2*t4)
           + 3.*C2*C7*t5 - 3.*C8*C9*t6)) + 
           amLam*C5*C6*t4*
           (2.*C3*C5*
           (-6.*t1 - 3.*t2 + 
           3.*(amk*amk)*C1*(C4*C4)*t3 + 
           (amk*amk)*(C5*C5)*C6*t4 + 
           amLam*(C5*C5)*C6*Ekaon*t4) + 
           Fm2*
           (-2.*amLam*C5*C6*Ekaon*Fm2*t4 + 
           (amk*amk)*
           (6.*C1*C4*Fm1*t3 - 2.*C5*C6*Fm2*t4)
           - 3.*C2*C7*t5 + 3.*C8*C9*t6))) + 
           am*((aml*aml)*
           (-18.*(amk*amk)*(C1*C1)*(C4*C4)*
           (C3*(C4*C4) - 1.*(Fm1*Fm1))*(t3*t3) + 
           C5*C6*t4*
           (C3*
           (6.*C5*(2.*t1 + t2) - 
           2.*(amk*amk)*(C5*C5*C5)*C6*t4) + 
           Fm2*
           (2.*(amk*amk)*C5*C6*Fm2*t4 + 
           3.*C2*C7*t5 - 3.*C8*C9*t6)) + 
           3.*C1*C4*t3*
           (2.*C3*C4*
           (6.*t1 + 3.*t2 - 
           2.*(amk*amk)*(C5*C5)*C6*t4) + 
           Fm1*
           (-4.*(amk*amk)*C5*C6*Fm2*t4 - 
           3.*C2*C7*t5 + 3.*C8*C9*t6))) + 
           3.*Elep*
           (C5*C6*t4*
           (-1.*(3. + C5)*
           (C2*C7*t5 - 1.*C8*C9*t6) + 
           Fm2*
           (4.*t1 - 4.*f*t1 + 
           2.*C2*C7*Ekaon*t5 - 
           2.*C8*C9*Ekaon*t6)) + 
           C1*C4*t3*
           (-2.*amSig*C5*C6*Ekaon*Fm2*t4 - 
           2.*amSig*C4*C5*C6*Ekaon*Fm2*t4 + 
           2.*amLam*(1. + C4)*C5*C6*Ekaon*Fm2*
           t4 + 3.*C2*C7*t5 - 3.*C2*C4*C7*t5 - 
           3.*C8*C9*t6 + 3.*C4*C8*C9*t6 + 
           2.*Fm1*
           (6.*(-1. + f)*t1 + 
           Ekaon*
           (amLam*(-3. + C5)*C5*C6*t4 - 
           1.*amSig*(-3. + C5)*C5*C6*t4 - 
           3.*C2*C7*t5 + 3.*C8*C9*t6)))))) - 
           1.*akpk*(6.*am*(Elep*Elep)*
           (amLam*C5*C6*
           ((amk*amk)*C1*C4*(3. + C5)*Fm1*t3 - 
           1.*Fm2*
           (4.*(1. + f)*t1 - 
           1.*(amk*amk)*C1*(-1. + C4)*C4*t3))*t4 + 
           amSig*C1*C4*t3*
           (-1.*(amk*amk)*(-1. + C4)*C5*C6*Fm2*t4 + 
           Fm1*
           (12.*(1. + f)*t1 - 
           1.*(amk*amk)*C5*(3. + C5)*C6*t4)) - 
           2.*am*
           (2.*C5*C6*(1. + f)*Fm2*t1*t4 + 
           C1*C4*t3*
           ((amLam - 1.*amSig)*(-1. + C4)*C5*C6*
           Ekaon*Fm2*t4 - 
           1.*Fm1*
           (6.*(1. + f)*t1 - 
           1.*(amLam - 1.*amSig)*C5*(3. + C5)*
           C6*Ekaon*t4)))) + 
           (aml*aml*aml*aml)*(-18.*(amSig*amSig)*(C1*C1)*C3*(C4*C4*C4*C4)*Ekaon*
           (t3*t3) + 
           2.*(am*am)*C3*Ekaon*
           ((3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4)*(3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4)) - 
           3.*amSig*C1*C4*t3*
           (2.*C3*C4*
           (-6.*t1 - 3.*t2 + 
           3.*(amk*amk)*C1*(C4*C4)*t3 + 
           (amk*amk)*(C5*C5)*C6*t4 + 
           2.*amLam*(C5*C5)*C6*Ekaon*t4) - 
           3.*C2*C7*Fm1*t5 + 3.*C8*C9*Fm1*t6) - 
           1.*amLam*C5*C6*t4*
           (2.*C3*C5*
           (-6.*t1 - 3.*t2 + 
           3.*(amk*amk)*C1*(C4*C4)*t3 + 
           (amk*amk)*(C5*C5)*C6*t4 + 
           amLam*(C5*C5)*C6*Ekaon*t4) + 
           3.*Fm2*(C2*C7*t5 - 1.*C8*C9*t6)) + 
           am*(-18.*(amk*amk)*(C1*C1)*C3*(C4*C4*C4*C4)*(t3*t3) + 
           C5*C6*t4*
           (C3*
           (6.*C5*(2.*t1 + t2) - 
           2.*(amk*amk)*(C5*C5*C5)*C6*t4) - 
           3.*C2*C7*Fm2*t5 + 3.*C8*C9*Fm2*t6) + 
           3.*C1*C4*t3*
           (2.*C3*C4*
           (6.*t1 + 3.*t2 - 
           2.*(amk*amk)*(C5*C5)*C6*t4) + 
           3.*Fm1*(C2*C7*t5 - 1.*C8*C9*t6)))) + 
           (aml*aml)*(72.*amSig*C1*C4*t1*t3 - 
           72.*amSig*C1*(C4*C4)*t1*t3 + 
           72.*amSig*C1*C4*f*t1*t3 - 
           72.*amSig*C1*(C4*C4)*f*t1*t3 - 
           72.*(amk*amk)*C1*C4*Fm1*t1*t3 - 
           72.*amSig*C1*C4*Ekaon*Fm1*t1*t3 + 
           36.*(amk*amk)*C1*C4*f*Fm1*t1*t3 + 
           72.*amSig*C1*C4*Ekaon*f*Fm1*t1*t3 + 
           18.*(amk*amk)*amSig*(C1*C1)*(C4*C4)*(t3*t3) - 
           36.*(amk*amk)*amSig*(C1*C1)*(C4*C4*C4)*(t3*t3) + 
           18.*(amk*amk)*amSig*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) + 
           18.*(amSig*amSig)*(C1*C1)*(C4*C4)*Ekaon*(t3*t3) - 
           36.*(amSig*amSig)*(C1*C1)*(C4*C4*C4)*Ekaon*(t3*t3) + 
           18.*(amSig*amSig)*(C1*C1)*(C4*C4*C4*C4)*Ekaon*(t3*t3) + 
           9.*(amk*amk*amk*amk)*(C1*C1)*(C4*C4)*Fm1*(t3*t3) + 
           9.*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4)*Fm1*(t3*t3) + 
           18.*(amk*amk*amk*amk)*(C1*C1)*(C4*C4*C4)*Fm1*(t3*t3) + 
           18.*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4*C4)*Fm1*
           (t3*t3) + 
           18.*(amk*amk)*amSig*(C1*C1)*(C4*C4)*Ekaon*Fm1*
           (t3*t3) + 
           36.*(amk*amk)*amSig*(C1*C1)*(C4*C4*C4)*Ekaon*Fm1*
           (t3*t3) - 
           18.*(amk*amk*amk*amk)*(C1*C1)*(C4*C4)*Enu*(Fm1*Fm1)*(t3*t3) + 
           18.*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4)*Enu*(Fm1*Fm1)*
           (t3*t3) - 72.*amLam*C5*C6*t1*t4 - 
           24.*amLam*(C5*C5)*C6*t1*t4 - 
           72.*amLam*C5*C6*f*t1*t4 - 
           24.*amLam*(C5*C5)*C6*f*t1*t4 + 
           24.*(amk*amk)*C5*C6*Fm2*t1*t4 + 
           24.*amLam*C5*C6*Ekaon*Fm2*t1*t4 - 
           12.*(amk*amk)*C5*C6*f*Fm2*t1*t4 - 
           24.*amLam*C5*C6*Ekaon*f*Fm2*t1*t4 - 
           18.*(amk*amk)*amLam*C1*C4*C5*C6*t3*t4 - 
           18.*(amk*amk)*amSig*C1*C4*C5*C6*t3*t4 + 
           18.*(amk*amk)*amLam*C1*(C4*C4)*C5*C6*t3*t4 + 
           18.*(amk*amk)*amSig*C1*(C4*C4)*C5*C6*t3*t4 - 
           6.*(amk*amk)*amLam*C1*C4*(C5*C5)*C6*t3*t4 - 
           6.*(amk*amk)*amSig*C1*C4*(C5*C5)*C6*t3*t4 + 
           6.*(amk*amk)*amLam*C1*(C4*C4)*(C5*C5)*C6*t3*t4 + 
           6.*(amk*amk)*amSig*C1*(C4*C4)*(C5*C5)*C6*t3*t4 - 
           36.*amLam*amSig*C1*C4*C5*C6*Ekaon*t3*t4 + 
           36.*amLam*amSig*C1*(C4*C4)*C5*C6*Ekaon*t3*
           t4 - 
           12.*amLam*amSig*C1*C4*(C5*C5)*C6*Ekaon*t3*
           t4 + 
           12.*amLam*amSig*C1*(C4*C4)*(C5*C5)*C6*Ekaon*t3*
           t4 - 9.*(amk*amk*amk*amk)*C1*C4*C5*C6*Fm1*t3*t4 - 
           9.*(amk*amk)*amLam*amSig*C1*C4*C5*C6*Fm1*t3*
           t4 + 6.*(amk*amk*amk*amk)*C1*C4*(C5*C5)*C6*Fm1*t3*t4 + 
           6.*(amk*amk)*amLam*amSig*C1*C4*(C5*C5)*C6*Fm1*
           t3*t4 - 
           18.*(amk*amk)*amSig*C1*C4*C5*C6*Ekaon*Fm1*t3*
           t4 + 
           6.*(amk*amk)*amLam*C1*C4*(C5*C5)*C6*Ekaon*Fm1*
           t3*t4 + 
           6.*(amk*amk)*amSig*C1*C4*(C5*C5)*C6*Ekaon*Fm1*
           t3*t4 - 3.*(amk*amk*amk*amk)*C1*C4*C5*C6*Fm2*t3*t4 - 
           3.*(amk*amk)*amLam*amSig*C1*C4*C5*C6*Fm2*t3*
           t4 - 6.*(amk*amk*amk*amk)*C1*(C4*C4)*C5*C6*Fm2*t3*t4 - 
           6.*(amk*amk)*amLam*amSig*C1*(C4*C4)*C5*C6*Fm2*
           t3*t4 - 
           6.*(amk*amk)*amLam*C1*C4*C5*C6*Ekaon*Fm2*t3*
           t4 - 
           6.*(amk*amk)*amLam*C1*(C4*C4)*C5*C6*Ekaon*Fm2*
           t3*t4 - 
           6.*(amk*amk)*amSig*C1*(C4*C4)*C5*C6*Ekaon*Fm2*
           t3*t4 + 
           12.*(amk*amk*amk*amk)*C1*C4*C5*C6*Enu*Fm1*Fm2*t3*t4 - 
           12.*(amk*amk)*amLam*amSig*C1*C4*C5*C6*Enu*Fm1*
           Fm2*t3*t4 + 
           18.*(amk*amk)*amLam*(C5*C5)*(C6*C6)*(t4*t4) + 
           12.*(amk*amk)*amLam*(C5*C5*C5)*(C6*C6)*(t4*t4) + 
           2.*(amk*amk)*amLam*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) + 
           18.*(amLam*amLam)*(C5*C5)*(C6*C6)*Ekaon*(t4*t4) + 
           12.*(amLam*amLam)*(C5*C5*C5)*(C6*C6)*Ekaon*(t4*t4) + 
           2.*(amLam*amLam)*(C5*C5*C5*C5)*(C6*C6)*Ekaon*(t4*t4) + 
           3.*(amk*amk*amk*amk)*(C5*C5)*(C6*C6)*Fm2*(t4*t4) + 
           3.*(amk*amk)*(amLam*amLam)*(C5*C5)*(C6*C6)*Fm2*(t4*t4) - 
           2.*(amk*amk*amk*amk)*(C5*C5*C5)*(C6*C6)*Fm2*(t4*t4) - 
           2.*(amk*amk)*(amLam*amLam)*(C5*C5*C5)*(C6*C6)*Fm2*(t4*t4) + 
           6.*(amk*amk)*amLam*(C5*C5)*(C6*C6)*Ekaon*Fm2*
           (t4*t4) - 
           4.*(amk*amk)*amLam*(C5*C5*C5)*(C6*C6)*Ekaon*Fm2*
           (t4*t4) - 
           2.*(amk*amk*amk*amk)*(C5*C5)*(C6*C6)*Enu*(Fm2*Fm2)*(t4*t4) + 
           2.*(amk*amk)*(amLam*amLam)*(C5*C5)*(C6*C6)*Enu*(Fm2*Fm2)*
           (t4*t4) - 
           2.*(am*am*am)*Ekaon*
           (3.*C1*C4*Fm1*t3 - 1.*C5*C6*Fm2*t4)*
           (3.*C1*C4*(-1. + 2.*C4 + 2.*Enu*Fm1)*t3 + 
           C5*C6*(3. + 2.*C5 - 2.*Enu*Fm2)*t4) - 
           1.*(am*am)*
           (9.*(C1*C1)*(C4*C4)*
           (2.*((-1.+C4)*(-1.+C4))*Ekaon + 
           8.*(Ekaon*Ekaon)*Enu*(Fm1*Fm1) + 
           (amk*amk)*Fm1*(1. - 6.*C4 - 6.*Enu*Fm1))
           *(t3*t3) + 
           C5*C6*t4*
           (2.*C5*((3. + C5)*(3. + C5))*C6*Ekaon*t4 + 
           2.*C5*C6*(-3.*(amk*amk) + 4.*(Ekaon*Ekaon))*
           Enu*(Fm2*Fm2)*t4 - 
           3.*Fm2*
           (8.*(1. + f)*t1 - 
           1.*(amk*amk)*C5*(1. + 2.*C5)*C6*t4))
           + 3.*C1*C4*t3*
           (C5*C6*
           ((amk*amk)*(-1. + 6.*C4)*Fm2 + 
           2.*Ekaon*
           (-6. - 2.*C5 + 2.*C4*(3. + C5) + 
           amLam*Fm2 - 1.*amSig*Fm2))*t4 + 
           Fm1*
           (24.*(1. + f)*t1 - 
           1.*C5*C6*
           (3.*(amk*amk)*
           (1. + 2.*C5 - 4.*Enu*Fm2) + 
           2.*Ekaon*
           (3.*amLam - 3.*amSig + 
           8.*Ekaon*Enu*Fm2))*t4))) + 
           2.*C3*
           ((3.*am*C1*(C4*C4)*t3 + 
           3.*amSig*C1*(C4*C4)*t3 + 
           am*(C5*C5)*C6*t4 + amLam*(C5*C5)*C6*t4)*
           (-3.*((amk*amk) - 2.*am*Ekaon)*t2 + 
           2.*apkk1*
           (6.*t1 + 3.*t2 - 
           3.*(amk*amk)*C1*(C4*C4)*t3 + 
           3.*am*C1*(C4*C4)*Ekaon*t3 - 
           3.*amSig*C1*(C4*C4)*Ekaon*t3 - 
           1.*(amk*amk)*(C5*C5)*C6*t4 + 
           am*(C5*C5)*C6*Ekaon*t4 - 
           1.*amLam*(C5*C5)*C6*Ekaon*t4)) + 
           Enu*
           (2.*(am*am*am)*Ekaon*
           ((3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4)*(3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4)) - 
           1.*(am*am)*
           (3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4)*
           (-6.*t2 + 
           (3.*(amk*amk) - 4.*(Ekaon*Ekaon))*
           (3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4)) + 
           (amk*amk)*
           (9.*((amk*amk) - 1.*(amSig*amSig))*(C1*C1)*
           (C4*C4*C4*C4)*(t3*t3) + 
           (C5*C5)*C6*t4*
           (-6.*t2 + 
           ((amk*amk) - 1.*(amLam*amLam))*(C5*C5)*C6*t4)
           - 6.*C1*(C4*C4)*t3*
           (3.*t2 + 
           (-1.*(amk*amk) + amLam*amSig)*(C5*C5)*
           C6*t4)) - 
           2.*am*
           (18.*(amk*amk)*(C1*C1)*(C4*C4*C4*C4)*Ekaon*
           (t3*t3) + 
           9.*(amSig*amSig)*(C1*C1)*(C4*C4*C4*C4)*Ekaon*
           (t3*t3) + 
           6.*C1*(C4*C4)*t3*
           (-3.*Ekaon*t2 + 
           (amk*amk)*amLam*(C5*C5)*C6*t4 + 
           2.*(amk*amk)*(C5*C5)*C6*Ekaon*t4) + 
           3.*amSig*C1*(C4*C4)*t3*
           (-3.*t2 + 6.*(amk*amk)*C1*(C4*C4)*t3 + 
           2.*(amk*amk)*(C5*C5)*C6*t4 + 
           2.*amLam*(C5*C5)*C6*Ekaon*t4) + 
           (C5*C5)*C6*t4*
           (-3.*amLam*t2 - 6.*Ekaon*t2 + 
           2.*(amk*amk)*amLam*(C5*C5)*C6*t4 + 
           2.*(amk*amk)*(C5*C5)*C6*Ekaon*t4 + 
           (amLam*amLam)*(C5*C5)*C6*Ekaon*t4)) + 
           12.*t1*
           (6.*t2 + 
           (am*am)*
           (3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4) - 
           1.*(amk*amk)*
           (3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4) + 
           am*
           (3.*amSig*C1*(C4*C4)*t3 + 
           6.*C1*(C4*C4)*Ekaon*t3 + 
           (C5*C5)*C6*(amLam + 2.*Ekaon)*t4))))
           - 108.*C2*C7*f*t1*t5 - 
           27.*(amk*amk)*C1*C2*C4*C7*t3*t5 - 
           36.*amSig*C1*C2*C4*C7*Ekaon*t3*t5 + 
           18.*amSig*C1*C2*(C4*C4)*C7*Ekaon*t3*t5 + 
           18.*amSig*apkk1*C1*C2*C4*C7*Fm1*t3*t5 + 
           18.*(amk*amk)*C1*C2*C4*C7*Enu*Fm1*t3*t5 + 
           27.*(amk*amk)*C2*C5*C6*C7*t4*t5 + 
           36.*amLam*C2*C5*C6*C7*Ekaon*t4*t5 + 
           6.*amLam*C2*(C5*C5)*C6*C7*Ekaon*t4*t5 - 
           6.*amLam*apkk1*C2*C5*C6*C7*Fm2*t4*t5 - 
           6.*(amk*amk)*C2*C5*C6*C7*Enu*Fm2*t4*t5 + 
           18.*(C2*C2)*(C7*C7)*Ekaon*(t5*t5) - 
           18.*(C2*C2)*(C7*C7)*Enu*(t5*t5) + 
           108.*C8*C9*f*t1*t6 + 
           27.*(amk*amk)*C1*C4*C8*C9*t3*t6 + 
           36.*amSig*C1*C4*C8*C9*Ekaon*t3*t6 - 
           18.*amSig*C1*(C4*C4)*C8*C9*Ekaon*t3*t6 - 
           18.*amSig*apkk1*C1*C4*C8*C9*Fm1*t3*t6 - 
           18.*(amk*amk)*C1*C4*C8*C9*Enu*Fm1*t3*t6 - 
           27.*(amk*amk)*C5*C6*C8*C9*t4*t6 - 
           36.*amLam*C5*C6*C8*C9*Ekaon*t4*t6 - 
           6.*amLam*(C5*C5)*C6*C8*C9*Ekaon*t4*t6 + 
           6.*amLam*apkk1*C5*C6*C8*C9*Fm2*t4*t6 + 
           6.*(amk*amk)*C5*C6*C8*C9*Enu*Fm2*t4*t6 - 
           36.*C2*C7*C8*C9*Ekaon*t5*t6 + 
           36.*C2*C7*C8*C9*Enu*t5*t6 + 
           18.*(C8*C8)*(C9*C9)*Ekaon*(t6*t6) - 
           18.*(C8*C8)*(C9*C9)*Enu*(t6*t6) + 
           am*(18.*(C1*C1)*(C4*C4)*
           (amSig*Ekaon*Fm1*
           (-2.*(Ekaon + 2.*C4*Ekaon) + 
           amSig*(-1. + 2.*C4 + 2.*Enu*Fm1))
           + (amk*amk)*
           (1. + (C4*C4) + 4.*amSig*Enu*(Fm1*Fm1) + 
           C4*
           (-2. + 4.*amSig*Fm1 - 2.*Ekaon*Fm1)
           + Ekaon*Fm1*(-1. + 4.*Enu*Fm1)))*
           (t3*t3) + 
           2.*C5*C6*t4*
           (-12.*
           (3. + C5 + 3.*f + C5*f - 
           1.*amLam*Fm2 + Ekaon*Fm2 - 
           1.*amLam*f*Fm2)*t1 - 
           3.*(amLam*amLam)*C5*C6*Ekaon*Fm2*t4 - 
           2.*(amLam*amLam)*(C5*C5)*C6*Ekaon*Fm2*t4 - 
           6.*amLam*C5*C6*(Ekaon*Ekaon)*Fm2*t4 + 
           4.*amLam*(C5*C5)*C6*(Ekaon*Ekaon)*Fm2*t4 + 
           2.*(amLam*amLam)*C5*C6*Ekaon*Enu*(Fm2*Fm2)*
           t4 + 
           (amk*amk)*C5*C6*
           (9. + (C5*C5) + 4.*amLam*Enu*(Fm2*Fm2) + 
           C5*
           (6. - 4.*amLam*Fm2 + 2.*Ekaon*Fm2)
           + Ekaon*Fm2*(-3. + 4.*Enu*Fm2))*t4
           - 9.*C2*C7*Ekaon*t5 + 
           3.*C2*C5*C7*Ekaon*t5 - 
           3.*apkk1*C2*C7*Fm2*t5 + 
           6.*C2*C7*Ekaon*Enu*Fm2*t5 + 
           9.*C8*C9*Ekaon*t6 - 
           3.*C5*C8*C9*Ekaon*t6 + 
           3.*apkk1*C8*C9*Fm2*t6 - 
           6.*C8*C9*Ekaon*Enu*Fm2*t6) - 
           3.*C1*C4*t3*
           (24.*
           (-1. + C4 - 1.*f + C4*f + 
           amSig*Fm1 - 1.*Ekaon*Fm1 + 
           amSig*f*Fm1)*t1 + 
           (amk*amk)*C5*C6*
           (12. + 3.*amLam*Fm1 - 
           3.*amSig*Fm1 - 6.*Ekaon*Fm1 - 
           4.*C5*
           (-1. + amLam*Fm1 + amSig*Fm1 - 
           1.*Ekaon*Fm1) - 1.*amLam*Fm2 + 
           amSig*Fm2 - 2.*Ekaon*Fm2 + 
           8.*amLam*Enu*Fm1*Fm2 + 
           8.*amSig*Enu*Fm1*Fm2 + 
           16.*Ekaon*Enu*Fm1*Fm2 - 
           4.*C4*
           (3. + C5 - 1.*amLam*Fm2 - 
           1.*amSig*Fm2 + Ekaon*Fm2))*t4 - 
           2.*
           (2.*amSig*C5*C6*(Ekaon*Ekaon)*
           (3.*Fm1 - 1.*C5*Fm1 + C4*Fm2)*t4 + 
           amLam*C5*C6*Ekaon*
           (2.*Ekaon*
           (-1.*C5*Fm1 + Fm2 + C4*Fm2) + 
           amSig*
           (Fm2 - 2.*C4*Fm2 + 
           Fm1*(3. + 2.*C5 - 4.*Enu*Fm2)))*t4
           + 3.*
           (Ekaon + C4*Ekaon + apkk1*Fm1 - 
           2.*Ekaon*Enu*Fm1)*
           (C2*C7*t5 - 1.*C8*C9*t6))))) + 
           Elep*(144.*((-1. + f)*(-1. + f))*(t1*t1) + 
           36.*(amk*amk)*(aml*aml)*C1*C3*(C4*C4)*t2*t3 + 
           9.*(amk*amk*amk*amk)*(C1*C1)*(C4*C4)*(t3*t3) - 
           9.*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4)*(t3*t3) + 
           18.*(amk*amk*amk*amk)*(C1*C1)*(C4*C4*C4)*(t3*t3) - 
           18.*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4*C4)*(t3*t3) + 
           9.*(amk*amk*amk*amk)*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) - 
           9.*(amk*amk)*(amSig*amSig)*(C1*C1)*(C4*C4*C4*C4)*(t3*t3) - 
           9.*(amk*amk*amk*amk)*(aml*aml)*(C1*C1)*C3*(C4*C4*C4*C4)*(t3*t3) + 
           9.*(amk*amk)*(aml*aml)*(amSig*amSig)*(C1*C1)*C3*(C4*C4*C4*C4)*
           (t3*t3) + 
           12.*(amk*amk)*(aml*aml)*C3*(C5*C5)*C6*t2*t4 - 
           18.*(amk*amk*amk*amk)*C1*C4*C5*C6*t3*t4 + 
           18.*(amk*amk)*amLam*amSig*C1*C4*C5*C6*t3*t4 - 
           18.*(amk*amk*amk*amk)*C1*(C4*C4)*C5*C6*t3*t4 + 
           18.*(amk*amk)*amLam*amSig*C1*(C4*C4)*C5*C6*t3*
           t4 + 6.*(amk*amk*amk*amk)*C1*C4*(C5*C5)*C6*t3*t4 - 
           6.*(amk*amk)*amLam*amSig*C1*C4*(C5*C5)*C6*t3*
           t4 + 6.*(amk*amk*amk*amk)*C1*(C4*C4)*(C5*C5)*C6*t3*t4 - 
           6.*(amk*amk)*amLam*amSig*C1*(C4*C4)*(C5*C5)*C6*t3*
           t4 - 
           6.*(amk*amk*amk*amk)*(aml*aml)*C1*C3*(C4*C4)*(C5*C5)*C6*t3*
           t4 + 
           6.*(amk*amk)*(aml*aml)*amLam*amSig*C1*C3*(C4*C4)*
           (C5*C5)*C6*t3*t4 - 
           9.*(amk*amk)*(aml*aml)*amLam*C1*C4*C5*C6*Fm1*t3*
           t4 + 
           9.*(amk*amk)*(aml*aml)*amSig*C1*C4*C5*C6*Fm1*t3*
           t4 - 
           18.*(amk*amk)*amLam*apkk1*C1*C4*C5*C6*Fm1*t3*
           t4 + 
           18.*(amk*amk)*amSig*apkk1*C1*C4*C5*C6*Fm1*t3*
           t4 - 
           6.*(amk*amk)*(aml*aml)*amLam*C1*C4*(C5*C5)*C6*Fm1*
           t3*t4 + 
           6.*(amk*amk)*(aml*aml)*amSig*C1*C4*(C5*C5)*C6*Fm1*
           t3*t4 - 
           6.*(amk*amk)*amLam*apkk1*C1*C4*(C5*C5)*C6*Fm1*
           t3*t4 + 
           6.*(amk*amk)*amSig*apkk1*C1*C4*(C5*C5)*C6*Fm1*
           t3*t4 + 
           3.*(amk*amk)*(aml*aml)*amLam*C1*C4*C5*C6*Fm2*t3*
           t4 - 
           3.*(amk*amk)*(aml*aml)*amSig*C1*C4*C5*C6*Fm2*t3*
           t4 + 
           6.*(amk*amk)*amLam*apkk1*C1*C4*C5*C6*Fm2*t3*
           t4 - 
           6.*(amk*amk)*amSig*apkk1*C1*C4*C5*C6*Fm2*t3*
           t4 - 
           6.*(amk*amk)*(aml*aml)*amLam*C1*(C4*C4)*C5*C6*Fm2*
           t3*t4 + 
           6.*(amk*amk)*(aml*aml)*amSig*C1*(C4*C4)*C5*C6*Fm2*
           t3*t4 - 
           6.*(amk*amk)*amLam*apkk1*C1*(C4*C4)*C5*C6*Fm2*
           t3*t4 + 
           6.*(amk*amk)*amSig*apkk1*C1*(C4*C4)*C5*C6*Fm2*
           t3*t4 + 9.*(amk*amk*amk*amk)*(C5*C5)*(C6*C6)*(t4*t4) - 
           9.*(amk*amk)*(amLam*amLam)*(C5*C5)*(C6*C6)*(t4*t4) - 
           6.*(amk*amk*amk*amk)*(C5*C5*C5)*(C6*C6)*(t4*t4) + 
           6.*(amk*amk)*(amLam*amLam)*(C5*C5*C5)*(C6*C6)*(t4*t4) + 
           (amk*amk*amk*amk)*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) - 
           1.*(amk*amk)*(amLam*amLam)*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) - 
           1.*(amk*amk*amk*amk)*(aml*aml)*C3*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) + 
           (amk*amk)*(aml*aml)*(amLam*amLam)*C3*(C5*C5*C5*C5)*(C6*C6)*
           (t4*t4) - 
           2.*(am*am*am)*Ekaon*
           (9.*(C1*C1)*(C4*C4)*
           (-1. + 2.*C4 + 
           (-1. + (aml*aml)*C3)*(C4*C4))*(t3*t3) + 
           6.*C1*C4*C5*
           (3. + C5 + 
           C4*(-3. - 1.*C5 + (aml*aml)*C3*C5))*C6*
           t3*t4 + 
           (C5*C5)*
           (-9. - 6.*C5 + 
           (-1. + (aml*aml)*C3)*(C5*C5))*(C6*C6)*(t4*t4))
           - 1.*(am*am)*
           (3.*(amk*amk)*
           (3.*(C1*C1)*(C4*C4)*
           (3. - 2.*C4 + 3.*(C4*C4))*(t3*t3) + 
           2.*C1*C4*C5*
           (-9. - 1.*C5 + 3.*C4*(1. + C5))*C6*
           t3*t4 + 
           (C5*C5)*(9. + 2.*C5 + (C5*C5))*(C6*C6)*
           (t4*t4)) - 
           4.*Ekaon*
           (9.*(C1*C1)*(C4*C4)*((1.+C4)*(1.+C4))*Ekaon*
           (t3*t3) + 
           3.*C1*C4*C5*C6*
           (2.*(1. + C4)*(-3. + C5)*Ekaon + 
           amSig*
           (3.*C4 + C5 + 3.*Enu*Fm1 - 
           1.*C5*Enu*Fm1 - 1.*Enu*Fm2 - 
           1.*C4*Enu*Fm2) + 
           amLam*
           (-1.*C5 - 3.*Enu*Fm1 + 
           C5*Enu*Fm1 + Enu*Fm2 + 
           C4*(-3. + Enu*Fm2)))*t3*t4 + 
           ((-3. + C5)*(-3. + C5))*(C5*C5)*(C6*C6)*Ekaon*(t4*t4)
           ) - 
           1.*(aml*aml)*C3*
           (3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4)*
           (-12.*t2 + 
           (3.*(amk*amk) - 4.*(Ekaon*Ekaon))*
           (3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4))) - 
           12.*t1*
           (-2.*
           ((amk*amk)*(-1. + f)*
           (3.*C1*C4*(1. + C4)*t3 + 
           (-3. + C5)*C5*C6*t4) - 
           1.*apkk1*(1. + f)*
           (3.*amSig*C1*C4*Fm1*t3 - 
           1.*amLam*C5*C6*Fm2*t4) + 
           (am*am)*
           (3.*C1*C4*
           (C4 - 1.*f + 2.*C4*f - 
           1.*Enu*Fm1 + Enu*f*Fm1)*t3 + 
           C5*C6*
           (C5 + 3.*f + 2.*C5*f + Enu*Fm2 - 
           1.*Enu*f*Fm2)*t4) + 
           am*
           (-3.*C1*C4*
           (2.*(1. + C4)*Ekaon*(-1. + f) + 
           apkk1*(1. + f)*Fm1)*t3 + 
           3.*amSig*C1*C4*
           (C4 + 2.*C4*f - 1.*Enu*Fm1 + 
           f*(-1. + Enu*Fm1))*t3 + 
           C5*C6*
           (-2.*(-3. + C5)*Ekaon*(-1. + f) + 
           apkk1*(1. + f)*Fm2 + 
           amLam*
           (C5 + 3.*f + 2.*C5*f + Enu*Fm2 - 
           1.*Enu*f*Fm2))*t4)) + 
           (aml*aml)*
           ((2. + f)*
           (3.*am*C1*C4*Fm1*t3 + 
           3.*amSig*C1*C4*Fm1*t3 - 
           1.*am*C5*C6*Fm2*t4 - 
           1.*amLam*C5*C6*Fm2*t4) + 
           C3*
           (12.*t2 - 
           1.*((amk*amk) - 2.*am*Ekaon)*
           (3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4))))
           + 18.*(amk*amk)*amSig*C1*C2*C4*C7*t3*t5 - 
           9.*(aml*aml)*amSig*C1*C2*C4*C7*t3*t5 - 
           18.*amSig*apkk1*C1*C2*C4*C7*t3*t5 - 
           18.*amSig*apkk1*C1*C2*(C4*C4)*C7*t3*t5 - 
           9.*(amk*amk)*(aml*aml)*C1*C2*C4*C7*Fm1*t3*t5 - 
           18.*(amk*amk)*apkk1*C1*C2*C4*C7*Fm1*t3*t5 - 
           18.*(amk*amk)*amLam*C2*C5*C6*C7*t4*t5 + 
           9.*(aml*aml)*amLam*C2*C5*C6*C7*t4*t5 + 
           18.*amLam*apkk1*C2*C5*C6*C7*t4*t5 - 
           6.*amLam*apkk1*C2*(C5*C5)*C6*C7*t4*t5 + 
           3.*(amk*amk)*(aml*aml)*C2*C5*C6*C7*Fm2*t4*t5 + 
           6.*(amk*amk)*apkk1*C2*C5*C6*C7*Fm2*t4*t5 + 
           18.*(aml*aml)*(C2*C2)*(C7*C7)*(t5*t5) + 
           36.*apkk1*(C2*C2)*(C7*C7)*(t5*t5) - 
           18.*(amk*amk)*amSig*C1*C4*C8*C9*t3*t6 + 
           9.*(aml*aml)*amSig*C1*C4*C8*C9*t3*t6 + 
           18.*amSig*apkk1*C1*C4*C8*C9*t3*t6 + 
           18.*amSig*apkk1*C1*(C4*C4)*C8*C9*t3*t6 + 
           9.*(amk*amk)*(aml*aml)*C1*C4*C8*C9*Fm1*t3*t6 + 
           18.*(amk*amk)*apkk1*C1*C4*C8*C9*Fm1*t3*t6 + 
           18.*(amk*amk)*amLam*C5*C6*C8*C9*t4*t6 - 
           9.*(aml*aml)*amLam*C5*C6*C8*C9*t4*t6 - 
           18.*amLam*apkk1*C5*C6*C8*C9*t4*t6 + 
           6.*amLam*apkk1*(C5*C5)*C6*C8*C9*t4*t6 - 
           3.*(amk*amk)*(aml*aml)*C5*C6*C8*C9*Fm2*t4*t6 - 
           6.*(amk*amk)*apkk1*C5*C6*C8*C9*Fm2*t4*t6 - 
           36.*(aml*aml)*C2*C7*C8*C9*t5*t6 - 
           72.*apkk1*C2*C7*C8*C9*t5*t6 + 
           18.*(aml*aml)*(C8*C8)*(C9*C9)*(t6*t6) + 
           36.*apkk1*(C8*C8)*(C9*C9)*(t6*t6) - 
           1.*am*
           (-1.*(aml*aml)*
           (36.*(amk*amk)*(C1*C1)*C3*(C4*C4*C4*C4)*Ekaon*
           (t3*t3) + 
           18.*(amSig*amSig)*(C1*C1)*C3*(C4*C4*C4*C4)*Ekaon*
           (t3*t3) + 
           6.*amSig*C1*C4*t3*
           (C5*C6*Ekaon*
           (-1.*(3. + 2.*C5)*Fm1 + Fm2 - 
           2.*C4*Fm2)*t4 + 
           2.*C3*C4*
           (-3.*t2 + 
           amLam*(C5*C5)*C6*Ekaon*t4 + 
           (amk*amk)*
           (3.*C1*(C4*C4)*t3 + (C5*C5)*C6*t4))) + 
           3.*C1*C4*t3*
           (2.*amLam*C5*C6*Ekaon*
           (3.*Fm1 + 2.*C5*Fm1 - 1.*Fm2 + 
           2.*C4*Fm2)*t4 + 
           4.*C3*C4*
           (-6.*Ekaon*t2 + 
           (amk*amk)*amLam*(C5*C5)*C6*t4 + 
           2.*(amk*amk)*(C5*C5)*C6*Ekaon*t4) + 
           3.*(-1. + 2.*Ekaon*Fm1)*
           (C2*C7*t5 - 1.*C8*C9*t6)) + 
           C5*C6*t4*
           (2.*(amLam*amLam)*C3*(C5*C5*C5)*C6*Ekaon*t4 + 
           4.*C3*C5*Ekaon*
           (-6.*t2 + (amk*amk)*(C5*C5)*C6*t4) + 
           4.*amLam*C3*C5*
           (-3.*t2 + (amk*amk)*(C5*C5)*C6*t4) - 
           3.*(-3. + 2.*Ekaon*Fm2)*
           (C2*C7*t5 - 1.*C8*C9*t6))) + 
           2.*
           (9.*(amSig*amSig)*(C1*C1)*(C4*C4)*Ekaon*(t3*t3) - 
           18.*(amSig*amSig)*(C1*C1)*(C4*C4*C4)*Ekaon*
           (t3*t3) + 
           9.*(amSig*amSig)*(C1*C1)*(C4*C4*C4*C4)*Ekaon*
           (t3*t3) - 
           18.*amLam*amSig*C1*C4*C5*C6*Ekaon*t3*
           t4 + 
           18.*amLam*amSig*C1*(C4*C4)*C5*C6*Ekaon*
           t3*t4 - 
           6.*amLam*amSig*C1*C4*(C5*C5)*C6*Ekaon*
           t3*t4 + 
           6.*amLam*amSig*C1*(C4*C4)*(C5*C5)*C6*
           Ekaon*t3*t4 - 
           18.*amLam*apkk1*C1*C4*C5*C6*Ekaon*
           Fm1*t3*t4 + 
           18.*amSig*apkk1*C1*C4*C5*C6*Ekaon*
           Fm1*t3*t4 - 
           6.*amLam*apkk1*C1*C4*(C5*C5)*C6*Ekaon*
           Fm1*t3*t4 + 
           6.*amSig*apkk1*C1*C4*(C5*C5)*C6*Ekaon*
           Fm1*t3*t4 + 
           6.*amLam*apkk1*C1*C4*C5*C6*Ekaon*Fm2*
           t3*t4 - 
           6.*amSig*apkk1*C1*C4*C5*C6*Ekaon*Fm2*
           t3*t4 - 
           6.*amLam*apkk1*C1*(C4*C4)*C5*C6*Ekaon*
           Fm2*t3*t4 + 
           6.*amSig*apkk1*C1*(C4*C4)*C5*C6*Ekaon*
           Fm2*t3*t4 + 
           9.*(amLam*amLam)*(C5*C5)*(C6*C6)*Ekaon*
           (t4*t4) + 
           6.*(amLam*amLam)*(C5*C5*C5)*(C6*C6)*Ekaon*
           (t4*t4) + 
           (amLam*amLam)*(C5*C5*C5*C5)*(C6*C6)*Ekaon*(t4*t4) + 
           9.*apkk1*C1*C2*C4*C7*t3*t5 + 
           9.*apkk1*C1*C2*(C4*C4)*C7*t3*t5 - 
           18.*apkk1*C1*C2*C4*C7*Ekaon*Fm1*t3*
           t5 - 9.*apkk1*C2*C5*C6*C7*t4*t5 + 
           3.*apkk1*C2*(C5*C5)*C6*C7*t4*t5 + 
           6.*apkk1*C2*C5*C6*C7*Ekaon*Fm2*t4*
           t5 - 9.*apkk1*C1*C4*C8*C9*t3*t6 - 
           9.*apkk1*C1*(C4*C4)*C8*C9*t3*t6 + 
           18.*apkk1*C1*C4*C8*C9*Ekaon*Fm1*t3*
           t6 + 9.*apkk1*C5*C6*C8*C9*t4*t6 - 
           3.*apkk1*(C5*C5)*C6*C8*C9*t4*t6 - 
           6.*apkk1*C5*C6*C8*C9*Ekaon*Fm2*t4*
           t6 + 
           (amk*amk)*
           (18.*(C1*C1)*(C4*C4)*((1.+C4)*(1.+C4))*Ekaon*
           (t3*t3) + 
           3.*amSig*C1*C4*t3*
           (6.*C1*C4*(1. + (C4*C4))*t3 + 
           C5*C6*
           (-6. + C5 + 3.*Enu*Fm1 - 
           1.*C5*Enu*Fm1 - 1.*Enu*Fm2 + 
           C4*(3. + 2.*C5 - 1.*Enu*Fm2))*t4)
           + C5*C6*t4*
           (2.*amLam*C5*(9. + (C5*C5))*C6*t4 + 
           18.*C5*C6*Ekaon*t4 - 
           12.*(C5*C5)*C6*Ekaon*t4 + 
           2.*(C5*C5*C5)*C6*Ekaon*t4 + 
           9.*C2*C7*t5 - 9.*C8*C9*t6) + 
           3.*C1*C4*t3*
           (-12.*(1. + C4)*C5*C6*Ekaon*t4 + 
           4.*(1. + C4)*(C5*C5)*C6*Ekaon*t4 + 
           amLam*C5*C6*
           (-6. - 1.*C5 - 3.*Enu*Fm1 + 
           C5*Enu*Fm1 + Enu*Fm2 + 
           C4*(-3. + 2.*C5 + Enu*Fm2))*t4 - 
           3.*C2*C7*t5 + 3.*C8*C9*t6))))) - 
           2.*apkk1*
           (-18.*(amSig*amSig)*(C1*C1)*(C4*C4)*(1. + (C4*C4))*
           Ekaon*(t3*t3) + 24.*amLam*(C5*C5)*C6*t1*t4 + 
           72.*amLam*C5*C6*f*t1*t4 - 
           12.*amLam*C5*C6*Enu*Fm2*t1*t4 + 
           12.*amLam*C5*C6*Enu*f*Fm2*t1*t4 + 
           18.*(amk*amk)*amLam*C1*C4*C5*C6*t3*t4 - 
           6.*(amk*amk)*amLam*C1*(C4*C4)*(C5*C5)*C6*t3*t4 - 
           9.*(amk*amk)*amLam*C1*C4*C5*C6*Enu*Fm1*t3*
           t4 + 
           3.*(amk*amk)*amLam*C1*C4*(C5*C5)*C6*Enu*Fm1*t3*
           t4 + 
           3.*(amk*amk)*amLam*C1*C4*C5*C6*Enu*Fm2*t3*
           t4 + 
           3.*(amk*amk)*amLam*C1*(C4*C4)*C5*C6*Enu*Fm2*t3*
           t4 - 18.*(amk*amk)*amLam*(C5*C5)*(C6*C6)*(t4*t4) - 
           2.*(amk*amk)*amLam*(C5*C5*C5*C5)*(C6*C6)*(t4*t4) - 
           18.*(amLam*amLam)*(C5*C5)*(C6*C6)*Ekaon*(t4*t4) - 
           2.*(amLam*amLam)*(C5*C5*C5*C5)*(C6*C6)*Ekaon*(t4*t4) + 
           2.*(am*am)*Ekaon*
           (9.*(C1*C1)*(C4*C4)*(1. + (C4*C4))*(t3*t3) + 
           6.*C1*C4*C5*(-3. + C4*C5)*C6*t3*t4 + 
           (C5*C5)*(9. + (C5*C5))*(C6*C6)*(t4*t4)) + 
           72.*C2*C7*f*t1*t5 + 
           18.*(amk*amk)*C1*C2*C4*C7*t3*t5 - 
           9.*(amk*amk)*C1*C2*C4*C7*Enu*Fm1*t3*t5 - 
           18.*(amk*amk)*C2*C5*C6*C7*t4*t5 - 
           36.*amLam*C2*C5*C6*C7*Ekaon*t4*t5 + 
           9.*amLam*C2*C5*C6*C7*Enu*t4*t5 + 
           3.*amLam*C2*(C5*C5)*C6*C7*Enu*t4*t5 + 
           3.*(amk*amk)*C2*C5*C6*C7*Enu*Fm2*t4*t5 - 
           18.*(C2*C2)*(C7*C7)*Ekaon*(t5*t5) + 
           18.*(C2*C2)*(C7*C7)*Enu*(t5*t5) - 
           72.*C8*C9*f*t1*t6 - 
           18.*(amk*amk)*C1*C4*C8*C9*t3*t6 + 
           9.*(amk*amk)*C1*C4*C8*C9*Enu*Fm1*t3*t6 + 
           18.*(amk*amk)*C5*C6*C8*C9*t4*t6 + 
           36.*amLam*C5*C6*C8*C9*Ekaon*t4*t6 - 
           9.*amLam*C5*C6*C8*C9*Enu*t4*t6 - 
           3.*amLam*(C5*C5)*C6*C8*C9*Enu*t4*t6 - 
           3.*(amk*amk)*C5*C6*C8*C9*Enu*Fm2*t4*t6 + 
           36.*C2*C7*C8*C9*Ekaon*t5*t6 - 
           36.*C2*C7*C8*C9*Enu*t5*t6 - 
           18.*(C8*C8)*(C9*C9)*Ekaon*(t6*t6) + 
           18.*(C8*C8)*(C9*C9)*Enu*(t6*t6) - 
           3.*amSig*C1*C4*t3*
           (24.*f*t1 - 12.*Enu*Fm1*t1 + 
           12.*Enu*f*Fm1*t1 + 
           6.*(amk*amk)*C1*(C4*C4*C4)*t3 - 
           6.*(amk*amk)*C5*C6*t4 - 
           12.*amLam*C5*C6*Ekaon*t4 - 
           3.*(amk*amk)*C5*C6*Enu*Fm1*t4 + 
           (amk*amk)*(C5*C5)*C6*Enu*Fm1*t4 + 
           (amk*amk)*C5*C6*Enu*Fm2*t4 - 
           12.*C2*C7*Ekaon*t5 + 3.*C2*C7*Enu*t5 + 
           12.*C8*C9*Ekaon*t6 - 3.*C8*C9*Enu*t6 + 
           C4*
           (-24.*t1 + 
           4.*amLam*(C5*C5)*C6*Ekaon*t4 + 
           (amk*amk)*
           (6.*C1*t3 + 
           C5*C6*(2.*C5 + Enu*Fm2)*t4) - 
           3.*C2*C7*Enu*t5 + 3.*C8*C9*Enu*t6))
           - 1.*am*
           (18.*(amk*amk)*(C1*C1)*(C4*C4)*(1. + (C4*C4))*
           (t3*t3) - 
           3.*C1*C4*t3*
           (12.*Enu*Fm1*t1 - 
           12.*f*(2. + Enu*Fm1)*t1 + 
           12.*(amk*amk)*C5*C6*t4 + 
           6.*amLam*C5*C6*Ekaon*Enu*Fm1*t4 - 
           6.*amSig*C5*C6*Ekaon*Enu*Fm1*t4 - 
           2.*amLam*(C5*C5)*C6*Ekaon*Enu*Fm1*t4 + 
           2.*amSig*(C5*C5)*C6*Ekaon*Enu*Fm1*t4 - 
           2.*amLam*C5*C6*Ekaon*Enu*Fm2*t4 + 
           2.*amSig*C5*C6*Ekaon*Enu*Fm2*t4 - 
           3.*C2*C7*Enu*t5 + 
           6.*C2*C7*Ekaon*Enu*Fm1*t5 + 
           3.*C8*C9*Enu*t6 - 
           6.*C8*C9*Ekaon*Enu*Fm1*t6 + 
           C4*
           (24.*t1 - 4.*(amk*amk)*(C5*C5)*C6*t4 + 
           Enu*
           (-2.*amLam*C5*C6*Ekaon*Fm2*t4 + 
           2.*amSig*C5*C6*Ekaon*Fm2*t4 + 
           3.*C2*C7*t5 - 3.*C8*C9*t6))) + 
           C5*C6*t4*
           (2.*(amk*amk)*(C5*C5*C5)*C6*t4 - 
           3.*C5*
           (8.*t1 - 6.*(amk*amk)*C6*t4 + 
           C2*C7*Enu*t5 - 1.*C8*C9*Enu*t6) + 
           3.*
           (-4.*f*(6. + Enu*Fm2)*t1 + 
           Enu*
           (-3.*C2*C7*t5 + 3.*C8*C9*t6 + 
           2.*Fm2*
           (2.*t1 + C2*C7*Ekaon*t5 - 
           1.*C8*C9*Ekaon*t6))))))))));
  
  if (sol <= 0.0) {
    LOG( "SKXSec", pERROR )
      << "Matrix element is 0 or negative! This is bad and really should throw an exception!";
  }  

  return sol;
}
//____________________________________________________________________________
