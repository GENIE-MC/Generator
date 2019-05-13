//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Joshua Berger <jberger \at physics.wisc.edu>
         University of Wisconsin-Madison

         Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 
*/
//____________________________________________________________________________

#include <TMath.h>

//#include "Framework/Algorithm/AlgConfigPool.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/Conventions/KineVar.h"
#include "Physics/BoostedDarkMatter/XSection/AhrensDMELPXSec.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Numerical/MathUtils.h"
#include "Physics/NuclearState/NuclearUtils.h"
#include "Framework/Utils/RunOpt.h"

using namespace genie;
using namespace genie::utils;
using namespace genie::constants;

//____________________________________________________________________________
AhrensDMELPXSec::AhrensDMELPXSec() :
XSecAlgorithmI("genie::AhrensDMELPXSec")
{

}
//____________________________________________________________________________
AhrensDMELPXSec::AhrensDMELPXSec(string config) :
XSecAlgorithmI("genie::AhrensDMELPXSec", config)
{

}
//____________________________________________________________________________
AhrensDMELPXSec::~AhrensDMELPXSec()
{

}
//____________________________________________________________________________
double AhrensDMELPXSec::XSec(
    const Interaction * interaction, KinePhaseSpace_t kps) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  const InitialState & init_state = interaction -> InitState();
  const Kinematics &   kinematics = interaction -> Kine();
  const Target &       target     = init_state.Tgt();

  LOG("AhrensDMEL", pDEBUG) << "Using v^" << fVelMode << " dependence";
  
  double E    = init_state.ProbeE(kRfHitNucRest);
  double ml   = init_state.GetProbeP4(kRfHitNucRest)->M();
  double Q2   = kinematics.Q2();
  double M    = target.HitNucMass();
  double M2   = TMath::Power(M, 2.);
  double E2   = TMath::Power(E, 2.);
  double ml2  = TMath::Power(ml,2.);
  double qma2 = TMath::Power(1 + Q2/fMa2, 2);

  //-- handle terms changing sign for antineutrinos and isospin rotations
  // int nusign  = 1; // comment-out unused variable to eliminate warnings
  int nucsign = 1;
  // int nupdgc  = init_state.ProbePdg(); // // comment-out unused variable to eliminate warnings
  int nucpdgc = target.HitNucPdg();
  if( pdg::IsNeutron(nucpdgc) ) nucsign = -1;

  //-- compute axial form factor terms
  // For dark matter scattering, the factor of 1/2 should be removed
  // This factor should go to ~ 1 in the limit q^2 -> 0
  // The 1/2 is an isospin piece that is irrelevant for DM
  double Ga1 = -fFa0 * (1 + (nucsign) * fEta) / qma2;  

  //-- compute form factors
  double Ga = (nucsign) * Ga1;
  double Ga2 = TMath::Power(Ga,2);

  //-- compute the free nucleon cross section
  double xsec = 0.;
  double tau   = 0.25 * Q2/M2;
  double fa    = 1-M*tau/E;
  double fa2   = TMath::Power(fa,2);
  double fb    = tau*(tau+1)*M2/E2;
  double fc    = (tau+2)*ml2/E2;
  LOG("AhrensDMEL", pDEBUG)
    << "Using a mediator mass of " << fMedMass;
  double MZ2   = TMath::Power(fMedMass,2);
  double fd    = 8*ml2*M2*tau*(2*M2*tau+MZ2) / (MZ2*MZ2*E2);
  // double gZp   = RunOpt::Instance()->ZpCoupling();
  double gZp   = fgZp;
  double gZp4  = TMath::Power(gZp,4);
  double prop  = 1. / (Q2 + MZ2);
  double prop2 = TMath::Power(prop,2);
  double B     = 0.;
  double xsec0 = 0.;
  switch (fVelMode) {
  // v^0 fermion cross-section
  case 0:
    B     = (Ga2) * (fa2+fb+fc+fd);
    xsec0 = 0.25*prop2*gZp4*E2/kPi/(E2-ml2);
    xsec  = xsec0 * (B);
    break;
  // v^2 scalar cross-section
  case 2:
    fc    = fc - ml2/E2;
    B     = (Ga2) * (fa2-fb-fc);
    xsec0 = 0.25*prop2*gZp4*E2/kPi/(E2-ml2);
    xsec  = xsec0 * (B);
    break;
  }

  LOG("AhrensDMEL", pDEBUG)
    << "dXSec[vN,El]/dQ2 [FreeN](Ev = "<< E<< ", Q2 = "<< Q2 << ") = "<< xsec;

  //-- The algorithm computes dxsec/dQ2
  //   Check whether variable tranformation is needed
  if(kps!=kPSQ2fE) {
    double J = utils::kinematics::Jacobian(interaction,kPSQ2fE,kps);
    xsec *= J;
  }

  //-- if requested return the free nucleon xsec even for input nuclear tgt
  if( interaction->TestBit(kIAssumeFreeNucleon) ) return xsec;

  //-- compute nuclear suppression factor
  //   (R(Q2) is adapted from NeuGEN - see comments therein)
  double R = nuclear::NuclQELXSecSuppression("Default", 0.5, interaction);

  //-- number of scattering centers in the target
  int NNucl = (pdg::IsProton(nucpdgc)) ? target.Z() : target.N();

  LOG("AhrensDMEL", pDEBUG)
       << "Nuclear suppression factor R(Q2) = " << R << ", NNucl = " << NNucl;

  //-- compute nuclear cross section
  xsec *= (R*NNucl); 

  return xsec;
}
//____________________________________________________________________________
double AhrensDMELPXSec::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool AhrensDMELPXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;
  return true;
}
//____________________________________________________________________________
void AhrensDMELPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void AhrensDMELPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void AhrensDMELPXSec::LoadConfig(void)
{
  // alpha and gamma
  double thw ;
  this->GetParam( "WeinbergAngle", thw ) ;
  double sin2thw = TMath::Power(TMath::Sin(thw), 2);
  fkAlpha = 1.-2.*sin2thw;
  fkGamma = -0.66666667*sin2thw;

  // eta and Fa(q2=0)
  this->GetParam( "EL-Axial-Eta", fEta ) ;
  this->GetParam( "QEL-FA0", fFa0 ) ;

  // axial and vector masses
  double ma, mv ;
  this->GetParam( "QEL-Ma", ma ) ;
  this->GetParam( "QEL-Mv", mv ) ;
  fMa2 = TMath::Power(ma,2);
  fMv2 = TMath::Power(mv,2);

  // anomalous magnetic moments
  this->GetParam( "AnomMagnMoment-P", fMuP ) ;
  this->GetParam( "AnomMagnMoment-N", fMuN ) ;

  // velocity dependence of interaction
  this->GetParamDef("velocity-mode", fVelMode, 0 );

  // mediator coupling
  this->GetParam("ZpCoupling", fgZp ) ;

  // mediator mass
  fMedMass = PDGLibrary::Instance()->Find(kPdgMediator)->Mass();

  // load XSec Integrator
  fXSecIntegrator =
      dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);
}
//____________________________________________________________________________
