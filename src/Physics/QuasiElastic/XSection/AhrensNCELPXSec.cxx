//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 19, 2009 - CA
   NuNucElasticPXSec -> AhrensNCELPXSec

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/Conventions/KineVar.h"
#include "Physics/QuasiElastic/XSection/AhrensNCELPXSec.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Numerical/MathUtils.h"
#include "Physics/NuclearState/NuclearUtils.h"

using namespace genie;
using namespace genie::utils;
using namespace genie::constants;

//____________________________________________________________________________
AhrensNCELPXSec::AhrensNCELPXSec() :
XSecAlgorithmI("genie::AhrensNCELPXSec")
{

}
//____________________________________________________________________________
AhrensNCELPXSec::AhrensNCELPXSec(string config) :
XSecAlgorithmI("genie::AhrensNCELPXSec", config)
{

}
//____________________________________________________________________________
AhrensNCELPXSec::~AhrensNCELPXSec()
{

}
//____________________________________________________________________________
double AhrensNCELPXSec::XSec(
    const Interaction * interaction, KinePhaseSpace_t kps) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  const InitialState & init_state = interaction -> InitState();
  const Kinematics &   kinematics = interaction -> Kine();
  const Target &       target     = init_state.Tgt();

  double E    = init_state.ProbeE(kRfHitNucRest);
  double Q2   = kinematics.Q2();
  double M    = target.HitNucMass();
  double M2   = TMath::Power(M, 2.);
  double E2   = TMath::Power(E, 2.);
  double qmv2 = TMath::Power(1 + Q2/fMv2, 2);
  double qma2 = TMath::Power(1 + Q2/fMa2, 2);

  //-- handle terms changing sign for antineutrinos and isospin rotations
  int nusign  = 1;
  int nucsign = 1;
  int nupdgc  = init_state.ProbePdg();
  int nucpdgc = target.HitNucPdg();
  if( pdg::IsAntiNeutrino(nupdgc) ) nusign  = -1;
  if( pdg::IsNeutron(nucpdgc)     ) nucsign = -1;

  //-- compute isoscalar form factor terms
  double Ge0 = 1.5 * fkGamma / qmv2;
  double Gm0 = 1.5 * fkGamma * (fMuP+fMuN) / qmv2;

  //-- compute isovector form factor terms
  double Ge1 = 0.5 * fkAlpha / qmv2;
  double Gm1 = 0.5 * fkAlpha * (fMuP-fMuN) / qmv2;
  double Ga1 = -0.5 * fFa0 * (1 + (nucsign) * fEta) / qma2;  

  //-- compute form factors
  double Ge  = Ge0 + (nucsign) * Ge1;
  double Gm  = Gm0 + (nucsign) * Gm1;
  double Ga  = (nucsign) * Ga1;
  double Ge2 = TMath::Power(Ge,2);
  double Gm2 = TMath::Power(Gm,2);
  double Ga2 = TMath::Power(Ga,2);

  //-- compute the free nucleon cross section
  double tau   = 0.25 * Q2/M2;
  double fa    = 1-M*tau/E;
  double fa2   = TMath::Power(fa,2);
  double fb    = tau*(tau+1)*M2/E2; 
  double A     = (Ge2/(1+tau))           * (fa2-fb);
  double B     = (Ga2 + tau*Gm2/(1+tau)) * (fa2+fb);
  double C     = 4*tau*(M/E)*Gm*Ga       * fa;
  double xsec0 = 0.5*kGF2/kPi;
  double xsec  = xsec0 * (A + B + (nusign)*C);

  LOG("AhrensNCEL", pDEBUG)
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

  LOG("AhrensNCEL", pDEBUG)
       << "Nuclear suppression factor R(Q2) = " << R << ", NNucl = " << NNucl;

  //-- compute nuclear cross section
  xsec *= (R*NNucl); 

  return xsec;
}
//____________________________________________________________________________
double AhrensNCELPXSec::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool AhrensNCELPXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;
  return true;
}
//____________________________________________________________________________
void AhrensNCELPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void AhrensNCELPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void AhrensNCELPXSec::LoadConfig(void)
{
  // alpha and gamma
  double thw ;
  GetParam( "WeinbergAngle", thw ) ;
  double sin2thw = TMath::Power(TMath::Sin(thw), 2);
  fkAlpha = 1.-2.*sin2thw;
  fkGamma = -0.66666667*sin2thw;

  // eta and Fa(q2=0)
  GetParam( "EL-Axial-Eta", fEta ) ;
  GetParam( "QEL-FA0", fFa0 ) ;

  // axial and vector masses
  double ma, mv ;
  GetParam( "QEL-Ma", ma ) ;
  GetParam( "QEL-Mv", mv ) ;
  fMa2 = TMath::Power(ma,2);
  fMv2 = TMath::Power(mv,2);

  // anomalous magnetic moments
  GetParam( "AnomMagnMoment-P", fMuP ) ;
  GetParam( "AnomMagnMoment-N", fMuN ) ;

  // load XSec Integrator
  fXSecIntegrator =
      dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);
}
//____________________________________________________________________________
