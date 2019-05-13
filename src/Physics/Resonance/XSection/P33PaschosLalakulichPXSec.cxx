//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: This class is based on code written by the model authors (Olga
         Lalakulich, 17.02.2005). The code was modified to fit into the
         GENIE framework by Costas Andreopoulos.

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/ParticleData/BaryonResonance.h"
#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/Utils/Range1.h"
#include "Physics/Resonance/XSection/P33PaschosLalakulichPXSec.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
P33PaschosLalakulichPXSec::P33PaschosLalakulichPXSec() :
XSecAlgorithmI("genie::P33PaschosLalakulichPXSec")
{

}
//____________________________________________________________________________
P33PaschosLalakulichPXSec::P33PaschosLalakulichPXSec(string config) :
XSecAlgorithmI("genie::P33PaschosLalakulichPXSec", config)
{

}
//____________________________________________________________________________
P33PaschosLalakulichPXSec::~P33PaschosLalakulichPXSec()
{

}
//____________________________________________________________________________
double P33PaschosLalakulichPXSec::XSec(
                  const Interaction * interaction, KinePhaseSpace_t kps) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  //-- Get initial state and kinematic variables
  const InitialState & init_state = interaction -> InitState();
  const Kinematics &   kinematics = interaction -> Kine();
  const Target &       target     = init_state.Tgt();

  double E    = init_state.ProbeE(kRfHitNucRest);
  double E2   = TMath::Power(E,2);
  double Q2   = kinematics.Q2();
  double W    = kinematics.W();
  double MN   = target.HitNucMass();
  double MN2  = TMath::Power(MN,2);
  double Mmu2 = kMuonMass2;
  double Mpi2 = kPionMass2;

  LOG("PaschLal", pDEBUG) << "Input kinematics: W = " << W << ", Q2 = " << Q2;

  //-- Retrieve P33(1232) information
  double Gamma_R0  = utils::res::Width (kP33_1232);
  double MR        = utils::res::Mass  (kP33_1232);
  double MR2      = TMath::Power(MR,2);
  double MR3      = TMath::Power(MR,3);

  //-- Auxiliary params

  const double kPlRes_f3_P1232_V =  1.95/MN;
  const double kPlRes_f4_P1232_V = -1.95/MN;
  const double kPlRes_f5_P1232_V =  0;
  const double kPlRes_f5_P1232_A =  1.2;
  const double kPlRes_f4_P1232_A = -0.3/MN2;
  const double kPlRes_f3_P1232_A =  0;
  const double kPlRes_f6_P1232_A =  kPlRes_f5_P1232_A;

  double MA2    = TMath::Power( fMa, 2 );
  double MV2    = TMath::Power( fMv, 2 );
  double ftmp1a = TMath::Power( 1 + Q2/MA2, 2 );
  double ftmp1v = TMath::Power( 1 + Q2/MA2, 2 );
  double ftmp2a = 1 + Q2/3/MA2;
  double ftmp2v = 1 + Q2/4/MV2;
  double f3A    = kPlRes_f3_P1232_A/ftmp1a/ftmp2a;
  double f4A    = kPlRes_f4_P1232_A/ftmp1a/ftmp2a;
  double f5A    = kPlRes_f5_P1232_A/ftmp1a/ftmp2a;
  double f6A    = kPlRes_f6_P1232_A/ftmp1a/ftmp2a/(Q2+Mpi2);
  double f3V    = kPlRes_f3_P1232_V/ftmp1v/ftmp2v;
  double f4V    = kPlRes_f4_P1232_V/ftmp1v/ftmp2v/W;
  double f5V    = kPlRes_f5_P1232_V/ftmp1v/ftmp2v;
  double f3V4A  = f3V*f4A;
  double f3V5A  = f3V*f5A;
  double f4V4A  = f4V*f4A;
  double f4V5A  = f4V*f5A;
  double f3A2   = TMath::Power( f3A, 2 );
  double f4A2   = TMath::Power( f4A, 2 );
  double f5A2   = TMath::Power( f5A, 2 );
  double f6A2   = TMath::Power( f6A, 2 );
  double f3V2   = TMath::Power( f3V, 2 );
  double f4V2   = TMath::Power( f4V, 2 );
  double f5V2   = TMath::Power( f5V, 2 );

  //-- move these running gamma definitions into GENIE's breit-wigner
  //   functions so that they can be directly used here

  // model of the running Gamma from Paschos [default]
  double Gamma_R=Gamma_R0*pow((this->PPiStar(W,MN)/this->PPiStar(MR,MN)),3);

  // check for other option
  if ( GetConfig().Exists("running-gamma") ) {

     string gamma_model = GetConfig().GetString("running-gamma");

     if ( gamma_model.find("Hagiwara") != string::npos )
     {
         // model of the running Gamma from Hagiwara et. al.
         Gamma_R = Gamma_R0*MR/W*pow((this->PPiStar(W,MN)/this->PPiStar(MR,MN)),1);
     } else
     if ( gamma_model.find("Galster") != string::npos )
     {
         // model of the running Gamma similar to Galster-1972
         double gtmp1 = TMath::Power( this->PPiStar(W,MN) / this->PPiStar(MR,MN), 3);
         double gtmp2 = TMath::Power( this->PPiStar(W,MN) / this->PPiStar(MR,MN), 2);
         Gamma_R = Gamma_R0*gtmp1/(1+gtmp2);
     }
  }
  double Breit_Wigner = TMath::Power(W*W-MR2,2) + MR2 * TMath::Power(Gamma_R,2);

  //-- Include Pauli suppression [if the option was turned on by the user]
  double pauli = 1.;
  if(fTurnOnPauliCorrection) pauli = this->Pauli(Q2,W,MN);

  //-- Kinematic variables
  double nu  = this->Nu(Q2,W,MN);
  double pq  = MN*nu;
  double qk  = -(Q2+Mmu2)/2.;
  double pk  = MN*E;
  double pq2 = TMath::Power(pq,2);
  double pq3 = TMath::Power(pq,3);
  double Q4  = TMath::Power(Q2,2);

  //-- Compute Wi's, i=1-5
  double W1=0, W2=0, W3=0, W4=0, W5=0;

  W1 = 3.*(2*f5A2*MN2*MR2+2*f5A2*MN*MR3+2*f3A*f5A*MN2*MR*pq+2*f5A2*MR2*pq
       +4*f3A*f5A*MN*MR2*pq+4*f4A*f5A*MN2*MR2*pq+2*f3A*f5A*MR3*pq
       +4*f4A*f5A*MN*MR3*pq+2*f3A2*MN2*pq2+2*f3V2*MN2*pq2+2*f3A*f5A*MR*pq2
       +2*f3A*f4A*MN2*MR*pq2+2*f3V*f4V*MN2*MR*pq2+2*f3V*f5V*MN2*MR*pq2
       +2*f3A2*MR2*pq2+2*f3V2*MR2*pq2+4*f4A*f5A*MR2*pq2+4*f3A*f4A*MN*MR2*pq2
       -4*f3V*f4V*MN*MR2*pq2-4*f3V*f5V*MN*MR2*pq2+2*f4A2*MN2*MR2*pq2+2*f4V2*MN2*MR2*pq2
       +4*f4V*f5V*MN2*MR2*pq2+2*f5V2*MN2*MR2*pq2+2*f3A*f4A*MR3*pq2+2*f3V*f4V*MR3*pq2
       +2*f3V*f5V*MR3*pq2+2*f4A2*MN*MR3*pq2-2*f4V2*MN*MR3*pq2-4*f4V*f5V*MN*MR3*pq2
       -2*f5V2*MN*MR3*pq2+2*f3A2*pq3+2*f3V2*pq3+2*f3A*f4A*MR*pq3+2*f3V*f4V*MR*pq3
       +2*f3V*f5V*MR*pq3+2*f4A2*MR2*pq3+2*f4V2*MR2*pq3+4*f4V*f5V*MR2*pq3+2*f5V2*MR2*pq3
       -2*f3A*f5A*MN2*MR*Q2-4*f3A*f5A*MN*MR2*Q2+2*f3A2*MN2*MR2*Q2
       +2*f3V2*MN2*MR2*Q2-4*f4A*f5A*MN2*MR2*Q2-2*f3A2*MN*MR3*Q2+2*f3V2*MN*MR3*Q2
       -4*f4A*f5A*MN*MR3*Q2-4*f3A2*MN2*pq*Q2-4*f3V2*MN2*pq*Q2-2*f3A*f5A*MR*pq*Q2
       -4*f3A*f4A*MN2*MR*pq*Q2-4*f3V*f4V*MN2*MR*pq*Q2-2*f3V*f5V*MN2*MR*pq*Q2
       -4*f4A*f5A*MR2*pq*Q2-8*f3A*f4A*MN*MR2*pq*Q2+8*f3V*f4V*MN*MR2*pq*Q2
       +4*f3V*f5V*MN*MR2*pq*Q2-4*f4A2*MN2*MR2*pq*Q2-4*f4V2*MN2*MR2*pq*Q2
       -4*f4V*f5V*MN2*MR2*pq*Q2-2*f3A*f4A*MR3*pq*Q2-2*f3V*f4V*MR3*pq*Q2
       -4*f4A2*MN*MR3*pq*Q2+4*f4V2*MN*MR3*pq*Q2+4*f4V*f5V*MN*MR3*pq*Q2
       -4*f3A2*pq2*Q2-4*f3V2*pq2*Q2-4*f3A*f4A*MR*pq2*Q2-4*f3V*f4V*MR*pq2*Q2
       -2*f3V*f5V*MR*pq2*Q2-4*f4A2*MR2*pq2*Q2-4*f4V2*MR2*pq2*Q2-4*f4V*f5V*MR2*pq2*Q2
       +2*f3A2*MN2*Q4+2*f3V2*MN2*Q4+2*f3A*f4A*MN2*MR*Q4+2*f3V*f4V*MN2*MR*Q4
       +4*f3A*f4A*MN*MR2*Q4-4*f3V*f4V*MN*MR2*Q4+2*f4A2*MN2*MR2*Q4+2*f4V2*MN2*MR2*Q4
       +2*f4A2*MN*MR3*Q4-2*f4V2*MN*MR3*Q4+2*f3A2*pq*Q4+2*f3V2*pq*Q4+2*f3A*f4A*MR*pq*Q4
       +2*f3V*f4V*MR*pq*Q4+2*f4A2*MR2*pq*Q4+2*f4V2*MR2*pq*Q4)/(3.*MR2);

  W2 = 3.*(2*(f5A2*MN2
       +f5A2*MN*MR+f5A2*pq+f3A2*MN2*Q2+f3V2*MN2*Q2+f3A*f5A*MR*Q2+f3A*f4A*MN2*MR*Q2
       +f3V*f4V*MN2*MR*Q2+f3V*f5V*MN2*MR*Q2+f3A2*MR2*Q2+f3V2*MR2*Q2
       +2*f3A*f4A*MN*MR2*Q2-2*f3V*f4V*MN*MR2*Q2-2*f3V*f5V*MN*MR2*Q2
       +f4A2*MN2*MR2*Q2+f4V2*MN2*MR2*Q2+2*f4V*f5V*MN2*MR2*Q2+f5V2*MN2*MR2*Q2+f3A*f4A*MR3*Q2
       +f3V*f4V*MR3*Q2+f3V*f5V*MR3*Q2+f4A2*MN*MR3*Q2-f4V2*MN*MR3*Q2
       -2*f4V*f5V*MN*MR3*Q2-f5V2*MN*MR3*Q2+f3A2*pq*Q2+f3V2*pq*Q2+f3A*f4A*MR*pq*Q2
       +f3V*f4V*MR*pq*Q2+f3V*f5V*MR*pq*Q2+f4A2*MR2*pq*Q2+f4V2*MR2*pq*Q2
       +2*f4V*f5V*MR2*pq*Q2+f5V2*MR2*pq*Q2+f5V2*MN2*Q4+f3V*f5V*MR*Q4
       -f5V2*MN*MR*Q4+f5V2*pq*Q4))/(3.*MR2);

  W3 = 3.*((f3V4A*(Q2-pq)-f3V5A)*(2*MR2+2*MN*MR+Q2-pq)*4./3./MR
       -(Q2-pq)*(f4V4A*(Q2-pq)-f4V5A)*4./3.);


  W4 = 3.*(2*(f5A2*MN2+f5A2*MN*MR+f3A*f5A*MN2*MR
       +2*f3A*f5A*MN*MR2-f3A2*MN2*MR2-f3V2*MN2*MR2+2*f4A*f5A*MN2*MR2-2*f5A*f6A*MN2*MR2
       +f3A2*MN*MR3-f3V2*MN*MR3+2*f4A*f5A*MN*MR3-2*f5A*f6A*MN*MR3+f5A2*pq
       +2*f3A2*MN2*pq+2*f3V2*MN2*pq+2*f5A*f6A*MN2*pq+2*f3A*f5A*MR*pq
       +2*f5A*f6A*MN*MR*pq+2*f3A*f4A*MN2*MR*pq+2*f3V*f4V*MN2*MR*pq
       +f3V*f5V*MN2*MR*pq-f3A*f6A*MN2*MR*pq+2*f4A*f5A*MR2*pq
       -2*f5A*f6A*MR2*pq+4*f3A*f4A*MN*MR2*pq-4*f3V*f4V*MN*MR2*pq
       -2*f3V*f5V*MN*MR2*pq-2*f3A*f6A*MN*MR2*pq+2*f4A2*MN2*MR2*pq
       +2*f4V2*MN2*MR2*pq+2*f4V*f5V*MN2*MR2*pq-2*f4A*f6A*MN2*MR2*pq+f3A*f4A*MR3*pq
       +f3V*f4V*MR3*pq-f3A*f6A*MR3*pq+2*f4A2*MN*MR3*pq-2*f4V2*MN*MR3*pq
       -2*f4V*f5V*MN*MR3*pq-2*f4A*f6A*MN*MR3*pq+2*f3A2*pq2+2*f3V2*pq2
       +2*f5A*f6A*pq2+f5V2*MN2*pq2+f6A2*MN2*pq2+2*f3A*f4A*MR*pq2+2*f3V*f4V*MR*pq2
       +2*f3V*f5V*MR*pq2-f5V2*MN*MR*pq2+f6A2*MN*MR*pq2+2*f4A2*MR2*pq2+2*f4V2*MR2*pq2
       +2*f4V*f5V*MR2*pq2-2*f4A*f6A*MR2*pq2+f5V2*pq3+f6A2*pq3-f3A2*MN2*Q2-f3V2*MN2*Q2
       -2*f5A*f6A*MN2*Q2-2*f5A*f6A*MN*MR*Q2-f3A*f4A*MN2*MR*Q2
       -f3V*f4V*MN2*MR*Q2-2*f3A*f4A*MN*MR2*Q2+2*f3V*f4V*MN*MR2*Q2
       -f4A2*MN2*MR2*Q2-f4V2*MN2*MR2*Q2+f6A2*MN2*MR2*Q2-f4A2*MN*MR3*Q2+f4V2*MN*MR3*Q2
       +f6A2*MN*MR3*Q2-f3A2*pq*Q2-f3V2*pq*Q2-2*f5A*f6A*pq*Q2-2*f6A2*MN2*pq*Q2
       -f3A*f4A*MR*pq*Q2-f3V*f4V*MR*pq*Q2-f3A*f6A*MR*pq*Q2
       -2*f6A2*MN*MR*pq*Q2-f4A2*MR2*pq*Q2-f4V2*MR2*pq*Q2+f6A2*MR2*pq*Q2
       -2*f6A2*pq2*Q2+f6A2*MN2*Q4+f6A2*MN*MR*Q4+f6A2*pq*Q4))/(3.*MR2);

  W5 = 3.*(2*f5A2*MN2
       +2*f5A2*MN*MR+f3A*f5A*MN2*MR+2*f3A*f5A*MN*MR2+2*f4A*f5A*MN2*MR2
       +f3A*f5A*MR3+2*f4A*f5A*MN*MR3+2*f5A2*pq+2*f3A2*MN2*pq+2*f3V2*MN2*pq
       +2*f5A*f6A*MN2*pq+2*f3A*f5A*MR*pq+2*f5A*f6A*MN*MR*pq
       +2*f3A*f4A*MN2*MR*pq+2*f3V*f4V*MN2*MR*pq+2*f3V*f5V*MN2*MR*pq
       +2*f3A2*MR2*pq+2*f3V2*MR2*pq+2*f4A*f5A*MR2*pq+4*f3A*f4A*MN*MR2*pq
       -4*f3V*f4V*MN*MR2*pq-4*f3V*f5V*MN*MR2*pq+2*f4A2*MN2*MR2*pq
       +2*f4V2*MN2*MR2*pq+4*f4V*f5V*MN2*MR2*pq+2*f5V2*MN2*MR2*pq+2*f3A*f4A*MR3*pq
       +2*f3V*f4V*MR3*pq+2*f3V*f5V*MR3*pq+2*f4A2*MN*MR3*pq-2*f4V2*MN*MR3*pq
       -4*f4V*f5V*MN*MR3*pq-2*f5V2*MN*MR3*pq+2*f3A2*pq2+2*f3V2*pq2+2*f5A*f6A*pq2
       +2*f3A*f4A*MR*pq2+2*f3V*f4V*MR*pq2+2*f3V*f5V*MR*pq2+2*f4A2*MR2*pq2
       +2*f4V2*MR2*pq2+4*f4V*f5V*MR2*pq2+2*f5V2*MR2*pq2-2*f5A*f6A*MN2*Q2+f3A*f5A*MR*Q2
       -2*f5A*f6A*MN*MR*Q2-f3A*f6A*MN2*MR*Q2-2*f3A*f6A*MN*MR2*Q2
       -2*f4A*f6A*MN2*MR2*Q2-f3A*f6A*MR3*Q2-2*f4A*f6A*MN*MR3*Q2
       -2*f5A*f6A*pq*Q2+2*f5V2*MN2*pq*Q2+2*f3V*f5V*MR*pq*Q2
       -2*f5V2*MN*MR*pq*Q2-2*f4A*f6A*MR2*pq*Q2+2*f5V2*pq2*Q2-f3A*f6A*MR*Q4)/(3.*MR2);

  double s1 =    W1 * (Q2+Mmu2)
               + W2 * (2*pk*pk-2*pq*pk+MN*qk)
               - W3 * (pq*qk+Q2*pk)
               + W4 * Mmu2*(Q2+Mmu2)/2.
               - W5 * 2*Mmu2*pk;

  double xsec = kGF2/4./kPi*fCos28c/MN2/E2*W*MR*Gamma_R/kPi/Breit_Wigner*pauli*s1;

  //-- The algorithm computes d^2xsec/dWdQ2
  //   Check whether variable tranformation is needed
  if(kps!=kPSWQ2fE) {
    double J = utils::kinematics::Jacobian(interaction,kPSWQ2fE,kps);
    xsec *= J;
  }

  //-- If requested return the free nucleon xsec even for input nuclear tgt
  if( interaction->TestBit(kIAssumeFreeNucleon) ) return xsec;

  //-- number of scattering centers in the target
  bool isp = pdg::IsProton(target.HitNucPdg());
  int NNucl = (isp) ? target.Z() : target.N();

  xsec*=NNucl; // nuclear xsec (no nuclear suppression factor)

  return xsec;
}
//____________________________________________________________________________
double P33PaschosLalakulichPXSec::Integral(
                                        const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool P33PaschosLalakulichPXSec::ValidProcess(
                                        const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;
  return true;
}
//____________________________________________________________________________
void P33PaschosLalakulichPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void P33PaschosLalakulichPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void P33PaschosLalakulichPXSec::LoadConfig(void)
{
  GetParam( "RES-Ma", fMa ) ;
  GetParam( "RES-Mv", fMv ) ;

  double thc ;
  GetParam( "CabibboAngle", thc ) ;
  fCos28c = TMath::Power( TMath::Cos(thc), 2 );

  GetParamDef( "TurnOnPauliSuppr", fTurnOnPauliCorrection, false ) ;

  //-- load the differential cross section integrator
  fXSecIntegrator =
      dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);
}
//____________________________________________________________________________
double P33PaschosLalakulichPXSec::Pauli(double Q2, double W, double MN) const
{
// Pauli suppression for deuterium with Fermi momentum 0.160 GeV

  // ---remove--- this value from here & ask GENIE for D2 Fermi momentum
  double qF=0.160; /* 0.160 deuterium */

  double Paulii = 0;

  double p_pi_star = this->PPiStar(W,MN);
  double nu_star   = this->NuStar(Q2,W,MN);

  double p_pi_star_2 = TMath::Power(p_pi_star,   2);
  double p_pi_star_4 = TMath::Power(p_pi_star_2, 2);
  double nu_star_2   = TMath::Power(nu_star,     2);

  double q3  = TMath::Sqrt(Q2+nu_star_2);
  double q6  = TMath::Power(q3,2);
  double q12 = TMath::Power(q6,2);

  double qF2 = TMath::Power(qF,2);
  double qF3 = TMath::Power(qF,3);

  if ( q3+p_pi_star < 2*qF )
  {
     Paulii = ( (3*q6 + p_pi_star_2)/2/qF
               -(5*q12+p_pi_star_4 + 10*q6*p_pi_star_2)/40/qF3
              )/2/q3;
  }
  if ( (q3+p_pi_star > 2*qF) && (q3-p_pi_star < 2*qF) )
  {
     double tmp1 = TMath::Power( q3+p_pi_star, 2 );
     double tmp2 = TMath::Power( q3-p_pi_star, 3 );
     double tmp3 = TMath::Power( q3-p_pi_star, 5 );

     Paulii = (tmp1-4.0*qF2/5.0 - tmp2/2/qF + tmp3/40/qF3)/4/p_pi_star/q3;
  }
  if ( q3-p_pi_star > 2*qF )
  {
     Paulii = 1;
  }

  return Paulii;
}
//____________________________________________________________________________
double P33PaschosLalakulichPXSec::Nu(double Q2, double W, double MN) const
{
  return (TMath::Power(W,2) - TMath::Power(MN,2) + Q2)/2/MN;
}
//____________________________________________________________________________
double P33PaschosLalakulichPXSec::PPiStar(double W, double MN) const
{
  double W2 = TMath::Power(W,2);
  double a  = TMath::Power(MN+kPionMass,2);
  double b  = TMath::Power(MN-kPionMass,2);

  return TMath::Sqrt( (W2-a)*(W2-b) )/2/W;
}
//____________________________________________________________________________
double P33PaschosLalakulichPXSec::NuStar(double Q2, double W, double MN) const
{
  return (TMath::Power(W,2) - TMath::Power(MN,2) - Q2)/2/W;
}
//____________________________________________________________________________


