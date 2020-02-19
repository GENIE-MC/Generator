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

#include "Framework/Algorithm/AlgConfigPool.h"
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
  LOG("AhrensDMEL", pNOTICE) << "Form factor masses are " << fMv2 << ", " << fMa2;
  double qmv2 = TMath::Power(1 + Q2/fMv2, 2);
  double qma2 = TMath::Power(1 + Q2/fMa2, 2);

  //-- handle terms changing sign for antineutrinos and isospin rotations
  int nusign  = 1; 
  int nucsign = 1;
  int nupdgc  = init_state.ProbePdg(); 
  int nucpdgc = target.HitNucPdg();
  if( pdg::IsAntiDarkMatter(nupdgc) ) nusign = -1;
  if( pdg::IsNeutron(nucpdgc)       ) nucsign = -1;

  LOG("AhrensDMEL", pNOTICE) << "Calculating for nuclear sign " << nucsign;

  //-- compute up quark form factor terms
  double Geu = fQuV * (1.5 + nucsign*0.5) / qmv2;
  double Gmu = fQuV * ((1.5 + nucsign*0.5) * fMuP + (1.5 - nucsign*0.5) * fMuN) / qmv2;
  double FAu = fQuA * (nucsign > 0 ? fDelu : fDeld) / qma2;
		
  //-- compute down quark form factor terms
  double Ged = fQdV * (1.5 - nucsign*0.5) / qmv2;
  double Gmd = fQdV * ((1.5 - nucsign*0.5) * fMuP + (1.5 + nucsign*0.5) * fMuN) / qmv2;
  double FAd = fQdA * (nucsign > 0 ? fDeld : fDelu) / qma2;

  //-- compute the induced pseudoscalar form factors
  double pole3 = 4.0*M2 / (Q2 + fMpi2);
  double pole0 = 4.0*M2 / (Q2 + fMeta2);
  double FPu = 0.5 * (pole3 * (FAu - FAd) + pole0 * (FAu + FAd));
  double FPd = 0.5 * (- pole3 * (FAu - FAd) + pole0 * (FAu + FAd));
  
  //-- compute strange quark form factor terms
  double Ges = 0.0;
  double Gms = 0.0;
  double FAs = fQsA * fDels / qma2;
  double FPs = 0.0; // fQsA * 2.0 * M2 * fDels / qmp2 / fMpi2;

  //-- compute form factors
  double Ge = Geu + Ged + Ges;
  double Gm = Gmu + Gmd + Gms;
  double FA = FAu + FAd + FAs;
  double FP = FPu + FPd + FPs;
  double tau = 0.25 * Q2/M2;
  double F1 = (Ge + tau * Gm) / (1.0 + tau);
  double F2 = (Gm - Ge) / (1.0 + tau);
  double F12 = TMath::Power(F1,2);
  double F22 = TMath::Power(F2,2);
  double FA2 = TMath::Power(FA,2);

  //-- compute the free nucleon cross section
  double xsec = 0.; 
  double del   = ml2 / M2;
  double AT_F1F1 = 0.;
  double AT_F2F2 = 0.;
  double AT_FAFA = 0.;
  double AT_F1F2 = 0.;
  double AL = 0.;
  double B = 0.;
  double C = 0.;
  if (fVelMode == 0) {
    double QchiV2 = TMath::Power(fQchiV,2); 
    double QchiA2 = TMath::Power(fQchiA,2);
    C = (QchiA2 + QchiV2) * (F12 + F22 * tau + FA2);
    B = 8. * fQchiA * fQchiV * tau * FA * (F1 + F2);
    AL = 16. * QchiA2 * del * TMath::Power(tau*(FA - 2.*FP*tau),2);
    AT_F1F1 = QchiA2*(tau-1.)*(del+tau) + QchiV2*tau*(-del+tau-1);
    AT_F2F2 = -tau*(QchiA2*(tau-1.)*(del+tau) + QchiV2*(del + (tau-1.)*tau));
    AT_FAFA = (1.+tau)*(QchiA2*(del+tau) + QchiV2*(tau-del));
    AT_F1F2 = 2.*tau*(2.*QchiA2*(del+tau) - QchiV2*(del-2.*tau));    
  }
  else if (fVelMode == 2) {
    double QchiS2 = TMath::Power(fQchiS,2);
    C = QchiS2 * (F12 + F22 * tau + FA2);
    AT_F1F1 = -QchiS2 * tau * (del + tau);
    AT_F2F2 = AT_F1F1;
    AT_FAFA = -QchiS2 * (tau + 1.) * (del + tau);
    AT_F1F2 = 2.*AT_F1F1;
  }
  double smu = E/M - tau;
  double MZ2   = TMath::Power(fMedMass,2);
  double lon = TMath::Power(M2 / MZ2 + 0.25/tau,2);
  LOG("AhrensDMEL", pDEBUG)
    << "Using a mediator mass of " << fMedMass;
  //  double fd    = 8*ml2*M2*tau*(2*M2*tau+MZ2) / (MZ2*MZ2*E2);
  double gZp   = fgZp;
  double gZp4  = TMath::Power(gZp,4);
  double prop  = 1. / (Q2 + MZ2);
  double prop2 = TMath::Power(prop,2);
  double xsec0 = gZp4 * M2 * prop2 / (4. * kPi * (E2 - ml2));
  xsec = xsec0 * (AL * lon + AT_F1F1 * F12 + AT_F2F2 * F22 + AT_FAFA * FA2 + AT_F1F2 * F1 * F2 + nusign * B * smu + C * smu * smu); 

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
  // dark matter couplings to mediator
  double QchiL, QchiR;
  this->GetParam( "DarkLeftCharge", QchiL ) ;
  this->GetParam( "DarkRightCharge", QchiR ) ;
  this->GetParam( "DarkScalarCharge", fQchiS ) ;
  fQchiV = 0.5*(QchiL + QchiR);
  fQchiA = 0.5*(- QchiL + QchiR);

  // quark couplings to mediator
  double QuL, QuR, QdL, QdR, QsL, QsR;
  this->GetParam( "UpLeftCharge", QuL ) ;
  this->GetParam( "UpRightCharge", QuR ) ;
  this->GetParam( "DownLeftCharge", QdL ) ;
  this->GetParam( "DownRightCharge", QdR ) ;
  this->GetParam( "StrangeLeftCharge", QsL ) ;
  this->GetParam( "StrangeRightCharge", QsR ) ;
  fQuV = 0.5*(QuL + QuR);
  fQuA = 0.5*(- QuL + QuR);
  fQdV = 0.5*(QdL + QdR);
  fQdA = 0.5*(- QdL + QdR);
  fQsV = 0.5*(QsL + QsR);
  fQsA = 0.5*(- QsL + QsR);

  // axial and vector masses
  double ma, mv, mp, mpi, meta ;
  this->GetParam( "QEL-Ma", ma ) ;
  this->GetParam( "QEL-Mv", mv ) ;
  this->GetParam( "DMEL-Mp", mp ) ;
  this->GetParam( "DMEL-Mpi", mpi ) ;
  this->GetParam( "DMEL-Meta", meta ) ;
  fMa2 = TMath::Power(ma,2);
  fMv2 = TMath::Power(mv,2);
  fMp2 = TMath::Power(mp,2);
  fMpi2 = TMath::Power(mpi,2);
  fMeta2 = TMath::Power(meta,2);

  // anomalous magnetic moments
  this->GetParam( "AnomMagnMoment-P", fMuP ) ;
  this->GetParam( "AnomMagnMoment-N", fMuN ) ;

  // Axial-vector spin charge
  // Since we have a more general axial dependence,
  // we need a more complex treatment than the usual model
  this->GetParam( "AxialVectorSpin-u", fDelu );
  this->GetParam( "AxialVectorSpin-d", fDeld );
  this->GetParam( "AxialVectorSpin-s", fDels );

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
