//____________________________________________________________________________
/*
 Copyright (c) 2003-2023, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Júlia Tena Vidal <jtenavidal \at tauex.tau.ac.il>
 Tel Aviv University

 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Physics/DeepInelastic/XSection/BYStrucFunc2021.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGUtils.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
BYStrucFunc2021::BYStrucFunc2021() :
QPMDISStrucFuncBase("genie::BYStrucFunc2021")
{
  this->Init();
}
//____________________________________________________________________________
BYStrucFunc2021::BYStrucFunc2021(string config):
QPMDISStrucFuncBase("genie::BYStrucFunc2021", config)
{
  this->Init();
}
//____________________________________________________________________________
BYStrucFunc2021::~BYStrucFunc2021()
{

}
//____________________________________________________________________________
void BYStrucFunc2021::Configure(const Registry & config)
{
// Overload Algorithm::Configure() to read the config. registry and set
// private data members.
// QPMDISStrucFuncBase::Configure() creates the owned PDF object that gets
// configured with the specified PDFModelI
// For the ReadBYParams() method see below

  QPMDISStrucFuncBase::Configure(config);
  this->ReadBYParams();
}
//____________________________________________________________________________
void BYStrucFunc2021::Configure(string param_set)
{
  QPMDISStrucFuncBase::Configure(param_set);
  this->ReadBYParams();
}
//____________________________________________________________________________
void BYStrucFunc2021::ReadBYParams(void)
{
// Get the Bodek-Yang model parameters A,B,Csea,Cv1,Cv2 from the config.
// registry and set some private data members so as not to accessing the
// registry at every calculation.
//
  GetParam( "BY-A", fA ) ;
  GetParam( "BY-B", fB ) ;
  GetParam( "BY-CsU", fCsU ) ;
  GetParam( "BY-CsD", fCsD ) ;
  GetParam( "BY-CvLW", fCvLW ) ;
  GetParam( "BY-Cv1U", fCv1U ) ;
  GetParam( "BY-Cv2U", fCv2U ) ;
  GetParam( "BY-Cv1D", fCv1D ) ;
  GetParam( "BY-Cv2D", fCv2D ) ;
  GetParam( "BY-CsS" , fCsS  ) ;
  GetParam( "BY-PsA" , fPsA  ) ;
  GetParam( "BY-PvA" , fPvA  ) ;
  GetParam( "BY-CsA" , fCsA  ) ;
  GetParam( "BY-CaLW-nubar" , fCaLW_nu  ) ;
  GetParam( "BY-CaLW-nubar" , fCaLW_nubar  ) ;
  GetParam( "BY-H0" , fH0  ) ;
  GetParam( "BY-H1" , fH1  ) ;
  GetParam( "BY-H2" , fH2  ) ;
  GetParam( "BY-H3" , fH3  ) ;
}
//____________________________________________________________________________
void BYStrucFunc2021::Init(void)
{
  fA    = 0;
  fB    = 0;
  fCsU  = 0;
  fCsD  = 0;
  fCvLW = 0;
  fCv1U = 0;
  fCv2U = 0;
  fCv1D = 0;
  fCv2D = 0;
  fCsS  = 0;
  fPsA  = 0;
  fPvA  = 0;
  fCsA  = 0;
  fCaLW_nu = 0;
  fCaLW_nubar = 0;
  fH0   = 0;
  fH1   = 0;
  fH2   = 0;
  fH3   = 0;
}
//____________________________________________________________________________
double BYStrucFunc2021::ScalingVar(const Interaction * interaction, double Mf ) const
{
// Overrides QPMDISStrucFuncBase::ScalingVar() to compute the BY scaling var

  const Kinematics & kine  = interaction->Kine();
  double x  = kine.x();
  double myQ2 = this->Q2(interaction);
  //myQ2 = TMath::Max(Q2,fQ2min);
  LOG("BodekYang", pDEBUG) << "Q2 at scaling var calculation = " << myQ2;

  double a  = TMath::Power( 2*kProtonMass*x, 2 ) / myQ2;
  double Mf2 = TMath::Power( Mf, 2 ) ; 
  double xw =  2*x*(myQ2+Mf2+fB) / (myQ2*(1.+TMath::Sqrt(1+a)) +  2*fA*x);
  return xw;
}
//____________________________________________________________________________
void BYStrucFunc2021::KFactors(const Interaction * interaction,
   double & kuv, double & kdv, double & kus, double & kds, double & kss ) const
{
// Overrides QPMDISStrucFuncBase::KFactors() to compute the BY K factors for
// u(valence), d(valence), u(sea), d(sea), s(sea);

  double myQ2  = this->Q2(interaction);
  double GD  = 1. / TMath::Power(1.+myQ2/0.71, 2); // p elastic form factor
  double GD2 = TMath::Power(GD,2);

  double q0 = this->q0(interaction);
  double q02 = TMath::Power(q0,2);
  double KLW = 1 ;
  // The K factor blows up at q0 = 0. A. Bodek recomends to use it until W = 1.1 GeV
  double W = interaction->Kine().W();      
  if ( W >= 1.1 ) KLW = ( q02 + fCvLW ) / q02;
  
  kuv = KLW * (1.-GD2)*(myQ2+fCv2U)/(myQ2+fCv1U); // K - u(valence)
  kdv = KLW * (1.-GD2)*(myQ2+fCv2D)/(myQ2+fCv1D); // K - d(valence)
  kus = myQ2/(myQ2+fCsU);    // K - u(sea)
  kds = myQ2/(myQ2+fCsD);    // K - d(sea)
  kss = myQ2/(myQ2+fCsS);    // K - s(sea)	
}
//____________________________________________________________________________
void BYStrucFunc2021::KAxialFactors (const Interaction * interaction, double & ksea, double & kvalance ) const {
  // https://arxiv.org/pdf/2108.09240 Sec 11.2 
  double myQ2  = this->Q2(interaction);
  ksea = ( myQ2 + fPsA*fCsA ) / ( myQ2 + fCsA ) ;
  kvalance = ( myQ2 + fPvA * 0.18 ) / ( myQ2 + 0.18 ) ; 

  // Compute Low-q0 modificaion factor for neutrinos and anti-neutrinos
  double q0 = this->q0(interaction);
  double q02 = TMath::Power(q0,2);
  const ProcessInfo &  proc_info  = interaction->ProcInfo();
  const InitialState & init_state = interaction->InitState();
  int  probe_pdgc  = init_state.ProbePdg();
  bool is_nu       = pdg::IsNeutrino     ( probe_pdgc  );
  bool is_nubar    = pdg::IsAntiNeutrino ( probe_pdgc  );
  double KLW = 1; 
  if( is_nu ) {
    KLW = ( q02 + fCaLW_nu ) / q02 ;
  } else if( is_nubar ) {
    KLW = ( q02 + fCaLW_nubar ) / q02 ;
  }
  // double check it only affects valance factors:
  kvalance *= KLW ; 
}
//____________________________________________________________________________
double BYStrucFunc2021::H(const Interaction * interaction) const {
  // Overrides QPMDISStrucFuncBase::H() function to compute the correction of the BY 2021 update
  // The correction is given by Eq. 34 of https://arxiv.org/pdf/2108.09240
  const Kinematics & kinematics = interaction->Kine();
  double bjx = kinematics.x();
  return fH0 + fH1 * bjx + fH2 * pow(bjx,2) + fH3 * pow(bjx,3);
}
//____________________________________________________________________________
double BYStrucFunc2021::R(const Interaction * interaction) const {

  // Evaluate correction for Q2 < 0.3 GeV2/c4 according to Sec. 7 of https://arxiv.org/pdf/2108.09240
  double Q2 = this->Q2(interaction);
  double R1998 = QPMDISStrucFuncBase::R(interaction);

  if( Q2 < 0.3 ) {
    Q2 = 0.3 ;
    double Q4 = pow(Q2,2);
    R1998 *= 3.633 * Q2 / ( Q4 + 1 ) ;
  }
  
  return R1998 ; 
}
//____________________________________________________________________________
double BYStrucFunc2021::NuclMod(const Interaction * interaction) const {
// Nuclear modification to Fi
// The scaling variable can be overwritten to include corrections

  if( interaction->TestBit(kIAssumeFreeNucleon)   ) return 1.0;
  if( interaction->TestBit(kINoNuclearCorrection) ) return 1.0;

  double f = 1.;
  if(fIncludeNuclMod) {
    const Target & tgt  = interaction->InitState().Tgt();
     const Kinematics & kinematics = interaction->Kine();
     double x = kinematics.x();
     int    A = tgt.A();

     double xv     = TMath::Min(0.75, TMath::Max(0.05, x));
     double xv2    = xv  * xv;
     double xv3    = xv2 * xv;
     double xv4    = xv3 * xv;
     double xv5    = xv4 * xv;

     double f = 1.;
     
     // first factor goes from free nucleons to deuterium
     if(A >= 2) {
       f= 0.985*(1.+0.422*xv - 2.745*xv2 + 7.570*xv3 - 10.335*xv4 + 5.422*xv5);
     }

     double chiTM = this->ScalingVar(interaction);
     // 2nd factor goes from deuterium to iso-scalar iron
     if(A > 2) {
       // Fermi motion is already accounted for at high chiTM. To avoid double-counting, we limit the chiTM range:
       if(! interaction->TestBit(kIAssumeFreeNucleon) && chiTM > 0.65 ) chiTM = 0.65;
       f *= (1.096 - 0.38*chiTM - 0.3 * TMath::Exp(-23*chiTM) + 8 * pow(chiTM,15) ) ;
     }

     if( A == 79 || A == 82 ) { // Gold nd Lead F2(Au,Pb)/F2(Fe)
       f *= (0.932 + 2.461 * chiTM - 24.23*pow(chiTM,2) + 101.03*pow(chiTM,3) - 203.47*pow(chiTM,4) + 193.85*pow(chiTM,5) - 69.82*pow(chiTM,6)); 
     }

     if( A == 12 ) { // F2(Fe)/F2(C)
       f /= ( 0.919 + 1.844*chiTM - 12.73*pow(chiTM,2) + 36.89*pow(chiTM,3) - 46.77*pow(chiTM,4) + 21.22*pow(chiTM,5));
     }
     
  return f;

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
     LOG("DISSF", pDEBUG) << "Nuclear factor for x of " << x << "  = " << f;
#endif
  }

  return f;  
}
