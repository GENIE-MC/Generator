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
#include "Physics/DeepInelastic/XSection/BY21StrucFunc.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGUtils.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
BY21StrucFunc::BY21StrucFunc() :
  QPMDISStrucFuncBase("genie::BY21StrucFunc")
{
  this->Init();
}
//____________________________________________________________________________
BY21StrucFunc::BY21StrucFunc(string config):
  QPMDISStrucFuncBase("genie::BY21StrucFunc", config)
{
  this->Init();
}
//____________________________________________________________________________
BY21StrucFunc::~BY21StrucFunc()
{

}
//____________________________________________________________________________
void BY21StrucFunc::Configure(const Registry & config)
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
void BY21StrucFunc::Configure(string param_set)
{
  QPMDISStrucFuncBase::Configure(param_set);
  this->ReadBYParams();
}
//____________________________________________________________________________
void BY21StrucFunc::ReadBYParams(void)
{
  // vector mass 
  GetParam( "EL-Mv",fMv ) ;
  fMv2 = TMath::Power(fMv,2);
  
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
  GetParam( "BY-RQ2min", fRQ2min);
  GetParamDef( "BY-IncludeH", fIncludeH, true );
  GetParamDef( "BY-IncludeKCharm", fIncludeKCharm, true );
  GetParamDef( "BY-IncludeAxial", fIncludeAxial, true );
}
//____________________________________________________________________________
void BY21StrucFunc::Init(void)
{
  fMv   = 0;
  fMv2  = 0;
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
  fIncludeH = true;
  fIncludeAxial = true;
  fIncludeKCharm = true;
  fH0   = 0;
  fH1   = 0;
  fH2   = 0;
  fH3   = 0;
  fRQ2min = 0;
}
//____________________________________________________________________________
double BY21StrucFunc::ScalingVar(const Interaction * interaction, double Mf ) const
{
  // Overrides QPMDISStrucFuncBase::ScalingVar() to compute the BY scaling var
  // Eq. 15, https://arxiv.org/abs/2108.09240 
  const Kinematics & kine  = interaction->Kine();
  double x  = kine.x();
  double myQ2 = this->Q2(interaction);
  //myQ2 = TMath::Max(Q2,fQ2min);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__ 
  LOG("BodekYang", pDEBUG) << "Q2 at scaling var calculation = " << myQ2;
#endif

  double a  = TMath::Power( 2*kProtonMass*x, 2 ) / myQ2;
  double Mf2 = TMath::Power( Mf, 2 ) ; // Accounts for slow rescaling
  double xw =  2*x*(myQ2+Mf2+fB) / (myQ2*(1.+TMath::Sqrt(1+a)) + 2*fA*x);

  return xw;
}
//____________________________________________________________________________
void BY21StrucFunc::KVectorFactors(const Interaction * interaction,
				     double & kuv, double & kdv, double & kus, double & kds, double & kss ) const
{
  // Overrides QPMDISStrucFuncBase::KFactors() to compute the BY K factors for
  // u(valence), d(valence), u(sea), d(sea), s(sea);

  double myQ2  = this->Q2(interaction);
  double GD  = 1. / TMath::Power(1.+myQ2/fMv2, 2); // p elastic form factor
  double GD2 = TMath::Power(GD,2);

  // We include the BY21 K factors. The arxiv.org/pdf/2108.09240 publicatio also accounts for a low energy transfer correction (KLW)
  // The KLW is only important in the SIS region (W<2GeV). As we scale it altogether with the RES region, we do not need this factor.
  // It also causes issues at low-W. After discussing with A. Bodek, we agreed to comment this part out. 
  kuv = (1.-GD2)*(myQ2+fCv2U)/(myQ2+fCv1U); // K - u(valence)
  kdv = (1.-GD2)*(myQ2+fCv2D)/(myQ2+fCv1D); // K - d(valence)
  kus = myQ2/(myQ2+fCsU);                   // K - u(sea)
  kds = myQ2/(myQ2+fCsD);                   // K - d(sea)
  kss = myQ2/(myQ2+fCsS);                   // K - s(sea)	
}
//____________________________________________________________________________
void BY21StrucFunc::KAxialFactors(const Interaction * interaction,
				    double & kuv, double & kdv, double & kus, double & kds, double & kss ) const {

  // If requested, assume axial is equal to vector;
  if( !fIncludeAxial ) {
    KVectorFactors( interaction, kuv, kdv, kus, kds, kss );
    return ;
  }
  
  // https://arxiv.org/pdf/2108.09240 Sec 11.2 
  // We apply a different K factor for axial contribtions, as given in Eq. 49-50.
  // There is a single correction for valance quarks, and a different single correction for sea quarks.
  double myQ2  = this->Q2(interaction);
  double ksea = ( myQ2 + fPsA*fCsA ) / ( myQ2 + fCsA ) ;
  double kvalance = ( myQ2 + fPvA * fMv2/4. ) / ( myQ2 + fMv2/4. ) ; 

  kuv = kvalance;
  kdv = kvalance;
  kus = ksea;
  kds = ksea;
  kss = ksea;

}

//____________________________________________________________________________
double BY21StrucFunc::H(const Interaction * interaction) const {
  if( !fIncludeH ) return 1;
  
  // Overrides QPMDISStrucFuncBase::H() function to compute the correction of the BY 2021 update
  // The correction is given by Eq. 34 of https://arxiv.org/pdf/2108.09240
  // It accounts for the difference in the QCD higher order corrections in F2 and xF3
  const Kinematics & kinematics = interaction->Kine();
  double bjx = kinematics.x();
  return fH0 + fH1 * bjx + fH2 * pow(bjx,2) + fH3 * pow(bjx,3);
}

double BY21StrucFunc::KCharm(const Interaction * interaction, double Mf) const {
  // Overrides QPMDISStrucFuncBase::KCharm() function to compute the correction of the BY 2021 update
  // For production of heavy quarks, the relation between the structure function changes and
  // instead of kinematics x, we need to use chi. See Eq. 5.1 of PhysRevD.14.1829 for details.
  // For a simple modification of Bjorken x, chi'=x(1+Mf2/Q2), the correction factor is equivalent
  // to KCharm, defined in Sec 8 of https://arxiv.org/pdf/2108.09240
  // In our case, to keep it general, we define it explicitly as x/chi_w
  if( !fIncludeKCharm ) return 1.;
  
  double K = interaction->Kine().x()/ScalingVar(interaction, fMc);

  return K;
  
}
//____________________________________________________________________________
double BY21StrucFunc::R(const Interaction * interaction) const {
  if(!fIncludeR) return 0;
  
  // Evaluate correction for Q2 < 0.3 GeV2/c4 according to Sec. 7 of https://arxiv.org/pdf/2108.09240
  double Q2 = this->Q2(interaction);
  double Q2_int = this->Q2(interaction) ;

  // We use the R1998 parameterization only at Q2>0.3 GeV2.
  // At lower Q2, we freeze the function at Q2 =0.3 GeV2
  // The parameter, fRQ2min, is setup to 0.3 in the configuration. The default value can be changed.

  if( Q2_int < fRQ2min ) {
    // Freeze R at Q2 = 0.3 GeV2
    Q2 = fRQ2min ;
  }
  double Q4 = pow(Q2,2);
  double Q8 = pow(Q4,2);

  const Kinematics & kinematics = interaction->Kine();
  double x = kinematics.x();
  double x2 = pow(x,2);
  double x3 = pow(x,3);
  
  double Theta = 1 + 12.0 * ( Q2 / (Q2+1.) ) * pow(0.125,2)/(pow(0.125,2)+x2);
  
  // R1998 is defined as the average of Ra, Rb and Rc, each parameterized to accomodate new data at low x
  // arXiv:hep-ex/9808028
  const double a1 = 0.0485;
  const double a2 = 0.5470;
  const double a3 = 2.0621;
  const double a4 = -0.3804;
  const double a5 = 0.5090;
  const double a6 = -0.0285;
  
  double Ra = (a1 * Theta / TMath::Log(Q2/0.04) ) ;
  Ra += a2 * ( 1 + a4 * x + a5 * x2 ) * pow( x, a6 ) / pow( Q8 + pow(a3,4), 1./4. ) ;

  const double b1 = 0.0481;
  const double b2 = 0.6114;
  const double b3 = -0.3509;
  const double b4 = -0.4611;
  const double b5 = 0.7172;
  const double b6 = -0.0317;

  double Rb = b1 * Theta / TMath::Log(Q2/0.04) ;
  Rb += (b2 / Q2 + b3 / (Q4 + 0.09) ) * (1 + b4*x + b5*x2) * pow(x,b6);
	 
  const double c1 = 0.0577;
  const double c2 = 0.4644;
  const double c3 = 1.8288;
  const double c4 = 12.3708;
  const double c5 = -43.1043;
  const double c6 = 41.7415;
  double Q2thr = c4*x + c5*x2 + c6*x3;
  double Rc = c1 * Theta / TMath::Log( Q2/0.04 ) ;
  Rc += c2 / sqrt( pow(Q2 - Q2thr,2) + pow(c3,2));
  
  double R1998 = (Ra + Rb + Rc) / 3. ;
      
  // At Q2 < 0.3 GeV2, we add a K factor multipying R(Q2=0.3)
  // for a smooth transtion down to Q2 = 0 (the photonproduction limit)
  if( Q2_int < fRQ2min ) {
    double Q4_int = pow(Q2_int,2);
    double R0 = (pow(fRQ2min,2) + 1 ) / fRQ2min ;
    R1998 *= R0  * Q2_int / ( Q4_int + 1 ) ;
  }
  
  return R1998 ; 
}

