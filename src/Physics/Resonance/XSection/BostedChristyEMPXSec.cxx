//_________________________________________________________________________
/*
 Copyright (c) 2003-2021, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 For the class documentation see the corresponding header file.
*/
//_________________________________________________________________________
#include <TMath.h>

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Physics/Resonance/XSection/BostedChristyEMPXSec.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::utils;

//____________________________________________________________________________
BostedChristyEMPXSec::BostedChristyEMPXSec() :
XSecAlgorithmI("genie::BostedChristyEMPXSec")
{
  
}
//____________________________________________________________________________
BostedChristyEMPXSec::BostedChristyEMPXSec(string config) :
XSecAlgorithmI("genie::BostedChristyEMPXSec", config)
{
  
}
//____________________________________________________________________________
BostedChristyEMPXSec::~BostedChristyEMPXSec()
{
  
}
//____________________________________________________________________________
// resonance cross section fit of electron-proton scattering data
// sf=0 \sigma_T
// sf=1 \sigma_L
double BostedChristyEMPXSec::sigmaR(int sf, double Q2, double W) const
{
  double W2 = W*W;
  // proton mass
  double Mp = 0.9382727;
  // pion mass
  double Mpi = 0.135;
  // eta-meson mass
  double Meta = 0.547;
  
  double Mp2 = Mp*Mp;
  double Mpi2 = Mpi*Mpi;
  double Meta2 = Meta*Meta;
    
  // Calculate kinematics needed for threshold Relativistic B-W
  
  // Eq. (10)
  double k   = (W2 - Mp2)/2./Mp;
  // Eq. (11)
  double kcm = (W2 - Mp2)/2./W;
  // mesons energy and momentim
  double Epicm = (W2 + Mpi2 - Mp2)/2./W;                                  // pion energy in CMS
  double ppicm = TMath::Sqrt(TMath::Max(0.0, Epicm*Epicm - Mpi2));        // pion momentum in CMS
  double Epi2cm = (W2 + 4*Mpi2 - Mp2)/2./W;                               // two pion energy in CMS
  double ppi2cm = TMath::Sqrt(TMath::Max(0.0, Epi2cm*Epi2cm - 4*Mpi2));        // two pion energi n CMS
  double Eetacm = (W2 + Meta2 - Mp2 )/2./W;                               // eta energy in CMS
  double petacm = TMath::Sqrt(TMath::Max(0.0, Eetacm*Eetacm - Meta2));         // eta energy in CMS

  double xsec = 0.0;
  
  // going through seven resonances
  for (int i=0;i<7;i++)
  {
    // resonance mass squared
    double fmass2 = fmass[i]*fmass[i];
    // resonance damping parameter
    double x0 = i!=0?0.215:0.14462;
    // Eq. (12)
    double kr = (fmass2-Mp2)/2./Mp;                                 
    // Eq. (13)
    double kcmr = (fmass2-Mp2)/2./fmass[i];                        
    
    // formulas analogous to the above with substitution W->MR_i
    double Epicmr = (fmass2 + Mpi2 - Mp2)/2./fmass[i];
    double ppicmr = TMath::Sqrt(TMath::Max(0.0, Epicmr*Epicmr - Mpi2));
    double Epi2cmr = (fmass2 + 4.*Mpi2 - Mp2)/2./fmass[i];
    double ppi2cmr = TMath::Sqrt(TMath::Max(0.0, Epi2cmr*Epi2cmr - 4.*Mpi2));
    double Eetacmr = (fmass2 + Meta2 - Mp2)/2./fmass[i];
    double petacmr =  TMath::Sqrt(TMath::Max(0.0, Eetacmr*Eetacmr - Meta2));

    //   Calculate partial widths 
    // Eq. (15) for single pion
    double pwid0 = fwidth[i]*TMath::Power(ppicm/ppicmr, 1.+2.*fang[i])*
                             TMath::Power((ppicmr*ppicmr + x0*x0)/(ppicm*ppicm+x0*x0), fang[i]);                        // 1-pion decay mode
    // Eq. (16) for two pions
    double pwid1 = W/fmass[i]*fwidth[i]*TMath::Power(ppi2cm/ppi2cmr, 4.+2.*fang[i])*
                                        TMath::Power((ppi2cmr*ppi2cmr + x0*x0)/(ppi2cm*ppi2cm + x0*x0), 2.+fang[i]);   //  2-pion decay mode
    double pwid2 = 0.;                                                                                                 //  eta decay mode

    // Eq. (15) for eta
    if(i==1 || i==4)
      pwid2 =  fwidth[i]*TMath::Power(petacm/petacmr, 1.+2.*fang[i])*TMath::Power((petacmr*petacmr + x0*x0)/(petacm*petacm + x0*x0), fang[i]);  // eta decay only for S11's 
    
    // Eq. (17)
    double pgam = fwidth[i]*(kcm/kcmr)*(kcm/kcmr)*(kcmr*kcmr+x0*x0)/(kcm*kcm+x0*x0);
    //Eq. (14)
    double width = fbr[i][0]*pwid0+fbr[i][1]*pwid1+fbr[i][2]*pwid2;
    
    //    Begin resonance Q^2 dependence calculations
    double A;
    
    if (sf==0)
       // Eq. (18)
       A = frescoefT[i][0]*(1.+frescoefT[i][1]*Q2/(1.+frescoefT[i][2]*Q2))/TMath::Power(1.+Q2/0.91, frescoefT[i][3]);
    else
       // Eq. (19)
       A = frescoefL[i][0]*Q2/(1.+frescoefL[i][1]*Q2)*TMath::Exp(-1.*frescoefL[i][2]*Q2);
    
    
    // Eq. (9)
    double BW = kr/k*kcmr/kcm/fwidth[i]*width*pgam/((W2 - fmass2)*(W2 - fmass2) + fmass2*width*width);
    
    // Eq. (8) divided by W
    xsec += BW*A*A;   
  }
  // restore factor W in Eq. (8)
  return W*xsec;
}
//____________________________________________________________________________
// nonresonance cross section fit of electron-proton scattering data
// sf=0 \sigma_T
// sf=1 \sigma_L
double BostedChristyEMPXSec::sigmaNR(int sf, double Q2, double W) const
{
  
  double W2 = W*W;
  double Mp = 0.9382727;
  double Mpi = 0.135;
  double Mp2 = Mp*Mp;
  
  double Wdif = W - (Mp + Mpi);
  
  double m0 = 4.2802; // GeV
  double Q20 = sf==0?0.05:0.125; 
  double xpr = 1./(1.+(W2-(Mp+Mpi)*(Mp+Mpi))/(Q2+Q20));      // Eq.(22)
  
  double xsec = 0.0;
  
      
  if (sf==0)
  {
    
    for (int i=0;i<2;i++)
    {
      double h_nr = fnrcoefT[i][0]/TMath::Power(Q2+fnrcoefT[i][1], fnrcoefT[i][2]+fnrcoefT[i][3]*Q2+fnrcoefT[i][4]*Q2*Q2);   // Eq. (21)
      xsec += h_nr*TMath::Power(Wdif, 1.5+i);
    }

    xsec *= xpr;
  }
  else
  {
     double xb = Q2/(Q2+W2-Mp2);
     double t = TMath::Log(TMath::Log((Q2+m0)/0.330/0.330)/TMath::Log(m0/0.330/0.330));                      // Eq. (24)
     xsec += fnrcoefL[0]*TMath::Power(1.-xpr, fnrcoefL[2]+fnrcoefL[1]*t)/(1.-xb)/(Q2+Q20)
                     *TMath::Power(Q2/(Q2+Q20), fnrcoefL[3])*TMath::Power(xpr, fnrcoefL[4]+fnrcoefL[5]*t);   // Eq. (23)
  }
                       
  return xsec;
}
//____________________________________________________________________________
double BostedChristyEMPXSec::XSec(
                 const Interaction * interaction, KinePhaseSpace_t kps) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;
  
  // Get kinematical parameters
  const Kinematics & kinematics = interaction -> Kine();
  const InitialState & init_state = interaction -> InitState();
  
  double E  = init_state.ProbeE(kRfHitNucRest);
  double W  = kinematics.W();
  double Q2 = kinematics.Q2();
  double W2 = W*W;
  double Mp = 0.9382727;
  double Mp2 = Mp*Mp;
  
  double sigmaT = sigmaR(0, Q2, W) + sigmaNR(0, Q2, W);
  double sigmaL = sigmaR(1, Q2, W) + sigmaNR(1, Q2, W);

  double nu = (W2-Mp2+Q2)/2./Mp;
  double nu2 = nu*nu;
  double Eprime = E - nu;
  double sin2theta_2 = Q2/4/E/Eprime;
  double cos2theta_2 = 1 - sin2theta_2;
  double tan2theta_2 = sin2theta_2/cos2theta_2;
  double eps = 1./(1. + 2*(1+nu2/Q2)*tan2theta_2);            // Eq. (4)
  double Gamma = kAem*Eprime*(W2-Mp2)/Q2/Mp/E/(1-eps)/4/kPi2;   // Eq. (5)
  double xsec = Gamma*(sigmaT+eps*sigmaL);                // Eq. (3) d2xsec/dOmegadEprime
  double jacobian = W*kPi/E/Eprime/Mp;  
  xsec*= jacobian;                                            // d2xsec/dOmegadEprime-> d2xsec/dWdQ2
  
  // The algorithm computes d^2xsec/dWdQ2
  // Check whether variable tranformation is needed
  if(kps!=kPSWQ2fE) {
    double J = utils::kinematics::Jacobian(interaction,kPSWQ2fE,kps);
    xsec *= J;
  }
  
  return xsec*units::ub;                              
  
}
//____________________________________________________________________________
double BostedChristyEMPXSec::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool BostedChristyEMPXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const InitialState & init_state = interaction->InitState();
  const ProcessInfo &  proc_info  = interaction->ProcInfo();

  if(!proc_info.IsDeepInelastic()) return false;
 

  int  hitnuc = init_state.Tgt().HitNucPdg();
  bool is_pn = (pdg::IsProton(hitnuc) || pdg::IsNeutron(hitnuc));

  if (!is_pn) return false;

  int  probe   = init_state.ProbePdg();
  bool is_em   = proc_info.IsEM();

  bool l_em    = (pdg::IsChargedLepton(probe) && is_em  );

  if (!l_em) return false;

  return true;
}

//____________________________________________________________________________
void BostedChristyEMPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void BostedChristyEMPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void BostedChristyEMPXSec::LoadConfig(void)
{  
  // single pion branching ratios
  fbr[0][0] = 1.0;       //  P33(1232)       
  fbr[1][0] = 0.45;      //  S11(1535)   
  fbr[2][0] = 0.65;      //  D13(1520)
  fbr[3][0] = 0.65;      //  F15(1680)
  fbr[4][0] = 0.4;       //  S11(1650)
  fbr[5][0] = 0.65;      //  P11(1440) roper 
  fbr[6][0] = 0.50;      //  F37(1950)

  //  eta branching ratios
  fbr[0][2] = 0.0;       //  P33(1232)
  fbr[1][2] = 0.45;      //  S11(1535) 
  fbr[2][2] = 0.0;       //  D13(1520)
  fbr[3][2] = 0.0;       //  F15(1680)
  fbr[4][2] = 0.1;       //  S11(1650)
  fbr[5][2] = 0.0;       //  P11(1440) roper   
  fbr[6][2] = 0.0;       //  F37(1950)

  //  2-pion branching ratios 
  for (int i=0;i<7;i++)
    fbr[i][1] = 1.-fbr[i][0]-fbr[i][2];
    
  
  // Meson angular momentum
  fang[0] = 1;          //  P33(1232)
  fang[1] = 0;          //  S11(1535)
  fang[2] = 2;          //  D13(1520)
  fang[3] = 3;          //  F15(1680)
  fang[4] = 0;          //  S15(1650)
  fang[5] = 1;          //  P11(1440) roper   
  fang[6] = 3;          //  F37(1950)
  
  // Meson mass
  fmass[0] = 1.2298;    //  P33(1232)
  fmass[1] = 1.5304;    //  S11(1535)
  fmass[2] = 1.5057;    //  D13(1520)
  fmass[3] = 1.6980;    //  F15(1680)
  fmass[4] = 1.6650;    //  S15(1650)
  fmass[5] = 1.4333;    //  P11(1440) roper   
  fmass[6] = 1.9341;    //  F37(1950)

  // Meson width        
  fwidth[0] = 0.135730; //  P33(1232)
  fwidth[1] = 0.220000; //  S11(1535)
  fwidth[2] = 0.082956; //  D13(1520)
  fwidth[3] = 0.095782; //  F15(1680)
  fwidth[4] = 0.109360; //  S15(1650)
  fwidth[5] = 0.379440; //  P11(1440) roper   
  fwidth[6] = 0.380000; //  F37(1950)
  
  //  P33(1232)
  frescoefT[0][0] = 7.7805;   // Table III, A^i_T(0)
  frescoefT[0][1] = 4.2291;   // Table III, a_i
  frescoefT[0][2] = 1.2598;   // Table III, b_i
  frescoefT[0][3] = 2.1242;   // Table III, c_i
  //  S11(1535)
  frescoefT[1][0] = 6.3351;   // Table III, A^i_T(0)
  frescoefT[1][1] = 6823.2;   // Table III, a_i
  frescoefT[1][2] = 33521.0;  // Table III, b_i
  frescoefT[1][3] = 2.5686;   // Table III, c_i
  //  D13(1520)
  frescoefT[2][0] = 0.60347;  // Table III, A^i_T(0)
  frescoefT[2][1] = 21.24;    // Table III, a_i
  frescoefT[2][2] = 0.055746; // Table III, b_i
  frescoefT[2][3] = 2.4886;   // Table III, c_i
   //  F15(1680)
  frescoefT[3][0] = 2.3305;   // Table III, A^i_T(0)
  frescoefT[3][1] = -0.28789; // Table III, a_i
  frescoefT[3][2] = 0.18607;  // Table III, b_i
  frescoefT[3][3] = 0.063534; // Table III, c_i
  //  S15(1650)
  frescoefT[4][0] = 1.979;    // Table III, A^i_T(0)
  frescoefT[4][1] = -0.56175; // Table III, a_i
  frescoefT[4][2] = 0.38964;  // Table III, b_i
  frescoefT[4][3] = 0.54883;  // Table III, c_i
  //  P11(1440) roper 
  frescoefT[5][0] = 0.022506; // Table III, A^i_T(0) 
  frescoefT[5][1] = 462.13;   // Table III, a_i
  frescoefT[5][2] = 0.19221;  // Table III, b_i
  frescoefT[5][3] = 1.9141;   // Table III, c_i
  //  F37(1950)
  frescoefT[6][0] = 3.4187;   // Table III, A^i_T(0)
  frescoefT[6][1] = 0.;       // Table III, a_i
  frescoefT[6][2] = 0.;       // Table III, b_i
  frescoefT[6][3] = 1.;       // Table III, c_i
  
  //  P33(1232)
  frescoefL[0][0] =  29.414;   // Table III, A^i_L(0)
  frescoefL[0][1] =  19.91;    // Table III, d_i
  frescoefL[0][2] =  0.22587;  // Table III, e_i
  //  S11(1535)                
  frescoefL[1][0] =  0.;       // Table III, A^i_L(0)
  frescoefL[1][1] =  3856.5;   // Table III, d_i
  frescoefL[1][2] =  0.65717;  // Table III, e_i
  //  D13(1520)                
  frescoefL[2][0] =  157.92;   // Table III, A^i_L(0)
  frescoefL[2][1] =  97.046;   // Table III, d_i
  frescoefL[2][2] =  0.31042;  // Table III, e_i
   //  F15(1680)               
  frescoefL[3][0] =  4.216;    // Table III, A^i_L(0)
  frescoefL[3][1] =  0.0382;   // Table III, d_i
  frescoefL[3][2] =  1.2182;   // Table III, e_i
  //  S15(1650)                
  frescoefL[4][0] =  13.764;   // Table III, A^i_L(0)
  frescoefL[4][1] =  0.31393;  // Table III, d_i
  frescoefL[4][2] =  2.9997;   // Table III, e_i    
  //  P11(1440) roper          
  frescoefL[5][0] =  5.5124;   // Table III, A^i_L(0)
  frescoefL[5][1] =  0.053743; // Table III, d_i  
  frescoefL[5][2] =  1.3091;   // Table III, e_i
  //  F37(1950)                
  frescoefL[6][0] =  11.0;     // Table III, A^i_L(0)
  frescoefL[6][1] =  1.8951;   // Table III, d_i
  frescoefL[6][2] =  0.51376;  // Table III, e_i
  
  fnrcoefT[0][0] =   246.06;          // Table IV, \sigma^NR,1_T(0)
  fnrcoefT[0][1] =   0.067469;        // Table IV, aT_1
  fnrcoefT[0][2] =   1.3501;          // Table IV, bT_1
  fnrcoefT[0][3] =   0.12054;         // Table IV, cT_1
  fnrcoefT[0][4] =   -.38495E-02;     // Table IV, dT_1
  fnrcoefT[1][0] =   -89.36;          // Table IV, \sigma^NR,1_T(0) 
  fnrcoefT[1][1] =   0.20977;         // Table IV, aT_2
  fnrcoefT[1][2] =   1.5715;          // Table IV, bT_2
  fnrcoefT[1][3] =   0.90736E-01;     // Table IV, cT_2
  fnrcoefT[1][4] =   1.0362e-02;      // Table IV, dT_2
  
  fnrcoefL[0] = 86.746;               // Table IV, \sigma^NR,1_L(0) 
  fnrcoefL[1] = 4.0864e-005;          // Table IV, aL_2
  fnrcoefL[2] = 4.0294;               // Table IV, bL_2
  fnrcoefL[3] = 3.1285;               // Table IV, cL_2
  fnrcoefL[4] = 0.33403;              // Table IV, dL_2
  fnrcoefL[5] = 4.9623;               // Table IV, eL_2
 

}
