//_________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 For the class documentation see the corresponding header file.
*/
//_________________________________________________________________________
#include <vector>
#include <string>
#include <sstream>

#include <TMath.h>
#include <TDecayChannel.h>

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/BaryonResUtils.h"
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
// resonance cross section fit of electron-proton and electron-deuterium  scattering data
// sf=0 \sigma_T
// sf=1 \sigma_L
double BostedChristyEMPXSec::sigmaR(int sf, double Q2, double W, bool isDeuterium=false) const
{
  const std::array<std::array<double, 3>, 7> &fBR = !isDeuterium?fBRp:fBRD;
  const std::array<std::array<double, 4>, 7> &fRescoefT = !isDeuterium?fRescoefTp:fRescoefTD;
  if (isDeuterium)
    sf=0;
  
  double W2 = W*W;
  // proton mass
  double Mp = fMP;
  // pion mass
  double Mpi = fMpi0;
  // eta-meson mass
  double Meta = fMeta;
  
  double Mp2 = Mp*Mp;
  double Mpi2 = Mpi*Mpi;
  double Meta2 = Meta*Meta;
    
  // Calculate kinematics needed for threshold Relativistic B-W
  
  // Ref.1, Eq. (10)
  double k   = (W2 - Mp2)/2./Mp;
  // Ref.1, Eq. (11)
  double kcm = (W2 - Mp2)/2./W;
  // mesons energy and momentim
  double Epicm = (W2 + Mpi2 - Mp2)/2./W;                                  // pion energy in CMS
  double ppicm = TMath::Sqrt(TMath::Max(0.0, Epicm*Epicm - Mpi2));        // pion momentum in CMS
  double Epi2cm = (W2 + 4*Mpi2 - Mp2)/2./W;                               // two pion energy in CMS
  double ppi2cm = TMath::Sqrt(TMath::Max(0.0, Epi2cm*Epi2cm - 4*Mpi2));   // two pion energi n CMS
  double Eetacm = (W2 + Meta2 - Mp2 )/2./W;                               // eta energy in CMS
  double petacm = TMath::Sqrt(TMath::Max(0.0, Eetacm*Eetacm - Meta2));    // eta energy in CMS

  double xsec = 0.0;
  
  // going through seven resonances
  for (int i=0;i<7;i++)
  {
    // resonance mass squared
    double MassRes2 = fMassRes[i]*fMassRes[i];
    // resonance damping parameter
    double x0 = i!=0?0.215:!isDeuterium?0.14462:0.1446;
    // Ref.1, Eq. (12)
    double kr = (MassRes2-Mp2)/2./Mp;                                 
    // Ref.1, Eq. (13)
    double kcmr = (MassRes2-Mp2)/2./fMassRes[i];                        
    
    // formulas analogous to the above with substitution W->MR_i
    double Epicmr = (MassRes2 + Mpi2 - Mp2)/2./fMassRes[i];
    double ppicmr = TMath::Sqrt(TMath::Max(0.0, Epicmr*Epicmr - Mpi2));
    double Epi2cmr = (MassRes2 + 4.*Mpi2 - Mp2)/2./fMassRes[i];
    double ppi2cmr = TMath::Sqrt(TMath::Max(0.0, Epi2cmr*Epi2cmr - 4.*Mpi2));
    double Eetacmr = (MassRes2 + Meta2 - Mp2)/2./fMassRes[i];
    double petacmr =  TMath::Sqrt(TMath::Max(0.0, Eetacmr*Eetacmr - Meta2));

    //   Calculate partial widths 
    // Ref.1, Eq. (15) for single pion
    double pwid0 = fWidthRes[i]*TMath::Power(ppicm/ppicmr, 1.+2.*fAngRes[i])*
                             TMath::Power((ppicmr*ppicmr + x0*x0)/(ppicm*ppicm+x0*x0), fAngRes[i]);                       // 1-pion decay mode
    // Ref.1, Eq. (16) for two pions
    double pwid1 = 0;
    if (!isDeuterium || (isDeuterium && i!=1))
      pwid1 = W/fMassRes[i]*fWidthRes[i]*TMath::Power(ppi2cm/ppi2cmr, 4.+2.*fAngRes[i])*
                                        TMath::Power((ppi2cmr*ppi2cmr + x0*x0)/(ppi2cm*ppi2cm + x0*x0), 2.+fAngRes[i]);   //  2-pion decay mode
    else
      pwid1 =  fWidthRes[i]*TMath::Power(petacm/petacmr, 1.+2.*fAngRes[i])*TMath::Power((ppi2cmr*ppi2cmr + x0*x0)/(ppi2cm*ppi2cm + x0*x0), fAngRes[i]);
    
    
    double pwid2 = 0.;                                                                                                    //  eta decay mode
    // Ref.1, Eq. (15) for eta
    if(!isDeuterium && (i==1 || i==4))
      pwid2 =  fWidthRes[i]*TMath::Power(petacm/petacmr, 1.+2.*fAngRes[i])*TMath::Power((petacmr*petacmr + x0*x0)/(petacm*petacm + x0*x0), fAngRes[i]);  // eta decay only for S11's 
    
    // Ref.1, Eq. (17)
    double pgam = fWidthRes[i]*(kcm/kcmr)*(kcm/kcmr)*(kcmr*kcmr+x0*x0)/(kcm*kcm+x0*x0);
    // Ref.1, Eq. (14)
    double width = fBR[i][0]*pwid0+fBR[i][1]*pwid1+fBR[i][2]*pwid2;
    
    //  Begin resonance Q^2 dependence calculations
    double A;
    
    if (sf==0)
       // Ref.1, Eq. (18)
       A = fRescoefT[i][0]*(1.+fRescoefT[i][1]*Q2/(1.+fRescoefT[i][2]*Q2))/TMath::Power(1.+Q2/0.91, fRescoefT[i][3]);
    else
       // Ref.1, Eq. (19)
       A = fRescoefL[i][0]*Q2/(1.+fRescoefL[i][1]*Q2)*TMath::Exp(-1.*fRescoefL[i][2]*Q2);
    
    
    // Ref.1, Eq. (9)
    double BW = kr/k*kcmr/kcm/fWidthRes[i]*width*pgam/((W2 - MassRes2)*(W2 - MassRes2) + MassRes2*width*width);
    
    // Ref.1, Eq. (8) divided by W
    xsec += BW*A*A;   
  }
  // restore factor W in Ref.1, Eq. (8)
  return W*xsec*units::ub;
}
//____________________________________________________________________________
// nonresonance cross section fit of electron-proton and electron-deuterium scattering data
// sf=0 \sigma_T
// sf=1 \sigma_L
double BostedChristyEMPXSec::sigmaNR(int sf, double Q2, double W, bool isDeuterium=false) const
{
  const std::array<std::array<double, 5>, 2> &fNRcoefT = !isDeuterium?fNRcoefTp:fNRcoefTD;
  if (isDeuterium)
    sf=0;
  double W2 = W*W;
  double Mp = fMP;
  double Mpi = fMpi0;
  double Mp2 = Mp*Mp;
  
  double Wdif = W - (Mp + Mpi);
  
  double m0 = (sf==0) ? 0.125 : 4.2802;                      //Ref.1, Eqs.(22, 24)
  
  double Q20 = (sf==0) ? 0.05 : 0.125;                       //Ref.1, Eqs.(22, 24)
  
  double xpr = 1./(1.+(W2-(Mp+Mpi)*(Mp+Mpi))/(Q2+Q20));      // Ref.1, Eq.(22)
  
  double xsec = 0.0;
  
      
  if (sf==0)
  {
    
    for (int i=0;i<2;i++)
    {
      double h_nr = fNRcoefT[i][0]/TMath::Power(Q2+fNRcoefT[i][1], fNRcoefT[i][2]+fNRcoefT[i][3]*Q2+fNRcoefT[i][4]*Q2*Q2);   // Ref.1, Eq. (21)
      xsec += h_nr*TMath::Power(Wdif, 1.5+i);
    }

    xsec *= xpr;
  }
  else
  {
     double xb = Q2/(Q2+W2-Mp2);
     double t = TMath::Log(TMath::Log((Q2+m0)/0.330/0.330)/TMath::Log(m0/0.330/0.330));                      // Ref.1, Eq. (24)
     xsec += fNRcoefL[0]*TMath::Power(1.-xpr, fNRcoefL[2]+fNRcoefL[1]*t)/(1.-xb)/(Q2+Q20)
                     *TMath::Power(Q2/(Q2+Q20), fNRcoefL[3])*TMath::Power(xpr, fNRcoefL[4]+fNRcoefL[5]*t);   // Ref.1, Eq. (23)
  }
                       
  return xsec*units::ub;
}
//___________________________________________________________________________
// Calculate proton and neutron with Fermi smearing of a nulei
void BostedChristyEMPXSec::FermiSmearingA(double Q2, double W, double pF, double Es, double & F1p, double & F1d, double & sigmaT, double & sigmaL) const
{
  // The numbers in arrays bellow were not supposed to change in the original 
  // fortran code and therefore are not configurable
  static constexpr std::array<double, 99> fyp
  {0.0272,0.0326,0.0390,0.0464,0.0551,0.0651,0.0766,0.0898,0.1049,
   0.1221,0.1416,0.1636,0.1883,0.2159,0.2466,0.2807,0.3182,0.3595,
   0.4045,0.4535,0.5066,0.5637,0.6249,0.6901,0.7593,0.8324,0.9090,
   0.9890,1.0720,1.1577,1.2454,1.3349,1.4254,1.5163,1.6070,1.6968,
   1.7849,1.8705,1.9529,2.0313,2.1049,2.1731,2.2350,2.2901,2.3379,
   2.3776,2.4090,2.4317,2.4454,2.4500,2.4454,2.4317,2.4090,2.3776,
   2.3379,2.2901,2.2350,2.1731,2.1049,2.0313,1.9529,1.8705,1.7849,
   1.6968,1.6070,1.5163,1.4254,1.3349,1.2454,1.1577,1.0720,0.9890,
   0.9090,0.8324,0.7593,0.6901,0.6249,0.5637,0.5066,0.4535,0.4045,
   0.3595,0.3182,0.2807,0.2466,0.2159,0.1883,0.1636,0.1416,0.1221,
   0.1049,0.0898,0.0766,0.0651,0.0551,0.0464,0.0390,0.0326,0.0272};
   
  static constexpr std::array<double, 99> xxp
  {-3.000,-2.939,-2.878,-2.816,-2.755,-2.694,-2.633,-2.571,-2.510,
   -2.449,-2.388,-2.327,-2.265,-2.204,-2.143,-2.082,-2.020,-1.959,
   -1.898,-1.837,-1.776,-1.714,-1.653,-1.592,-1.531,-1.469,-1.408,
   -1.347,-1.286,-1.224,-1.163,-1.102,-1.041,-0.980,-0.918,-0.857,
   -0.796,-0.735,-0.673,-0.612,-0.551,-0.490,-0.429,-0.367,-0.306,
   -0.245,-0.184,-0.122,-0.061, 0.000, 0.061, 0.122, 0.184, 0.245,
    0.306, 0.367, 0.429, 0.490, 0.551, 0.612, 0.673, 0.735, 0.796,
    0.857, 0.918, 0.980, 1.041, 1.102, 1.163, 1.224, 1.286, 1.347,
    1.408, 1.469, 1.531, 1.592, 1.653, 1.714, 1.776, 1.837, 1.898,
    1.959, 2.020, 2.082, 2.143, 2.204, 2.265, 2.327, 2.388, 2.449,
    2.510, 2.571, 2.633, 2.694, 2.755, 2.816, 2.878, 2.939, 3.000};
   
  double MN = fPM;
  double MN2 = MN*MN;
  double Mp = fMP;
  double Mp2 = Mp*Mp;
  double W2 = W*W;
  
  double nu = (W2 - MN2 + Q2)/2./MN;
  double qv = TMath::Sqrt(nu*nu + Q2);
  // assume this is 2*pf*qv
  double dW2dpF = 2.*qv;
  double dW2dEs = 2.*(nu + MN); 
  // switched to using 99 bins!
  F1p = 0;
  F1d = 0;
  sigmaT = 0;
  sigmaL = 0;
  for (int i=0; i<99; i++)
  {
    double fyuse = fyp[i]/100.;
    double W2p = W2 + xxp[i]*pF*dW2dpF - Es*dW2dEs;
    if(W2p>1.159)
    {
       //proton 
       double Wp = TMath::Sqrt(W2p);
       double sigmaTp = sigmaR(0, Q2, Wp) + sigmaNR(0, Q2, Wp);
       double sigmaLp = sigmaR(1, Q2, Wp) + sigmaNR(1, Q2, Wp); 
       double F1pp = sigmaTp*(W2p-Mp2)/8./kPi2/kAem;
       //neutron
       double sigmaTd = sigmaR(0, Q2, Wp, true) + sigmaNR(0, Q2, Wp, true);
       double F1dp = sigmaTd*(W2p-Mp2)/8./kPi2/kAem;
       F1d += F1dp*fyuse;
       F1p += F1pp*fyuse;
       sigmaT += sigmaTp*fyuse;
       sigmaL += sigmaLp*fyuse;
    }
  }

}
//___________________________________________________________________________
// Calculate proton and neutron with Fermi smearing of a deuteron 
void BostedChristyEMPXSec::FermiSmearingD(double Q2, double W, double & F1, double & R, double & sigmaT, double & sigmaL, bool isDeuterium=false) const
{
  // The numbers in arrays bellow were not supposed to change in the original 
  // fortran code and therefore are not configurable
  static constexpr std::array<double, 20> fyd
  {0.4965, 0.4988, 0.4958, 0.5008, 0.5027, 0.5041, 0.5029, 0.5034, 
   0.4993, 0.5147, 0.5140, 0.4975, 0.5007, 0.4992, 0.4994, 0.4977, 
   0.5023, 0.4964, 0.4966, 0.4767};
   
  static constexpr std::array<double, 20> avpz
  {-0.1820,-0.0829,-0.0590,-0.0448,-0.0345,-0.0264,-0.0195, -0.0135,
   -0.0079,-0.0025, 0.0029, 0.0083, 0.0139, 0.0199, 0.0268,  0.0349, 
    0.0453, 0.0598, 0.0844, 0.1853};
  
  static constexpr std::array<double, 20> avp2
  {0.0938, 0.0219, 0.0137, 0.0101, 0.0081, 0.0068, 0.0060, 0.0054, 
   0.0051, 0.0049, 0.0050, 0.0051, 0.0055, 0.0060, 0.0069, 0.0081, 
   0.0102, 0.0140, 0.0225, 0.0964};
  
  // Look up tables for deuteron in fine bins for sub threshold
  static constexpr std::array<double, 200> fydf
  {0.00001,0.00002,0.00003,0.00005,0.00006,0.00009,0.00010,0.00013,
   0.00015,0.00019,0.00021,0.00026,0.00029,0.00034,0.00038,0.00044,
   0.00049,0.00057,0.00062,0.00071,0.00078,0.00089,0.00097,0.00109,
   0.00119,0.00134,0.00146,0.00161,0.00176,0.00195,0.00211,0.00232,
   0.00252,0.00276,0.00299,0.00326,0.00352,0.00383,0.00412,0.00447,
   0.00482,0.00521,0.00560,0.00603,0.00648,0.00698,0.00747,0.00803,
   0.00859,0.00921,0.00985,0.01056,0.01126,0.01205,0.01286,0.01376,
   0.01467,0.01569,0.01671,0.01793,0.01912,0.02049,0.02196,0.02356,
   0.02525,0.02723,0.02939,0.03179,0.03453,0.03764,0.04116,0.04533,
   0.05004,0.05565,0.06232,0.07015,0.07965,0.09093,0.10486,0.12185,
   0.14268,0.16860,0.20074,0.24129,0.29201,0.35713,0.44012,0.54757,
   0.68665,0.86965,1.11199,1.43242,1.86532,2.44703,3.22681,4.24972,
   5.54382,7.04016,8.48123,9.40627,9.40627,8.48123,7.04016,5.54382,
   4.24972,3.22681,2.44703,1.86532,1.43242,1.11199,0.86965,0.68665,
   0.54757,0.44012,0.35713,0.29201,0.24129,0.20074,0.16860,0.14268,
   0.12185,0.10486,0.09093,0.07965,0.07015,0.06232,0.05565,0.05004,
   0.04533,0.04116,0.03764,0.03453,0.03179,0.02939,0.02723,0.02525,
   0.02356,0.02196,0.02049,0.01912,0.01793,0.01671,0.01569,0.01467,
   0.01376,0.01286,0.01205,0.01126,0.01056,0.00985,0.00921,0.00859,
   0.00803,0.00747,0.00698,0.00648,0.00603,0.00560,0.00521,0.00482,
   0.00447,0.00412,0.00383,0.00352,0.00326,0.00299,0.00276,0.00252,
   0.00232,0.00211,0.00195,0.00176,0.00161,0.00146,0.00134,0.00119,
   0.00109,0.00097,0.00089,0.00078,0.00071,0.00062,0.00057,0.00049,
   0.00044,0.00038,0.00034,0.00029,0.00026,0.00021,0.00019,0.00015,
   0.00013,0.00010,0.00009,0.00006,0.00005,0.00003,0.00002,0.00001};
  
  static constexpr std::array<double, 200> avp2f
  {1.0,0.98974,0.96975,0.96768,0.94782,0.94450,0.92494,0.92047,
   0.90090,0.89563,0.87644,0.87018,0.85145,0.84434,0.82593,0.81841,
   0.80021,0.79212,0.77444,0.76553,0.74866,0.73945,0.72264,0.71343,
   0.69703,0.68740,0.67149,0.66182,0.64631,0.63630,0.62125,0.61154,
   0.59671,0.58686,0.57241,0.56283,0.54866,0.53889,0.52528,0.51581,
   0.50236,0.49291,0.47997,0.47063,0.45803,0.44867,0.43665,0.42744,
   0.41554,0.40656,0.39511,0.38589,0.37488,0.36611,0.35516,0.34647,
   0.33571,0.32704,0.31656,0.30783,0.29741,0.28870,0.27820,0.26945,
   0.25898,0.25010,0.23945,0.23023,0.21943,0.20999,0.19891,0.18911,
   0.17795,0.16793,0.15669,0.14667,0.13553,0.12569,0.11504,0.10550,
   0.09557,0.08674,0.07774,0.06974,0.06184,0.05484,0.04802,0.04203,
   0.03629,0.03129,0.02654,0.02247,0.01867,0.01545,0.01251,0.01015,
   0.00810,0.00664,0.00541,0.00512,0.00512,0.00541,0.00664,0.00810,
   0.01015,0.01251,0.01545,0.01867,0.02247,0.02654,0.03129,0.03629,
   0.04203,0.04802,0.05484,0.06184,0.06974,0.07774,0.08674,0.09557,
   0.10550,0.11504,0.12569,0.13553,0.14667,0.15669,0.16793,0.17795,
   0.18911,0.19891,0.20999,0.21943,0.23023,0.23945,0.25010,0.25898,
   0.26945,0.27820,0.28870,0.29741,0.30783,0.31656,0.32704,0.33571,
   0.34647,0.35516,0.36611,0.37488,0.38589,0.39511,0.40656,0.41554,
   0.42744,0.43665,0.44867,0.45803,0.47063,0.47997,0.49291,0.50236,
   0.51581,0.52528,0.53889,0.54866,0.56283,0.57241,0.58686,0.59671,
   0.61154,0.62125,0.63630,0.64631,0.66182,0.67149,0.68740,0.69703,
   0.71343,0.72264,0.73945,0.74866,0.76553,0.77444,0.79212,0.80021,
   0.81841,0.82593,0.84434,0.85145,0.87018,0.87644,0.89563,0.90090,
   0.92047,0.92494,0.94450,0.94782,0.96768,0.96975,0.98974,1.0};
  
  double W2=W*W;
  double MN = fAM;
  double MN2 = MN*MN;
  double MD = fMD;
  double Mp = fMP;
  double Mp2 = Mp*Mp;
  double nu = (W2 - MN2 + Q2)/2./MN;
  double qv = TMath::Sqrt(nu*nu + Q2);
  F1 = 0.;
  R = 0.;
  sigmaT = 0.;
  sigmaL = 0.;
  // Do fast 20 bins if abvoe threshold
  if(W2>1.30)
  {
    for (int ism = 0; ism<20; ism++)
    {
      double W2p = TMath::Power(MD + nu - sqrt(MN2 + avp2[ism]),2) - qv*qv + 2.*qv*avpz[ism] - avp2[ism];
      if(W2p>1.155)
      {
         double Wp = TMath::Sqrt(W2p);
         double sigtp = sigmaR(0, Q2, Wp, isDeuterium) + sigmaNR(0, Q2, Wp, isDeuterium);
         double F1p = sigtp*(W2p-Mp2)/8./kPi2/kAem;
         F1 += F1p*fyd[ism]/10.;
         if (!isDeuterium)
         {
            double siglp = sigmaR(1, Q2, Wp) + sigmaNR(1, Q2, Wp); 
            sigmaL += siglp*fyd[ism]/10.;
            sigmaT += sigtp*fyd[ism]/10.;
         }
      }
    }
  }
  else
  {
     for (int ism = 0;ism<200;ism++)
     {
       double pz = -1. + 0.01*(ism + 0.5);
       // Need avp2f term to get right behavior x>1! 
       double W2p = TMath::Power(MD + nu - sqrt(MN2 + avp2f[ism]),2) - qv*qv + 2.*qv*pz - avp2f[ism];
       if(W2p>1.155)
       {
          double Wp = TMath::Sqrt(W2p);
          double sigtp = sigmaR(0, Q2, Wp, isDeuterium) + sigmaNR(0, Q2, Wp, isDeuterium);
          double F1p = sigtp*(W2p-Mp2)/8./kPi2/kAem;
          F1 += F1p*fydf[ism]/100.;
          if (!isDeuterium)
          {
            double siglp = sigmaR(1, Q2, Wp) + sigmaNR(1, Q2, Wp); 
            sigmaT += sigtp*fydf[ism]/100.;
            sigmaL += siglp*fydf[ism]/100.;
          }
       }
     }
  }
  if (isDeuterium && fUseMEC)
     // Ref.2, Eq. (20)
     F1 += fMECcoef[0]*TMath::Exp(-(W - fMECcoef[1])*(W - fMECcoef[1])/fMECcoef[2])/
             TMath::Power(1. + TMath::Max(0.3,Q2)/fMECcoef[3],fMECcoef[4])*TMath::Power(nu, fMECcoef[5]);
  if(!isDeuterium && sigmaT!=0.) 
    R = sigmaL/sigmaT;
  
}
//____________________________________________________________________________
double BostedChristyEMPXSec::MEC2009(int A, double Q2, double W) const
{
  double F1 = 0.0;
  double W2 = W*W;
  double Mp = fAM;
  double Mp2 = Mp*Mp;
  if(W2<=0.0)
    return F1;
  double nu = (W2 - Mp2 + Q2)/2./Mp;
  double x  = Q2/2.0/Mp/nu;
  
  if(A<=2) 
    return F1;
    
  double p18;
  for (const auto& kv : fMEC2009p18) 
  {
    p18 = kv.second;
    if (A<=kv.first)
      break;
  }
  
  F1 = fMEC2009coef[0]*TMath::Exp(-(W - fMEC2009coef[1])*(W - fMEC2009coef[1])/fMEC2009coef[2])/
             TMath::Power(1. + TMath::Max(0.3,Q2)/fMEC2009coef[3],fMEC2009coef[4])*TMath::Power(nu, fMEC2009coef[5])*(1.0 + 
                                                                     p18*TMath::Power(A, fMEC2009coef[6] + fMEC2009coef[7]*x));
                                                                     
  if(F1<=1.0E-9) 
    F1 = 0.0;
    
  return F1;
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
  const Target & target = init_state.Tgt();
  int A = target.A();
  int Z = target.Z();
  double E  = init_state.ProbeE(kRfHitNucRest);
  double W  = kinematics.W();
  double Q2 = kinematics.Q2();
  double Wsq = W*W;
  // Cross section for proton or neutron
  
  double Mp = fMP;
  double Mp2 = Mp*Mp;
  double MN = fPM;
  double MN2 = MN*MN;
  
  double nu = (Wsq - MN2 + Q2)/2./MN;
  double x  = Q2/2./MN/nu;

  double sigmaT, sigmaL, F1p, R, W1;
  // Cross section for proton or neutron
  if (A<2 && Wsq>1.155)
  {
     double xb = Q2/(Wsq+Q2-Mp2);
     sigmaT = sigmaR(0, Q2, W) + sigmaNR(0, Q2, W);
     sigmaL = sigmaR(1, Q2, W) + sigmaNR(1, Q2, W); 
     F1p = sigmaT*(Wsq-Mp2)/8./kPi2/kAem;
     R = sigmaL/sigmaT;
     // If neutron, subtract proton from deuteron. Factor of two to
     // convert from per nucleon to per deuteron
     if(Z==0)
     {
          sigmaT = sigmaR(0, Q2, W, true) + sigmaNR(0, Q2, W, true);
          double F1d =  sigmaT*(Wsq-Mp2)/8./kPi2/kAem;
          F1p = 2.*F1d - F1p;
     }
     W1 = F1p/MN;
  }
  
  // For deuteron
  if(A==2)
  {
     double Rd, F1c, F1d;
     //get Fermi-smeared R from Erics proton fit
     FermiSmearingD(Q2, W, F1c, R, sigmaT, sigmaL);
     //get fit to F1 in deuteron, per nucleon
     FermiSmearingD(Q2, W, F1d, Rd, sigmaT, sigmaL, true);
     //convert to W1 per deuteron
     W1 = F1d/MN*2.0;
  }
  
  //For nuclei
  if (A>2)
  {
    // Modifed to use Superscaling from Ref. 3
    double Es, pF, kF;
    for (const auto& kv : fNucRmvE) 
    {
      Es = kv.second;
      if (A<=kv.first)
        break;
    }
    for (const auto& kv : fKFTable) 
    {
      kF = kv.second;
      if (A<=kv.first)
        break;
    }
    // adjust pf to give right width based on kf
    pF = 0.5*kF;
    double F1d;
    FermiSmearingA(Q2, W, pF, Es, F1p, F1d, sigmaT, sigmaL);
    R = 0.;
    if(sigmaT>0.) 
      R = sigmaL/sigmaT;
    W1 = (2.*Z*F1d + (A - 2.*Z)*(2.*F1d - F1p))/MN; 

    W1 *= (fAfitcoef[0] + x*(fAfitcoef[1] + x*(fAfitcoef[2] + x*(fAfitcoef[3] + x*(fAfitcoef[4] + x*fAfitcoef[5])))));

    if(W>0.) 
      W1 *= TMath::Power(fAfitcoef[6] + (fAfitcoef[7]*W + fAfitcoef[8]*Wsq)/(fAfitcoef[9] + fAfitcoef[10]*Q2),2);
      
    double F1M = MEC2009(A, Q2, W);

    W1 += F1M;
    if(Wsq>0.) 
      R *= (fAfitcoef[11] + fAfitcoef[12]*A);
  }
  
  double emcfac = FitEMC(x, A);

  W1 *= emcfac;
    
  double nu2 = nu*nu;
  double Eprime = E - nu;
  double sin2theta_2 = Q2/4/E/Eprime;
  double cos2theta_2 = 1 - sin2theta_2;
  double W2 = W1*(1 + R)/ (1+nu2/Q2);
  double xsec = 4*Eprime*Eprime*kAem2/Q2/Q2*(2*W1*sin2theta_2 + W2*cos2theta_2);   // d2xsec/dOmegadEprime
  double jacobian = W*kPi/E/Eprime/MN;  
  xsec*= jacobian;                                                                 // d2xsec/dOmegadEprime-> d2xsec/dWdQ2
  
  // The algorithm computes d^2xsec/dWdQ2
  // Check whether variable tranformation is needed
  if(kps!=kPSWQ2fE) {
    double J = utils::kinematics::Jacobian(interaction,kPSWQ2fE,kps);
    xsec *= J;
  }
  
  return xsec;                              
  
}
//____________________________________________________________________________
// Fit to EMC effect.  Steve Rock 8/3/94                                 
// Funciton returns value of sigma(A)/sigma(d) 
// with no isoscalerity correction
// A = atomic number                                                      
// x = Bjorken x.                                                        
//                                                                       
// Fit of sigma(A)/sigma(d) to form C*A**alpha where A is atomic number  
// First data at each x was fit to form C*A**alpha.  The point A=2(d)    
//  was includded with a value of 1 and an error of about 2%.             
// For x>=.125 Javier Gomez fit of 7/93 to E139 data was used.           
// For .09 >=x>=.0085 NMC data from Amaudruz et al Z. Phys C. 51,387(91) 
//  Steve did the fit for alpha and C to the He/d. C/d and Ca/d NMC data.
// Alpha(x) was fit to a 9 term polynomial a0 +a1*x +a2*x**2 + a3*x**3 ..
// C(x) was fit to a 3 term polynomial in natural logs as                
//  Ln(C) = c0 + c1*Ln(x) + c2*[Ln(x)]**2.                               
//
// 6/2/98 *****  Bug (which set x= .00885 if x was out of range) fixed
//                    also gave value at x=.0085 if x>.88
// 11/05 PYB PEB modified to use value at x=0.7 if x>0.7, because
//    beyond that assume Fermi motion is giving the rise, and we
//    already are taking that into account with the y-smearing of
//    the inelastic
//____________________________________________________________________________
double BostedChristyEMPXSec::FitEMC(double x, int A) const
{                                                                     
  double fitemc = 1.;
  if(A<=2) 
    return fitemc;
  
  double x_u;
  if (x>0.70 || x<0.0085)   
  //Out of range of fit   
  {
   if(x<0.0085) 
     x_u = .0085;
   if(x>0.70) 
     x_u = 0.70;
  }
  else
   x_u = x;
  
  double ln_c = fEMCc[0];               
  for (int i = 1; i<=2; i++)
     ln_c += fEMCc[i]*TMath::Power(TMath::Log(x_u), i); 
  double c = TMath::Exp(ln_c);                  
  
  double alpha = fEMCalpha[0];          
  for (int i = 1; i<=8; i++)                                 
   alpha += fEMCalpha[i]*TMath::Power(x_u, i);
                           
  fitemc = c*TMath::Power(A, alpha); 
  return fitemc;
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

  if(!proc_info.IsResonant() ) return false;
 

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
bool BostedChristyEMPXSec::ValidKinematics(const Interaction * interaction) const
{
   const Kinematics & kinematics = interaction -> Kine();
   double W  = kinematics.W();
   double Q2 = kinematics.Q2();
   if (W<fWmin || W>fWmax)
     return false;
   if (Q2<fQ2min || Q2>fQ2max)
     return false;
     
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
void BostedChristyEMPXSec::BranchingRatios(int respdg, double & brpi, double & breta) const
{
  brpi  = 0.;
  breta = 0.;
  PDGLibrary * pdglib = PDGLibrary::Instance();
  TParticlePDG * res_pdg = pdglib->Find(respdg);
  if (res_pdg != 0)
  {
    for (int nch = 0; nch < res_pdg->NDecayChannels(); nch++)
    {
      TDecayChannel * ch = res_pdg->DecayChannel(nch);
      if (ch->NDaughters() == 2)
      {
        int first_daughter_pdg  = ch->DaughterPdgCode (0);
        int second_daughter_pdg = ch->DaughterPdgCode (1);
        if ((genie::pdg::IsNucleon(first_daughter_pdg ) && genie::pdg::IsPion(second_daughter_pdg)) ||
        (genie::pdg::IsNucleon(second_daughter_pdg) && genie::pdg::IsPion(first_daughter_pdg )))
           brpi += ch->BranchingRatio();
        if (first_daughter_pdg == kPdgEta || second_daughter_pdg == kPdgEta)    
          breta += ch->BranchingRatio();
      }
    }
  }
}
//____________________________________________________________________________
void BostedChristyEMPXSec::LoadConfig(void)
{  
  
  PDGLibrary * pdglib = PDGLibrary::Instance();
  GetParamDef("BostedChristyFitEM-PM",     fPM,     pdglib->Find(kPdgProton)->Mass());
  GetParamDef("BostedChristyFitEM-MP",     fMP,     pdglib->Find(kPdgProton)->Mass());
  GetParamDef("BostedChristyFitEM-AM",     fAM,     pdglib->Find(kPdgProton)->Mass());
  GetParamDef("BostedChristyFitEM-MD",     fMD,     pdglib->Find(kPdgTgtDeuterium)->Mass());
  GetParamDef("BostedChristyFitEM-Mpi0",   fMpi0,   pdglib->Find(kPdgPi0)->Mass());
  GetParamDef("BostedChristyFitEM-Meta",   fMeta,   pdglib->Find(kPdgEta)->Mass());
  GetParamDef("BostedChristyFitEM-Wmin",   fWmin,   0.0);
  GetParamDef("BostedChristyFitEM-Wmax",   fWmax,   3.0);
  GetParamDef("BostedChristyFitEM-Q2min",  fQ2min,  0.0);
  GetParamDef("BostedChristyFitEM-Q2max",  fQ2max,  10.0);
  GetParamDef("BostedChristyFitEM-UseMEC", fUseMEC, true);
  
  double BRpi, BReta;
  double brpi, breta;
    
  std::vector<double> vBRpi1;
  std::vector<double> vBRpi2;
  std::vector<double> vBReta;
  
  // load braching ratios for pi
  bool useBRpi1Default = (GetParamVect("BostedChristyFitEM-PionBRp", vBRpi1, false)<7);
  bool useBRpi2Default = (GetParamVect("BostedChristyFitEM-PionBRD", vBRpi2, false)<7);
  // load braching ratios for eta
  bool useBRetaDefault = (GetParamVect("BostedChristyFitEM-EtaBR", vBReta, false)<7);
  
  if (useBRpi1Default || useBRpi2Default || useBRetaDefault)
  {
     // use default branching ratios from PDG table
     //  P33(1232)
     BRpi  = 0.;
     BReta = 0.;
     BranchingRatios(kPdgP33m1232_DeltaM, brpi, breta);
     BRpi  += brpi;
     BReta += breta;
     BranchingRatios(kPdgP33m1232_Delta0, brpi, breta);
     BRpi  += brpi;
     BReta += breta;
     BranchingRatios(kPdgP33m1232_DeltaP, brpi, breta);
     BRpi  += brpi;
     BReta += breta;
     BranchingRatios(kPdgP33m1232_DeltaPP, brpi, breta);
     BRpi  += brpi;
     BReta += breta;
     BRpi  /= 4.;
     BReta /= 4.;
     fBRp[0][0] = BRpi;
     fBRp[0][2] = BReta;
     fBRD[0][0] = BRpi;
     
     //  S11(1535)
     BRpi  = 0.;
     BReta = 0.;
     BranchingRatios(kPdgS11m1535_N0, brpi, breta);
     BRpi  += brpi;
     BReta += breta;
     BranchingRatios(kPdgS11m1535_NP, brpi, breta);
     BRpi  += brpi;
     BReta += breta;
     BRpi  /= 2.;
     BReta /= 2.;
     fBRp[1][0] = BRpi;
     fBRp[1][2] = BReta;
     fBRD[1][0] = BRpi;
     
     //  D13(1520)
     BRpi  = 0.;
     BReta = 0.;
     BranchingRatios(kPdgD13m1520_N0, brpi, breta);
     BRpi  += brpi;
     BReta += breta;
     BranchingRatios(kPdgD13m1520_NP, brpi, breta);
     BRpi  += brpi;
     BReta += breta;
     BRpi  /= 2.;
     BReta /= 2.;
     fBRp[2][0] = BRpi;
     fBRp[2][2] = BReta;
     fBRD[2][0] = BRpi;
     
     //  F15(1680)
     BRpi  = 0.;
     BReta = 0.;
     BranchingRatios(kPdgF15m1680_N0, brpi, breta);
     BRpi  += brpi;
     BReta += breta;
     BranchingRatios(kPdgF15m1680_NP, brpi, breta);
     BRpi  += brpi;
     BReta += breta;
     BRpi  /= 2.;
     BReta /= 2.;
     fBRp[3][0] = BRpi;
     fBRp[3][2] = BReta;
     fBRD[3][0] = BRpi;
     
     //  S11(1650)
     BRpi  = 0.;
     BReta = 0.;
     BranchingRatios(kPdgS11m1650_N0, brpi, breta);
     BRpi  += brpi;
     BReta += breta;
     BranchingRatios(kPdgS11m1650_NP, brpi, breta);
     BRpi  += brpi;
     BReta += breta;
     BRpi  /= 2.;
     BReta /= 2.;
     fBRp[4][0] = BRpi;
     fBRp[4][2] = BReta;
     fBRD[4][0] = BRpi;
     
     //  P11(1440)
     BRpi  = 0.;
     BReta = 0.;
     BranchingRatios(kPdgP11m1440_N0, brpi, breta);
     BRpi  += brpi;
     BReta += breta;
     BranchingRatios(kPdgP11m1440_NP, brpi, breta);
     BRpi  += brpi;
     BReta += breta;
     BRpi  /= 2.;
     BReta /= 2.;
     fBRp[5][0] = BRpi;
     fBRp[5][2] = BReta;
     fBRD[5][0] = BRpi;
     
     //  F37(1950)
     BRpi  = 0.;
     BReta = 0.;
     BranchingRatios(kPdgF37m1950_DeltaM , brpi, breta);
     BRpi  += brpi;
     BReta += breta;
     BranchingRatios(kPdgF37m1950_Delta0, brpi, breta);
     BRpi  += brpi;
     BReta += breta;
     BranchingRatios(kPdgF37m1950_DeltaP, brpi, breta);
     BRpi  += brpi;
     BReta += breta;
     BranchingRatios(kPdgF37m1950_DeltaPP, brpi, breta);
     BRpi  += brpi;
     BReta += breta;
     BRpi  /= 4.;
     BReta /= 4.;
     fBRp[6][0] = BRpi;
     fBRp[6][2] = BReta;
     fBRD[6][0] = BRpi;
  }
  if (!useBRpi1Default)
    // single pion branching ratios from config file
    for (int i=0; i<7; i++)
      fBRp[i][0] =  vBRpi1[i];
  if (!useBRpi2Default)
    // single pion branching ratios from config file
    for (int i=0; i<7; i++)
      fBRD[i][0] =  vBRpi2[i];
  if (!useBRetaDefault) 
    //  eta branching ratios from config file
    for (int i=0; i<7; i++)
      fBRp[i][2] =  vBReta[i];
      
  for (int i=0; i<7; i++)
    fBRD[i][2] =  0.;
  
  if (useBRpi1Default || useBRpi2Default)
     LOG("BostedChristyEMPXSec", pALERT)  << "*** Use branching ratios for pion decay from PDG table";
     
  if (useBRetaDefault)
     LOG("BostedChristyEMPXSec", pALERT)  << "*** Use branching ratios for eta decay from PDG table";

  //  2-pion branching ratios 
  for (int i=0;i<7;i++)
  {
    fBRp[i][1] = 1.-fBRp[i][0]-fBRp[i][2];
    fBRD[i][1] = 1.-fBRD[i][0]-fBRD[i][2];
  }

  // Meson angular momentum
  fAngRes[0] = res::OrbitalAngularMom(res::FromPdgCode(kPdgP33m1232_Delta0));      //  P33(1232)
  fAngRes[1] = res::OrbitalAngularMom(res::FromPdgCode(kPdgS11m1535_N0));          //  S11(1535)
  fAngRes[2] = res::OrbitalAngularMom(res::FromPdgCode(kPdgD13m1520_N0));          //  D13(1520)
  fAngRes[3] = res::OrbitalAngularMom(res::FromPdgCode(kPdgF15m1680_N0));          //  F15(1680)
  fAngRes[4] = res::OrbitalAngularMom(res::FromPdgCode(kPdgS11m1650_N0));          //  S15(1650)
  fAngRes[5] = res::OrbitalAngularMom(res::FromPdgCode(kPdgP11m1440_N0));          //  P11(1440) roper   
  fAngRes[6] = res::OrbitalAngularMom(res::FromPdgCode(kPdgF37m1950_Delta0));      //  F37(1950)
  
  std::vector<double> vResMass;
  // load resonance masses
  bool useResMassDefault = (GetParamVect("BostedChristyFitEM-ResMass", vResMass, false)<7);
  
  if (useResMassDefault)
  {
    LOG("BostedChristyEMPXSec", pALERT)  << "*** Use resonance masses from PDG table";
    // Resonance mass
    fMassRes[0] = res::Mass(res::FromPdgCode(kPdgP33m1232_Delta0));    //  P33(1232)
    fMassRes[1] = res::Mass(res::FromPdgCode(kPdgS11m1535_N0));        //  S11(1535)
    fMassRes[2] = res::Mass(res::FromPdgCode(kPdgD13m1520_N0));        //  D13(1520)
    fMassRes[3] = res::Mass(res::FromPdgCode(kPdgF15m1680_N0));        //  F15(1680)
    fMassRes[4] = res::Mass(res::FromPdgCode(kPdgS11m1650_N0));        //  S15(1650)
    fMassRes[5] = res::Mass(res::FromPdgCode(kPdgP11m1440_N0));        //  P11(1440) roper   
    fMassRes[6] = res::Mass(res::FromPdgCode(kPdgF37m1950_Delta0));    //  F37(1950)
  }
  else
    //  eta branching ratios from config file
    for (int i=0; i<7; i++)
      fMassRes[i] =  vResMass[i];
        
  std::vector<double> vResWidth;
  // load resonance masses
  bool useResWidthDefault = (GetParamVect("BostedChristyFitEM-ResWidth", vResWidth, false)<7);
  
  if (useResWidthDefault)
  {
    LOG("BostedChristyEMPXSec", pALERT)  << "*** Use resonance widths from PDG table";
    // Resonance width
    fWidthRes[0] = res::Width(res::FromPdgCode(kPdgP33m1232_Delta0));    //  P33(1232)
    fWidthRes[1] = res::Width(res::FromPdgCode(kPdgS11m1535_N0));        //  S11(1535)
    fWidthRes[2] = res::Width(res::FromPdgCode(kPdgD13m1520_N0));        //  D13(1520)
    fWidthRes[3] = res::Width(res::FromPdgCode(kPdgF15m1680_N0));        //  F15(1680)
    fWidthRes[4] = res::Width(res::FromPdgCode(kPdgS11m1650_N0));        //  S15(1650)
    fWidthRes[5] = res::Width(res::FromPdgCode(kPdgP11m1440_N0));        //  P11(1440) roper   
    fWidthRes[6] = res::Width(res::FromPdgCode(kPdgF37m1950_Delta0));    //  F37(1950)
  }
  else
    //  eta branching ratios from config file
    for (int i=0; i<7; i++)
      fWidthRes[i] =  vResWidth[i];

  int length;
  
  std::vector<double> vRescoef;
  length = 7;
  bool isOk = (GetParamVect("BostedChristyFitEM-ResAT0p", vRescoef)>=length);
  if (!isOk)
  {
    LOG("BostedChristyEMPXSec", pFATAL)  << "*** Can't find enough proton AT(0)-parameters for xsec^R_T in the config file!";
    exit(1);
  }
  // Ref.1, Table III, AT(0)
  for (int i=0;i<length;i++)
    fRescoefTp[i][0] = vRescoef[i];
    
  length = 7;
  isOk = (GetParamVect("BostedChristyFitEM-Resap", vRescoef)>=length);
  if (!isOk)
  {
    LOG("BostedChristyEMPXSec", pFATAL)  << "*** Can't find enough proton a-parameters for xsec^R_T  in the config file!";
    exit(1);
  }
  // Ref.1, Table III, a
  for (int i=0;i<length;i++)
    fRescoefTp[i][1] = vRescoef[i];
    
  length = 7;
  isOk = (GetParamVect("BostedChristyFitEM-Resbp", vRescoef)>=length);
  if (!isOk)
  {
    LOG("BostedChristyEMPXSec", pFATAL)  << "*** Can't find enough proton b-parameters parameters for xsec^R_T  in the config file!";
    exit(1);
  }
  // Ref.1, Table III, b
  for (int i=0;i<length;i++)
    fRescoefTp[i][2] = vRescoef[i];
    
  length = 7;
  isOk = (GetParamVect("BostedChristyFitEM-Rescp", vRescoef)>=length);
  if (!isOk)
  {
    LOG("BostedChristyEMPXSec", pFATAL)  << "*** Can't find enough proton c-parameters parameters for xsec^R_T  in the config file!";
    exit(1);
  }
  // Ref.1, Table III, c
  for (int i=0;i<length;i++)
    fRescoefTp[i][3] = vRescoef[i];

  length = 7;
  isOk = (GetParamVect("BostedChristyFitEM-ResAT0D", vRescoef)>=length);
  if (!isOk)
  {
    LOG("BostedChristyEMPXSec", pFATAL)  << "*** Can't find enough deuterium AT(0)-parameters for xsec^R_T in the config file!";
    exit(1);
  }
  // Ref.2, Table III, AT(0)
  for (int i=0;i<length;i++)
    fRescoefTD[i][0] = vRescoef[i];
    
  length = 7;
  isOk = (GetParamVect("BostedChristyFitEM-ResaD", vRescoef)>=length);
  if (!isOk)
  {
    LOG("BostedChristyEMPXSec", pFATAL)  << "*** Can't find enough deuterium a-parameters for xsec^R_T  in the config file!";
    exit(1);
  }
  // Ref.2, Table III, a
  for (int i=0;i<length;i++)
    fRescoefTD[i][1] = vRescoef[i];
    
  length = 7;
  isOk = (GetParamVect("BostedChristyFitEM-ResbD", vRescoef)>=length);
  if (!isOk)
  {
    LOG("BostedChristyEMPXSec", pFATAL)  << "*** Can't find enough deuterium b-parameters parameters for xsec^R_T  in the config file!";
    exit(1);
  }
  // Ref.2, Table III, b
  for (int i=0;i<length;i++)
    fRescoefTD[i][2] = vRescoef[i];
    
  length = 7;
  isOk = (GetParamVect("BostedChristyFitEM-RescD", vRescoef)>=length);
  if (!isOk)
  {
    LOG("BostedChristyEMPXSec", pFATAL)  << "*** Can't find enough deuterium c-parameters parameters for xsec^R_T  in the config file!";
    exit(1);
  }
  // Ref.2, Table III, c
  for (int i=0;i<length;i++)
    fRescoefTD[i][3] = vRescoef[i];
    
  length = 7;
  isOk = (GetParamVect("BostedChristyFitEM-ResAL0", vRescoef)>=length);
  if (!isOk)
  {
    LOG("BostedChristyEMPXSec", pFATAL)  << "*** Can't find enough proton AL0-parameters parameters for xsec^R_T  in the config file!";
    exit(1);
  }
  // Ref.1, Table III, AL(0)
  for (int i=0;i<length;i++)
    fRescoefL[i][0] = vRescoef[i];
      
  length = 7;
  isOk = (GetParamVect("BostedChristyFitEM-Resd", vRescoef)>=length);
  if (!isOk)
  {
    LOG("BostedChristyEMPXSec", pFATAL)  << "*** Can't find enough proton d-parameters parameters for xsec^R_L  in the config file!";
    exit(1);
  }
  // Ref.1, Table III, d
  for (int i=0;i<length;i++)
    fRescoefL[i][1] = vRescoef[i];
 
  length = 7; 
  isOk = (GetParamVect("BostedChristyFitEM-Rese", vRescoef)>=length);
  if (!isOk)
  {
    LOG("BostedChristyEMPXSec", pFATAL)  << "*** Can't find enough proton e-parameters parameters for xsec^R_L  in the config file!";
    exit(1);
  }
  // Ref.1, Table III, e
  for (int i=0;i<length;i++)
    fRescoefL[i][2] = vRescoef[i];
    
    
  std::vector<double> vNRcoef;
  length = 5; 
  isOk = (GetParamVect("BostedChristyFitEM-NRXSecT1p", vNRcoef)>=length);
  if (!isOk)
  {
    LOG("BostedChristyEMPXSec", pFATAL)  << "*** Can't find enough proton bkg parameters for xsec^NR_T in the config file!";
    exit(1);
  }
  // Ref.1, Table IV: \sigma^NR,1_T(0), aT_1, bT_1, cT_1, dT_1
  for (int i=0;i<length;i++)
    fNRcoefTp[0][i] = vNRcoef[i];
    
  length = 5; 
  isOk = (GetParamVect("BostedChristyFitEM-NRXSecT2p", vNRcoef)>=length);
  if (!isOk)
  {
    LOG("BostedChristyEMPXSec", pFATAL)  << "*** Can't find enough proton bkg parameters for xsec^NR_T in the config file!";
    exit(1);
  }
  // Ref.1, Table IV: \sigma^NR,2_T(0), aT_2, bT_2, cT_2, dT_2
  for (int i=0;i<length;i++)
    fNRcoefTp[1][i] = vNRcoef[i];
    
  length = 5; 
  isOk = (GetParamVect("BostedChristyFitEM-NRXSecT1D", vNRcoef)>=length);
  if (!isOk)
  {
    LOG("BostedChristyEMPXSec", pFATAL)  << "*** Can't find enough deuterium bkg parameters for xsec^NR_T in the config file!";
    exit(1);
  }
  // Ref.2, Table IV: \sigma^NR,1_T(0), aT_1, bT_1, cT_1, dT_1
  for (int i=0;i<length;i++)
    fNRcoefTD[0][i] = vNRcoef[i];
    
  length = 5; 
  isOk = (GetParamVect("BostedChristyFitEM-NRXSecT2D", vNRcoef)>=length);
  if (!isOk)
  {
    LOG("BostedChristyEMPXSec", pFATAL)  << "*** Can't find enough deuterium bkg parameters for xsec^NR_T in the config file!";
    exit(1);
  }
  // Ref.2, Table IV: \sigma^NR,2_T(0), aT_2, bT_2, cT_2, dT_2
  for (int i=0;i<length;i++)
    fNRcoefTD[1][i] = vNRcoef[i];
    
  length = 6; 
  isOk = (GetParamVect("BostedChristyFitEM-NRXSecL", vNRcoef)>=length);
  if (!isOk)
  {
    LOG("BostedChristyEMPXSec", pFATAL)  << "*** Can't find enough proton bkg parameters for xsec^NR_L in the config file!";
    exit(1);
  }
  // Ref.1, Table IV: \sigma^NR_L, aL, bL, cL, dL, eL
  for (int i=0;i<length;i++)
    fNRcoefL[i] = vNRcoef[i];
    
  length = 6; 
  isOk = (GetParamVect("BostedChristyFitEM-MEC", vNRcoef)>=length);
  if (!isOk)
  {
    LOG("BostedChristyEMPXSec", pFATAL)  << "*** Can't find enough parameters for MEC in the config file!";
    exit(1);
  }
  for (int i=0;i<length;i++)
    fMECcoef[i] = vNRcoef[i];
       
  length = 8; 
  isOk = (GetParamVect("BostedChristyFitEM-MEC2009", vNRcoef)>=length);
  if (!isOk)
  {
    LOG("BostedChristyEMPXSec", pFATAL)  << "*** Can't find enough parameters for MEC2009 in the config file!";
    exit(1);
  }
  for (int i=0;i<length;i++)
    fMEC2009coef[i] = vNRcoef[i];
    
  length = 13; 
  isOk = (GetParamVect("BostedChristyFitEM-Afit", vNRcoef)>=length);
  if (!isOk)
  {
    LOG("BostedChristyEMPXSec", pFATAL)  << "*** Can't find enough parameters for nuclei fit (A-fit) in the config file!";
    exit(1);
  }
  for (int i=0;i<length;i++)
    fAfitcoef[i] = vNRcoef[i]; 
    
    
  length = 9; 
  isOk = (GetParamVect("BostedChristyFitEM-EMCalpha", vNRcoef)>=length);
  if (!isOk)
  {
    LOG("BostedChristyEMPXSec", pFATAL)  << "*** Can't find enough alpha coefficients for EMC correction in the config file!";
    exit(1);
  }
  for (int i=0;i<length;i++)
    fEMCalpha[i] = vNRcoef[i];
    
  length = 3; 
  isOk = (GetParamVect("BostedChristyFitEM-EMCc", vNRcoef)>=length);
  if (!isOk)
  {
    LOG("BostedChristyEMPXSec", pFATAL)  << "*** Can't find enough c coefficients for EMC correction in the config file!";
    exit(1);
  }
  for (int i=0;i<length;i++)
    fEMCc[i] = vNRcoef[i]; 
    
    
  std::string keyStart = "BostedChristy-SeparationE@Pdg=";
  RgIMap entries = GetConfig().GetItemMap();
  for(RgIMap::const_iterator it = entries.begin(); it != entries.end(); ++it)
  {
    const std::string& key = it->first;
    int pdg = 0;
    int A = 0;
    if (0 == key.compare(0, keyStart.size(), keyStart.c_str()))
    {
      pdg = atoi(key.c_str() + keyStart.size());
      A = pdg::IonPdgCodeToA(pdg);
    }
    if (0 != pdg && A != 0) 
    {
      std::ostringstream key_ss ;
      key_ss << keyStart << pdg;
      RgKey rgkey   = key_ss.str();
      double eb;
      GetParam( rgkey, eb) ;
      eb = TMath::Max(eb, 0.);
      fNucRmvE.insert(map<int,double>::value_type(A,eb));
    }
  }
  
  keyStart = "BostedChristy-FermiMomentum@Pdg=";
  for(RgIMap::const_iterator it = entries.begin(); it != entries.end(); ++it)
  {
    const std::string& key = it->first;
    int pdg = 0;
    int A = 0;
    if (0 == key.compare(0, keyStart.size(), keyStart.c_str()))
    {
      pdg = atoi(key.c_str() + keyStart.size());
      A = pdg::IonPdgCodeToA(pdg);
    }
    if (0 != pdg && A != 0) 
    {
      std::ostringstream key_ss ;
      key_ss << keyStart << pdg;
      RgKey rgkey   = key_ss.str();
      double pf;
      GetParam( rgkey, pf) ;
      pf = TMath::Max(pf, 0.);
      fKFTable.insert(map<int,double>::value_type(A,pf));
    }
  }
  
  keyStart = "BostedChristy-p18@Pdg=";
  for(RgIMap::const_iterator it = entries.begin(); it != entries.end(); ++it)
  {
    const std::string& key = it->first;
    int pdg = 0;
    int A = 0;
    if (0 == key.compare(0, keyStart.size(), keyStart.c_str()))
    {
      pdg = atoi(key.c_str() + keyStart.size());
      A = pdg::IonPdgCodeToA(pdg);
    }
    if (0 != pdg && A != 0) 
    {
      std::ostringstream key_ss ;
      key_ss << keyStart << pdg;
      RgKey rgkey   = key_ss.str();
      double p18;
      GetParam( rgkey, p18) ;
      fMEC2009p18.insert(map<int,double>::value_type(A,p18));
    }
  }
  
}
