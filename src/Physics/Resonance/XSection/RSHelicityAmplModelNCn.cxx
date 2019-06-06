//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Oct 05, 2009 - CA
   Compute() now returns a `const RSHelicityAmpl &' and avoids creating a new
   RSHelicityAmpl at each call.                      

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/Resonance/XSection/RSHelicityAmplModelNCn.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
RSHelicityAmplModelNCn::RSHelicityAmplModelNCn() :
RSHelicityAmplModelI("genie::RSHelicityAmplModelNCn")
{

}
//____________________________________________________________________________
RSHelicityAmplModelNCn::RSHelicityAmplModelNCn(string config) :
RSHelicityAmplModelI("genie::RSHelicityAmplModelNCn", config)
{

}
//____________________________________________________________________________
RSHelicityAmplModelNCn::~RSHelicityAmplModelNCn()
{

}
//____________________________________________________________________________
const RSHelicityAmpl & 
  RSHelicityAmplModelNCn::Compute(
      Resonance_t res, const FKR & fkr) const
{
  double xi = fSin28w;

  switch(res) {

   case (kP33_1232) :
   {
     double rx     = 2 * xi * fkr.R;
     double Rm2xiR = fkr.Rminus + rx;
     double Rp2xiR = fkr.Rplus  + rx;

     fAmpl.fMinus1 =  -kSqrt2 * Rm2xiR;
     fAmpl.fPlus1  =   kSqrt2 * Rp2xiR;
     fAmpl.fMinus3 =  -kSqrt6 * Rm2xiR;
     fAmpl.fPlus3  =   kSqrt6 * Rp2xiR;
     fAmpl.f0Minus = 2*kSqrt2 * fkr.C;
     fAmpl.f0Plus  =   fAmpl.f0Minus;
     break;
   }
   case (kS11_1535) :
   {
     double xt      = 2*xi*fkr.T;
     double xr      = xi*fkr.R;
     double Tm2xiT  = fkr.Tminus + xt;
     double Tp2xiT  = fkr.Tplus  + xt;
     double LRmxiR  = fkr.Lamda * (fkr.Rminus + xr);
     double LRpxiR  = fkr.Lamda * (fkr.Rplus  + xr);
     double a       = kSqrt3_2 * (1-2*xi) * fkr.Lamda * fkr.S;
     double b       = kSqrt2_3 * (fkr.Lamda * fkr.C - 3*fkr.B);

     fAmpl.fMinus1 = -1*kSqrt3 * Tm2xiT - kSqrt2_3 * LRmxiR;
     fAmpl.fPlus1  =    kSqrt3 * Tp2xiT + kSqrt2_3 * LRpxiR;
     fAmpl.fMinus3 =  0.;
     fAmpl.fPlus3  =  0.;
     fAmpl.f0Minus =  a-b;
     fAmpl.f0Plus  = -a-b;
     break;
   }
   case (kD13_1520) :
   {
     double xt      = 2* xi * fkr.T;
     double xr      =    xi * fkr.R;
     double Tm2xiT  = fkr.Tminus + xt;
     double Tp2xiT  = fkr.Tplus  + xt;
     double LRmxiR  = fkr.Lamda  * (fkr.Rminus + xr);
     double LRpxiR  = fkr.Lamda  * (fkr.Rplus  + xr);
     double a       = kSqrt3 * (1-2*xi) * fkr.Lamda * fkr.S;
     double b       = k2_Sqrt3 * fkr.Lamda * fkr.C;

     fAmpl.fMinus1 = -kSqrt3_2 * Tm2xiT + k2_Sqrt3 * LRmxiR;
     fAmpl.fPlus1  = -kSqrt3_2 * Tp2xiT + k2_Sqrt3 * LRpxiR;
     fAmpl.fMinus3 = -k3_Sqrt2 * Tm2xiT;
     fAmpl.fPlus3  = -k3_Sqrt2 * Tp2xiT;
     fAmpl.f0Minus =  a - b;
     fAmpl.f0Plus  =  a + b;
     break;
   }
   case (kS11_1650) :
   {
     double xr      = 4* xi * fkr.R;
     double LRm4xiR = fkr.Lamda * (fkr.Rminus + xr);
     double LRp4xiR = fkr.Lamda * (fkr.Rplus  + xr);

     fAmpl.fMinus1 = -k1_Sqrt24 * LRm4xiR;
     fAmpl.fPlus1  =  k1_Sqrt24 * LRp4xiR;
     fAmpl.fMinus3 =  0.;
     fAmpl.fPlus3  =  0.;
     fAmpl.f0Minus =  k1_Sqrt6 * (fkr.Lamda * fkr.C - 3*fkr.B);
     fAmpl.f0Plus  =  fAmpl.f0Minus;
     break;
   }
   case (kD13_1700) :
   {
     double xr      = 4* xi * fkr.R;
     double LRm4xiR = fkr.Lamda * (fkr.Rminus + xr);
     double LRp4xiR = fkr.Lamda * (fkr.Rplus  + xr);

     fAmpl.fMinus1 = -k1_Sqrt120 * LRm4xiR;
     fAmpl.fPlus1  = -k1_Sqrt120 * LRp4xiR;
     fAmpl.fMinus3 = -k3_Sqrt40  * LRm4xiR;
     fAmpl.fPlus3  = -k3_Sqrt40  * LRp4xiR;
     fAmpl.f0Minus = -k1_Sqrt30  * fkr.Lamda * fkr.C;
     fAmpl.f0Plus  =  -1.* fAmpl.f0Minus;
     break;
   }
   case (kD15_1675) :
   {
     double xr      = 4* xi * fkr.R;
     double LRm4xiR = fkr.Lamda * (fkr.Rminus + xr);
     double LRp4xiR = fkr.Lamda * (fkr.Rplus  + xr);

     fAmpl.fMinus1 =  kSqrt3_40 * LRm4xiR;
     fAmpl.fPlus1  = -kSqrt3_40 * LRp4xiR;
     fAmpl.fMinus3 =  kSqrt3_20 * LRm4xiR;
     fAmpl.fPlus3  = -kSqrt3_20 * LRp4xiR;
     fAmpl.f0Minus = -kSqrt3_10 * (fkr.Lamda * fkr.C);
     fAmpl.f0Plus  =  fAmpl.f0Minus;
     break;
   }
   case (kS31_1620) :
   {
     double xt      = 2*xi*fkr.T;
     double xr      = 2*xi*fkr.R;
     double Tm2xiT  = fkr.Tminus + xt;
     double Tp2xiT  = fkr.Tplus  + xt;
     double LRm2xiR = fkr.Lamda * (fkr.Rminus + xr);
     double LRp2xiR = fkr.Lamda * (fkr.Rplus  + xr);
     double a       = kSqrt3_2 * (1-2*xi) * fkr.Lamda * fkr.S;
     double b       = k1_Sqrt6 * (fkr.Lamda * fkr.C - 3*fkr.B);

     fAmpl.fMinus1 =  kSqrt3 * Tm2xiT - k1_Sqrt6 * LRm2xiR;
     fAmpl.fPlus1  = -kSqrt3 * Tp2xiT + k1_Sqrt6 * LRp2xiR;
     fAmpl.fMinus3 =  0.;
     fAmpl.fPlus3  =  0.;
     fAmpl.f0Minus = -a-b;
     fAmpl.f0Plus  =  a-b;
     break;
   }
   case (kD33_1700) :
   {
     double xt      = 2*xi*fkr.T;
     double xr      = 2*xi*fkr.R;
     double Tm2xiT  = fkr.Tminus + xt;
     double Tp2xiT  = fkr.Tplus  + xt;
     double LRm2xiR = fkr.Lamda * (fkr.Rminus + xr);
     double LRp2xiR = fkr.Lamda * (fkr.Rplus  + xr);
     double a       = kSqrt3 * (1-2*xi) * fkr.Lamda * fkr.S;
     double b       = k1_Sqrt3 * fkr.Lamda * fkr.C;

     fAmpl.fMinus1 = kSqrt3_2 * Tm2xiT + k1_Sqrt3 * LRm2xiR;
     fAmpl.fPlus1  = kSqrt3_2 * Tp2xiT + k1_Sqrt3 * LRp2xiR;
     fAmpl.fMinus3 = k3_Sqrt2 * Tm2xiT;
     fAmpl.fPlus3  = k3_Sqrt2 * Tp2xiT;
     fAmpl.f0Minus = -a-b;
     fAmpl.f0Plus  = -a+b;
     break;
   }
   case (kP11_1440) :
   {
     double c       = (5./12.)*kSqrt3;
     double xr      = (8./5.)* xi * fkr.R;
     double L2      = TMath::Power(fkr.Lamda, 2);
     double L2RmxiR = L2 * (fkr.Rminus + xr);
     double L2RpxiR = L2 * (fkr.Rplus  + xr);
     double a       = 0.25*kSqrt3 * L2 * fkr.S;
     double b       = c * (L2 * fkr.C - 2 * fkr.Lamda * fkr.B);

     fAmpl.fMinus1 = c * L2RmxiR;
     fAmpl.fPlus1  = c * L2RpxiR;
     fAmpl.fMinus3 = 0.;
     fAmpl.fPlus3  = 0.;
     fAmpl.f0Minus = a - b;
     fAmpl.f0Plus  = a + b;
     break;
   }
   case (kP33_1600) :
   {
     double xr      = 2*xi*fkr.R;
     double L2      = TMath::Power(fkr.Lamda, 2);
     double L2RmxiR = L2 * (fkr.Rminus + xr);
     double L2RpxiR = L2 * (fkr.Rplus  + xr);

     fAmpl.fMinus1 =  k1_Sqrt6 * L2RmxiR;
     fAmpl.fPlus1  = -k1_Sqrt6 * L2RpxiR;
     fAmpl.fMinus3 =  k1_Sqrt2 * L2RmxiR;
     fAmpl.fPlus3  = -k1_Sqrt2 * L2RpxiR;
     fAmpl.f0Minus = -kSqrt2_3 * (L2 * fkr.C - 2 * fkr.Lamda * fkr.B);
     fAmpl.f0Plus  =  fAmpl.f0Minus;
     break;
   }
   case (kP13_1720) :
   {
     double L2      = TMath::Power(fkr.Lamda, 2);
     double xr      = (8./5.)*xi*fkr.R;
     double LTm     = fkr.Lamda * fkr.Tminus;
     double LTp     = fkr.Lamda * fkr.Tplus;
     double L2RmxiR = L2 * (fkr.Rminus + xr);
     double L2RpxiR = L2 * (fkr.Rplus  + xr);
     double a       = kSqrt3_20 * L2 * fkr.S;
     double b       = kSqrt5_12 * (L2 * fkr.C - 5 * fkr.Lamda * fkr.B);

     fAmpl.fMinus1 =  kSqrt27_40 * LTm + kSqrt5_12 * L2RmxiR;
     fAmpl.fPlus1  = -kSqrt27_40 * LTp - kSqrt5_12 * L2RpxiR;
     fAmpl.fMinus3 = -kSqrt9_40 * LTm;
     fAmpl.fPlus3  =  kSqrt9_40 * LTp;
     fAmpl.f0Minus = -a+b;
     fAmpl.f0Plus  =  a+b;
     break;
   }
   case (kF15_1680) :
   {
     double L2      = TMath::Power(fkr.Lamda, 2);
     double xr      = (8./5.)*xi*fkr.R;
     double LTm     = fkr.Lamda * fkr.Tminus;
     double LTp     = fkr.Lamda * fkr.Tplus;
     double L2RmxiR = L2 * (fkr.Rminus + xr);
     double L2RpxiR = L2 * (fkr.Rplus  + xr);
     double a       = k3_Sqrt40 * L2 * fkr.S;
     double b       = kSqrt5_8  * L2 * fkr.C;

     fAmpl.fMinus1 =  k3_Sqrt20 * LTm - kSqrt5_8 * L2RmxiR;
     fAmpl.fPlus1  =  k3_Sqrt20 * LTp - kSqrt5_8 * L2RpxiR;
     fAmpl.fMinus3 =  kSqrt18_20 * LTm;
     fAmpl.fPlus3  =  kSqrt18_20 * LTp;
     fAmpl.f0Minus =  -a+b;
     fAmpl.f0Plus  =  -a-b;
     break;
   }
   case (kP31_1910) :
   {
     double xr       = 2*xi*fkr.R;
     double L2       = TMath::Power(fkr.Lamda, 2);

     fAmpl.fMinus1 = -k1_Sqrt15 * L2 * (fkr.Rminus + xr);
     fAmpl.fPlus1  = -k1_Sqrt15 * L2 * (fkr.Rplus  + xr);
     fAmpl.fMinus3 =  0.;
     fAmpl.fPlus3  =  0.;
     fAmpl.f0Minus = -kSqrt4_15 * (L2 * fkr.C - 5 * fkr.Lamda * fkr.B);
     fAmpl.f0Plus  = -1.* fAmpl.f0Minus;
     break;
   }
   case (kP33_1920) :
   {
     double xr       = 2*xi*fkr.R;
     double L2       = TMath::Power(fkr.Lamda, 2);
     double L2Rm2xiR = L2 * (fkr.Rminus + xr);
     double L2Rp2xiR = L2 * (fkr.Rplus  + xr);

     fAmpl.fMinus1 =  k1_Sqrt15 * L2Rm2xiR;
     fAmpl.fPlus1  = -k1_Sqrt15 * L2Rp2xiR;
     fAmpl.fMinus3 = -k1_Sqrt5  * L2Rm2xiR;
     fAmpl.fPlus3  =  k1_Sqrt5  * L2Rp2xiR;
     fAmpl.f0Minus = -(2./kSqrt15) * (L2 * fkr.C - 5 * fkr.Lamda * fkr.B);
     fAmpl.f0Plus  =  fAmpl.f0Minus;
     break;
   }
   case (kF35_1905) :
   {
     double xr       = 2*xi*fkr.R;
     double L2       = TMath::Power(fkr.Lamda, 2);
     double L2Rm2xiR = L2 * (fkr.Rminus + xr);
     double L2Rp2xiR = L2 * (fkr.Rplus  + xr);

     fAmpl.fMinus1 =  k1_Sqrt35  * L2Rm2xiR;
     fAmpl.fPlus1  =  k1_Sqrt35  * L2Rp2xiR;
     fAmpl.fMinus3 =  kSqrt18_35 * L2Rm2xiR;
     fAmpl.fPlus3  =  kSqrt18_35 * L2Rp2xiR;
     fAmpl.f0Minus =  k2_Sqrt35  * L2 * fkr.C;
     fAmpl.f0Plus  =  -1. * fAmpl.f0Minus;
     break;
   }
   case (kF37_1950) :
   {
     double xr       = 2*xi*fkr.R;
     double L2       = TMath::Power(fkr.Lamda, 2);
     double L2Rm2xiR = L2 * (fkr.Rminus + xr);
     double L2Rp2xiR = L2 * (fkr.Rplus  + xr);

     fAmpl.fMinus1 =  -kSqrt6_35 * L2Rm2xiR;
     fAmpl.fPlus1  =   kSqrt6_35 * L2Rp2xiR;
     fAmpl.fMinus3 =  -kSqrt2_7  * L2Rm2xiR;
     fAmpl.fPlus3  =   kSqrt2_7  * L2Rp2xiR;
     fAmpl.f0Minus = 2*kSqrt6_35 * L2 * fkr.C;
     fAmpl.f0Plus  =  fAmpl.f0Minus;
     break;
   }
   case (kP11_1710) :
   {
     double L2      = TMath::Power(fkr.Lamda, 2);
     double xr      = xi*fkr.R;
     double L2RmxiR = L2 * (fkr.Rminus + xr);
     double L2RpxiR = L2 * (fkr.Rplus  + xr);
     double a       = kSqrt3_8 * (1-2*xi) * L2 * fkr.S;
     double b       = k1_Sqrt6 * (L2 * fkr.C - 2 * fkr.Lamda * fkr.B);

     fAmpl.fMinus1 = -k1_Sqrt6 * L2RmxiR;
     fAmpl.fPlus1  = -k1_Sqrt6 * L2RpxiR;
     fAmpl.fMinus3 = 0.;
     fAmpl.fPlus3  = 0.;
     fAmpl.f0Minus = -a+b;
     fAmpl.f0Plus  = -a-b;
     break;
   }
   case (kF17_1970) :
   {
     fAmpl.fMinus1 = 0.;
     fAmpl.fPlus1  = 0.;
     fAmpl.fMinus3 = 0.;
     fAmpl.fPlus3  = 0.;
     fAmpl.f0Minus = 0.;
     fAmpl.f0Plus  = 0.;
     break;
   }
   default:
   {
     LOG("RSHAmpl", pWARN) << "*** UNRECOGNIZED RESONANCE!";
     fAmpl.fMinus1 = 0.;
     fAmpl.fPlus1  = 0.;
     fAmpl.fMinus3 = 0.;
     fAmpl.fPlus3  = 0.;
     fAmpl.f0Minus = 0.;
     fAmpl.f0Plus  = 0.;
     break;
   }

  }//switch

  return fAmpl;
}
//____________________________________________________________________________
void RSHelicityAmplModelNCn::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void RSHelicityAmplModelNCn::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void RSHelicityAmplModelNCn::LoadConfig(void)
{
  double thw ;
  GetParam( "WeinbergAngle", thw ) ;
  fSin28w = TMath::Power( TMath::Sin(thw), 2 );
}
//____________________________________________________________________________

