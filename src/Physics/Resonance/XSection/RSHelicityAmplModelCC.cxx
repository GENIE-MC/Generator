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

#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/Resonance/XSection/RSHelicityAmplModelCC.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
RSHelicityAmplModelCC::RSHelicityAmplModelCC() :
RSHelicityAmplModelI("genie::RSHelicityAmplModelCC")
{

}
//____________________________________________________________________________
RSHelicityAmplModelCC::RSHelicityAmplModelCC(string config) :
RSHelicityAmplModelI("genie::RSHelicityAmplModelCC", config)
{

}
//____________________________________________________________________________
RSHelicityAmplModelCC::~RSHelicityAmplModelCC()
{

}
//____________________________________________________________________________
const RSHelicityAmpl & 
  RSHelicityAmplModelCC::Compute(
      Resonance_t res, const FKR & fkr) const
{
  switch(res) {

   case (kP33_1232) :
   {
     fAmpl.fMinus1 =    kSqrt2 * fkr.Rminus;
     fAmpl.fPlus1  =   -kSqrt2 * fkr.Rplus;
     fAmpl.fMinus3 =    kSqrt6 * fkr.Rminus;
     fAmpl.fPlus3  =   -kSqrt6 * fkr.Rplus;
     fAmpl.f0Minus = -2*kSqrt2 * fkr.C;
     fAmpl.f0Plus  =    fAmpl.f0Minus;
     break;
   }
   case (kS11_1535) :
   {
     double c = 4./kSqrt6;
     double d = 2.*kSqrt3;
     double a = kSqrt6 * fkr.Lamda * fkr.S;
     double b = 2 * kSqrt2_3 * (fkr.Lamda * fkr.C - 3.* fkr.B);
     
     fAmpl.fMinus1 =  d * fkr.Tminus + c * fkr.Lamda * fkr.Rminus;
     fAmpl.fPlus1  = -d * fkr.Tplus  - c * fkr.Lamda * fkr.Rplus;
     fAmpl.fMinus3 =  0;
     fAmpl.fPlus3  =  0;
     fAmpl.f0Minus = -a+b;
     fAmpl.f0Plus  =  a+b;
     break;
   }
   case (kD13_1520) :
   {
     double c = 4./kSqrt3;
     double d = 6./kSqrt2;
     double a = 2.* kSqrt3 * fkr.Lamda * fkr.S;
     double b = (4./kSqrt3)* fkr.Lamda * fkr.C;

     fAmpl.fMinus1 =  kSqrt6 * fkr.Tminus - c * fkr.Lamda * fkr.Rminus;
     fAmpl.fPlus1  =  kSqrt6 * fkr.Tplus  - c * fkr.Lamda * fkr.Rplus;
     fAmpl.fMinus3 =  d * fkr.Tminus;
     fAmpl.fPlus3  =  d * fkr.Tplus;
     fAmpl.f0Minus =  -a+b;
     fAmpl.f0Plus  =  -a-b;
     break;
   }
   case (kS11_1650) :
   {
     fAmpl.fMinus1 =  k1_Sqrt6 * fkr.Lamda * fkr.Rminus;
     fAmpl.fPlus1  = -k1_Sqrt6 * fkr.Lamda * fkr.Rplus;
     fAmpl.fMinus3 =  0;
     fAmpl.fPlus3  =  0;
     fAmpl.f0Minus = -kSqrt2_3 * (fkr.Lamda * fkr.C - 3.* fkr.B);
     fAmpl.f0Plus  =  fAmpl.f0Minus;
     break;
   }
   case (kD13_1700) :
   {
     double LRm = fkr.Lamda * fkr.Rminus;
     double LRp = fkr.Lamda * fkr.Rplus;

     fAmpl.fMinus1 =  k1_Sqrt30 * LRm;
     fAmpl.fPlus1  =  k1_Sqrt30 * LRp;
     fAmpl.fMinus3 =  k3_Sqrt10 * LRm;
     fAmpl.fPlus3  =  k3_Sqrt10 * LRp;
     fAmpl.f0Minus =  kSqrt2_15 * fkr.Lamda * fkr.C;
     fAmpl.f0Plus  =  -1. * fAmpl.f0Minus;
     break;
   }
   case (kD15_1675) :
   {
     double LRm = fkr.Lamda * fkr.Rminus;
     double LRp = fkr.Lamda * fkr.Rplus;

     fAmpl.fMinus1 = -kSqrt3_10 * LRm;
     fAmpl.fPlus1  =  kSqrt3_10 * LRp;
     fAmpl.fMinus3 = -kSqrt3_5  * LRm;
     fAmpl.fPlus3  =  kSqrt3_5  * LRp;
     fAmpl.f0Minus =  kSqrt6_5  * fkr.Lamda * fkr.C;
     fAmpl.f0Plus  =  fAmpl.f0Minus;
     break;
   }
   case (kS31_1620) :
   {
     double a = kSqrt3_2 * fkr.Lamda * fkr.S;
     double b = k1_Sqrt6 * (fkr.Lamda * fkr.C - 3.* fkr.B);

     fAmpl.fMinus1 = -kSqrt3 * fkr.Tminus + k1_Sqrt6 * fkr.Lamda * fkr.Rminus;
     fAmpl.fPlus1  =  kSqrt3 * fkr.Tplus  - k1_Sqrt6 * fkr.Lamda * fkr.Rplus;
     fAmpl.fMinus3 =  0;
     fAmpl.fPlus3  =  0;
     fAmpl.f0Minus =  a+b;
     fAmpl.f0Plus  = -a+b;
     break;
   }
   case (kD33_1700) :
   {
     double a = kSqrt3   * fkr.Lamda * fkr.S;
     double b = k1_Sqrt3 * fkr.Lamda * fkr.C;

     fAmpl.fMinus1 = -kSqrt3_2 * fkr.Tminus - k1_Sqrt3 * fkr.Lamda * fkr.Rminus;
     fAmpl.fPlus1  = -kSqrt3_2 * fkr.Tplus  - k1_Sqrt3 * fkr.Lamda * fkr.Rplus;
     fAmpl.fMinus3 = -k3_Sqrt2 * fkr.Tminus;
     fAmpl.fPlus3  = -k3_Sqrt2 * fkr.Tplus;
     fAmpl.f0Minus =  a + b;
     fAmpl.f0Plus  =  a - b;
     break;
   }
   case (kP11_1440) :
   {
     double L2 = TMath::Power(fkr.Lamda, 2);
     double c  = 5.*kSqrt3/6.;
     double a  = kSqrt3_4 * L2 * fkr.S;
     double b  = c * (L2 * fkr.C - 2 * fkr.Lamda * fkr.B);

     fAmpl.fMinus1 =  -c * L2 * fkr.Rminus;
     fAmpl.fPlus1  =  -c * L2 * fkr.Rplus;
     fAmpl.fMinus3 =   0;
     fAmpl.fPlus3  =   0;
     fAmpl.f0Minus =  -a+b;
     fAmpl.f0Plus  =  -a-b;
     break;
   }
   case (kP33_1600) :
   {
     double L2      = TMath::Power(fkr.Lamda, 2);
     double L2Rm    = L2 * fkr.Rminus;
     double L2Rp    = L2 * fkr.Rplus;

     fAmpl.fMinus1 = -k1_Sqrt6 * L2Rm;
     fAmpl.fPlus1  =  k1_Sqrt6 * L2Rp;
     fAmpl.fMinus3 = -k1_Sqrt2 * L2Rm;
     fAmpl.fPlus3  =  k1_Sqrt2 * L2Rp;
     fAmpl.f0Minus =  kSqrt2_3 * (L2 * fkr.C - 2 * fkr.Lamda * fkr.B);
     fAmpl.f0Plus  =  fAmpl.f0Minus;
     break;
   }
   case (kP13_1720) :
   {
     double L2      = TMath::Power(fkr.Lamda, 2);
     double L2Rm    = L2 * fkr.Rminus;
     double L2Rp    = L2 * fkr.Rplus;
     double LTm     = fkr.Lamda * fkr.Tminus;
     double LTp     = fkr.Lamda * fkr.Tplus;
     double a       = kSqrt3_5 * L2 * fkr.S;
     double b       = kSqrt5_3 * (L2 * fkr.C - 5 * fkr.Lamda * fkr.B);

     fAmpl.fMinus1 =  -kSqrt27_10 * LTm - kSqrt5_3 * L2Rm;
     fAmpl.fPlus1  =   kSqrt27_10 * LTp + kSqrt5_3 * L2Rp;
     fAmpl.fMinus3 =   k3_Sqrt10 * LTm;
     fAmpl.fPlus3  =  -k3_Sqrt10 * LTp;
     fAmpl.f0Minus =   a-b;
     fAmpl.f0Plus  =  -a-b;
     break;
   }
   case (kF15_1680) :
   {
     double L2  = TMath::Power(fkr.Lamda, 2);
     double LTm = fkr.Lamda * fkr.Tminus;
     double LTp = fkr.Lamda * fkr.Tplus;
     double a   = kSqrt9_10 * L2 * fkr.S;
     double b   = kSqrt5_2  * L2 * fkr.C;

     fAmpl.fMinus1 = -k3_Sqrt5  * LTm + kSqrt5_2 * L2 * fkr.Rminus;
     fAmpl.fPlus1  = -k3_Sqrt5  * LTp + kSqrt5_2 * L2 * fkr.Rplus;
     fAmpl.fMinus3 = -kSqrt18_5 * LTm;
     fAmpl.fPlus3  = -kSqrt18_5 * LTp;
     fAmpl.f0Minus =  a - b;
     fAmpl.f0Plus  =  a + b;
     break;
   }
   case (kP31_1910) :
   {
     double L2 = TMath::Power(fkr.Lamda, 2);

     fAmpl.fMinus1 =  k1_Sqrt15 * L2 * fkr.Rminus;
     fAmpl.fPlus1  =  k1_Sqrt15 * L2 * fkr.Rplus;
     fAmpl.fMinus3 =  0;
     fAmpl.fPlus3  =  0;
     fAmpl.f0Minus =  k2_Sqrt15 * (L2 * fkr.C - 5 * fkr.Lamda * fkr.B);
     fAmpl.f0Plus  = -1.* fAmpl.f0Minus;
     break;
   }
   case (kP33_1920) :
   {
     double L2   = TMath::Power(fkr.Lamda, 2);
     double L2Rm = L2 * fkr.Rminus;
     double L2Rp = L2 * fkr.Rplus;

     fAmpl.fMinus1 = -k1_Sqrt15 * L2Rm;
     fAmpl.fPlus1  =  k1_Sqrt15 * L2Rp;
     fAmpl.fMinus3 =  k1_Sqrt5  * L2Rm;
     fAmpl.fPlus3  = -k1_Sqrt5  * L2Rp;
     fAmpl.f0Minus =  k2_Sqrt15 * (L2 * fkr.C - 5 * fkr.Lamda * fkr.B);
     fAmpl.f0Plus  =  fAmpl.f0Minus;
     break;
   }
   case (kF35_1905) :
   {
     double L2      = TMath::Power(fkr.Lamda, 2);
     double L2Rm    = L2 * fkr.Rminus;
     double L2Rp    = L2 * fkr.Rplus;

     fAmpl.fMinus1 =  -k1_Sqrt35  * L2Rm;
     fAmpl.fPlus1  =  -k1_Sqrt35  * L2Rp;
     fAmpl.fMinus3 =  -kSqrt18_35 * L2Rm;
     fAmpl.fPlus3  =  -kSqrt18_35 * L2Rp;
     fAmpl.f0Minus =  -k2_Sqrt35  * L2 * fkr.C;
     fAmpl.f0Plus  =  -1.* fAmpl.f0Minus;
     break;
   }
   case (kF37_1950) :
   {
     double L2      = TMath::Power(fkr.Lamda, 2);
     double L2Rm    = L2 * fkr.Rminus;
     double L2Rp    = L2 * fkr.Rplus;

     fAmpl.fMinus1 =  kSqrt6_35  * L2Rm;
     fAmpl.fPlus1  = -kSqrt6_35  * L2Rp;
     fAmpl.fMinus3 =  kSqrt2_7   * L2Rm;
     fAmpl.fPlus3  = -kSqrt2_7   * L2Rp;
     fAmpl.f0Minus = -kSqrt24_35 * L2 * fkr.C;
     fAmpl.f0Plus  =  fAmpl.f0Minus;
     break;
   }
   case (kP11_1710) :
   {
     double L2 = TMath::Power(fkr.Lamda, 2);
     double a  = kSqrt3_2 * L2 * fkr.S;
     double b  = kSqrt2_3 * (L2 * fkr.C - 2 * fkr.Lamda * fkr.B);

     fAmpl.fMinus1 = kSqrt2_3 * L2 * fkr.Rminus;
     fAmpl.fPlus1  = kSqrt2_3 * L2 * fkr.Rplus;
     fAmpl.fMinus3 = 0;
     fAmpl.fPlus3  = 0;
     fAmpl.f0Minus = a - b;
     fAmpl.f0Plus  = a + b;
     break;
   }
   case (kF17_1970) :
   {
     double L2   = TMath::Power(fkr.Lamda, 2);
     double L2Rm = L2 * fkr.Rminus;
     double L2Rp = L2 * fkr.Rplus;

     fAmpl.fMinus1 =  -kSqrt3_35 * L2Rm;
     fAmpl.fPlus1  =   kSqrt3_35 * L2Rp;
     fAmpl.fMinus3 =  -k1_Sqrt7  * L2Rm;
     fAmpl.fPlus3  =   k1_Sqrt7  * L2Rp;
     fAmpl.f0Minus =   kSqrt6_35 * L2 * fkr.C;
     fAmpl.f0Plus  =   fAmpl.f0Minus;
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


