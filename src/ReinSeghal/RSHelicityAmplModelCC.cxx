//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - May 03, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "BaryonResonance/BaryonResUtils.h"
#include "Conventions/Constants.h"
#include "ReinSeghal/RSHelicityAmplModelCC.h"
#include "ReinSeghal/RSHelicityAmpl.h"
#include "Messenger/Messenger.h"

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
RSHelicityAmpl * RSHelicityAmplModelCC::Compute(
                                       Resonance_t res, const FKR & fkr) const
{
  RSHelicityAmpl * hampl = new RSHelicityAmpl;

  switch(res) {

   case (kP33_1232) :
   {
     hampl->fMinus1 =    kSqrt2 * fkr.Rminus;
     hampl->fPlus1  =   -kSqrt2 * fkr.Rplus;
     hampl->fMinus3 =    kSqrt6 * fkr.Rminus;
     hampl->fPlus3  =   -kSqrt6 * fkr.Rplus;
     hampl->f0Minus = -2*kSqrt2 * fkr.C;
     hampl->f0Plus  =    hampl->f0Minus;
     break;
   }
   case (kS11_1535) :
   {
     double c = 4./kSqrt6;
     double d = 2.*kSqrt3;
     double a = kSqrt6 * fkr.Lamda * fkr.S;
     double b = 2 * kSqrt2_3 * (fkr.Lamda * fkr.C - 3.* fkr.B);
     
     hampl->fMinus1 =  d * fkr.Tminus + c * fkr.Lamda * fkr.Rminus;
     hampl->fPlus1  = -d * fkr.Tplus  - c * fkr.Lamda * fkr.Rplus;
     hampl->fMinus3 =  0;
     hampl->fPlus3  =  0;
     hampl->f0Minus = -a+b;
     hampl->f0Plus  =  a+b;
     break;
   }
   case (kD13_1520) :
   {
     double c = 4./kSqrt3;
     double d = 6./kSqrt2;
     double a = 2.* kSqrt3 * fkr.Lamda * fkr.S;
     double b = (4./kSqrt3)* fkr.Lamda * fkr.C;

     hampl->fMinus1 =  kSqrt6 * fkr.Tminus - c * fkr.Lamda * fkr.Rminus;
     hampl->fPlus1  =  kSqrt6 * fkr.Tplus  - c * fkr.Lamda * fkr.Rplus;
     hampl->fMinus3 =  d * fkr.Tminus;
     hampl->fPlus3  =  d * fkr.Tplus;
     hampl->f0Minus =  -a+b;
     hampl->f0Plus  =  -a-b;
     break;
   }
   case (kS11_1650) :
   {
     hampl->fMinus1 =  k1_Sqrt6 * fkr.Lamda * fkr.Rminus;
     hampl->fPlus1  = -k1_Sqrt6 * fkr.Lamda * fkr.Rplus;
     hampl->fMinus3 =  0;
     hampl->fPlus3  =  0;
     hampl->f0Minus = -kSqrt2_3 * (fkr.Lamda * fkr.C - 3.* fkr.B);
     hampl->f0Plus  =  hampl->f0Minus;
     break;
   }
   case (kD13_1700) :
   {
     double LRm = fkr.Lamda * fkr.Rminus;
     double LRp = fkr.Lamda * fkr.Rplus;

     hampl->fMinus1 =  k1_Sqrt30 * LRm;
     hampl->fPlus1  =  k1_Sqrt30 * LRp;
     hampl->fMinus3 =  k3_Sqrt10 * LRm;
     hampl->fPlus3  =  k3_Sqrt10 * LRp;
     hampl->f0Minus =  kSqrt2_15 * fkr.Lamda * fkr.C;
     hampl->f0Plus  =  -1. * hampl->f0Minus;
     break;
   }
   case (kD15_1675) :
   {
     double LRm = fkr.Lamda * fkr.Rminus;
     double LRp = fkr.Lamda * fkr.Rplus;

     hampl->fMinus1 = -kSqrt3_10 * LRm;
     hampl->fPlus1  =  kSqrt3_10 * LRp;
     hampl->fMinus3 = -kSqrt3_5  * LRm;
     hampl->fPlus3  =  kSqrt3_5  * LRp;
     hampl->f0Minus =  kSqrt6_5  * fkr.Lamda * fkr.C;
     hampl->f0Plus  =  hampl->f0Minus;
     break;
   }
   case (kS31_1620) :
   {
     double a = kSqrt3_2 * fkr.Lamda * fkr.S;
     double b = k1_Sqrt6 * (fkr.Lamda * fkr.C - 3.* fkr.B);

     hampl->fMinus1 = -kSqrt3 * fkr.Tminus + k1_Sqrt6 * fkr.Lamda * fkr.Rminus;
     hampl->fPlus1  =  kSqrt3 * fkr.Tplus  - k1_Sqrt6 * fkr.Lamda * fkr.Rplus;
     hampl->fMinus3 =  0;
     hampl->fPlus3  =  0;
     hampl->f0Minus =  a+b;
     hampl->f0Plus  = -a+b;
     break;
   }
   case (kD33_1700) :
   {
     double a = kSqrt3   * fkr.Lamda * fkr.S;
     double b = k1_Sqrt3 * fkr.Lamda * fkr.C;

     hampl->fMinus1 = -kSqrt3_2 * fkr.Tminus - k1_Sqrt3 * fkr.Lamda * fkr.Rminus;
     hampl->fPlus1  = -kSqrt3_2 * fkr.Tplus  - k1_Sqrt3 * fkr.Lamda * fkr.Rplus;
     hampl->fMinus3 = -k3_Sqrt2 * fkr.Tminus;
     hampl->fPlus3  = -k3_Sqrt2 * fkr.Tplus;
     hampl->f0Minus =  a + b;
     hampl->f0Plus  =  a - b;
     break;
   }
   case (kP11_1440) :
   {
     double L2 = TMath::Power(fkr.Lamda, 2);
     double c  = 5.*kSqrt3/6.;
     double a  = kSqrt3_4 * L2 * fkr.S;
     double b  = c * (L2 * fkr.C - 2 * fkr.Lamda * fkr.B);

     hampl->fMinus1 =  -c * L2 * fkr.Rminus;
     hampl->fPlus1  =  -c * L2 * fkr.Rplus;
     hampl->fMinus3 =   0;
     hampl->fPlus3  =   0;
     hampl->f0Minus =  -a+b;
     hampl->f0Plus  =  -a-b;
     break;
   }
   case (kP33_1600) :
   {
     double L2      = TMath::Power(fkr.Lamda, 2);
     double L2Rm    = L2 * fkr.Rminus;
     double L2Rp    = L2 * fkr.Rplus;

     hampl->fMinus1 = -k1_Sqrt6 * L2Rm;
     hampl->fPlus1  =  k1_Sqrt6 * L2Rp;
     hampl->fMinus3 = -k1_Sqrt2 * L2Rm;
     hampl->fPlus3  =  k1_Sqrt2 * L2Rp;
     hampl->f0Minus =  kSqrt2_3 * (L2 * fkr.C - 2 * fkr.Lamda * fkr.B);
     hampl->f0Plus  =  hampl->f0Minus;
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

     hampl->fMinus1 =  -kSqrt27_10 * LTm - kSqrt5_3 * L2Rm;
     hampl->fPlus1  =   kSqrt27_10 * LTm + kSqrt5_3 * L2Rp;
     hampl->fMinus3 =   k3_Sqrt10 * LTm;
     hampl->fPlus3  =  -k3_Sqrt10 * LTp;
     hampl->f0Minus =   a-b;
     hampl->f0Plus  =  -a-b;
     break;
   }
   case (kF15_1680) :
   {
     double L2  = TMath::Power(fkr.Lamda, 2);
     double LTm = fkr.Lamda * fkr.Tminus;
     double LTp = fkr.Lamda * fkr.Tplus;
     double a   = kSqrt9_10 * L2 * fkr.S;
     double b   = kSqrt5_2  * L2 * fkr.C;

     hampl->fMinus1 = -k3_Sqrt5  * LTm + kSqrt5_2 * L2 * fkr.Rminus;
     hampl->fPlus1  = -k3_Sqrt5  * LTp + kSqrt5_2 * L2 * fkr.Rplus;
     hampl->fMinus3 = -kSqrt18_5 * LTm;
     hampl->fPlus3  = -kSqrt18_5 * LTp;
     hampl->f0Minus =  a - b;
     hampl->f0Plus  =  a + b;
     break;
   }
   case (kP31_1910) :
   {
     double L2 = TMath::Power(fkr.Lamda, 2);

     hampl->fMinus1 =  k1_Sqrt15 * L2 * fkr.Rminus;
     hampl->fPlus1  =  k1_Sqrt15 * L2 * fkr.Rplus;
     hampl->fMinus3 =  0;
     hampl->fPlus3  =  0;
     hampl->f0Minus =  k2_Sqrt15 * (L2 * fkr.C - 5 * fkr.Lamda * fkr.B);
     hampl->f0Plus  = -1.* hampl->f0Minus;
     break;
   }
   case (kP33_1920) :
   {
     double L2   = TMath::Power(fkr.Lamda, 2);
     double L2Rm = L2 * fkr.Rminus;
     double L2Rp = L2 * fkr.Rplus;

     hampl->fMinus1 = -k1_Sqrt15 * L2Rm;
     hampl->fPlus1  =  k1_Sqrt15 * L2Rp;
     hampl->fMinus3 =  k1_Sqrt5  * L2Rm;
     hampl->fPlus3  = -k1_Sqrt5  * L2Rp;
     hampl->f0Minus =  k2_Sqrt15 * (L2 * fkr.C - 5 * fkr.Lamda * fkr.B);
     hampl->f0Plus  =  hampl->f0Minus;
     break;
   }
   case (kF35_1905) :
   {
     double L2      = TMath::Power(fkr.Lamda, 2);
     double L2Rm    = L2 * fkr.Rminus;
     double L2Rp    = L2 * fkr.Rplus;

     hampl->fMinus1 =  -k1_Sqrt35  * L2Rm;
     hampl->fPlus1  =  -k1_Sqrt35  * L2Rp;
     hampl->fMinus3 =  -kSqrt18_35 * L2Rm;
     hampl->fPlus3  =  -kSqrt18_35 * L2Rp;
     hampl->f0Minus =  -k2_Sqrt35  * L2 * fkr.C;
     hampl->f0Plus  =  -1.* hampl->f0Minus;
     break;
   }
   case (kF37_1950) :
   {
     double L2      = TMath::Power(fkr.Lamda, 2);
     double L2Rm    = L2 * fkr.Rminus;
     double L2Rp    = L2 * fkr.Rplus;

     hampl->fMinus1 =  kSqrt6_35  * L2Rm;
     hampl->fPlus1  = -kSqrt6_35  * L2Rp;
     hampl->fMinus3 =  kSqrt2_7   * L2Rm;
     hampl->fPlus3  = -kSqrt2_7   * L2Rp;
     hampl->f0Minus = -kSqrt24_35 * L2 * fkr.C;
     hampl->f0Plus  =  hampl->f0Minus;
     break;
   }
   case (kP11_1710) :
   {
     double L2 = TMath::Power(fkr.Lamda, 2);
     double a  = kSqrt3_2 * L2 * fkr.S;
     double b  = kSqrt2_3 * (L2 * fkr.C - 2 * fkr.Lamda * fkr.B);

     hampl->fMinus1 = kSqrt2_3 * L2 * fkr.Rminus;
     hampl->fPlus1  = kSqrt2_3 * L2 * fkr.Rplus;
     hampl->fMinus3 = 0;
     hampl->fPlus3  = 0;
     hampl->f0Minus = a - b;
     hampl->f0Plus  = a + b;
     break;
   }
   case (kF17_1970) :
   {
     double L2   = TMath::Power(fkr.Lamda, 2);
     double L2Rm = L2 * fkr.Rminus;
     double L2Rp = L2 * fkr.Rplus;

     hampl->fMinus1 =  -kSqrt3_35 * L2Rm;
     hampl->fPlus1  =   kSqrt3_35 * L2Rp;
     hampl->fMinus3 =  -k1_Sqrt7  * L2Rm;
     hampl->fPlus3  =   k1_Sqrt7  * L2Rp;
     hampl->f0Minus =   kSqrt6_35 * L2 * fkr.C;
     hampl->f0Plus  =   hampl->f0Minus;
     break;
   }
   default:
     LOG("RSHAmpl", pWARN) << "*** UNRECOGNIZED RESONANCE!";
     delete hampl;
     hampl = 0;
     break;
   }
  return hampl;
}
//____________________________________________________________________________


