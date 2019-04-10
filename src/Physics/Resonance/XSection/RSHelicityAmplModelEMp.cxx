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
#include "Physics/Resonance/XSection/RSHelicityAmplModelEMp.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
RSHelicityAmplModelEMp::RSHelicityAmplModelEMp() :
RSHelicityAmplModelI("genie::RSHelicityAmplModelEMp")
{

}
//____________________________________________________________________________
RSHelicityAmplModelEMp::RSHelicityAmplModelEMp(string config) :
RSHelicityAmplModelI("genie::RSHelicityAmplModelEMp", config)
{

}
//____________________________________________________________________________
RSHelicityAmplModelEMp::~RSHelicityAmplModelEMp()
{

}
//____________________________________________________________________________
const RSHelicityAmpl & 
    RSHelicityAmplModelEMp::Compute(
          Resonance_t res, const FKR & fkr) const
{
  switch(res) {

   case (kP33_1232) :
   {
     fAmpl.fPlus1  =  kSqrt2 * fkr.R;
     fAmpl.fPlus3  =  kSqrt6 * fkr.R;
     fAmpl.fMinus1 = -1 * fAmpl.fPlus1;
     fAmpl.fMinus3 = -1 * fAmpl.fPlus3;
     fAmpl.f0Minus =  0.;
     fAmpl.f0Plus  =  0.;
     break;
   }
   case (kS11_1535) :
   {
     fAmpl.fMinus1 =  kSqrt3 * fkr.T + kSqrt3_2 * fkr.Lamda * fkr.R;
     fAmpl.f0Minus = -kSqrt3_2 * fkr.Lamda * fkr.S;
     fAmpl.fPlus1  = -1. * fAmpl.fMinus1;
     fAmpl.f0Plus  = -1. * fAmpl.f0Minus;
     fAmpl.fMinus3 =  0.;
     fAmpl.fPlus3  =  0.;
     break;
   }
   case (kD13_1520) :
   {
     fAmpl.fMinus1 =  kSqrt3_2 * fkr.T - kSqrt3 * fkr.Lamda * fkr.R;
     fAmpl.fMinus3 =  k3_Sqrt2 * fkr.T;
     fAmpl.f0Minus = -kSqrt3 * fkr.Lamda * fkr.S;
     fAmpl.fPlus1  =  fAmpl.fMinus1;
     fAmpl.fPlus3  =  fAmpl.fMinus3;
     fAmpl.f0Plus  =  fAmpl.f0Minus;
     break;
   }
   case (kS11_1650) :
   {
     fAmpl.fMinus1 = 0.;
     fAmpl.fPlus1  = 0.;
     fAmpl.fMinus3 = 0.;
     fAmpl.fPlus3  = 0.;
     fAmpl.f0Minus = 0.;
     fAmpl.f0Plus  = 0.;
     break;
   }
   case (kD13_1700) :
   {
     fAmpl.fMinus1 = 0.;
     fAmpl.fPlus1  = 0.;
     fAmpl.fMinus3 = 0.;
     fAmpl.fPlus3  = 0.;
     fAmpl.f0Minus = 0.;
     fAmpl.f0Plus  = 0.;
     break;
   }
   case (kD15_1675) :
   {
     fAmpl.fMinus1 = 0.;
     fAmpl.fPlus1  = 0.;
     fAmpl.fMinus3 = 0.;
     fAmpl.fPlus3  = 0.;
     fAmpl.f0Minus = 0.;
     fAmpl.f0Plus  = 0.;
     break;
   }
   case (kS31_1620) :
   {
     fAmpl.fMinus1 =  kSqrt3 * fkr.T - k1_Sqrt6 * fkr.Lamda * fkr.R;
     fAmpl.f0Minus = -kSqrt3_2 * fkr.Lamda * fkr.S;
     fAmpl.fPlus1  = -1. * fAmpl.fMinus1;
     fAmpl.f0Plus  = -1. * fAmpl.f0Minus;
     fAmpl.fMinus3 = 0.;
     fAmpl.fPlus3  = 0.;
     break;
   }
   case (kD33_1700) :
   {
     fAmpl.fMinus1 =  kSqrt3_2 * fkr.T + k1_Sqrt3 * fkr.Lamda * fkr.R;
     fAmpl.fMinus3 =  k3_Sqrt2 * fkr.T;
     fAmpl.f0Minus = -kSqrt3 * fkr.Lamda * fkr.S;
     fAmpl.fPlus1  =  fAmpl.fMinus1;
     fAmpl.fPlus3  =  fAmpl.fMinus3;
     fAmpl.f0Plus  =  fAmpl.f0Minus;
     break;
   }
   case (kP11_1440) :
   {
     double L2  = TMath::Power(fkr.Lamda, 2);

     fAmpl.fMinus1 = -0.5*kSqrt3 * L2 * fkr.R;
     fAmpl.fPlus1  =  fAmpl.fMinus1;
     fAmpl.fMinus3 =  0.;
     fAmpl.fPlus3  =  0.;
     fAmpl.f0Minus = -0.5*kSqrt3 * L2 * fkr.S;
     fAmpl.f0Plus  =  fAmpl.f0Minus;
     break;
   }
   case (kP33_1600) :
   {
     double L2R  = TMath::Power(fkr.Lamda, 2) * fkr.R;

     fAmpl.fMinus1 = k1_Sqrt6 * L2R;
     fAmpl.fMinus3 = k1_Sqrt2 * L2R;
     fAmpl.fPlus1  = -1. * fAmpl.fMinus1;
     fAmpl.fPlus3  = -1. * fAmpl.fMinus3;
     fAmpl.f0Minus = 0.;
     fAmpl.f0Plus  = 0.;
     break;
   }
   case (kP13_1720) :
   {
     double L2  = TMath::Power(fkr.Lamda, 2);
     double LT  = fkr.Lamda * fkr.T;

     fAmpl.fMinus1 = -kSqrt27_10 * LT - kSqrt3_5 * L2 * fkr.R;
     fAmpl.fMinus3 =  k3_Sqrt10 * LT;
     fAmpl.f0Minus =  kSqrt3_5  * L2 * fkr.S;
     fAmpl.fPlus1  = -1. * fAmpl.fMinus1;
     fAmpl.fPlus3  = -1. * fAmpl.fMinus3;
     fAmpl.f0Plus  = -1. * fAmpl.f0Minus;
     break;
   }
   case (kF15_1680) :
   {
     double L2  = TMath::Power(fkr.Lamda, 2);
     double LT  = fkr.Lamda * fkr.T;

     fAmpl.fMinus1 =  -k3_Sqrt5  * LT + k3_Sqrt10 * L2 * fkr.R;
     fAmpl.fMinus3 =  -kSqrt18_5 * LT;
     fAmpl.f0Minus =   k3_Sqrt10 * L2 * fkr.S;
     fAmpl.fPlus1  =  fAmpl.fMinus1;
     fAmpl.fPlus3  =  fAmpl.fMinus3;
     fAmpl.f0Plus  =  fAmpl.f0Minus;
     break;
   }
   case (kP31_1910) :
   {
     fAmpl.fMinus1 = -k1_Sqrt15 * TMath::Power(fkr.Lamda, 2) * fkr.R;
     fAmpl.fPlus1  = fAmpl.fMinus1;
     fAmpl.fMinus3 = 0.;
     fAmpl.fPlus3  = 0.;
     fAmpl.f0Minus = 0.;
     fAmpl.f0Plus  = 0.;
     break;
   }
   case (kP33_1920) :
   {
     double L2R  = TMath::Power(fkr.Lamda, 2) * fkr.R;

     fAmpl.fMinus1 =  k1_Sqrt15 * L2R;
     fAmpl.fMinus3 = -k1_Sqrt5  * L2R;
     fAmpl.fPlus1  = -1.* fAmpl.fMinus1;
     fAmpl.fPlus3  = -1.* fAmpl.fMinus3;
     fAmpl.f0Minus =  0.;
     fAmpl.f0Plus  =  0.;
     break;
   }
   case (kF35_1905) :
   {
     double L2R  = TMath::Power(fkr.Lamda, 2) * fkr.R;

     fAmpl.fMinus1 = k1_Sqrt35  * L2R;
     fAmpl.fMinus3 = kSqrt18_35 * L2R;
     fAmpl.fPlus1  = fAmpl.fMinus1;
     fAmpl.fPlus3  = fAmpl.fMinus3;
     fAmpl.f0Minus = 0.;
     fAmpl.f0Plus  = 0.;
     break;
   }
   case (kF37_1950) :
   {
     double L2R  = TMath::Power(fkr.Lamda, 2) * fkr.R;

     fAmpl.fMinus1 = -kSqrt6_35 * L2R;
     fAmpl.fMinus3 = -kSqrt2_7  * L2R;
     fAmpl.fPlus1  = -1. * fAmpl.fMinus1;
     fAmpl.fPlus3  = -1. * fAmpl.fMinus3;
     fAmpl.f0Minus = 0.;
     fAmpl.f0Plus  = 0.;
     break;
   }
   case (kP11_1710) :
   {
     double L2  = TMath::Power(fkr.Lamda, 2);

     fAmpl.fMinus1 = kSqrt3_8 * L2 * fkr.R;
     fAmpl.f0Minus = kSqrt3_8 * L2 * fkr.S;
     fAmpl.fPlus1  = fAmpl.fMinus1;
     fAmpl.f0Plus  = fAmpl.f0Minus;
     fAmpl.fMinus3 = 0.;
     fAmpl.fPlus3  = 0.;
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


