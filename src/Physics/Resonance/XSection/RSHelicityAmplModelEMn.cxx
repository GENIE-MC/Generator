//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool
*/
//____________________________________________________________________________

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/Resonance/XSection/RSHelicityAmplModelEMn.h"
#include "Physics/Resonance/XSection/RSHelicityAmpl.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
RSHelicityAmplModelEMn::RSHelicityAmplModelEMn() :
RSHelicityAmplModelI("genie::RSHelicityAmplModelEMn")
{

}
//____________________________________________________________________________
RSHelicityAmplModelEMn::RSHelicityAmplModelEMn(string config) :
RSHelicityAmplModelI("genie::RSHelicityAmplModelEMn", config)
{

}
//____________________________________________________________________________
RSHelicityAmplModelEMn::~RSHelicityAmplModelEMn()
{

}
//____________________________________________________________________________
const RSHelicityAmpl &
   RSHelicityAmplModelEMn::Compute(
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
     fAmpl.fPlus1  =  kSqrt3   * fkr.T + k1_Sqrt6 * fkr.Lamda * fkr.R;
     fAmpl.f0Minus =  kSqrt3_2 * fkr.Lamda * fkr.S;
     fAmpl.fMinus1 = -1 * fAmpl.fPlus1;
     fAmpl.f0Plus  = -1 * fAmpl.f0Minus;
     fAmpl.fMinus3 =  0.;
     fAmpl.fPlus3  =  0.;

     break;
   }
   case (kD13_1520) :
   {
     fAmpl.fMinus1 = -kSqrt3_2 * fkr.T + k1_Sqrt3 * fkr.Lamda * fkr.R;
     fAmpl.fMinus3 = -k3_Sqrt2 * fkr.T;
     fAmpl.f0Minus =  kSqrt3 * fkr.Lamda * fkr.S;
     fAmpl.fPlus1  =  fAmpl.fMinus1;
     fAmpl.fPlus3  =  fAmpl.fMinus3;
     fAmpl.f0Plus  =  fAmpl.f0Minus;
     break;
   }
   case (kS11_1650) :
   {
     fAmpl.fPlus1  =  k1_Sqrt6 * fkr.Lamda * fkr.R;
     fAmpl.fMinus1 = -1 * fAmpl.fPlus1;
     fAmpl.fMinus3 =  0.;
     fAmpl.fPlus3  =  0.;
     fAmpl.f0Minus =  0.;
     fAmpl.f0Plus  =  0.;
     break;
   }
   case (kD13_1700) :
   {
     double LR = fkr.Lamda * fkr.R;

     fAmpl.fMinus1 = -(1./kSqrt30) * LR;
     fAmpl.fMinus3 = -(3./kSqrt10) * LR;
     fAmpl.fPlus1  =  fAmpl.fMinus1;
     fAmpl.fPlus3  =  fAmpl.fMinus3;
     fAmpl.f0Minus =  0.;
     fAmpl.f0Plus  =  0.;
     break;
   }
   case (kD15_1675) :
   {
     double LR = fkr.Lamda * fkr.R;

     fAmpl.fMinus1 = kSqrt3_10 * LR;
     fAmpl.fMinus3 = kSqrt3_5  * LR;
     fAmpl.fPlus1  = -1 * fAmpl.fMinus1;
     fAmpl.fPlus3  = -1 * fAmpl.fMinus3;
     fAmpl.f0Minus =  0.;
     fAmpl.f0Plus  =  0.;
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
     fAmpl.fMinus1 = k1_Sqrt3 * TMath::Power(fkr.Lamda, 2) * fkr.R;
     fAmpl.fPlus1  = fAmpl.fMinus1;
     fAmpl.fMinus3 =  0.;
     fAmpl.fPlus3  =  0.;
     fAmpl.f0Minus =  0.;
     fAmpl.f0Plus  =  0.;
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
     fAmpl.fMinus1 = k2_Sqrt15 * TMath::Power(fkr.Lamda, 2) * fkr.R;
     fAmpl.fPlus1  = -1 * fAmpl.fMinus1;
     fAmpl.fMinus3 =  0.;
     fAmpl.fPlus3  =  0.;
     fAmpl.f0Minus =  0.;
     fAmpl.f0Plus  =  0.;
     break;
   }
   case (kF15_1680) :
   {
     fAmpl.fMinus1 =  -kSqrt2_5 * TMath::Power(fkr.Lamda, 2) * fkr.R;
     fAmpl.fPlus1  =  fAmpl.fMinus1;
     fAmpl.fMinus3 =  0.;
     fAmpl.fPlus3  =  0.;
     fAmpl.f0Minus =  0.;
     fAmpl.f0Plus  =  0.;
     break;
   }
   case (kP31_1910) :
   {
     fAmpl.fMinus1 =  -k1_Sqrt15 * TMath::Power(fkr.Lamda, 2) * fkr.R;
     fAmpl.fPlus1  =  fAmpl.fMinus1;
     fAmpl.fMinus3 =  0.;
     fAmpl.fPlus3  =  0.;
     fAmpl.f0Minus =  0.;
     fAmpl.f0Plus  =  0.;
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
     double L2 = TMath::Power(fkr.Lamda, 2);

     fAmpl.fMinus1 = -k1_Sqrt24 * L2 * fkr.R;
     fAmpl.f0Minus = -kSqrt3_8  * L2 * fkr.S;
     fAmpl.fPlus1  = fAmpl.fMinus1;
     fAmpl.f0Plus  = fAmpl.f0Minus;
     fAmpl.fMinus3 = 0.;
     fAmpl.fPlus3  = 0.;

     break;
   }
   case (kF17_1970) :
   {
     double L2R = TMath::Power(fkr.Lamda, 2) * fkr.R;

     fAmpl.fMinus1 = kSqrt3_35 * L2R;
     fAmpl.fPlus1  = -1 * fAmpl.fMinus1;
     fAmpl.fMinus3 = k1_Sqrt7  * L2R;
     fAmpl.fPlus3  = -1 * fAmpl.fMinus3;
     fAmpl.f0Minus =  0.;
     fAmpl.f0Plus  =  0.;
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
