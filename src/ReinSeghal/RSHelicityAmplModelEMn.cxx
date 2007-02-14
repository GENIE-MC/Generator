//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - March 30, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "Algorithm/AlgConfigPool.h"
#include "BaryonResonance/BaryonResUtils.h"
#include "Conventions/Constants.h"
#include "ReinSeghal/RSHelicityAmplModelEMn.h"
#include "ReinSeghal/RSHelicityAmpl.h"
#include "Messenger/Messenger.h"

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
RSHelicityAmpl * RSHelicityAmplModelEMn::Compute(
                                       Resonance_t res, const FKR & fkr) const
{
  RSHelicityAmpl * hampl = new RSHelicityAmpl;

  switch(res) {

   case (kP33_1232) :
   {
     hampl->fPlus1  =  kSqrt2 * fkr.R;
     hampl->fPlus3  =  kSqrt6 * fkr.R;
     hampl->fMinus1 = -1 * hampl->fPlus1;
     hampl->fMinus3 = -1 * hampl->fPlus3;
     hampl->f0Minus =  0.;
     hampl->f0Plus  =  0.;
     break;
   }
   case (kS11_1535) :
   {
     hampl->fPlus1  =  kSqrt3   * fkr.T + k1_Sqrt6 * fkr.Lamda * fkr.R;
     hampl->f0Minus =  kSqrt3_2 * fkr.Lamda * fkr.S;
     hampl->fMinus1 = -1 * hampl->fPlus1;
     hampl->f0Plus  = -1 * hampl->f0Minus;
     hampl->fMinus3 =  0.;
     hampl->fPlus3  =  0.;

     break;
   }
   case (kD13_1520) :
   {
     hampl->fMinus1 = -kSqrt3_2 * fkr.T + k1_Sqrt3 * fkr.Lamda * fkr.R;
     hampl->fMinus3 = -k3_Sqrt2 * fkr.T;
     hampl->f0Minus =  kSqrt3 * fkr.Lamda * fkr.S;
     hampl->fPlus1  =  hampl->fMinus1;
     hampl->fPlus3  =  hampl->fMinus3;
     hampl->f0Plus  =  hampl->f0Minus;
     break;
   }
   case (kS11_1650) :
   {
     hampl->fPlus1  =  k1_Sqrt6 * fkr.Lamda * fkr.R;
     hampl->fMinus1 = -1 * hampl->fPlus1;
     hampl->fMinus3 =  0.;
     hampl->fPlus3  =  0.;
     hampl->f0Minus =  0.;
     hampl->f0Plus  =  0.;
     break;
   }
   case (kD13_1700) :
   {
     double LR = fkr.Lamda * fkr.R;

     hampl->fMinus1 = -(1./kSqrt30) * LR;
     hampl->fMinus3 = -(3./kSqrt10) * LR;
     hampl->fPlus1  =  hampl->fMinus1;
     hampl->fPlus3  =  hampl->fMinus3;
     hampl->f0Minus =  0.;
     hampl->f0Plus  =  0.;
     break;
   }
   case (kD15_1675) :
   {
     double LR = fkr.Lamda * fkr.R;

     hampl->fMinus1 = kSqrt3_10 * LR;
     hampl->fMinus3 = kSqrt3_5  * LR;
     hampl->fPlus1  = -1 * hampl->fMinus1;
     hampl->fPlus3  = -1 * hampl->fMinus3;
     hampl->f0Minus =  0.;
     hampl->f0Plus  =  0.;
     break;
   }
   case (kS31_1620) :
   {
     hampl->fMinus1 =  kSqrt3 * fkr.T - k1_Sqrt6 * fkr.Lamda * fkr.R;
     hampl->f0Minus = -kSqrt3_2 * fkr.Lamda * fkr.S;
     hampl->fPlus1  = -1. * hampl->fMinus1;
     hampl->f0Plus  = -1. * hampl->f0Minus;
     hampl->fMinus3 = 0.;
     hampl->fPlus3  = 0.;
     break;
   }
   case (kD33_1700) :
   {
     hampl->fMinus1 =  kSqrt3_2 * fkr.T + k1_Sqrt3 * fkr.Lamda * fkr.R;
     hampl->fMinus3 =  k3_Sqrt2 * fkr.T;
     hampl->f0Minus = -kSqrt3 * fkr.Lamda * fkr.S;
     hampl->fPlus1  =  hampl->fMinus1;
     hampl->fPlus3  =  hampl->fMinus3;
     hampl->f0Plus  =  hampl->f0Minus;
     break;
   }
   case (kP11_1440) :
   {
     hampl->fMinus1 = k1_Sqrt3 * TMath::Power(fkr.Lamda, 2) * fkr.R;
     hampl->fPlus1  = hampl->fMinus1;
     hampl->fMinus3 =  0.;
     hampl->fPlus3  =  0.;
     hampl->f0Minus =  0.;
     hampl->f0Plus  =  0.;
     break;
   }
   case (kP33_1600) :
   {
     double L2R  = TMath::Power(fkr.Lamda, 2) * fkr.R;

     hampl->fMinus1 = k1_Sqrt6 * L2R;
     hampl->fMinus3 = k1_Sqrt2 * L2R;
     hampl->fPlus1  = -1. * hampl->fMinus1;
     hampl->fPlus3  = -1. * hampl->fMinus3;
     hampl->f0Minus = 0.;
     hampl->f0Plus  = 0.;
     break;
   }
   case (kP13_1720) :
   {
     hampl->fMinus1 = k2_Sqrt15 * TMath::Power(fkr.Lamda, 2) * fkr.R;
     hampl->fPlus1  = -1 * hampl->fMinus1;
     hampl->fMinus3 =  0.;
     hampl->fPlus3  =  0.;
     hampl->f0Minus =  0.;
     hampl->f0Plus  =  0.;
     break;
   }
   case (kF15_1680) :
   {
     hampl->fMinus1 =  -kSqrt2_5 * TMath::Power(fkr.Lamda, 2) * fkr.R;
     hampl->fPlus1  =  hampl->fMinus1;
     hampl->fMinus3 =  0.;
     hampl->fPlus3  =  0.;
     hampl->f0Minus =  0.;
     hampl->f0Plus  =  0.;
     break;
   }
   case (kP31_1910) :
   {
     hampl->fMinus1 =  -k1_Sqrt15 * TMath::Power(fkr.Lamda, 2) * fkr.R;
     hampl->fPlus1  =  hampl->fMinus1;
     hampl->fMinus3 =  0.;
     hampl->fPlus3  =  0.;
     hampl->f0Minus =  0.;
     hampl->f0Plus  =  0.;
     break;
   }
   case (kP33_1920) :
   {
     double L2R  = TMath::Power(fkr.Lamda, 2) * fkr.R;

     hampl->fMinus1 =  k1_Sqrt15 * L2R;
     hampl->fMinus3 = -k1_Sqrt5  * L2R;
     hampl->fPlus1  = -1.* hampl->fMinus1;
     hampl->fPlus3  = -1.* hampl->fMinus3;
     hampl->f0Minus =  0.;
     hampl->f0Plus  =  0.;
     break;
   }
   case (kF35_1905) :
   {
     double L2R  = TMath::Power(fkr.Lamda, 2) * fkr.R;

     hampl->fMinus1 = k1_Sqrt35  * L2R;
     hampl->fMinus3 = kSqrt18_35 * L2R;
     hampl->fPlus1  = hampl->fMinus1;
     hampl->fPlus3  = hampl->fMinus3;
     hampl->f0Minus = 0.;
     hampl->f0Plus  = 0.;
     break;
   }
   case (kF37_1950) :
   {
     double L2R  = TMath::Power(fkr.Lamda, 2) * fkr.R;

     hampl->fMinus1 = -kSqrt6_35 * L2R;
     hampl->fMinus3 = -kSqrt2_7  * L2R;
     hampl->fPlus1  = -1. * hampl->fMinus1;
     hampl->fPlus3  = -1. * hampl->fMinus3;
     hampl->f0Minus = 0.;
     hampl->f0Plus  = 0.;
     break;
   }
   case (kP11_1710) :
   {
     double L2 = TMath::Power(fkr.Lamda, 2);

     hampl->fMinus1 = -k1_Sqrt24 * L2 * fkr.R;
     hampl->f0Minus = -kSqrt3_8  * L2 * fkr.S;
     hampl->fPlus1  = hampl->fMinus1;
     hampl->f0Plus  = hampl->f0Minus;
     hampl->fMinus3 = 0.;
     hampl->fPlus3  = 0.;

     break;
   }
   case (kF17_1970) :
   {
     double L2R = TMath::Power(fkr.Lamda, 2) * fkr.R;

     hampl->fMinus1 = kSqrt3_35 * L2R;
     hampl->fPlus1  = -1 * hampl->fMinus1;
     hampl->fMinus3 = k1_Sqrt7  * L2R;
     hampl->fPlus3  = -1 * hampl->fMinus3;
     hampl->f0Minus =  0.;
     hampl->f0Plus  =  0.;
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
