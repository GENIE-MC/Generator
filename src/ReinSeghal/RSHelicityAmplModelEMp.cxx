//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory - March 30, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "BaryonResonance/BaryonResUtils.h"
#include "Conventions/Constants.h"
#include "ReinSeghal/RSHelicityAmplModelEMp.h"
#include "ReinSeghal/RSHelicityAmpl.h"
#include "Messenger/Messenger.h"

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
RSHelicityAmpl * RSHelicityAmplModelEMp::Compute(
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
     hampl->fMinus1 =  kSqrt3 * fkr.T + kSqrt3_2 * fkr.Lamda * fkr.R;
     hampl->f0Minus = -kSqrt3_2 * fkr.Lamda * fkr.S;
     hampl->fPlus1  = -1. * hampl->fMinus1;
     hampl->f0Plus  = -1. * hampl->f0Minus;
     hampl->fMinus3 =  0.;
     hampl->fPlus3  =  0.;
     break;
   }
   case (kD13_1520) :
   {
     hampl->fMinus1 =  kSqrt3_2 * fkr.T - kSqrt3 * fkr.Lamda * fkr.R;
     hampl->fMinus3 =  k3_Sqrt2 * fkr.T;
     hampl->f0Minus = -kSqrt3 * fkr.Lamda * fkr.S;
     hampl->fPlus1  =  hampl->fMinus1;
     hampl->fPlus3  =  hampl->fMinus3;
     hampl->f0Plus  =  hampl->f0Minus;
     break;
   }
   case (kS11_1650) :
   {
     hampl->fMinus1 = 0.;
     hampl->fPlus1  = 0.;
     hampl->fMinus3 = 0.;
     hampl->fPlus3  = 0.;
     hampl->f0Minus = 0.;
     hampl->f0Plus  = 0.;
     break;
   }
   case (kD13_1700) :
   {
     hampl->fMinus1 = 0.;
     hampl->fPlus1  = 0.;
     hampl->fMinus3 = 0.;
     hampl->fPlus3  = 0.;
     hampl->f0Minus = 0.;
     hampl->f0Plus  = 0.;
     break;
   }
   case (kD15_1675) :
   {
     hampl->fMinus1 = 0.;
     hampl->fPlus1  = 0.;
     hampl->fMinus3 = 0.;
     hampl->fPlus3  = 0.;
     hampl->f0Minus = 0.;
     hampl->f0Plus  = 0.;
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
     double L2  = TMath::Power(fkr.Lamda, 2);

     hampl->fMinus1 = -0.5*kSqrt3 * L2 * fkr.R;
     hampl->fPlus1  =  hampl->fMinus1;
     hampl->fMinus3 =  0.;
     hampl->fPlus3  =  0.;
     hampl->f0Minus = -0.5*kSqrt3 * L2 * fkr.S;
     hampl->f0Plus  =  hampl->f0Minus;
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
     double L2  = TMath::Power(fkr.Lamda, 2);
     double LT  = fkr.Lamda * fkr.T;

     hampl->fMinus1 = -kSqrt27_10 * LT - kSqrt3_5 * L2 * fkr.R;
     hampl->fMinus3 =  k3_Sqrt10 * LT;
     hampl->f0Minus =  kSqrt3_5  * L2 * fkr.S;
     hampl->fPlus1  = -1. * hampl->fMinus1;
     hampl->fPlus3  = -1. * hampl->fMinus3;
     hampl->f0Plus  = -1. * hampl->f0Minus;
     break;
   }
   case (kF15_1680) :
   {
     double L2  = TMath::Power(fkr.Lamda, 2);
     double LT  = fkr.Lamda * fkr.T;

     hampl->fMinus1 =  -k3_Sqrt5  * LT + k3_Sqrt10 * L2 * fkr.R;
     hampl->fMinus3 =  -kSqrt18_5 * LT;
     hampl->f0Minus =   k3_Sqrt10 * L2 * fkr.S;
     hampl->fPlus1  =  hampl->fMinus1;
     hampl->fPlus3  =  hampl->fMinus3;
     hampl->f0Plus  =  hampl->f0Minus;
     break;
   }
   case (kP31_1910) :
   {
     hampl->fMinus1 = -k1_Sqrt15 * TMath::Power(fkr.Lamda, 2) * fkr.R;
     hampl->fPlus1  = hampl->fMinus1;
     hampl->fMinus3 = 0.;
     hampl->fPlus3  = 0.;
     hampl->f0Minus = 0.;
     hampl->f0Plus  = 0.;
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
     double L2  = TMath::Power(fkr.Lamda, 2);

     hampl->fMinus1 = kSqrt3_8 * L2 * fkr.R;
     hampl->f0Minus = kSqrt3_8 * L2 * fkr.S;
     hampl->fPlus1  = hampl->fMinus1;
     hampl->f0Plus  = hampl->f0Minus;
     hampl->fMinus3 = 0.;
     hampl->fPlus3  = 0.;
     break;
   }
   case (kF17_1970) :
   {
     hampl->fMinus1 = 0.;
     hampl->fPlus1  = 0.;
     hampl->fMinus3 = 0.;
     hampl->fPlus3  = 0.;
     hampl->f0Minus = 0.;
     hampl->f0Plus  = 0.;
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


