//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - March 30, 2005

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
     hampl->fMinus1 = -kSqrt2 * fkr.R();
     hampl->fPlus1  = -1. * hampl->fMinus1;
     hampl->fMinus3 = -kSqrt6 * fkr.R();
     hampl->fPlus3  = -1. * hampl->fMinus3;
     hampl->f0Minus =  0.;
     hampl->f0Plus  =  0.;
     break;

   case (kS11_1535) :
     hampl->fMinus1 = kSqrt3* fkr.T() + (kSqrt3/kSqrt2) * fkr.LR();
     hampl->fPlus1  = -1. * hampl->fMinus1;
     hampl->fMinus3 =  0.;
     hampl->fPlus3  =  0.;
     hampl->f0Minus = -(kSqrt3/kSqrt2) * fkr.LS();
     hampl->f0Plus  = -1. * hampl->f0Minus;
     break;

   case (kD13_1520) :
     hampl->fMinus1 = (kSqrt3/kSqrt2) * fkr.T() - kSqrt3 * fkr.LR();
     hampl->fPlus1  = hampl->fMinus1;
     hampl->fMinus3 = (3./kSqrt2) * fkr.T();
     hampl->fPlus3  = hampl->fMinus3;
     hampl->f0Minus = -kSqrt3 * fkr.LS();
     hampl->f0Plus  = hampl->f0Minus;
     break;

   case (kS11_1650) :
     hampl->fMinus1 = 0.;
     hampl->fPlus1  = 0.;
     hampl->fMinus3 = 0.;
     hampl->fPlus3  = 0.;
     hampl->f0Minus = 0.;
     hampl->f0Plus  = 0.;
     break;

   case (kD13_1700) :
     hampl->fMinus1 = 0.;
     hampl->fPlus1  = 0.;
     hampl->fMinus3 = 0.;
     hampl->fPlus3  = 0.;
     hampl->f0Minus = 0.;
     hampl->f0Plus  = 0.;
     break;

   case (kD15_1675) :
     hampl->fMinus1 = 0.;
     hampl->fPlus1  = 0.;
     hampl->fMinus3 = 0.;
     hampl->fPlus3  = 0.;
     hampl->f0Minus = 0.;
     hampl->f0Plus  = 0.;
     break;

   case (kS31_1620) :
     hampl->fMinus1 =  kSqrt3* fkr.T() - (1./kSqrt6)* fkr.LR();
     hampl->fPlus1  = -1. * hampl->fMinus1;
     hampl->fMinus3 = 0.;
     hampl->fPlus3  = 0.;
     hampl->f0Minus = -(kSqrt3/kSqrt2)* fkr.LS();
     hampl->f0Plus  = -1. * hampl->f0Minus;
     break;

   case (kD33_1700) :
     hampl->fMinus1 = (kSqrt3/kSqrt2)* fkr.T() + (1./kSqrt3) * fkr.LR();
     hampl->fPlus1  = hampl->fMinus1;
     hampl->fMinus3 = (3./kSqrt2) * fkr.T();
     hampl->fPlus3  = hampl->fMinus3;
     hampl->f0Minus = -kSqrt3 * fkr.LS();
     hampl->f0Plus  = hampl->f0Minus;
     break;

   case (kP11_1440) :
     hampl->fMinus1 = -0.5*kSqrt3 * fkr.L2R();
     hampl->fPlus1  = -1 * hampl->fMinus1;
     hampl->fMinus3 = 0.;
     hampl->fPlus3  = 0.;
     hampl->f0Minus = 0.5*kSqrt3 * fkr.L2S();
     hampl->f0Plus  = -1. * hampl->f0Minus;
     break;

   case (kP33_1600) :
     hampl->fMinus1 = (1./kSqrt6) * fkr.L2R();
     hampl->fPlus1  = -1. * hampl->fMinus1;
     hampl->fMinus3 = (1./kSqrt2) * fkr.L2R();
     hampl->fPlus3  = -1. * hampl->fMinus3;
     hampl->f0Minus = 0.;
     hampl->f0Plus  = 0.;
     break;

   case (kP13_1720) :
     hampl->fMinus1 = -(kSqrt27/kSqrt10) * fkr.LT() - (kSqrt3/kSqrt5) * fkr.L2R();
     hampl->fPlus1  =  -1. * hampl->fMinus1;
     hampl->fMinus3 = (3./kSqrt10) * fkr.LT();
     hampl->fPlus3  = -1. * hampl->fMinus3;
     hampl->f0Minus = (kSqrt3/kSqrt5) * fkr.L2S();
     hampl->f0Plus  = -1. * hampl->f0Minus;
     break;

   case (kF15_1680) :
     hampl->fMinus1 = -(3./kSqrt5) * fkr.LT() + (kSqrt3/kSqrt10) * fkr.L2R();
     hampl->fPlus1  = -(3./kSqrt5) * fkr.LT() + (3./kSqrt10) * fkr.L2R();
     hampl->fMinus3 = -1*(kSqrt18/kSqrt5) * fkr.LT();
     hampl->fPlus3  = hampl->fMinus3;
     hampl->f0Minus = (3./kSqrt10) * fkr.L2S();
     hampl->f0Plus  = hampl->f0Minus;
     break;

   case (kP31_1910) :
     hampl->fMinus1 = -(1./kSqrt15) * fkr.L2R();
     hampl->fPlus1  = hampl->fMinus1;
     hampl->fMinus3 = 0.;
     hampl->fPlus3  = 0.;
     hampl->f0Minus = 0.;
     hampl->f0Plus  = 0.;
     break;

   case (kP33_1920) :
     hampl->fMinus1 =  (1./kSqrt15) * fkr.L2R();
     hampl->fPlus1  = -1. * hampl->fMinus1;
     hampl->fMinus3 = -(1./kSqrt5 ) * fkr.L2R();
     hampl->fPlus3  = -1. * hampl->fMinus3;
     hampl->f0Minus = 0.;
     hampl->f0Plus  = 0.;
     break;

   case (kF35_1905) :
     hampl->fMinus1 = (1./kSqrt35) * fkr.L2R();
     hampl->fPlus1  = hampl->fMinus1;
     hampl->fMinus3 = (kSqrt18/kSqrt35) * fkr.L2R();
     hampl->fPlus3  = hampl->fMinus3;
     hampl->f0Minus = 0.;
     hampl->f0Plus  = 0.;
     break;

   case (kF37_1950) :
     hampl->fMinus1 = -(kSqrt6/kSqrt35) * fkr.L2R();
     hampl->fPlus1  = -1. * hampl->fMinus1;
     hampl->fMinus3 = -(kSqrt2/kSqrt7) * fkr.L2R();
     hampl->fPlus3  = -1. * hampl->fMinus3;
     hampl->f0Minus = 0.;
     hampl->f0Plus  = 0.;
     break;

   case (kP11_1710) :
     hampl->fMinus1 = (kSqrt3/kSqrt8) * fkr.L2R();
     hampl->fPlus1  = -1. * hampl->fMinus1;
     hampl->fMinus3 = 0.;
     hampl->fPlus3  = 0.;
     hampl->f0Minus = (kSqrt3/kSqrt8) * fkr.L2S();
     hampl->f0Plus  = hampl->f0Minus;
     break;

   case (kF17_1970) :
     hampl->fMinus1 = 0.;
     hampl->fPlus1  = 0.;
     hampl->fMinus3 = 0.;
     hampl->fPlus3  = 0.;
     hampl->f0Minus = 0.;
     hampl->f0Plus  = 0.;
     break;

   default:
     LOG("RSHAmpl", pWARN) << "*** UNRECOGNIZED RESONANCE!";
     delete hampl;
     hampl = 0;
     break;
   }
  return hampl;
}
//____________________________________________________________________________


