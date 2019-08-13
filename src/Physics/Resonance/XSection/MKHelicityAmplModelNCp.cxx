//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Authors: Igor Kakorin <kakorin@jinr.ru> Joint Institute for Nuclear Research
          based on code by 
          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/Resonance/XSection/MKHelicityAmplModelNCp.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
MKHelicityAmplModelNCp::MKHelicityAmplModelNCp() :
MKHelicityAmplModelI("genie::MKHelicityAmplModelNCp")
{

}
//____________________________________________________________________________
MKHelicityAmplModelNCp::MKHelicityAmplModelNCp(string config) :
MKHelicityAmplModelI("genie::MKHelicityAmplModelNCp", config)
{

}
//____________________________________________________________________________
MKHelicityAmplModelNCp::~MKHelicityAmplModelNCp()
{

}
//____________________________________________________________________________
const MKHelicityAmpl & MKHelicityAmplModelNCp::Compute(Resonance_t res, const FKR_MK & fkr) const
{
  switch(res) {

   case (kP33_1232) :
   {
     fAmpl.fVMinus3    =-kSqrt6*fkr.Rv*(-1. + 2.*fSin28w);
     fAmpl.fVMinus1    =-kSqrt2*fkr.Rv*(-1. + 2.*fSin28w);
     fAmpl.fVPlus1     =-fAmpl.fVMinus1;
     fAmpl.fVPlus3     =-fAmpl.fVMinus3;
     fAmpl.fV0LPlus    = 0.;
     fAmpl.fV0LMinus   = 0.;
     fAmpl.fV0RPlus    = 0.;
     fAmpl.fV0RMinus   = 0.;
     
     fAmpl.fAMinus3    = kSqrt6*fkr.Ra;
     fAmpl.fAMinus1    = kSqrt2*fkr.Ra;
     fAmpl.fAPlus1     = fAmpl.fAMinus1;
     fAmpl.fAPlus3     = fAmpl.fAMinus3; 
     fAmpl.fA0LPlus    =-2.*kSqrt2*fkr.Cminus;
     fAmpl.fA0LMinus   = fAmpl.fA0LPlus;
     fAmpl.fA0RPlus    = 0.;
     fAmpl.fA0RMinus   = 0.;
     break;
   }
   case (kS11_1535) :
   {
     fAmpl.fVMinus3    = 0.;
     fAmpl.fVMinus1    = kSqrt3*fkr.Tv*(-1. + 2.*fSin28w) + kSqrt2_3*fkr.Lamda*fkr.Rv*(-1. + 3.*fSin28w);
     fAmpl.fVPlus1     =-fAmpl.fVMinus1;
     fAmpl.fVPlus3     = 0.;
     fAmpl.fV0LPlus    = kSqrt3_2*fkr.Lamda*fkr.Sminus*(1. - 2.*fSin28w);
     fAmpl.fV0LMinus   =-fAmpl.fV0LPlus;
     fAmpl.fV0RPlus    = 0.;
     fAmpl.fV0RMinus   = 0.;
     
     fAmpl.fAMinus3    = 0.;
     fAmpl.fAMinus1    =-kSqrt3*fkr.Ta - kSqrt2_3*fkr.Lamda*fkr.Ra;
     fAmpl.fAPlus1     = fAmpl.fAMinus1;
     fAmpl.fAPlus3     = 0.;
     fAmpl.fA0LPlus    =-kSqrt2_3*(fkr.Lamda*fkr.Cminus - 3.*fkr.Bminus);
     fAmpl.fA0LMinus   = fAmpl.fA0LPlus;
     fAmpl.fA0RPlus    = 0.;
     fAmpl.fA0RMinus   = 0.;
     break;
   }
   case (kD13_1520) :
   {
     fAmpl.fVMinus3    = k3_Sqrt2*fkr.Tv*(-1. + 2.*fSin28w);
     fAmpl.fVMinus1    = kSqrt3_2*fkr.Tv*(-1. + 2.*fSin28w) - 2./kSqrt3*fkr.Lamda*fkr.Rv*(-1. + 3.*fSin28w);
     fAmpl.fVPlus1     = fAmpl.fVMinus1;
     fAmpl.fVPlus3     = fAmpl.fVMinus3;
     fAmpl.fV0LPlus    =-kSqrt3*fkr.Lamda*fkr.Sminus*(1. - 2.*fSin28w);
     fAmpl.fV0LMinus   = fAmpl.fV0LPlus;
     fAmpl.fV0RPlus    = 0.;
     fAmpl.fV0RMinus   = 0.;
     
     fAmpl.fAMinus3    =-k3_Sqrt2*fkr.Ta;
     fAmpl.fAMinus1    =-kSqrt3_2*fkr.Ta + 2./kSqrt3*fkr.Lamda*fkr.Ra;
     fAmpl.fAPlus1     =-fAmpl.fAMinus1;
     fAmpl.fAPlus3     =-fAmpl.fAMinus3;
     fAmpl.fA0LPlus    = 2./kSqrt3*fkr.Lamda*fkr.Cminus;
     fAmpl.fA0LMinus   =-fAmpl.fA0LPlus;
     fAmpl.fA0RPlus    = 0.;
     fAmpl.fA0RMinus   = 0.;
     break;
   }
   case (kS11_1650) :
   {
     fAmpl.fVMinus3    = 0.;
     fAmpl.fVMinus1    =-k1_Sqrt24*fkr.Lamda*fkr.Rv;
     fAmpl.fVPlus1     =-fAmpl.fVMinus1;
     fAmpl.fVPlus3     = 0.;
     fAmpl.fV0LPlus    = 0.;
     fAmpl.fV0LMinus   = 0.;
     fAmpl.fV0RPlus    = 0.;
     fAmpl.fV0RMinus   = 0.;
     
     fAmpl.fAMinus3    = 0.;
     fAmpl.fAMinus1    =-k1_Sqrt24*fkr.Lamda*fkr.Ra;
     fAmpl.fAPlus1     = fAmpl.fAMinus1;
     fAmpl.fAPlus3     = 0.;
     fAmpl.fA0LPlus    = k1_Sqrt6*(fkr.Lamda*fkr.Cminus - 3.*fkr.Bminus);
     fAmpl.fA0LMinus   = fAmpl.fA0LPlus;
     fAmpl.fA0RPlus    = 0.;
     fAmpl.fA0RMinus   = 0.;
     break;
   }
   case (kD13_1700) :
   {
     fAmpl.fVMinus3    =-k3_Sqrt40*fkr.Lamda*fkr.Rv;
     fAmpl.fVMinus1    =-k1_Sqrt120*fkr.Lamda*fkr.Rv;
     fAmpl.fVPlus1     = fAmpl.fVMinus1;
     fAmpl.fVPlus3     = fAmpl.fVMinus3;
     fAmpl.fV0LPlus    = 0.;
     fAmpl.fV0LMinus   = 0.;
     fAmpl.fV0RPlus    = 0.;
     fAmpl.fV0RMinus   = 0.;
     
     fAmpl.fAMinus3    =-k3_Sqrt40*fkr.Lamda*fkr.Ra;
     fAmpl.fAMinus1    =-k1_Sqrt120*fkr.Lamda*fkr.Ra;
     fAmpl.fAPlus1     =-fAmpl.fAMinus1;
     fAmpl.fAPlus3     =-fAmpl.fAMinus3;
     fAmpl.fA0LPlus    = k1_Sqrt30*fkr.Lamda*fkr.Cminus;
     fAmpl.fA0LMinus   =-fAmpl.fA0LPlus;
     fAmpl.fA0RPlus    = 0.;
     fAmpl.fA0RMinus   = 0.;
     break;
   }
   case (kD15_1675) :
   {
     fAmpl.fVMinus3    = kSqrt3_20*fkr.Lamda*fkr.Rv;
     fAmpl.fVMinus1    = kSqrt3_40*fkr.Lamda*fkr.Rv;
     fAmpl.fVPlus1     =-fAmpl.fVMinus1;
     fAmpl.fVPlus3     =-fAmpl.fVMinus3;
     fAmpl.fV0LPlus    = 0.;
     fAmpl.fV0LMinus   = 0.;
     fAmpl.fV0RPlus    = 0.;
     fAmpl.fV0RMinus   = 0.;
                       
     fAmpl.fAMinus3    = kSqrt3_20*fkr.Lamda*fkr.Ra;
     fAmpl.fAMinus1    = kSqrt3_40*fkr.Lamda*fkr.Ra;
     fAmpl.fAPlus1     = fAmpl.fAMinus1;
     fAmpl.fAPlus3     = fAmpl.fAMinus3;
     fAmpl.fA0LPlus    =-kSqrt3_10*fkr.Lamda*fkr.Cminus;
     fAmpl.fA0LMinus   = fAmpl.fA0LPlus;
     fAmpl.fA0RPlus    = 0.;
     fAmpl.fA0RMinus   = 0.;
     break;
   }
   case (kS31_1620) :
   {
     fAmpl.fVMinus3    = 0.;
     fAmpl.fVMinus1    = (kSqrt3*fkr.Tv - k1_Sqrt6*fkr.Rv*fkr.Lamda)*(-1. + 2.*fSin28w);
     fAmpl.fVPlus1     =-fAmpl.fVMinus1;
     fAmpl.fVPlus3     = 0.;
     fAmpl.fV0LPlus    =-kSqrt3_2*fkr.Lamda*fkr.Sminus*(-1. + 2.*fSin28w);
     fAmpl.fV0LMinus   =-fAmpl.fV0LPlus;
     fAmpl.fV0RPlus    = 0.;
     fAmpl.fV0RMinus   = 0.;
                       
     fAmpl.fAMinus3    = 0.;
     fAmpl.fAMinus1    =-kSqrt3*fkr.Ta + k1_Sqrt6*fkr.Lamda*fkr.Ra;
     fAmpl.fAPlus1     = fAmpl.fAMinus1;
     fAmpl.fAPlus3     = 0.;
     fAmpl.fA0LPlus    = k1_Sqrt6*(fkr.Lamda*fkr.Cminus - 3.*fkr.Bminus);
     fAmpl.fA0LMinus   = fAmpl.fA0LPlus;
     fAmpl.fA0RPlus    = 0.;
     fAmpl.fA0RMinus   = 0.;
     break;
   }
   case (kD33_1700) :
   {
     fAmpl.fVMinus3    = k3_Sqrt2*fkr.Tv*(-1. + 2.*fSin28w);
     fAmpl.fVMinus1    = (kSqrt3_2*fkr.Tv + k1_Sqrt3*fkr.Lamda*fkr.Rv)*(-1. + 2.*fSin28w);
     fAmpl.fVPlus1     = fAmpl.fVMinus1;
     fAmpl.fVPlus3     = fAmpl.fVMinus3;
     fAmpl.fV0LPlus    = kSqrt3*fkr.Lamda*fkr.Sminus*(-1. + 2.*fSin28w);
     fAmpl.fV0LMinus   = fAmpl.fV0LPlus;
     fAmpl.fV0RPlus    = 0.;
     fAmpl.fV0RMinus   = 0.;
                       
     fAmpl.fAMinus3    =-k3_Sqrt2*fkr.Ta;
     fAmpl.fAMinus1    =-kSqrt3_2*fkr.Ta - k1_Sqrt3*fkr.Lamda*fkr.Ra;
     fAmpl.fAPlus1     =-fAmpl.fAMinus1;
     fAmpl.fAPlus3     =-fAmpl.fAMinus3;
     fAmpl.fA0LPlus    =-k1_Sqrt3*fkr.Lamda*fkr.Cminus;
     fAmpl.fA0LMinus   =-fAmpl.fA0LPlus;
     fAmpl.fA0RPlus    = 0.;
     fAmpl.fA0RMinus   = 0.;
     break;
   }
   case (kP11_1440) :
   {
     fAmpl.fVMinus3    = 0.;
     fAmpl.fVMinus1    =-kSqrt3*fkr.Lamda*fkr.Lamda*fkr.Rv*(-5./12. + fSin28w);
     fAmpl.fVPlus1     = fAmpl.fVMinus1;
     fAmpl.fVPlus3     = 0.;
     fAmpl.fV0LPlus    =-kSqrt3*fkr.Lamda*fkr.Lamda*fkr.Sminus*(0.25 - fSin28w);
     fAmpl.fV0LMinus   = fAmpl.fV0LPlus;
     fAmpl.fV0RPlus    = 0.;
     fAmpl.fV0RMinus   = 0.;
                       
     fAmpl.fAMinus3    = 0.;
     fAmpl.fAMinus1    = 5./12.*kSqrt3*fkr.Lamda*fkr.Lamda*fkr.Ra;
     fAmpl.fAPlus1     =-fAmpl.fAMinus1;
     fAmpl.fAPlus3     = 0.;
     fAmpl.fA0LPlus    = 5./12.*kSqrt3*fkr.Lamda*(fkr.Lamda*fkr.Cminus - 2.*fkr.Bminus);
     fAmpl.fA0LMinus   =-fAmpl.fA0LPlus; 
     fAmpl.fA0RPlus    = 0.;
     fAmpl.fA0RMinus   = 0.;
     break;
   }
   case (kP33_1600) :
   {
     fAmpl.fVMinus3    = k1_Sqrt2*fkr.Lamda*fkr.Lamda*fkr.Rv*(-1. + 2.*fSin28w);
     fAmpl.fVMinus1    = k1_Sqrt6*fkr.Lamda*fkr.Lamda*fkr.Rv*(-1. + 2.*fSin28w);
     fAmpl.fVPlus1     =-fAmpl.fVMinus1;
     fAmpl.fVPlus3     =-fAmpl.fVMinus3;
     fAmpl.fV0LPlus    = 0.;
     fAmpl.fV0LMinus   = 0.;
     fAmpl.fV0RPlus    = 0.;
     fAmpl.fV0RMinus   = 0.;
                        
     fAmpl.fAMinus3    =-k1_Sqrt2*fkr.Lamda*fkr.Lamda*fkr.Ra;
     fAmpl.fAMinus1    =-k1_Sqrt6*fkr.Lamda*fkr.Lamda*fkr.Ra;
     fAmpl.fAPlus1     = fAmpl.fAMinus1;
     fAmpl.fAPlus3     = fAmpl.fAMinus3;
     fAmpl.fA0LPlus    = kSqrt2_3*fkr.Lamda*(fkr.Lamda*fkr.Cminus - 2.*fkr.Bminus);
     fAmpl.fA0LMinus   = fAmpl.fA0LPlus;
     fAmpl.fA0RPlus    = 0.;
     fAmpl.fA0RMinus   = 0.;
     break;
   }
   case (kP13_1720) :
   {     
     fAmpl.fVMinus3    = k3_Sqrt40*fkr.Lamda*fkr.Tv*(-1. + 4.*fSin28w);
     fAmpl.fVMinus1    =-kSqrt27_40*fkr.Lamda*fkr.Tv*(-1. + 4.*fSin28w) - kSqrt5_12*fkr.Lamda*fkr.Lamda*fkr.Rv*(-1. + 12./5.*fSin28w);
     fAmpl.fVPlus1     =-fAmpl.fVMinus1; 
     fAmpl.fVPlus3     =-fAmpl.fVMinus3;
     fAmpl.fV0LPlus    = kSqrt3_20*fkr.Lamda*fkr.Lamda*fkr.Sminus*(-1. + 4.*fSin28w);
     fAmpl.fV0LMinus   =-fAmpl.fV0LPlus;
     fAmpl.fV0RPlus    = 0.;
     fAmpl.fV0RMinus   = 0.;
                        
     fAmpl.fAMinus3    =-k3_Sqrt40*fkr.Lamda*fkr.Ta;
     fAmpl.fAMinus1    = kSqrt27_40*fkr.Lamda*fkr.Ta + kSqrt5_12*fkr.Lamda*fkr.Lamda*fkr.Ra;
     fAmpl.fAPlus1     = fAmpl.fAMinus1; 
     fAmpl.fAPlus3     = fAmpl.fAMinus3;
     fAmpl.fA0LPlus    = kSqrt5_12*fkr.Lamda*(fkr.Lamda*fkr.Cminus - 5.*fkr.Bminus);
     fAmpl.fA0LMinus   = fAmpl.fA0LPlus;
     fAmpl.fA0RPlus    = 0.;
     fAmpl.fA0RMinus   = 0.;
     break;
   }
   case (kF15_1680) :
   {
     fAmpl.fVMinus3   =-kSqrt18_20*fkr.Lamda*fkr.Tv*(-1. + 4.*fSin28w);
     fAmpl.fVMinus1   =-k3_Sqrt20*fkr.Lamda*fkr.Tv*(-1. + 4.*fSin28w) + kSqrt5_8*fkr.Lamda*fkr.Lamda*fkr.Rv*(-1. + 12./5.*fSin28w);
     fAmpl.fVPlus1    = fAmpl.fVMinus1; 
     fAmpl.fVPlus3    = fAmpl.fVMinus3;
     fAmpl.fV0LPlus   = k3_Sqrt40*fkr.Lamda*fkr.Lamda*fkr.Sminus*(1. - 4.*fSin28w);
     fAmpl.fV0LMinus  = fAmpl.fV0LPlus;
     fAmpl.fV0RPlus   = 0.;
     fAmpl.fV0RMinus  = 0.;
                       
     fAmpl.fAMinus3   = kSqrt18_20*fkr.Lamda*fkr.Ta;
     fAmpl.fAMinus1   = k3_Sqrt20*fkr.Lamda*fkr.Ta - kSqrt5_8*fkr.Lamda*fkr.Lamda*fkr.Ra;
     fAmpl.fAPlus1    =-fAmpl.fAMinus1; 
     fAmpl.fAPlus3    =-fAmpl.fAMinus3;
     fAmpl.fA0LPlus   =-kSqrt5_8*fkr.Lamda*fkr.Lamda*fkr.Cminus;
     fAmpl.fA0LMinus  =-fAmpl.fA0LPlus;
     fAmpl.fA0RPlus   = 0.;
     fAmpl.fA0RMinus  = 0.;
     break;
   }
   case (kP31_1910) :
   {     
     fAmpl.fVMinus3   = 0.;
     fAmpl.fVMinus1   =-k1_Sqrt15*fkr.Lamda*fkr.Lamda*fkr.Rv*(-1. + 2.*fSin28w);
     fAmpl.fVPlus1    = fAmpl.fVMinus1;
     fAmpl.fVPlus3    = 0.;
     fAmpl.fV0LPlus   = 0.;
     fAmpl.fV0LMinus  = 0.;
     fAmpl.fV0RPlus   = 0.;
     fAmpl.fV0RMinus  = 0.;
                       
     fAmpl.fAMinus3   = 0.;
     fAmpl.fAMinus1   = k1_Sqrt15*fkr.Lamda*fkr.Lamda*fkr.Ra;
     fAmpl.fAPlus1    =-fAmpl.fAMinus1;
     fAmpl.fAPlus3    = 0.;
     fAmpl.fA0LPlus   =-kSqrt4_15*fkr.Lamda*(fkr.Lamda*fkr.Cminus - 5.*fkr.Bminus);
     fAmpl.fA0LMinus  =-fAmpl.fA0LPlus;
     fAmpl.fA0RPlus   = 0.;
     fAmpl.fA0RMinus  = 0.;
     break;
   }
   case (kP33_1920) :
   {
     fAmpl.fVMinus3   =-k1_Sqrt5*fkr.Lamda*fkr.Lamda*fkr.Rv*(-1. + 2.*fSin28w);
     fAmpl.fVMinus1   = k1_Sqrt15*fkr.Lamda*fkr.Lamda*fkr.Rv*(-1. + 2.*fSin28w);
     fAmpl.fVPlus1    =-fAmpl.fVMinus1;
     fAmpl.fVPlus3    =-fAmpl.fVMinus3;
     fAmpl.fV0LPlus   = 0.;
     fAmpl.fV0LMinus  = 0.;
     fAmpl.fV0RPlus   = 0.;
     fAmpl.fV0RMinus  = 0.;
                       
     fAmpl.fAMinus3   = k1_Sqrt5*fkr.Lamda*fkr.Lamda*fkr.Ra;
     fAmpl.fAMinus1   =-k1_Sqrt15*fkr.Lamda*fkr.Lamda*fkr.Ra;
     fAmpl.fAPlus1    = fAmpl.fAMinus1;
     fAmpl.fAPlus3    = fAmpl.fAMinus3;
     fAmpl.fA0LPlus   = (2./kSqrt15)*fkr.Lamda*(fkr.Lamda*fkr.Cminus - 5.*fkr.Bminus);
     fAmpl.fA0LMinus  = fAmpl.fA0LPlus;
     fAmpl.fA0RPlus   = 0.;
     fAmpl.fA0RMinus  = 0.;
     break;
   }
   case (kF35_1905) :
   {
     fAmpl.fVMinus3   = kSqrt18_35*fkr.Lamda*fkr.Lamda*fkr.Rv*(-1. +2.*fSin28w);
     fAmpl.fVMinus1   = k1_Sqrt35*fkr.Lamda*fkr.Lamda*fkr.Rv*(-1. +2.*fSin28w);
     fAmpl.fVPlus1    = fAmpl.fVMinus1;
     fAmpl.fVPlus3    = fAmpl.fVMinus3;
     fAmpl.fV0LPlus   = 0.;
     fAmpl.fV0LMinus  = 0.;
     fAmpl.fV0RPlus   = 0.;
     fAmpl.fV0RMinus  = 0.;
                       
     fAmpl.fAMinus3   =-kSqrt18_35*fkr.Lamda*fkr.Lamda*fkr.Ra;
     fAmpl.fAMinus1   =-k1_Sqrt35*fkr.Lamda*fkr.Lamda*fkr.Ra;
     fAmpl.fAPlus1    =-fAmpl.fAMinus1;
     fAmpl.fAPlus3    =-fAmpl.fAMinus3;
     fAmpl.fA0LPlus   = k2_Sqrt35*fkr.Lamda*fkr.Lamda*fkr.Cminus;
     fAmpl.fA0LMinus  =-fAmpl.fA0LPlus;
     fAmpl.fA0RPlus   = 0.;
     fAmpl.fA0RMinus  = 0.;
     break;
   }
   case (kF37_1950) :
   {
     fAmpl.fVMinus3   =-kSqrt2_7*fkr.Lamda*fkr.Lamda*fkr.Rv*(-1. + 2.*fSin28w);
     fAmpl.fVMinus1   =-kSqrt6_35*fkr.Lamda*fkr.Lamda*fkr.Rv*(-1. + 2.*fSin28w);
     fAmpl.fVPlus1    =-fAmpl.fVMinus1;
     fAmpl.fVPlus3    =-fAmpl.fVMinus3;
     fAmpl.fV0LPlus   = 0.;
     fAmpl.fV0LMinus  = 0.;
     fAmpl.fV0RPlus   = 0.;
     fAmpl.fV0RMinus  = 0.;
                       
     fAmpl.fAMinus3   = kSqrt2_7*fkr.Lamda*fkr.Lamda*fkr.Ra;
     fAmpl.fAMinus1   = kSqrt6_35*fkr.Lamda*fkr.Lamda*fkr.Ra;
     fAmpl.fAPlus1    = fAmpl.fAMinus1;
     fAmpl.fAPlus3    = fAmpl.fAMinus3;
     fAmpl.fA0LPlus   = -2*kSqrt6_35*fkr.Lamda*fkr.Lamda*fkr.Cminus;
     fAmpl.fA0LMinus  = fAmpl.fA0LPlus;
     fAmpl.fA0RPlus   = 0.;
     fAmpl.fA0RMinus  = 0.;
     break;
   }
   case (kP11_1710) :
   {
     fAmpl.fVMinus3   = 0.;
     fAmpl.fVMinus1   = k1_Sqrt6*fkr.Lamda*fkr.Lamda*fkr.Rv*(-1. + 3.*fSin28w);
     fAmpl.fVPlus1    = fAmpl.fVMinus1;
     fAmpl.fVPlus3    = 0.;
     fAmpl.fV0LPlus   = kSqrt3_8*fkr.Lamda*fkr.Lamda*fkr.Sminus*(1. - 2.*fSin28w);
     fAmpl.fV0LMinus  = fAmpl.fV0LPlus;
     fAmpl.fV0RPlus   = 0.;
     fAmpl.fV0RMinus  = 0.;
                       
     fAmpl.fAMinus3   = 0.;
     fAmpl.fAMinus1   =-k1_Sqrt6*fkr.Lamda*fkr.Lamda*fkr.Ra;
     fAmpl.fAPlus1    =-fAmpl.fAMinus1;
     fAmpl.fAPlus3    = 0.;
     fAmpl.fA0LPlus   =-k1_Sqrt6*fkr.Lamda*(fkr.Lamda*fkr.Cminus - 2.*fkr.Bminus);
     fAmpl.fA0LMinus  =-fAmpl.fA0LPlus;
     fAmpl.fA0RPlus   = 0.;
     fAmpl.fA0RMinus  = 0.;
     break;
   }
   case (kF17_1970) :
   {
     fAmpl.fVMinus3   = k1_Sqrt7/2*fkr.Lamda*fkr.Lamda*fkr.Rv;
     fAmpl.fVMinus1   = kSqrt3_35/2*fkr.Lamda*fkr.Lamda*fkr.Rv;
     fAmpl.fVPlus1    =-fAmpl.fVMinus1;
     fAmpl.fVPlus3    =-fAmpl.fVMinus3;
     fAmpl.fV0LPlus   = 0.;
     fAmpl.fV0LMinus  = 0.;
     fAmpl.fV0RPlus   = 0.;
     fAmpl.fV0RMinus  = 0.;
                         
     fAmpl.fAMinus3   = k1_Sqrt7/2*fkr.Lamda*fkr.Lamda*fkr.Ra;
     fAmpl.fAMinus1   = kSqrt3_35/2*fkr.Lamda*fkr.Lamda*fkr.Ra;
     fAmpl.fAPlus1    = fAmpl.fAMinus1;
     fAmpl.fAPlus3    = fAmpl.fAMinus3; 
     fAmpl.fA0LPlus   = kSqrt6_35/2*fkr.Lamda*fkr.Lamda*fkr.Cminus;
     fAmpl.fA0LMinus  = fAmpl.fA0LPlus;
     fAmpl.fA0RPlus   = 0.;
     fAmpl.fA0RMinus  = 0.;
     break;
   }
   default:
   {
     LOG("MKHAmpl", pWARN) << "*** UNRECOGNIZED RESONANCE!";
     fAmpl.fVMinus3   = 0.;
     fAmpl.fVMinus1   = 0.;
     fAmpl.fVPlus1    = 0.;
     fAmpl.fVPlus3    = 0.;
     fAmpl.fV0LPlus   = 0.;
     fAmpl.fV0LMinus  = 0.;
     fAmpl.fV0RPlus   = 0.;
     fAmpl.fV0RMinus  = 0.;
                    
     fAmpl.fAMinus3   = 0.;
     fAmpl.fAMinus1   = 0.;
     fAmpl.fAPlus1    = 0.;
     fAmpl.fAPlus3    = 0.;
     fAmpl.fA0LPlus   = 0.;
     fAmpl.fA0LMinus  = 0.;
     fAmpl.fA0RPlus   = 0.;
     fAmpl.fA0RMinus  = 0.;
     break;
   }

  }//switch

  return fAmpl;
}
//____________________________________________________________________________
void MKHelicityAmplModelNCp::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void MKHelicityAmplModelNCp::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void MKHelicityAmplModelNCp::LoadConfig(void)
{
  double thw;
  this->GetParam( "WeinbergAngle", thw );
  fSin28w = TMath::Power( TMath::Sin(thw), 2 );
}
//____________________________________________________________________________


