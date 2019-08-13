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
#include "Physics/Resonance/XSection/MKHelicityAmplModelCC.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
MKHelicityAmplModelCC::MKHelicityAmplModelCC() :
MKHelicityAmplModelI("genie::MKHelicityAmplModelCC")
{

}
//____________________________________________________________________________
MKHelicityAmplModelCC::MKHelicityAmplModelCC(string config) :
MKHelicityAmplModelI("genie::MKHelicityAmplModelCC", config)
{

}
//____________________________________________________________________________
MKHelicityAmplModelCC::~MKHelicityAmplModelCC()
{

}
//____________________________________________________________________________
const MKHelicityAmpl & MKHelicityAmplModelCC::Compute(Resonance_t res, const FKR_MK & fkr) const
{
  switch(res) {

   case (kP33_1232) :
   {
     fAmpl.fVMinus3    =-kSqrt6*fkr.Rv;
     fAmpl.fVMinus1    =-kSqrt2*fkr.Rv;
     fAmpl.fVPlus1     =-fAmpl.fVMinus1;
     fAmpl.fVPlus3     =-fAmpl.fVMinus3;
     fAmpl.fV0LPlus    = 0.;
     fAmpl.fV0LMinus   = 0.;
     fAmpl.fV0RPlus    = 0.;
     fAmpl.fV0RMinus   = 0.;
     
     fAmpl.fAMinus3    =-kSqrt6*fkr.Ra;
     fAmpl.fAMinus1    =-kSqrt2*fkr.Ra;
     fAmpl.fAPlus1     = fAmpl.fAMinus1;
     fAmpl.fAPlus3     = fAmpl.fAMinus3; 
     fAmpl.fA0LPlus    = 2.*kSqrt2*fkr.Cminus;
     fAmpl.fA0LMinus   = fAmpl.fA0LPlus;
     fAmpl.fA0RPlus    = 2.*kSqrt2*fkr.Cplus;
     fAmpl.fA0RMinus   = fAmpl.fA0RPlus;
     break;
   }
   case (kS11_1535) :
   {
     fAmpl.fVMinus3    = 0.;
     fAmpl.fVMinus1    =-2.*kSqrt3*fkr.Tv - 4./kSqrt6*fkr.Lamda*fkr.Rv;
     fAmpl.fVPlus1     =-fAmpl.fVMinus1;
     fAmpl.fVPlus3     = 0.;
     fAmpl.fV0LPlus    = kSqrt6*fkr.Lamda*fkr.Sminus;
     fAmpl.fV0LMinus   =-fAmpl.fV0LPlus;
     fAmpl.fV0RPlus    = kSqrt6*fkr.Lamda*fkr.Splus;
     fAmpl.fV0RMinus   =-fAmpl.fV0RPlus;
     
     fAmpl.fAMinus3    = 0.;
     fAmpl.fAMinus1    =-2.*kSqrt3*fkr.Ta - 4./kSqrt6*fkr.Lamda*fkr.Ra;
     fAmpl.fAPlus1     = fAmpl.fAMinus1;
     fAmpl.fAPlus3     = 0.;
     fAmpl.fA0LPlus    =-2.*kSqrt2_3*(fkr.Lamda*fkr.Cminus - 3.*fkr.Bminus);
     fAmpl.fA0LMinus   = fAmpl.fA0LPlus;
     fAmpl.fA0RPlus    =-2.*kSqrt2_3*(fkr.Lamda*fkr.Cplus  - 3.*fkr.Bplus);
     fAmpl.fA0RMinus   = fAmpl.fA0RPlus;
     break;
   }
   case (kD13_1520) :
   {
     fAmpl.fVMinus3    =-6./kSqrt2*fkr.Tv;
     fAmpl.fVMinus1    =-kSqrt6*fkr.Tv + 4./kSqrt3*fkr.Lamda*fkr.Rv;
     fAmpl.fVPlus1     = fAmpl.fVMinus1;
     fAmpl.fVPlus3     = fAmpl.fVMinus3;
     fAmpl.fV0LPlus    =-2.*kSqrt3*fkr.Lamda*fkr.Sminus;
     fAmpl.fV0LMinus   = fAmpl.fV0LPlus;
     fAmpl.fV0RPlus    =-2.*kSqrt3*fkr.Lamda*fkr.Splus;
     fAmpl.fV0RMinus   = fAmpl.fV0RPlus;
     
     fAmpl.fAMinus3    =-6./kSqrt2*fkr.Ta;
     fAmpl.fAMinus1    =-kSqrt6*fkr.Ta + 4./kSqrt3*fkr.Lamda*fkr.Ra;
     fAmpl.fAPlus1     =-fAmpl.fAMinus1;
     fAmpl.fAPlus3     =-fAmpl.fAMinus3;
     fAmpl.fA0LPlus    = 4./kSqrt3*fkr.Lamda*fkr.Cminus;
     fAmpl.fA0LMinus   =-fAmpl.fA0LPlus;
     fAmpl.fA0RPlus    = 4./kSqrt3*fkr.Lamda*fkr.Cplus;
     fAmpl.fA0RMinus   =-fAmpl.fA0RPlus;
     break;
   }
   case (kS11_1650) :
   {
     fAmpl.fVMinus3    = 0.;
     fAmpl.fVMinus1    =-k1_Sqrt6*fkr.Lamda*fkr.Rv;
     fAmpl.fVPlus1     =-fAmpl.fVMinus1;
     fAmpl.fVPlus3     = 0.;
     fAmpl.fV0LPlus    = 0.;
     fAmpl.fV0LMinus   = 0.;
     fAmpl.fV0RPlus    = 0.;
     fAmpl.fV0RMinus   = 0.;
     
     fAmpl.fAMinus3    = 0.;
     fAmpl.fAMinus1    =-k1_Sqrt6*fkr.Lamda*fkr.Ra;
     fAmpl.fAPlus1     = fAmpl.fAMinus1;
     fAmpl.fAPlus3     = 0.;
     fAmpl.fA0LPlus    = kSqrt2_3*(fkr.Lamda*fkr.Cminus - 3.*fkr.Bminus);
     fAmpl.fA0LMinus   = fAmpl.fA0LPlus;
     fAmpl.fA0RPlus    = kSqrt2_3*(fkr.Lamda*fkr.Cplus  - 3.*fkr.Bplus);
     fAmpl.fA0RMinus   = fAmpl.fA0RPlus;
     break;
   }
   case (kD13_1700) :
   {
     fAmpl.fVMinus3    =-k3_Sqrt10*fkr.Lamda*fkr.Rv;
     fAmpl.fVMinus1    =-k1_Sqrt30*fkr.Lamda*fkr.Rv;
     fAmpl.fVPlus1     = fAmpl.fVMinus1;
     fAmpl.fVPlus3     = fAmpl.fVMinus3;
     fAmpl.fV0LPlus    = 0.;
     fAmpl.fV0LMinus   = 0.;
     fAmpl.fV0RPlus    = 0.;
     fAmpl.fV0RMinus   = 0.;
     
     fAmpl.fAMinus3    =-k3_Sqrt10*fkr.Lamda*fkr.Ra;
     fAmpl.fAMinus1    =-k1_Sqrt30*fkr.Lamda*fkr.Ra;
     fAmpl.fAPlus1     =-fAmpl.fAMinus1;
     fAmpl.fAPlus3     =-fAmpl.fAMinus3;
     fAmpl.fA0LPlus    = kSqrt2_15*fkr.Lamda*fkr.Cminus;
     fAmpl.fA0LMinus   =-fAmpl.fA0LPlus;
     fAmpl.fA0RPlus    = kSqrt2_15*fkr.Lamda*fkr.Cplus;
     fAmpl.fA0RMinus   =-fAmpl.fA0RPlus;
     break;
   }
   case (kD15_1675) :
   {
     fAmpl.fVMinus3    = kSqrt3_5 *fkr.Lamda*fkr.Rv;
     fAmpl.fVMinus1    = kSqrt3_10*fkr.Lamda*fkr.Rv;
     fAmpl.fVPlus1     =-fAmpl.fVMinus1;
     fAmpl.fVPlus3     =-fAmpl.fVMinus3;
     fAmpl.fV0LPlus    = 0.;
     fAmpl.fV0LMinus   = 0.;
     fAmpl.fV0RPlus    = 0.;
     fAmpl.fV0RMinus   = 0.;
                       
     fAmpl.fAMinus3    = kSqrt3_5 *fkr.Lamda*fkr.Ra;
     fAmpl.fAMinus1    = kSqrt3_10*fkr.Lamda*fkr.Ra;
     fAmpl.fAPlus1     = fAmpl.fAMinus1;
     fAmpl.fAPlus3     = fAmpl.fAMinus3;
     fAmpl.fA0LPlus    =-kSqrt6_5 *fkr.Lamda*fkr.Cminus;
     fAmpl.fA0LMinus   = fAmpl.fA0LPlus;
     fAmpl.fA0RPlus    =-kSqrt6_5 *fkr.Lamda*fkr.Cplus;
     fAmpl.fA0RMinus   = fAmpl.fA0RPlus;
     break;
   }
   case (kS31_1620) :
   {
     fAmpl.fVMinus3    = 0.;
     fAmpl.fVMinus1    = kSqrt3*fkr.Tv - k1_Sqrt6*fkr.Lamda*fkr.Rv;
     fAmpl.fVPlus1     =-fAmpl.fVMinus1;
     fAmpl.fVPlus3     = 0.;
     fAmpl.fV0LPlus    =-kSqrt3_2*fkr.Lamda*fkr.Sminus;
     fAmpl.fV0LMinus   =-fAmpl.fV0LPlus;
     fAmpl.fV0RPlus    =-kSqrt3_2*fkr.Lamda*fkr.Splus;
     fAmpl.fV0RMinus   =-fAmpl.fV0RPlus;
                       
     fAmpl.fAMinus3    = 0.;
     fAmpl.fAMinus1    = kSqrt3*fkr.Ta - k1_Sqrt6*fkr.Lamda*fkr.Ra;
     fAmpl.fAPlus1     = fAmpl.fAMinus1;
     fAmpl.fAPlus3     = 0.;
     fAmpl.fA0LPlus    =-k1_Sqrt6*(fkr.Lamda*fkr.Cminus - 3.*fkr.Bminus);
     fAmpl.fA0LMinus   = fAmpl.fA0LPlus;
     fAmpl.fA0RPlus    =-k1_Sqrt6*(fkr.Lamda*fkr.Cplus  - 3.*fkr.Bplus);
     fAmpl.fA0RMinus   = fAmpl.fA0RPlus;
     break;
   }
   case (kD33_1700) :
   {
     fAmpl.fVMinus3    = k3_Sqrt2*fkr.Tv;
     fAmpl.fVMinus1    = kSqrt3_2*fkr.Tv + k1_Sqrt3*fkr.Lamda*fkr.Rv;
     fAmpl.fVPlus1     = fAmpl.fVMinus1;
     fAmpl.fVPlus3     = fAmpl.fVMinus3;
     fAmpl.fV0LPlus    = kSqrt3*fkr.Lamda*fkr.Sminus;
     fAmpl.fV0LMinus   = fAmpl.fV0LPlus;
     fAmpl.fV0RPlus    = kSqrt3*fkr.Lamda*fkr.Splus;
     fAmpl.fV0RMinus   = fAmpl.fV0RPlus;
                       
     fAmpl.fAMinus3    = k3_Sqrt2*fkr.Ta;
     fAmpl.fAMinus1    = kSqrt3_2*fkr.Ta + k1_Sqrt3*fkr.Lamda*fkr.Ra;
     fAmpl.fAPlus1     =-fAmpl.fAMinus1;
     fAmpl.fAPlus3     =-fAmpl.fAMinus3;
     fAmpl.fA0LPlus    = k1_Sqrt3*fkr.Lamda*fkr.Cminus;
     fAmpl.fA0LMinus   =-fAmpl.fA0LPlus;
     fAmpl.fA0RPlus    = k1_Sqrt3*fkr.Lamda*fkr.Cplus;
     fAmpl.fA0RMinus   =-fAmpl.fA0RPlus;
     break;
   }
   case (kP11_1440) :
   {
     fAmpl.fVMinus3    = 0.;
     fAmpl.fVMinus1    = 5.*kSqrt3/6.*fkr.Lamda*fkr.Lamda*fkr.Rv;
     fAmpl.fVPlus1     = fAmpl.fVMinus1;
     fAmpl.fVPlus3     = 0.;
     fAmpl.fV0LPlus    =-kSqrt3_4*fkr.Lamda*fkr.Lamda*fkr.Sminus;
     fAmpl.fV0LMinus   = fAmpl.fV0LPlus;
     fAmpl.fV0RPlus    =-kSqrt3_4*fkr.Lamda*fkr.Lamda*fkr.Splus;
     fAmpl.fV0RMinus   = fAmpl.fV0RPlus;
                       
     fAmpl.fAMinus3    = 0.;
     fAmpl.fAMinus1    = 5.*kSqrt3/6.*fkr.Lamda*fkr.Lamda*fkr.Ra;
     fAmpl.fAPlus1     =-fAmpl.fAMinus1;
     fAmpl.fAPlus3     = 0.;
     fAmpl.fA0LPlus    = 5.*kSqrt3/6.*fkr.Lamda*(fkr.Lamda*fkr.Cminus - 2.*fkr.Bminus);
     fAmpl.fA0LMinus   =-fAmpl.fA0LPlus; 
     fAmpl.fA0RPlus    = 5.*kSqrt3/6.*fkr.Lamda*(fkr.Lamda*fkr.Cplus  - 2.*fkr.Bplus);
     fAmpl.fA0RMinus   =-fAmpl.fA0RPlus;
     break;
   }
   case (kP33_1600) :
   {
     fAmpl.fVMinus3    = k1_Sqrt2*fkr.Lamda*fkr.Lamda*fkr.Rv;
     fAmpl.fVMinus1    = k1_Sqrt6*fkr.Lamda*fkr.Lamda*fkr.Rv;
     fAmpl.fVPlus1     =-fAmpl.fVMinus1;
     fAmpl.fVPlus3     =-fAmpl.fVMinus3;
     fAmpl.fV0LPlus    = 0.;
     fAmpl.fV0LMinus   = 0.;
     fAmpl.fV0RPlus    = 0.;
     fAmpl.fV0RMinus   = 0.;
                        
     fAmpl.fAMinus3    = k1_Sqrt2*fkr.Lamda*fkr.Lamda*fkr.Ra;
     fAmpl.fAMinus1    = k1_Sqrt6*fkr.Lamda*fkr.Lamda*fkr.Ra;
     fAmpl.fAPlus1     = fAmpl.fAMinus1;
     fAmpl.fAPlus3     = fAmpl.fAMinus3;
     fAmpl.fA0LPlus    =-kSqrt2_3*fkr.Lamda*(fkr.Lamda*fkr.Cminus - 2.*fkr.Bminus);
     fAmpl.fA0LMinus   = fAmpl.fA0LPlus;
     fAmpl.fA0RPlus    =-kSqrt2_3*fkr.Lamda*(fkr.Lamda*fkr.Cplus  -2.*fkr.Bplus);
     fAmpl.fA0RMinus   = fAmpl.fA0RPlus;
     break;
   }
   case (kP13_1720) :
   {     
     fAmpl.fVMinus3    =-k3_Sqrt10*fkr.Lamda*fkr.Tv;
     fAmpl.fVMinus1    = kSqrt27_10*fkr.Lamda*fkr.Tv + kSqrt5_3*fkr.Lamda*fkr.Lamda*fkr.Rv;
     fAmpl.fVPlus1     =-fAmpl.fVMinus1; 
     fAmpl.fVPlus3     =-fAmpl.fVMinus3;
     fAmpl.fV0LPlus    =-kSqrt3_5*fkr.Lamda*fkr.Lamda*fkr.Sminus;
     fAmpl.fV0LMinus   =-fAmpl.fV0LPlus; 
     fAmpl.fV0RPlus    =-kSqrt3_5*fkr.Lamda*fkr.Lamda*fkr.Splus;
     fAmpl.fV0RMinus   =-fAmpl.fV0RPlus;
                        
     fAmpl.fAMinus3    =-k3_Sqrt10*fkr.Lamda*fkr.Ta;
     fAmpl.fAMinus1    = kSqrt27_10*fkr.Lamda*fkr.Ta + kSqrt5_3*fkr.Lamda*fkr.Lamda*fkr.Ra;
     fAmpl.fAPlus1     = fAmpl.fAMinus1; 
     fAmpl.fAPlus3     = fAmpl.fAMinus3;
     fAmpl.fA0LPlus    = kSqrt5_3*fkr.Lamda*(fkr.Lamda*fkr.Cminus - 5.*fkr.Bminus);
     fAmpl.fA0LMinus   = fAmpl.fA0LPlus;
     fAmpl.fA0RPlus    = kSqrt5_3*fkr.Lamda*(fkr.Lamda*fkr.Cplus  - 5.*fkr.Bplus);
     fAmpl.fA0RMinus   = fAmpl.fA0RPlus;
     break;
   }
   case (kF15_1680) :
   {
     fAmpl.fVMinus3   = kSqrt18_5*fkr.Lamda*fkr.Tv;
     fAmpl.fVMinus1   = k3_Sqrt5*fkr.Lamda*fkr.Tv - kSqrt5_2*fkr.Lamda*fkr.Lamda*fkr.Rv;
     fAmpl.fVPlus1    = fAmpl.fVMinus1; 
     fAmpl.fVPlus3    = fAmpl.fVMinus3;
     fAmpl.fV0LPlus   = kSqrt9_10*fkr.Lamda*fkr.Lamda*fkr.Sminus;
     fAmpl.fV0LMinus  = fAmpl.fV0LPlus;
     fAmpl.fV0RPlus   = kSqrt9_10*fkr.Lamda*fkr.Lamda*fkr.Splus;
     fAmpl.fV0RMinus  = fAmpl.fV0RPlus;
                       
     fAmpl.fAMinus3   = kSqrt18_5*fkr.Lamda*fkr.Ta;
     fAmpl.fAMinus1   = k3_Sqrt5*fkr.Lamda*fkr.Ta - kSqrt5_2*fkr.Lamda*fkr.Lamda*fkr.Ra;
     fAmpl.fAPlus1    =-fAmpl.fAMinus1; 
     fAmpl.fAPlus3    =-fAmpl.fAMinus3;
     fAmpl.fA0LPlus   =-kSqrt5_2*fkr.Lamda*fkr.Lamda*fkr.Cminus;
     fAmpl.fA0LMinus  =-fAmpl.fA0LPlus;
     fAmpl.fA0RPlus   =-kSqrt5_2*fkr.Lamda*fkr.Lamda*fkr.Cplus;
     fAmpl.fA0RMinus  =-fAmpl.fA0RPlus;
     break;
   }
   case (kP31_1910) :
   {     
     fAmpl.fVMinus3   = 0.;
     fAmpl.fVMinus1   =-k1_Sqrt15*fkr.Lamda*fkr.Lamda*fkr.Rv;
     fAmpl.fVPlus1    = fAmpl.fVMinus1;
     fAmpl.fVPlus3    = 0.;
     fAmpl.fV0LPlus   = 0.;
     fAmpl.fV0LMinus  = 0.;
     fAmpl.fV0RPlus   = 0.;
     fAmpl.fV0RMinus  = 0.;
                       
     fAmpl.fAMinus3   = 0.;
     fAmpl.fAMinus1   =-k1_Sqrt15*fkr.Lamda*fkr.Lamda*fkr.Ra;
     fAmpl.fAPlus1    =-fAmpl.fAMinus1;
     fAmpl.fAPlus3    = 0.;
     fAmpl.fA0LPlus   = k2_Sqrt15*fkr.Lamda*(fkr.Lamda*fkr.Cminus - 5.*fkr.Bminus);
     fAmpl.fA0LMinus  =-fAmpl.fA0LPlus;
     fAmpl.fA0RPlus   = k2_Sqrt15*fkr.Lamda*(fkr.Lamda*fkr.Cplus  - 5.*fkr.Bplus);
     fAmpl.fA0RMinus  =-fAmpl.fA0RPlus;
     break;
   }
   case (kP33_1920) :
   {
     fAmpl.fVMinus3   =-k1_Sqrt5*fkr.Lamda*fkr.Lamda*fkr.Rv;
     fAmpl.fVMinus1   = k1_Sqrt15*fkr.Lamda*fkr.Lamda*fkr.Rv;
     fAmpl.fVPlus1    =-fAmpl.fVMinus1;
     fAmpl.fVPlus3    =-fAmpl.fVMinus3;
     fAmpl.fV0LPlus   = 0.;
     fAmpl.fV0LMinus  = 0.;
     fAmpl.fV0RPlus   = 0.;
     fAmpl.fV0RMinus  = 0.;
                       
     fAmpl.fAMinus3   =-k1_Sqrt5*fkr.Lamda*fkr.Lamda*fkr.Ra;
     fAmpl.fAMinus1   = k1_Sqrt15*fkr.Lamda*fkr.Lamda*fkr.Ra;
     fAmpl.fAPlus1    = fAmpl.fAMinus1;
     fAmpl.fAPlus3    = fAmpl.fAMinus3;
     fAmpl.fA0LPlus   =-k2_Sqrt15*fkr.Lamda*(fkr.Lamda*fkr.Cminus - 5.*fkr.Bminus);
     fAmpl.fA0LMinus  = fAmpl.fA0LPlus;
     fAmpl.fA0RPlus   =-k2_Sqrt15*fkr.Lamda*(fkr.Lamda*fkr.Cplus  - 5.*fkr.Bplus);
     fAmpl.fA0RMinus  = fAmpl.fA0RPlus;
     break;
   }
   case (kF35_1905) :
   {
     fAmpl.fVMinus3   = kSqrt18_35*fkr.Lamda*fkr.Lamda*fkr.Rv;
     fAmpl.fVMinus1   = k1_Sqrt35*fkr.Lamda*fkr.Lamda*fkr.Rv;
     fAmpl.fVPlus1    = fAmpl.fVMinus1;
     fAmpl.fVPlus3    = fAmpl.fVMinus3;
     fAmpl.fV0LPlus   = 0.;
     fAmpl.fV0LMinus  = 0.;
     fAmpl.fV0RPlus   = 0.;
     fAmpl.fV0RMinus  = 0.;
                       
     fAmpl.fAMinus3   = kSqrt18_35*fkr.Lamda*fkr.Lamda*fkr.Ra;
     fAmpl.fAMinus1   = k1_Sqrt35*fkr.Lamda*fkr.Lamda*fkr.Ra;
     fAmpl.fAPlus1    =-fAmpl.fAMinus1;
     fAmpl.fAPlus3    =-fAmpl.fAMinus3;
     fAmpl.fA0LPlus   =-k2_Sqrt35*fkr.Lamda*fkr.Lamda*fkr.Cminus;
     fAmpl.fA0LMinus  =-fAmpl.fA0LPlus;
     fAmpl.fA0RPlus   =-k2_Sqrt35*fkr.Lamda*fkr.Lamda*fkr.Cplus;
     fAmpl.fA0RMinus  =-fAmpl.fA0RPlus;
     break;
   }
   case (kF37_1950) :
   {
     fAmpl.fVMinus3   =-kSqrt2_7*fkr.Lamda*fkr.Lamda*fkr.Rv;
     fAmpl.fVMinus1   =-kSqrt6_35*fkr.Lamda*fkr.Lamda*fkr.Rv;
     fAmpl.fVPlus1    =-fAmpl.fVMinus1;
     fAmpl.fVPlus3    =-fAmpl.fVMinus3;
     fAmpl.fV0LPlus   = 0.;
     fAmpl.fV0LMinus  = 0.;
     fAmpl.fV0RPlus   = 0.;
     fAmpl.fV0RMinus  = 0.;
                       
     fAmpl.fAMinus3   =-kSqrt2_7*fkr.Lamda*fkr.Lamda*fkr.Ra;
     fAmpl.fAMinus1   =-kSqrt6_35*fkr.Lamda*fkr.Lamda*fkr.Ra;
     fAmpl.fAPlus1    = fAmpl.fAMinus1;
     fAmpl.fAPlus3    = fAmpl.fAMinus3;
     fAmpl.fA0LPlus   = 2.*kSqrt6_35*fkr.Lamda*fkr.Lamda*fkr.Cminus;
     fAmpl.fA0LMinus  = fAmpl.fA0LPlus;
     fAmpl.fA0RPlus   = 2.*kSqrt6_35*fkr.Lamda*fkr.Lamda*fkr.Cplus;
     fAmpl.fA0RMinus  = fAmpl.fA0RPlus;
     break;
   }
   case (kP11_1710) :
   {
     fAmpl.fVMinus3   = 0.;
     fAmpl.fVMinus1   =-kSqrt2_3*fkr.Lamda*fkr.Lamda*fkr.Rv;
     fAmpl.fVPlus1    = fAmpl.fVMinus1;
     fAmpl.fVPlus3    = 0.;
     fAmpl.fV0LPlus   = kSqrt3_2*fkr.Lamda*fkr.Lamda*fkr.Sminus;
     fAmpl.fV0LMinus  = fAmpl.fV0LPlus;
     fAmpl.fV0RPlus   = kSqrt3_2*fkr.Lamda*fkr.Lamda*fkr.Splus;
     fAmpl.fV0RMinus  = fAmpl.fV0RPlus;
                       
     fAmpl.fAMinus3   = 0.;
     fAmpl.fAMinus1   =-kSqrt2_3*fkr.Lamda*fkr.Lamda*fkr.Ra;
     fAmpl.fAPlus1    =-fAmpl.fAMinus1;
     fAmpl.fAPlus3    = 0.;
     fAmpl.fA0LPlus   =-kSqrt2_3*fkr.Lamda*(fkr.Lamda*fkr.Cminus - 2.*fkr.Bminus);
     fAmpl.fA0LMinus  =-fAmpl.fA0LPlus;
     fAmpl.fA0RPlus   =-kSqrt2_3*fkr.Lamda*(fkr.Lamda*fkr.Cplus  - 2.*fkr.Bplus);
     fAmpl.fA0RMinus  =-fAmpl.fA0RPlus;
     break;
   }
   case (kF17_1970) :
   {
     fAmpl.fVMinus3   = k1_Sqrt7*fkr.Lamda*fkr.Lamda*fkr.Rv;
     fAmpl.fVMinus1   = kSqrt3_35 *fkr.Lamda*fkr.Lamda*fkr.Rv;
     fAmpl.fVPlus1    =-fAmpl.fVMinus1;
     fAmpl.fVPlus3    =-fAmpl.fVMinus3;
     fAmpl.fV0LPlus   = 0.;
     fAmpl.fV0LMinus  = 0.;
     fAmpl.fV0RPlus   = 0.;
     fAmpl.fV0RMinus  = 0.;
                        
     fAmpl.fAMinus3   = k1_Sqrt7*fkr.Lamda*fkr.Lamda*fkr.Ra;
     fAmpl.fAMinus1   = kSqrt3_35*fkr.Lamda*fkr.Lamda*fkr.Ra;
     fAmpl.fAPlus1    = fAmpl.fAMinus1;
     fAmpl.fAPlus3    = fAmpl.fAMinus3; 
     fAmpl.fA0LPlus   =-kSqrt6_35*fkr.Lamda*fkr.Lamda*fkr.Cminus;
     fAmpl.fA0LMinus  = fAmpl.fA0LPlus;
     fAmpl.fA0RPlus   =-kSqrt6_35*fkr.Lamda*fkr.Lamda*fkr.Cplus;
     fAmpl.fA0RMinus  = fAmpl.fA0RPlus;
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


