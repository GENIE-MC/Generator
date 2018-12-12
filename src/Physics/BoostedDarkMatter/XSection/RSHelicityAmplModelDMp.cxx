//____________________________________________________________________________
/*
 Copyright (c) 2003-2018, The GENIE Collaboration
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
#include "Physics/BoostedDarkMatter/XSection/RSHelicityAmplModelDMp.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
RSHelicityAmplModelDMp::RSHelicityAmplModelDMp() :
RSHelicityAmplModelI("genie::RSHelicityAmplModelDMp")
{

}
//____________________________________________________________________________
RSHelicityAmplModelDMp::RSHelicityAmplModelDMp(string config) :
RSHelicityAmplModelI("genie::RSHelicityAmplModelDMp", config)
{

}
//____________________________________________________________________________
RSHelicityAmplModelDMp::~RSHelicityAmplModelDMp()
{

}
//____________________________________________________________________________
const RSHelicityAmpl & 
   RSHelicityAmplModelDMp::Compute(
        Resonance_t res, const FKR & fkr) const
{
  switch(res) {

   case (kP33_1232) :
   {
     double RVQVud = (fQVu - fQVd) * fkr.Rv;
     double RAQAud = (fQAu - fQAd) * fkr.Ra;

     fAmpl.fMinus1 =   kSqrt2 * (RVQVud + RAQAud);
     fAmpl.fPlus1  =   kSqrt2 * (-RVQVud + RAQAud);
     fAmpl.fMinus3 =   kSqrt6 * (RVQVud + RAQAud);
     fAmpl.fPlus3  =   kSqrt6 * (-RVQVud + RAQAud);
     fAmpl.f0Minus = 2*kSqrt2 * fkr.C * (fQAd - fQVd);
     fAmpl.f0Plus  =   fAmpl.f0Minus;
     break;
   }
   case (kS11_1535) :
   {
     double TVQVud = fkr.Tv * (fQVu - fQVd);
     double TAQAud = fkr.Ta * (fQAu - fQAd);
     double RVQV5ud = fkr.Rv * (5 * fQVu + fQVd);
     double RAQA5ud = fkr.Ra * (5 * fQAu + fQAd);
     double QA5ud = 5 * fQAu + fQAd;
     double QVud =  fQVu - fQVd;
     double a = kSqrt3_2 * QVud * fkr.Lamda * fkr.S;
     double b = - QA5ud * (fkr.Lamda * fkr.C - 3 * fkr.B) / kSqrt6;

     fAmpl.fMinus1 =    -kSqrt3 * ( TVQVud + TAQAud + ( RVQV5ud + RAQA5ud) / 3. / kSqrt2);
     fAmpl.fPlus1  =    -kSqrt3 * (-TVQVud + TAQAud + (-RVQV5ud + RAQA5ud) / 3. / kSqrt2);
     fAmpl.fMinus3 =  0.;
     fAmpl.fPlus3  =  0.;
     fAmpl.f0Minus = -a + b;
     fAmpl.f0Plus  =  a + b;
     break;
   }
   case (kD13_1520) :
   {
     double TVQVud = fkr.Tv * (fQVu - fQVd);
     double TAQAud = fkr.Ta * (fQAu - fQAd);
     double RVQV5ud = fkr.Rv * (5 * fQVu + fQVd);
     double RAQA5ud = fkr.Ra * (5 * fQAu + fQAd);
     double QA5ud = 5 * fQAu + fQAd;
     double QVud =  fQVu - fQVd;
     double a = kSqrt3 * QVud * fkr.Lamda * fkr.S;
     double b = - QA5ud * fkr.Lamda * fkr.C / kSqrt3;

     fAmpl.fMinus1 = -kSqrt3_2 * ( TVQVud + TAQAud) + k1_Sqrt3 * (RVQV5ud + RAQA5ud);
     fAmpl.fPlus1  = -kSqrt3_2 * (-TVQVud + TAQAud) + k1_Sqrt3 * (RVQV5ud - RAQA5ud);
     fAmpl.fMinus3 = -3 * k1_Sqrt2 * ( TAQAud + QVud * fkr.Tv);
     fAmpl.fPlus3  = -3 * k1_Sqrt2 * (-TAQAud + QVud * fkr.Tv);
     fAmpl.f0Minus = -a + b;
     fAmpl.f0Plus  = -a - b;
     break;
   }
   case (kS11_1650) :
   {
     double LRVQVu2d = fkr.Rv * (fQAu + 2 * fQAd) * fkr.Lamda;
     double LRAQAu2d = fkr.Ra * (fQAu + 2 * fQAd) * fkr.Lamda;
     
     fAmpl.fMinus1 =  ( LRVQVu2d + LRAQAu2d) / kSqrt6;
     fAmpl.fPlus1  =  (-LRVQVu2d + LRAQAu2d) / kSqrt6;
     fAmpl.fMinus3 =  0.;
     fAmpl.fPlus3  =  0.;
     fAmpl.f0Minus =  kSqrt2_3 * (2 * fQAd + fQAu) * (3 * fkr.B - fkr.Lamda * fkr.C);
     fAmpl.f0Plus  =  fAmpl.f0Minus;
     break;
   }
   case (kD13_1700) :
   {
     double LRVQVu2d = fkr.Rv * (fQAu + 2 * fQAd) * fkr.Lamda;
     double LRAQAu2d = fkr.Ra * (fQAu + 2 * fQAd) * fkr.Lamda;

     fAmpl.fMinus1 =  k1_Sqrt30 * ( LRVQVu2d + LRAQAu2d);
     fAmpl.fPlus1  =  k1_Sqrt30 * ( LRVQVu2d - LRAQAu2d);
     fAmpl.fMinus3 =  k3_Sqrt10  * ( LRVQVu2d + LRAQAu2d);
     fAmpl.fPlus3  =  k3_Sqrt10  * ( LRVQVu2d - LRAQAu2d);
     fAmpl.f0Minus =  kSqrt2_15 * fkr.Lamda * fkr.C * (fQAu + 2 * fQAd);
     fAmpl.f0Plus  =  -1.* fAmpl.f0Minus;
     break;
   }
   case (kD15_1675) :
   {
     double LRm     = fkr.Lamda * ((2 * fQAd + fQAu) * fkr.Ra + (2 * fQVd + fQVu) * fkr.Rv);
     double LRp     = fkr.Lamda * ((2 * fQAd + fQAu) * fkr.Ra - (2 * fQVd + fQVu) * fkr.Rv);

     fAmpl.fMinus1 = -kSqrt3_10 * LRm;
     fAmpl.fPlus1  = -kSqrt3_10 * LRp;
     fAmpl.fMinus3 = -kSqrt3_5  * LRm;
     fAmpl.fPlus3  = -kSqrt3_5  * LRp;
     fAmpl.f0Minus =  kSqrt6_5 * fkr.Lamda * fkr.C * (2 * fQAd + fQAu);
     fAmpl.f0Plus  =  fAmpl.f0Minus;
     break;
   }
   case (kS31_1620) :
   {
     double Tm2xiT  =  (fQVd - fQVu) * fkr.Tv + (fQAd - fQAu) * fkr.Ta;
     double Tp2xiT  = -(fQVd - fQVu) * fkr.Tv + (fQAd - fQAu) * fkr.Ta;
     double LRm2xiR = fkr.Lamda * ( (fQVu - fQVd) * fkr.Rv + (fQAu - fQAd) * fkr.Ra);
     double LRp2xiR = fkr.Lamda * (-(fQVu - fQVd) * fkr.Rv + (fQAu - fQAd) * fkr.Ra);
     double a       = kSqrt3_2 * fkr.Lamda * fkr.S * (fQVd - fQVu);
     double b       = k1_Sqrt6 * (fkr.Lamda * fkr.C * (-fQAd + fQAu) + 3 * fkr.B * (fQAd - fQAu));

     fAmpl.fMinus1 =  kSqrt3 * Tm2xiT + k1_Sqrt6 * LRm2xiR;
     fAmpl.fPlus1  =  kSqrt3 * Tp2xiT + k1_Sqrt6 * LRp2xiR;
     fAmpl.fMinus3 =  0.;
     fAmpl.fPlus3  =  0.;
     fAmpl.f0Minus =  a+b;
     fAmpl.f0Plus  = -a+b;
     break;
   }
   case (kD33_1700) :
   {
     double Tm2xiT  = (fQVd - fQVu) * fkr.Tv + (fQAd - fQAu) * fkr.Ta;
     double Tp2xiT  = (fQVd - fQVu) * fkr.Tv - (fQAd - fQAu) * fkr.Ta;
     double LRm2xiR = fkr.Lamda * ((fQVd - fQVu) * fkr.Rv + (fQAd - fQAu) * fkr.Ra);
     double LRp2xiR = fkr.Lamda * ((fQVd - fQVu) * fkr.Rv - (fQAd - fQAu) * fkr.Ra);
     double a       = kSqrt3 * fkr.Lamda * fkr.S * (fQVd - fQVu);
     double b       = k1_Sqrt3 * fkr.Lamda * fkr.C * (fQAd - fQAu);

     fAmpl.fMinus1 = kSqrt3_2 * Tm2xiT + k1_Sqrt3 * LRm2xiR;
     fAmpl.fPlus1  = kSqrt3_2 * Tp2xiT + k1_Sqrt3 * LRp2xiR;
     fAmpl.fMinus3 = k3_Sqrt2 * Tm2xiT;
     fAmpl.fPlus3  = k3_Sqrt2 * Tp2xiT;
     fAmpl.f0Minus = a-b;
     fAmpl.f0Plus  = a+b;
     break;
   }
   case (kP11_1440) :
   {
     double c       = 0.5 * k1_Sqrt3;
     double L2      = TMath::Power(fkr.Lamda, 2);
     double L2RmxiR = L2 * ((4 * fQVu - fQVd) * fkr.Rv + (4 * fQAu - fQAd) * fkr.Ra);
     double L2RpxiR = L2 * ((4 * fQVu - fQVd) * fkr.Rv - (4 * fQAu - fQAd) * fkr.Ra);
     double a       = 3 * c * L2 * fkr.S * (2 * fQVu + fQVd);
     double b       = c * (L2 * fkr.C - 2 * fkr.Lamda * fkr.B) * (4 * fQAu - fQAd);

     fAmpl.fMinus1 = c * L2RmxiR;
     fAmpl.fPlus1  = c * L2RpxiR;
     fAmpl.fMinus3 = 0.;
     fAmpl.fPlus3  = 0.;
     fAmpl.f0Minus = -a-b;
     fAmpl.f0Plus  = -a+b;
     break;
   }
   case (kP33_1600) :
   {
     double L2      = TMath::Power(fkr.Lamda, 2);
     double L2RmxiR = L2 * ( fkr.Rv * (fQVd - fQVu) + fkr.Ra * (fQAd - fQAu));
     double L2RpxiR = L2 * (-fkr.Rv * (fQVd - fQVu) + fkr.Ra * (fQAd - fQAu));

     fAmpl.fMinus1 =  k1_Sqrt6 * L2RmxiR;
     fAmpl.fPlus1  =  k1_Sqrt6 * L2RmxiR;
     fAmpl.fMinus3 =  k1_Sqrt2 * L2RmxiR;
     fAmpl.fPlus3  =  k1_Sqrt2 * L2RpxiR;
     fAmpl.f0Minus =  kSqrt2_3 * (L2 * fkr.C - 2 * fkr.Lamda * fkr.B) * (fQAd - fQAu);
     fAmpl.f0Plus  =  fAmpl.f0Minus;
     break;
   }
   case (kP13_1720) :
   {
     double L2      = TMath::Power(fkr.Lamda, 2);
     double LTm4xiT = fkr.Lamda * ( fkr.Tv * (2 * fQVu + fQVd) + fkr.Ta * (2 * fQAu + fQAd));
     double LTp4xiT = fkr.Lamda * (-fkr.Tv * (2 * fQVu + fQVd) + fkr.Ta * (2 * fQAu + fQAd));
     double L2RmxiR = L2        * ( fkr.Rv * (4 * fQVu - fQVd) + fkr.Ra * (4 * fQAu + fQAd));
     double L2RpxiR = L2        * (-fkr.Rv * (4 * fQVu - fQVd) + fkr.Ra * (4 * fQAu + fQAd));
     double a       = kSqrt3_5 * L2 * fkr.S * (2 * fQVu + fQVd);
     double b       = k1_Sqrt15 * (L2 * fkr.C - 5 * fkr.Lamda * fkr.B) * (4 * fQAu - fQAd);

     fAmpl.fMinus1 =  kSqrt27_10 * LTm4xiT + k1_Sqrt15 * L2RmxiR;
     fAmpl.fPlus1  =  kSqrt27_10 * LTp4xiT + k1_Sqrt15 * L2RpxiR;
     fAmpl.fMinus3 = -k3_Sqrt10  * LTm4xiT;
     fAmpl.fPlus3  = -k3_Sqrt10  * LTp4xiT;
     fAmpl.f0Minus =  a+b;
     fAmpl.f0Plus  = -a+b;
     break;
   }
   case (kF15_1680) :
   {
     double L2      = TMath::Power(fkr.Lamda, 2);
     double LTm4xiT = fkr.Lamda * (fkr.Tv * (2 * fQVu + fQVd) + fkr.Ta * (2 * fQAu + fQAd));
     double LTp4xiT = fkr.Lamda * (fkr.Tv * (2 * fQVu + fQVd) - fkr.Ta * (2 * fQAu + fQAd));
     double L2RmxiR = L2 * (fkr.Rv * (4 * fQVu - fQVd) + fkr.Ra * (4 * fQAu - fQAd));
     double L2RpxiR = L2 * (fkr.Rv * (4 * fQVu - fQVd) - fkr.Ra * (4 * fQAu - fQAd));
     double a       = k3_Sqrt10 * L2 * fkr.S * (2 * fQVu + fQVd);
     double b       = k1_Sqrt10 * L2 * fkr.C * (4 * fQAu - fQAd);

     fAmpl.fMinus1 =  k3_Sqrt5 * LTm4xiT - k1_Sqrt10 * L2RmxiR;
     fAmpl.fPlus1  =  k3_Sqrt5 * LTp4xiT - k1_Sqrt10 * L2RpxiR;
     fAmpl.fMinus3 =  kSqrt18_5 * LTm4xiT;
     fAmpl.fPlus3  =  kSqrt18_5 * LTp4xiT;
     fAmpl.f0Minus =  a + b;
     fAmpl.f0Plus  =  a - b;
     break;
   }
   case (kP31_1910) :
   {
     double L2       = TMath::Power(fkr.Lamda, 2);

     fAmpl.fMinus1 = -k1_Sqrt15 * L2 * (fkr.Rv * (fQVd - fQVu) + fkr.Ra * (fQAd - fQAu));
     fAmpl.fPlus1  = -k1_Sqrt15 * L2 * (fkr.Rv * (fQVd - fQVu) - fkr.Ra * (fQAd - fQAu));
     fAmpl.fMinus3 =  0.;
     fAmpl.fPlus3  =  0.;
     fAmpl.f0Minus = -kSqrt4_15 * (L2 * fkr.C - 5 * fkr.Lamda * fkr.B) * (fQAd - fQAu);
     fAmpl.f0Plus  = -1.* fAmpl.f0Minus;
     break;
   }
   case (kP33_1920) :
   {
     double L2       = TMath::Power(fkr.Lamda, 2);
     double L2Rm2xiR = L2 * (fkr.Rv * (fQVd - fQVu) + fkr.Ra * (fQAd - fQAu));
     double L2Rp2xiR = L2 * (fkr.Rv * (fQVd - fQVu) - fkr.Ra * (fQAd - fQAu));

     fAmpl.fMinus1 =  k1_Sqrt15 * L2Rm2xiR;
     fAmpl.fPlus1  = -k1_Sqrt15 * L2Rp2xiR;
     fAmpl.fMinus3 = -k1_Sqrt5  * L2Rm2xiR;
     fAmpl.fPlus3  =  k1_Sqrt5  * L2Rp2xiR;
     fAmpl.f0Minus = -(2./kSqrt15) * (L2 * fkr.C - 5 * fkr.Lamda * fkr.B) * (fQAd - fQAu);
     fAmpl.f0Plus  =  fAmpl.f0Minus;
     break;
   }
   case (kF35_1905) :
   {
     double L2       = TMath::Power(fkr.Lamda, 2);
     double L2Rm2xiR = L2 * (fkr.Rv * (fQVd - fQVu) + fkr.Ra * (fQAd - fQAu));
     double L2Rp2xiR = L2 * (fkr.Rv * (fQVd - fQVu) - fkr.Ra * (fQAd - fQAu));

     fAmpl.fMinus1 =  k1_Sqrt35  * L2Rm2xiR;
     fAmpl.fPlus1  =  k1_Sqrt35  * L2Rp2xiR;
     fAmpl.fMinus3 =  kSqrt18_35 * L2Rm2xiR;
     fAmpl.fPlus3  =  kSqrt18_35 * L2Rp2xiR;
     fAmpl.f0Minus =  k2_Sqrt35  * L2 * fkr.C * (fQAd - fQAu);
     fAmpl.f0Plus  =  -1. * fAmpl.f0Minus;
     break;
   }
   case (kF37_1950) :
   {
     double L2       = TMath::Power(fkr.Lamda, 2);
     double L2Rm2xiR = L2 * (fkr.Rv * (fQVd - fQVu) + fkr.Ra * (fQAd - fQAu));
     double L2Rp2xiR = L2 * (fkr.Rv * (fQVd - fQVu) - fkr.Ra * (fQAd - fQAu));

     fAmpl.fMinus1 =  -kSqrt6_35 * L2Rm2xiR;
     fAmpl.fPlus1  =   kSqrt6_35 * L2Rp2xiR;
     fAmpl.fMinus3 =  -kSqrt2_7  * L2Rm2xiR;
     fAmpl.fPlus3  =   kSqrt2_7  * L2Rp2xiR;
     fAmpl.f0Minus = 2*kSqrt6_35 * L2 * fkr.C * (fQAd - fQAu);
     fAmpl.f0Plus  =  fAmpl.f0Minus;
     break;
   }
   case (kP11_1710) :
   {
     double L2      = TMath::Power(fkr.Lamda, 2);
     double Rm3xiR  = fkr.Rv * (5 * fQVu + fQVd) + fkr.Ra * (5 * fQAu + fQAd);
     double Rp3xiR  = fkr.Rv * (5 * fQVu + fQVd) - fkr.Ra * (5 * fQAu + fQAd);
     double a       = kSqrt3_8 * L2 * fkr.S * (fQVu - fQVd);
     double b       = 0.5 * k1_Sqrt6 * (L2 * fkr.C - 2 * fkr.Lamda * fkr.B) * (5 * fQAu + fQAd);

     fAmpl.fMinus1 =  0.5 * k1_Sqrt6 * L2 * Rm3xiR;
     fAmpl.fPlus1  =  0.5 * k1_Sqrt6 * L2 * Rp3xiR;
     fAmpl.fMinus3 =  0.;
     fAmpl.fPlus3  =  0.;
     fAmpl.f0Minus =  a+b;
     fAmpl.f0Plus  =  a-b;
     break;
   }
   case (kF17_1970) :
   {
     double L2      = TMath::Power(fkr.Lamda, 2);
     double Rm3xiR  = fkr.Rv * (fQVu + 2 * fQVd) + fkr.Ra * (fQAu + 2 * fQAd);
     double Rp3xiR  = fkr.Rv * (fQVu + 2 * fQVd) - fkr.Ra * (fQAu + 2 * fQAd);

     fAmpl.fMinus1 = -kSqrt3_35 * k1_Sqrt2 * L2 * Rm3xiR;
     fAmpl.fPlus1  =  kSqrt3_35 * k1_Sqrt2 * L2 * Rp3xiR;
     fAmpl.fMinus3 = -k1_Sqrt7 * k1_Sqrt2 * L2 * Rm3xiR;
     fAmpl.fPlus3  =  k1_Sqrt7 * k1_Sqrt2 * L2 * Rp3xiR;
     fAmpl.f0Minus =  kSqrt6_35 * L2 * fkr.C * (fQAu + 2 * fQAd);
     fAmpl.f0Plus  =  fAmpl.f0Minus;
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
void RSHelicityAmplModelDMp::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void RSHelicityAmplModelDMp::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void RSHelicityAmplModelDMp::LoadConfig(void)
{
  double QLu, QRu, QLd, QRd ;
  GetParam( "UpLeftCharge", QLu ) ;
  GetParam( "UpRightCharge", QRu ) ;
  GetParam( "DownLeftCharge", QLd ) ;
  GetParam( "DownRightCharge", QRd ) ;
  fQVu = 0.5 * (QLu + QRu);
  fQAu = 0.5 * (-QLu + QRu);
  fQVd = 0.5 * (QLd + QRd);
  fQAd = 0.5 * (-QLd + QRd);
}
//____________________________________________________________________________

