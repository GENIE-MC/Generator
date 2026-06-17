//____________________________________________________________________________
/*
Dark Matter Resonant Scattering Helicity Amplitudes for proton Interactions


Zachary W. Orr, Colorado State University
*/
//____________________________________________________________________________

#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/BoostedDarkMatter/XSection/RSHelicityAmplModelDMp.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
RSHelicityAmplModelDMp::RSHelicityAmplModelDMp() :
RSHelicityAmplModelDMI("genie::RSHelicityAmplModelDMp")
{

}
//____________________________________________________________________________
RSHelicityAmplModelDMp::RSHelicityAmplModelDMp(string config) :
RSHelicityAmplModelDMI("genie::RSHelicityAmplModelDMp", config)
{

}
//____________________________________________________________________________
RSHelicityAmplModelDMp::~RSHelicityAmplModelDMp()
{

}
//____________________________________________________________________________
const RSHelicityAmplDM &
    RSHelicityAmplModelDMp::Compute(
          Resonance_t res, const FKRDM & fkrdm) const
{
  switch(res) {


   case (kP33_1232) :
   {
     fAmpl.fMinus3 = kSqrt6 * (fkrdm.Ra * (fQdA - fQuA) + fkrdm.Rv * (fQuV - fQdV));
     fAmpl.fMinus1 = kSqrt2 * (fkrdm.Ra * (fQdA - fQuA) + fkrdm.Rv * (fQuV - fQdV));
     fAmpl.fPlus1 = kSqrt2 * (fkrdm.Ra * (fQuA - fQdA) + fkrdm.Rv * (fQuV - fQdV));
     fAmpl.fPlus3 = kSqrt6 * (fkrdm.Ra * (fQuA - fQdA) + fkrdm.Rv * (fQuV - fQdV));
     fAmpl.f0Plus = 2.0 * kSqrt2 * fkrdm.Cs * (fQuA - fQdA);
     fAmpl.f0Minus = 2.0 * kSqrt2 * fkrdm.Cs * (fQuA - fQdA);
     fAmpl.fz0Plus = 2.0 * kSqrt2 * fkrdm.Cz * (fQuA - fQdA);
     fAmpl.fz0Minus = 2.0 * kSqrt2 * fkrdm.Cz * (fQuA - fQdA);

     break;
   }
   case (kS11_1535) :
   {
   fAmpl.fMinus3 = 0;
   fAmpl.fMinus1 = k1_Sqrt6 * fkrdm.Lamda * (fkrdm.Ra * (fQdA + 5 * fQuA) - fkrdm.Rv * (fQdV + 5 * fQuV)) + kSqrt3 * fkrdm.Ta * (fQuA - fQdA) + kSqrt3 * fkrdm.Tv * (fQdV - fQuV);
   fAmpl.fPlus1 = -k1_Sqrt6 * fkrdm.Lamda * (fkrdm.Ra * (fQdA + 5 * fQuA) + fkrdm.Rv * (fQdV + 5 * fQuV)) + kSqrt3 * fkrdm.Ta * (fQdA - fQuA) + kSqrt3 * fkrdm.Tv * (fQdV - fQuV);
   fAmpl.fPlus3 = 0;
   fAmpl.f0Plus = k1_Sqrt6 * (3 * fkrdm.Lamda * fkrdm.S * (fQuV - fQdV) - 3 * fkrdm.Bs * (fQdA + 5 * fQuA) + fkrdm.Lamda * fkrdm.Cs * (fQdA + 5 * fQuA));
   fAmpl.f0Minus = k1_Sqrt6 * (3 * fkrdm.Lamda * fkrdm.S * (fQdV - fQuV) - 3 * fkrdm.Bs * (fQdA + 5 * fQuA) + fkrdm.Lamda * fkrdm.Cs * (fQdA + 5 * fQuA));
   fAmpl.fz0Plus = -k1_Sqrt6 * ((fQdA + 5 * fQuA) * (3 * fkrdm.Bz - fkrdm.Lamda * fkrdm.Cz));
   fAmpl.fz0Minus = -k1_Sqrt6 * ((fQdA + 5 * fQuA) * (3 * fkrdm.Bz - fkrdm.Lamda * fkrdm.Cz));

     break;
   }
   case (kD13_1520) :
   {
   fAmpl.fMinus3 = -k3_Sqrt2 * (fkrdm.Tv * (fQuV - fQdV) + fkrdm.Ta * (fQdA - fQuA));
   fAmpl.fMinus1 = 0.5 * k1_Sqrt3 * (2.0 * fkrdm.Lamda * (fkrdm.Rv * (fQdV + 5 * fQuV) - 5 * fkrdm.Ra * fQuA) - fQdA * (2.0 * fkrdm.Lamda * fkrdm.Ra + 3 * kSqrt2 * fkrdm.Ta) + 3 * kSqrt2 * fkrdm.Tv * (fQdV - fQuV) + 3 * kSqrt2 * fQuA * fkrdm.Ta);
   fAmpl.fPlus1 = 0.5 * k1_Sqrt3 * (-2.0 * fkrdm.Lamda * (fkrdm.Rv * (fQdV + 5 * fQuV) + 5 * fkrdm.Ra * fQuA) - fQdA * (2.0 * fkrdm.Lamda * fkrdm.Ra + 3 * kSqrt2 * fkrdm.Ta) + 3 * kSqrt2 * fkrdm.Tv * (fQuV - fQdV) + 3 * kSqrt2 * fQuA * fkrdm.Ta);
   fAmpl.fPlus3 = -k3_Sqrt2 * (fkrdm.Tv * (fQdV - fQuV) + fkrdm.Ta * (fQdA - fQuA));
   fAmpl.f0Plus = -fkrdm.Lamda * k1_Sqrt3 * (3 * fkrdm.S * (fQuV - fQdV) + fkrdm.Cs * (fQdA + 5 * fQuA));
   fAmpl.f0Minus = fkrdm.Lamda * k1_Sqrt3 * (3 * fkrdm.S * (fQdV - fQuV) + fkrdm.Cs * (fQdA + 5 * fQuA));
   fAmpl.fz0Plus = -fkrdm.Lamda * fkrdm.Cz * k1_Sqrt3 * (fQdA + 5 * fQuA);
   fAmpl.fz0Minus = fkrdm.Lamda * fkrdm.Cz * k1_Sqrt3 * (fQdA + 5 * fQuA);

     break;
   }
   case (kS11_1650) :
   {
   fAmpl.fMinus3 = 0;
   fAmpl.fMinus1 = k1_Sqrt6 * fkrdm.Lamda *(fkrdm.Rv * (2.0 * fQdV + fQuV) - fkrdm.Ra * (2.0 * fQdA + fQuA));
   fAmpl.fPlus1 = k1_Sqrt6* fkrdm.Lamda *(fkrdm.Rv * (2.0 * fQdV + fQuV) + fkrdm.Ra * (2.0 * fQdA + fQuA));
   fAmpl.fPlus3 = 0;
   fAmpl.f0Plus = -kSqrt2_3 * (2.0 * fQdA + fQuA) * (3 * fkrdm.Bs - fkrdm.Lamda * fkrdm.Cs);
   fAmpl.f0Minus = -kSqrt2_3 * (2.0 * fQdA + fQuA) * (3 * fkrdm.Bs - fkrdm.Lamda * fkrdm.Cs);
   fAmpl.fz0Plus = -kSqrt2_3 * (2.0 * fQdA + fQuA) * (3 * fkrdm.Bz - fkrdm.Lamda * fkrdm.Cz);
   fAmpl.fz0Minus = -kSqrt2_3 * (2.0 * fQdA + fQuA) * (3 * fkrdm.Bz - fkrdm.Lamda * fkrdm.Cz);

     break;
   }
   case (kD13_1700) :
   {
   fAmpl.fMinus3 = (1.0/kSqrt10) * 3 * fkrdm.Lamda * (fkrdm.Rv * (2.0 * fQdV + fQuV) - fkrdm.Ra * (2.0 * fQdA + fQuA));
   fAmpl.fMinus1 = k1_Sqrt30 * fkrdm.Lamda * (fkrdm.Rv * (2.0 * fQdV + fQuV) - fkrdm.Ra * (2.0 * fQdA + fQuA));
   fAmpl.fPlus1 = -k1_Sqrt30 * fkrdm.Lamda * (fkrdm.Rv * (2.0 * fQdV + fQuV) + fkrdm.Ra * (2.0 * fQdA + fQuA));
   fAmpl.fPlus3 = -(1.0/kSqrt10) * 3 * fkrdm.Lamda * (fkrdm.Rv * (2.0 * fQdV + fQuV) + fkrdm.Ra * (2.0 * fQdA + fQuA));
   fAmpl.f0Plus = (kSqrt2 * k1_Sqrt15) * fkrdm.Lamda * fkrdm.Cs * (2.0 * fQdA + fQuA);
   fAmpl.f0Minus = -(kSqrt2 * k1_Sqrt15) * fkrdm.Lamda * fkrdm.Cs * (2.0 * fQdA + fQuA);
   fAmpl.fz0Plus = (kSqrt2 * k1_Sqrt15) * fkrdm.Lamda * fkrdm.Cz * (2.0 * fQdA + fQuA);
   fAmpl.fz0Minus = -(kSqrt2 * k1_Sqrt15) * fkrdm.Lamda * fkrdm.Cz * (2.0 * fQdA + fQuA);

     break;
   }
   case (kD15_1675) :
   {
   fAmpl.fMinus3 = kSqrt3_5 * fkrdm.Lamda * (fkrdm.Ra * (2.0 * fQdA + fQuA) - fkrdm.Rv * (2.0 * fQdV + fQuV));
   fAmpl.fMinus1 = kSqrt3_10 * fkrdm.Lamda * (fkrdm.Ra * (2.0 * fQdA + fQuA) - fkrdm.Rv * (2.0 * fQdV + fQuV));
   fAmpl.fPlus1 = -kSqrt3_10 * fkrdm.Lamda * (fkrdm.Ra * (2.0 * fQdA + fQuA) + fkrdm.Rv * (2.0 * fQdV + fQuV));
   fAmpl.fPlus3 = -kSqrt3_5 * fkrdm.Lamda * (fkrdm.Ra * (2.0 * fQdA + fQuA) + fkrdm.Rv * (2.0 * fQdV + fQuV));
   fAmpl.f0Plus = -(kSqrt3 * kSqrt2_5) * fkrdm.Lamda * fkrdm.Cs * (2.0 * fQdA + fQuA);
   fAmpl.f0Minus = -(kSqrt3 * kSqrt2_5) * fkrdm.Lamda * fkrdm.Cs * (2.0 * fQdA + fQuA);
   fAmpl.fz0Plus = -(kSqrt3 * kSqrt2_5) * fkrdm.Lamda * fkrdm.Cz * (2.0 * fQdA + fQuA);
   fAmpl.fz0Minus = -(kSqrt3 * kSqrt2_5) * fkrdm.Lamda * fkrdm.Cz * (2.0 * fQdA + fQuA);

     break;
   }
   case (kS31_1620) :
   {
   fAmpl.fMinus3 = 0;
   fAmpl.fMinus1 = k1_Sqrt6 * fkrdm.Lamda * (fkrdm.Ra * (fQdA - fQuA) + fkrdm.Rv * (fQuV - fQdV)) + kSqrt3 * fkrdm.Ta * (fQuA - fQdA) + kSqrt3 * fkrdm.Tv * (fQdV - fQuV);
   fAmpl.fPlus1 = k1_Sqrt6 * fkrdm.Lamda * (fkrdm.Ra * (fQuA - fQdA) + fkrdm.Rv * (fQuV - fQdV)) + kSqrt3 * fkrdm.Ta * (fQdA - fQuA) + kSqrt3 * fkrdm.Tv * (fQdV - fQuV);
   fAmpl.fPlus3 = 0;
   fAmpl.f0Plus = k1_Sqrt6 * (3 * fkrdm.Lamda * fkrdm.S * (fQuV - fQdV) + 3 * fkrdm.Bs * (fQuA - fQdA) + fkrdm.Lamda * fkrdm.Cs * (fQdA - fQuA));
   fAmpl.f0Minus = k1_Sqrt6 * (3 * fkrdm.Lamda * fkrdm.S * (fQdV - fQuV) + 3 * fkrdm.Bs * (fQuA - fQdA) + fkrdm.Lamda * fkrdm.Cs * (fQdA - fQuA));
   fAmpl.fz0Plus = -k1_Sqrt6 * (fQdA - fQuA) * (3 * fkrdm.Bz - fkrdm.Lamda * fkrdm.Cz);
   fAmpl.fz0Minus = -k1_Sqrt6 * (fQdA - fQuA) * (3 * fkrdm.Bz - fkrdm.Lamda * fkrdm.Cz);

     break;
   }
   case (kD33_1700) :
   {
   fAmpl.fMinus3 = -k3_Sqrt2 * (fkrdm.Tv * (fQuV - fQdV) + fkrdm.Ta * (fQdA - fQuA));
   fAmpl.fMinus1 = 0.5 * k1_Sqrt3 * ((fQuA - fQdA) * (2.0 * fkrdm.Lamda * fkrdm.Ra + 3 * kSqrt2 * fkrdm.Ta) + (fQdV - fQuV) * (2.0 * fkrdm.Lamda * fkrdm.Rv + 3 * kSqrt2 * fkrdm.Tv));
   fAmpl.fPlus1 = 0.5 * k1_Sqrt3 * ((fQuA - fQdA) * (2.0 * fkrdm.Lamda * fkrdm.Ra + 3 * kSqrt2 * fkrdm.Ta) - (fQdV - fQuV) * (2.0 * fkrdm.Lamda * fkrdm.Rv + 3 * kSqrt2 * fkrdm.Tv));
   fAmpl.fPlus3 = -k3_Sqrt2 * (fkrdm.Tv * (fQdV - fQuV) + fkrdm.Ta * (fQdA - fQuA));
   fAmpl.f0Plus = fkrdm.Lamda * k1_Sqrt3 * (3 * fkrdm.S * (fQdV - fQuV) + fkrdm.Cs * (fQuA - fQdA));
   fAmpl.f0Minus = fkrdm.Lamda * k1_Sqrt3 * (3 * fkrdm.S * (fQdV - fQuV) + fkrdm.Cs * (fQdA - fQuA));
   fAmpl.fz0Plus = fkrdm.Lamda * fkrdm.Cz * k1_Sqrt3 * (fQuA - fQdA);
   fAmpl.fz0Minus = fkrdm.Lamda * fkrdm.Cz * k1_Sqrt3 * (fQdA - fQuA);

     break;
   }
   case (kP11_1440) :
   {
     double L2  = TMath::Power(fkrdm.Lamda, 2);

     fAmpl.fMinus3 = 0;
     fAmpl.fMinus1 = 0.5 * k1_Sqrt3 * L2 * (fkrdm.Ra * (fQdA - 4 * fQuA) - fkrdm.Rv * (fQdV - 4 * fQuV));
     fAmpl.fPlus1 = 0.5 * k1_Sqrt3 * L2 * (fkrdm.Ra * (fQdA - 4 * fQuA) + fkrdm.Rv * (fQdV - 4 * fQuV));
     fAmpl.fPlus3 = 0;
     fAmpl.f0Plus = (fkrdm.Lamda/2.0) * k1_Sqrt3 * (-3 * fkrdm.Lamda * fkrdm.S * (fQdV + 2.0 * fQuV) - 2.0 * fkrdm.Bs * (fQdA - 4 * fQuA) + fkrdm.Lamda * fkrdm.Cs * (fQdA - 4 * fQuA));
     fAmpl.f0Minus = -(fkrdm.Lamda/2.0) * k1_Sqrt3 * (3 * fkrdm.Lamda * fkrdm.S * (fQdV + 2.0 * fQuV) - 2.0 * fkrdm.Bs * (fQdA - 4 * fQuA) + fkrdm.Lamda * fkrdm.Cs * (fQdA - 4 * fQuA));
     fAmpl.fz0Plus = (fkrdm.Lamda/2.0) * k1_Sqrt3 * ((fQdA - 4 * fQuA) * (fkrdm.Lamda * fkrdm.Cz - 2.0 * fkrdm.Bz));
     fAmpl.fz0Minus = -(fkrdm.Lamda/2.0) * k1_Sqrt3 * ((fQdA - 4 * fQuA) * (fkrdm.Lamda * fkrdm.Cz - 2.0 * fkrdm.Bz));

     break;
   }
   case (kP33_1600) :
   {
     double L2  = TMath::Power(fkrdm.Lamda, 2);

   fAmpl.fMinus3 = k1_Sqrt2 * L2 * (fkrdm.Ra * (fQuA - fQdA) + fkrdm.Rv * (fQdV - fQuV));
   fAmpl.fMinus1 = k1_Sqrt6 * L2 * (fkrdm.Ra * (fQuA - fQdA) + fkrdm.Rv * (fQdV - fQuV));
   fAmpl.fPlus1 = k1_Sqrt6 * L2 * (fkrdm.Ra * (fQdA - fQuA) + fkrdm.Rv * (fQdV - fQuV));
   fAmpl.fPlus3 = k1_Sqrt2 * L2 * (fkrdm.Ra * (fQdA - fQuA) + fkrdm.Rv * (fQdV - fQuV));
   fAmpl.f0Plus = kSqrt2_3 * fkrdm.Lamda * (fQdA - fQuA) * (fkrdm.Lamda * fkrdm.Cs - 2.0 * fkrdm.Bs);
   fAmpl.f0Minus = kSqrt2_3 * fkrdm.Lamda * (fQdA - fQuA) * (fkrdm.Lamda * fkrdm.Cs - 2.0 * fkrdm.Bs);
   fAmpl.fz0Plus = kSqrt2_3 * fkrdm.Lamda * (fQdA - fQuA) * (fkrdm.Lamda * fkrdm.Cz - 2.0 * fkrdm.Bz);
   fAmpl.fz0Minus = kSqrt2_3 * fkrdm.Lamda * (fQdA - fQuA) * (fkrdm.Lamda * fkrdm.Cz - 2.0 * fkrdm.Bz);

     break;
   }
   case (kP13_1720) :
   {
   fAmpl.fMinus3 = (1.0/kSqrt10) * 3 * fkrdm.Lamda * (fkrdm.Ta * (fQdA + 2.0 * fQuA) - fkrdm.Tv * (fQdV + 2.0 * fQuV));
   fAmpl.fMinus1 = k1_Sqrt15 * (fkrdm.Lamda/2.0) * (2.0 * fkrdm.Lamda * (fkrdm.Ra * (fQdA - 4 * fQuA) - fkrdm.Rv * (fQdV - 4 * fQuV)) + 9 * kSqrt2 * (fkrdm.Tv * (fQdV + 2.0 * fQuV) - fkrdm.Ta * (fQdA + 2.0 * fQuA)));
   fAmpl.fPlus1 = k1_Sqrt15 * (fkrdm.Lamda/2.0) * (9 * kSqrt2 * (fkrdm.Tv * (fQdV + 2.0 * fQuV) + fkrdm.Ta * (fQdA + 2.0 * fQuA)) - 2.0 * fkrdm.Lamda * (fkrdm.Ra * (fQdA - 4 * fQuA) + fkrdm.Rv * (fQdV - 4 * fQuV)));
   fAmpl.fPlus3 = -(1.0/kSqrt10) * 3 * fkrdm.Lamda * (fkrdm.Ta * (fQdA + 2.0 * fQuA) + fkrdm.Tv * (fQdV + 2.0 * fQuV));
   fAmpl.f0Plus = fkrdm.Lamda * k1_Sqrt15 * (-3 * fkrdm.Lamda * fkrdm.S * (fQdV + 2.0 * fQuV) - 5 * fkrdm.Bs * (fQdA - 4 * fQuA) + fkrdm.Lamda * fkrdm.Cs * (fQdA - 4 * fQuA));
   fAmpl.f0Minus = fkrdm.Lamda * k1_Sqrt15 * (3 * fkrdm.Lamda * fkrdm.S * (fQdV + 2.0 * fQuV) - 5 * fkrdm.Bs * (fQdA - 4 * fQuA) + fkrdm.Lamda * fkrdm.Cs * (fQdA - 4 * fQuA));
   fAmpl.fz0Plus = fkrdm.Lamda * k1_Sqrt15 * ((fQdA - 4 * fQuA) * (fkrdm.Lamda * fkrdm.Cz - 5 * fkrdm.Bz));
   fAmpl.fz0Minus = fkrdm.Lamda * k1_Sqrt15 * ((fQdA - 4 * fQuA) * (fkrdm.Lamda * fkrdm.Cz - 5 * fkrdm.Bz));

     break;
   }
   case (kF15_1680) :
   {
     double L2  = TMath::Power(fkrdm.Lamda, 2);

   fAmpl.fMinus3 = 3 * fkrdm.Lamda * kSqrt2_5 * (fkrdm.Tv * (fQdV + 2.0 * fQuV) - fkrdm.Ta * (fQdA + 2.0 * fQuA));
   fAmpl.fMinus1 = k1_Sqrt5 * (fkrdm.Lamda/2.0) * (kSqrt2 * fkrdm.Lamda * (fkrdm.Rv * (fQdV - 4 * fQuV) - fkrdm.Ra * (fQdA - 4 * fQuA)) + 6 * (fkrdm.Tv * (fQdV + 2.0 * fQuV) - fkrdm.Ta * (fQdA + 2.0 * fQuA)));
   fAmpl.fPlus1 = k1_Sqrt5 * (fkrdm.Lamda/2.0) * (-kSqrt2 * fkrdm.Lamda * (fkrdm.Rv * (fQdV - 4 * fQuV) + fkrdm.Ra * (fQdA - 4 * fQuA)) - 6 * (fkrdm.Tv * (fQdV + 2.0 * fQuV) + fkrdm.Ta * (fQdA + 2.0 * fQuA)));
   fAmpl.fPlus3 = -3 * fkrdm.Lamda * kSqrt2_5 * (fkrdm.Tv * (fQdV + 2.0 * fQuV) + fkrdm.Ta * (fQdA + 2.0 * fQuA));
   fAmpl.f0Plus = (1.0/kSqrt10) * L2 * (3 * fkrdm.S * (fQdV + 2.0 * fQuV) - fkrdm.Cs * (fQdA - 4 * fQuA));
   fAmpl.f0Minus = (1.0/kSqrt10) * L2 * (3 * fkrdm.S * (fQdV + 2.0 * fQuV) + fkrdm.Cs * (fQdA - 4 * fQuA));
   fAmpl.fz0Plus = -(1.0/kSqrt10) * L2 * fkrdm.Cz * (fQdA - 4 * fQuA);
   fAmpl.fz0Minus = (1.0/kSqrt10) * L2 * fkrdm.Cz * (fQdA - 4 * fQuA);

     break;
   }
   case (kP31_1910) :
   {
     double L2  = TMath::Power(fkrdm.Lamda, 2);

   fAmpl.fMinus3 = 0;
   fAmpl.fMinus1 = k1_Sqrt15 * L2 * (fkrdm.Ra * (fQdA - fQuA) + fkrdm.Rv * (fQuV - fQdV));
   fAmpl.fPlus1 = k1_Sqrt15 * L2 * (fkrdm.Ra * (fQdA - fQuA) + fkrdm.Rv * (fQdV - fQuV));
   fAmpl.fPlus3 = 0;
   fAmpl.f0Plus = -k1_Sqrt15 * 2.0 * fkrdm.Lamda * (fQdA - fQuA) * (fkrdm.Lamda * fkrdm.Cs - 5 * fkrdm.Bs);
   fAmpl.f0Minus = k1_Sqrt15 * 2.0 * fkrdm.Lamda * (fQdA - fQuA) * (fkrdm.Lamda * fkrdm.Cs - 5 * fkrdm.Bs);
   fAmpl.fz0Plus = -k1_Sqrt15 * 2.0 * fkrdm.Lamda * (fQdA - fQuA) * (fkrdm.Lamda * fkrdm.Cz - 5 * fkrdm.Bz);
   fAmpl.fz0Minus = k1_Sqrt15 * 2.0 * fkrdm.Lamda * (fQdA - fQuA) * (fkrdm.Lamda * fkrdm.Cz - 5 * fkrdm.Bz);

     break;
   }
   case (kP33_1920) :
   {
     double L2  = TMath::Power(fkrdm.Lamda, 2);

   fAmpl.fMinus3 = k1_Sqrt5 * L2 * (fkrdm.Ra * (fQdA - fQuA) + fkrdm.Rv * (fQuV - fQdV));
   fAmpl.fMinus1 = k1_Sqrt15 * L2 * (fkrdm.Ra * (fQuA - fQdA) + fkrdm.Rv * (fQdV - fQuV));
   fAmpl.fPlus1 = k1_Sqrt15 * L2 * (fkrdm.Ra * (fQdA - fQuA) + fkrdm.Rv * (fQdV - fQuV));
   fAmpl.fPlus3 = k1_Sqrt5 * L2 * (fkrdm.Ra * (fQuA - fQdA) + fkrdm.Rv * (fQuV - fQdV));
   fAmpl.f0Plus = k1_Sqrt15 * 2.0 * fkrdm.Lamda * (fQdA - fQuA) * (fkrdm.Lamda * fkrdm.Cs - 5 * fkrdm.Bs);
   fAmpl.f0Minus = k1_Sqrt15 * 2.0 * fkrdm.Lamda * (fQdA - fQuA) * (fkrdm.Lamda * fkrdm.Cs - 5 * fkrdm.Bs);
   fAmpl.fz0Plus = k1_Sqrt15 * 2.0 * fkrdm.Lamda * (fQdA - fQuA) * (fkrdm.Lamda * fkrdm.Cz - 5 * fkrdm.Bz);
   fAmpl.fz0Minus = k1_Sqrt15 * 2.0 * fkrdm.Lamda * (fQdA - fQuA) * (fkrdm.Lamda * fkrdm.Cz - 5 * fkrdm.Bz);

     break;
   }
   case (kF35_1905) :
   {
     double L2  = TMath::Power(fkrdm.Lamda, 2);

   fAmpl.fMinus3 = -3 * L2 * (kSqrt2 * k1_Sqrt35) * (fkrdm.Ra * (fQdA - fQuA) + fkrdm.Rv * (fQuV - fQdV));
   fAmpl.fMinus1 = L2 * k1_Sqrt35 * (fkrdm.Ra * (fQuA - fQdA) + fkrdm.Rv * (fQdV - fQuV));
   fAmpl.fPlus1 = L2 * k1_Sqrt35 * (fkrdm.Ra * (fQuA - fQdA) + fkrdm.Rv * (fQuV - fQdV));
   fAmpl.fPlus3 = -3 * L2 * (kSqrt2 * k1_Sqrt35) * (fkrdm.Ra * (fQdA - fQuA) + fkrdm.Rv * (fQdV - fQuV));
   fAmpl.f0Plus = 2.0 * L2 * k1_Sqrt35 * fkrdm.Cs * (fQdA - fQuA);
   fAmpl.f0Minus = 2.0 * L2 * k1_Sqrt35 * fkrdm.Cs * (fQuA - fQdA);
   fAmpl.fz0Plus = 2.0 * L2 * k1_Sqrt35 * fkrdm.Cz * (fQdA - fQuA);
   fAmpl.fz0Minus = 2.0 * L2 * k1_Sqrt35 * fkrdm.Cz * (fQuA - fQdA);

     break;
   }
   case (kF37_1950) :
   {
     double L2  = TMath::Power(fkrdm.Lamda, 2);

   fAmpl.fMinus3 = kSqrt2_7 * L2 * (fkrdm.Ra * (fQdA - fQuA) + fkrdm.Rv * (fQuV - fQdV));
   fAmpl.fMinus1 = kSqrt6_35 * L2 * (fkrdm.Ra * (fQdA - fQuA) + fkrdm.Rv * (fQuV - fQdV));
   fAmpl.fPlus1 = kSqrt6_35 * L2 * (fkrdm.Ra * (fQuA - fQdA) + fkrdm.Rv * (fQuV - fQdV));
   fAmpl.fPlus3 = kSqrt2_7 * L2 * (fkrdm.Ra * (fQuA - fQdA) + fkrdm.Rv * (fQuV - fQdV));
   fAmpl.f0Plus = 2.0 * kSqrt6_35 * L2 * fkrdm.Cs * (fQuA - fQdA);
   fAmpl.f0Minus = 2.0 * kSqrt6_35 * L2 * fkrdm.Cs * (fQuA - fQdA);
   fAmpl.fz0Plus = 2.0 * kSqrt6_35 * L2 * fkrdm.Cz * (fQuA - fQdA);
   fAmpl.fz0Minus = 2.0 * kSqrt6_35 * L2 * fkrdm.Cz * (fQuA - fQdA);

     break;
   }
   case (kP11_1710) :
   {
     double L2  = TMath::Power(fkrdm.Lamda, 2);

   fAmpl.fMinus3 = 0;
   fAmpl.fMinus1 = 0.5 * k1_Sqrt6 * L2 * (fkrdm.Ra * (fQdA + 5 * fQuA) - fkrdm.Rv * (fQdV + 5 * fQuV));
   fAmpl.fPlus1 = 0.5 * k1_Sqrt6 * L2 * (fkrdm.Ra * (fQdA + 5 * fQuA) + fkrdm.Rv * (fQdV + 5 * fQuV));
   fAmpl.fPlus3 = 0;
   fAmpl.f0Plus = fkrdm.Lamda * 0.5 * k1_Sqrt6 * (3 * fkrdm.Lamda * fkrdm.S * (fQuV - fQdV) + (fQdA + 5 * fQuA) * (fkrdm.Lamda * fkrdm.Cs - 2.0 * fkrdm.Bs));
   fAmpl.f0Minus = -fkrdm.Lamda * 0.5 * k1_Sqrt6 * (3 * fkrdm.Lamda * fkrdm.S * (fQdV - fQuV) + (fQdA + 5 * fQuA) * (fkrdm.Lamda * fkrdm.Cs - 2.0 * fkrdm.Bs));
   fAmpl.fz0Plus = fkrdm.Lamda * 0.5 * k1_Sqrt6 * (fQdA + 5 * fQuA) * (fkrdm.Lamda * fkrdm.Cz - 2.0 * fkrdm.Bz);
   fAmpl.fz0Minus = -fkrdm.Lamda * 0.5 * k1_Sqrt6 * (fQdA + 5 * fQuA) * (fkrdm.Lamda * fkrdm.Cz - 2.0 * fkrdm.Bz);

     break;
   }
   case (kF17_1970) :
   {
     double L2  = TMath::Power(fkrdm.Lamda, 2);

   fAmpl.fMinus3 = (k1_Sqrt2 * k1_Sqrt7) * L2 * (fkrdm.Ra * (2.0 * fQdA + fQuA) - fkrdm.Rv * (2.0 * fQdV + fQuV));
   fAmpl.fMinus1 = (kSqrt3_2 * k1_Sqrt35) * L2 * (fkrdm.Ra * (2.0 * fQdA + fQuA) - fkrdm.Rv * (2.0 * fQdV + fQuV));
   fAmpl.fPlus1 = -(kSqrt3_2 * k1_Sqrt35) * L2 * (fkrdm.Ra * (2.0 * fQdA + fQuA) + fkrdm.Rv * (2.0 * fQdV + fQuV));
   fAmpl.fPlus3 = -(k1_Sqrt2 * k1_Sqrt7) * L2 * (fkrdm.Ra * (2.0 * fQdA + fQuA) + fkrdm.Rv * (2.0 * fQdV + fQuV));
   fAmpl.f0Plus = -kSqrt6_35 * L2 * fkrdm.Cs * (2.0 * fQdA + fQuA);
   fAmpl.f0Minus = -kSqrt6_35 * L2 * fkrdm.Cs * (2.0 * fQdA + fQuA);
   fAmpl.fz0Plus = -kSqrt6_35 * L2 * fkrdm.Cz * (2.0 * fQdA + fQuA);
   fAmpl.fz0Minus = -kSqrt6_35 * L2 * fkrdm.Cz * (2.0 * fQdA + fQuA);

     break;
   }
   default:
   {
     LOG("DMHAmpl", pWARN) << "*** UNRECOGNIZED RESONANCE!";
     fAmpl.fMinus1 = 0.;
     fAmpl.fPlus1  = 0.;
     fAmpl.fMinus3 = 0.;
     fAmpl.fPlus3  = 0.;
     fAmpl.f0Minus = 0.;
     fAmpl.f0Plus  = 0.;
     //New terms specific to DM
     fAmpl.fz0Minus = 0.;
     fAmpl.fz0Plus  = 0.;
     break;
   }

  }//switch

  return fAmpl;
  }
  //___________________________________________________________________________
