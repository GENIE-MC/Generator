//____________________________________________________________________________
/*!

\class    genie::RSHelicityAmplModelEMn

\brief    The Helicity Amplitudes, for all baryon resonances, for Electro-
          Magnetic (EM) interactions on free neutrons, as computed in the
          Rein-Seghal's paper.

          Concrete implementation of the SPPHelicityAmplModelI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  March 30, 2005

*/
//____________________________________________________________________________

#include "BaryonResonance/BaryonResUtils.h"
#include "Conventions/Constants.h"
#include "ReinSeghal/RSHelicityAmplModelEMn.h"
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
double RSHelicityAmplModelEMn::AmpMinus1(
                       const Interaction * interaction, const FKR & fkr) const
{
  double ampl = 0;
  Resonance_t res = interaction->GetExclusiveTag().Resonance();

  switch(res) {

  case kP33_1232 :
     ampl = -kSqrt2 * fkr.R();
     break;

  case kS11_1535 :
     ampl = -kSqrt3 * fkr.T() - (1./kSqrt6) * fkr.LR();
     break;

  case kD13_1520 :
     ampl = -(kSqrt3/kSqrt2) * fkr.T() + (1./kSqrt3) * fkr.LR();
     break;

  case kS11_1650 :
     ampl = -(1./kSqrt6) * fkr.LR();
     break;

  case kD13_1700 :
     ampl = -(1./kSqrt30) * fkr.LR();
     break;

  case kD15_1675 :
     ampl = (kSqrt3/kSqrt10) * fkr.LR();
     break;

  case kS31_1620 :
     ampl = kSqrt3 * fkr.T() - (1./kSqrt6) * fkr.LR();
     break;

  case kD33_1700 :
     ampl = (kSqrt3/kSqrt2) * fkr.T() + (1./kSqrt3) * fkr.LR();
     break;

  case kP11_1440 :
     ampl = (1./kSqrt3) * fkr.L2R();
     break;

  case kP33_1600 :
     ampl = (1./kSqrt6) * fkr.L2R();
     break;

  case kP13_1720 :
     ampl = (2./kSqrt15) * fkr.L2R();
     break;

  case kF15_1680 :
     ampl = -(kSqrt2/kSqrt5) * fkr.L2R();
     break;

  case kP31_1910 :
     ampl = -(1./kSqrt15) * fkr.L2R();
     break;

  case kP33_1920 :
     ampl =  (1./kSqrt15) * fkr.L2R();
     break;

  case kF35_1905 :
     ampl =  (1./kSqrt35) * fkr.L2R();
     break;

  case kF37_1950 :
     ampl = -(kSqrt6/kSqrt35) * fkr.L2R();
     break;

  case kP11_1710 :
     ampl = -(1./kSqrt24) * fkr.L2R();
     break;

  case kF17_1970 :
     ampl = (kSqrt3/kSqrt35) * fkr.L2R();
     break;

  default:
     LOG("RSHAmpl", pWARN) << "*** UNRECOGNIZED RESONANCE!";
     ampl = 0;
  }

  LOG("RSHAmpl", pDEBUG)
             << "A(-1) for RES: " << utils::res::AsString(res) << " = " << ampl;

  return ampl;
}
//____________________________________________________________________________
double RSHelicityAmplModelEMn::AmpPlus1(
                       const Interaction * interaction, const FKR & fkr) const
{
  double ampl = 0;
  Resonance_t res = interaction->GetExclusiveTag().Resonance();

  switch(res) {

  case kP33_1232 :
     ampl = kSqrt2 * fkr.R();
     break;

  case kS11_1535 :
     ampl = kSqrt3 * fkr.T() + (1./kSqrt6) * fkr.LR();
     break;

  case kD13_1520 :
     ampl = -(kSqrt3/kSqrt2) * fkr.T() + (1./kSqrt3) * fkr.LR();
     break;

  case kS11_1650 :
     ampl = (1./kSqrt6) * fkr.LR();
     break;

  case kD13_1700 :
     ampl = -(1./kSqrt30) * fkr.LR();
     break;

  case kD15_1675 :
     ampl = -(kSqrt3/kSqrt10) * fkr.LR();
     break;

  case kS31_1620 :
     ampl = -kSqrt3 * fkr.T() + (1./kSqrt6) * fkr.LR();
     break;

  case kD33_1700 :
     ampl = (kSqrt3/kSqrt2)* fkr.T() + (1./kSqrt3) * fkr.LR();
     break;

  case kP11_1440 :
     ampl = (1./kSqrt3)*fkr.L2R();
     break;

  case kP33_1600 :
     ampl = -(1./kSqrt6)*fkr.L2R();
     break;

  case kP13_1720 :
     ampl = -(2./kSqrt15) * fkr.L2R();
     break;

  case kF15_1680 :
     ampl = -(kSqrt2/kSqrt5) * fkr.L2R();
     break;

  case kP31_1910 :
     ampl = -(1./kSqrt15) * fkr.L2R();
     break;

  case kP33_1920 :
     ampl = -(1./kSqrt15) * fkr.L2R();
     break;

  case kF35_1905 :
     ampl = (1./kSqrt35)  * fkr.L2R();
     break;

  case kF37_1950 :
     ampl = (kSqrt6/kSqrt35) * fkr.L2() * fkr.Rplus_2wR();
     break;

  case kP11_1710 :
     ampl = -(1./kSqrt24) * fkr.L2R();
     break;

  case kF17_1970 :
     ampl = -(kSqrt3/kSqrt35) * fkr.L2R();
     break;

  default:
     LOG("RSHAmpl", pWARN) << "*** UNRECOGNIZED RESONANCE!";
     ampl = 0;
  }

  LOG("RSHAmpl", pDEBUG)
           << "A(+1) for RES: " << utils::res::AsString(res) << " = " << ampl;

  return ampl;
}
//____________________________________________________________________________
double RSHelicityAmplModelEMn::AmpMinus3(
                       const Interaction * interaction, const FKR & fkr) const
{
  double ampl = 0;
  Resonance_t res = interaction->GetExclusiveTag().Resonance();

  switch(res) {

  case kP33_1232 :
     ampl = -1 * kSqrt6 * fkr.R();
     break;

  case kS11_1535 :
     ampl = 0;
     break;

  case kD13_1520 :
     ampl = (3./kSqrt2) * fkr.T();
     break;

  case kS11_1650 :
     ampl = 0;
     break;

  case kD13_1700 :
     ampl = -(3./kSqrt10) * fkr.LR();
     break;

  case kD15_1675 :
     ampl = (kSqrt3/kSqrt5) * fkr.LR();
     break;

  case kS31_1620 :
     ampl = 0;
     break;

  case kD33_1700 :
     ampl = (3./kSqrt2) * fkr.T();
     break;

  case kP11_1440 :
     ampl = 0;
     break;

  case kP33_1600 :
     ampl = (1./kSqrt2) * fkr.L2R();
     break;

  case kP13_1720 :
     ampl = 0;
     break;

  case kF15_1680 :
     ampl = 0;
     break;

  case kP31_1910 :
     ampl = 0;
     break;

  case kP33_1920 :
     ampl = -(1./kSqrt5) * fkr.L2R();
     break;

  case kF35_1905 :
     ampl = (kSqrt18/kSqrt35) * fkr.L2R();
     break;

  case kF37_1950 :
     ampl = -(kSqrt2/kSqrt7) * fkr.L2R();
     break;

  case kP11_1710 :
     ampl = 0;
     break;

  case kF17_1970 :
     ampl = (1./kSqrt7) * fkr.L2R();
     break;

  default:
     LOG("RSHAmpl", pWARN) << "*** UNRECOGNIZED RESONANCE!";
     ampl = 0;
  }

  LOG("RSHAmpl", pDEBUG)
            << "A(-3) for RES: " << utils::res::AsString(res) << " = " << ampl;

  return ampl;
}
//____________________________________________________________________________
double RSHelicityAmplModelEMn::AmpPlus3(
                       const Interaction * interaction, const FKR & fkr) const
{
  double ampl = 0;
  Resonance_t res = interaction->GetExclusiveTag().Resonance();

  switch(res) {

  case kP33_1232 :
     ampl = kSqrt6 * fkr.R();
     break;

  case kS11_1535 :
     ampl = 0;
     break;

  case kD13_1520 :
     ampl = -(3./kSqrt2) * fkr.T();
     break;

  case kS11_1650 :
     ampl = 0;
     break;

  case kD13_1700 :
     ampl = -(3./kSqrt10) * fkr.LR();
     break;

  case kD15_1675 :
     ampl = -(kSqrt3/kSqrt5) * fkr.LR();
     break;

  case kS31_1620 :
     ampl = 0;
     break;

  case kD33_1700 :
     ampl = (3./kSqrt2) * fkr.T();
     break;

  case kP11_1440 :
     ampl = 0;
     break;

  case kP33_1600 :
     ampl = -(1./kSqrt2) * fkr.L2R();
     break;

  case kP13_1720 :
     ampl = 0;
     break;

  case kF15_1680 :
     ampl = 0;
     break;

  case kP31_1910 :
     ampl = 0;
     break;

  case kP33_1920 :
     ampl = (1./kSqrt5) * fkr.L2R();
     break;

  case kF35_1905 :
     ampl = (kSqrt18/kSqrt35) * fkr.L2R();
     break;

  case kF37_1950 :
     ampl = (kSqrt2/kSqrt7) * fkr.L2R();
     break;

  case kP11_1710 :
     ampl = 0;
     break;

  case kF17_1970 :
     ampl = 0;
     break;

  default:
     LOG("RSHAmpl", pWARN) << "*** UNRECOGNIZED RESONANCE!";
     ampl = 0;
  }

  LOG("RSHAmpl", pDEBUG)
            << "A(+3) for RES: " << utils::res::AsString(res) << " = " << ampl;

  return ampl;
}
//____________________________________________________________________________
double RSHelicityAmplModelEMn::Amp0Minus(
                       const Interaction * interaction, const FKR & fkr) const
{
  double ampl = 0;
  Resonance_t res = interaction->GetExclusiveTag().Resonance();

  switch(res) {

  case kP33_1232 :
     ampl = 0;
     break;

  case kS11_1535 :
     ampl = (kSqrt3/kSqrt2) * fkr.LS();
     break;

  case kD13_1520 :
     ampl = (kSqrt3) * fkr.LS();
     break;

  case kS11_1650 :
     ampl = 0;
     break;

  case kD13_1700 :
     ampl = 0;
     break;

  case kD15_1675 :
     ampl = 0;
     break;

  case kS31_1620 :
     ampl = -(kSqrt3/kSqrt2) * fkr.LS();
     break;

  case kD33_1700 :
     ampl = -kSqrt3 * fkr.LS();
     break;

  case kP11_1440 :
     ampl = 0;
     break;

  case kP33_1600 :
     ampl = 0;
     break;

  case kP13_1720 :
     ampl = 0;
     break;

  case kF15_1680 :
     ampl = 0;
     break;

  case kP31_1910 :
     ampl = 0;
     break;

  case kP33_1920 :
     ampl = 0;
     break;

  case kF35_1905 :
     ampl = 0;
     break;

  case kF37_1950 :
     ampl = 0;
     break;

  case kP11_1710 :
     ampl = -(kSqrt3/kSqrt8) * fkr.L2S();
     break;

  case kF17_1970 :
     ampl = 0;
     break;

  default:
     LOG("RSHAmpl", pWARN) << "*** UNRECOGNIZED RESONANCE!";
     ampl = 0;
  }

  LOG("RSHAmpl", pDEBUG)
            << "A(0-) for RES: " << utils::res::AsString(res) << " = " << ampl;

  return ampl;
}
//____________________________________________________________________________
double RSHelicityAmplModelEMn::Amp0Plus(
                       const Interaction * interaction, const FKR & fkr) const
{
  double ampl = 0;
  Resonance_t res = interaction->GetExclusiveTag().Resonance();

  switch(res) {

  case kP33_1232 :
     ampl = 0;
     break;

  case kS11_1535 :
     ampl = -(kSqrt3/kSqrt2) * fkr.LS();
     break;

  case kD13_1520 :
     ampl = (kSqrt3) * fkr.LS();
     break;

  case kS11_1650 :
     ampl = 0;
     break;

  case kD13_1700 :
     ampl = 0;
     break;

  case kD15_1675 :
     ampl = 0;
     break;

  case kS31_1620 :
     ampl = (kSqrt3/kSqrt2) * fkr.LS();
     break;

  case kD33_1700 :
     ampl = -kSqrt3 * fkr.L2S();
     break;

  case kP11_1440 :
     ampl = 0;
     break;

  case kP33_1600 :
     ampl = 0;
     break;

  case kP13_1720 :
     ampl = 0;
     break;

  case kF15_1680 :
     ampl = 0;
     break;

  case kP31_1910 :
     ampl = 0;
     break;

  case kP33_1920 :
     ampl = 0;
     break;

  case kF35_1905 :
     ampl = 0;
     break;

  case kF37_1950 :
     ampl = 0;
     break;

  case kP11_1710 :
     ampl = -(kSqrt3/kSqrt8) * fkr.L2S();
     break;

  case kF17_1970 :
     ampl = 0;
     break;

  default:
     LOG("RSHAmpl", pWARN) << "*** UNRECOGNIZED RESONANCE!";
     ampl = 0;
  }

  LOG("RSHAmpl", pDEBUG)
            << "A(0+) for RES: " << utils::res::AsString(res) << " = " << ampl;

  return ampl;
}
//____________________________________________________________________________
