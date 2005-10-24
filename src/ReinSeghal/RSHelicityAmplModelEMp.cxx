//____________________________________________________________________________
/*!

\class    genie::RSHelicityAmplModelEMp

\brief    The Helicity Amplitudes, for all baryon resonances, for Electro-
          Magnetic (EM) interactions on free protons, as computed in the
          Rein-Seghal's paper.

          Concrete implementation of the SPPHelicityAmplModelI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  March 30, 2005

*/
//____________________________________________________________________________

#include "BaryonResonance/BaryonResUtils.h"
#include "Conventions/Constants.h"
#include "ReinSeghal/RSHelicityAmplModelEMp.h"
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
double RSHelicityAmplModelEMp::AmpMinus1(
                       const Interaction * interaction, const FKR & fkr) const
{
  double ampl = 0;

  Resonance_t res = utils::res::FromInteraction(interaction);

  switch(res) {

  case kP33_1232 :
     ampl = -1 * kSqrt2 * fkr.R();
     break;

  case kS11_1535 :
     ampl = kSqrt3* fkr.T() + (kSqrt3/kSqrt2) * fkr.LR();
     break;

  case kD13_1520 :
     ampl = (kSqrt3/kSqrt2) * fkr.T() - kSqrt3 * fkr.LR();
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
     ampl = kSqrt3* fkr.T() - (1./kSqrt6)* fkr.LR();
     break;

  case kD33_1700 :
     ampl = (kSqrt3/kSqrt2)* fkr.T() + (1./kSqrt3) * fkr.LR();
     break;

  case kP11_1440 :
     ampl = -0.5*kSqrt3*fkr.L2R();
     break;

  case kP33_1600 :
     ampl = (1./kSqrt6) * fkr.L2R();
     break;

  case kP13_1720 :
     ampl = -(kSqrt27/kSqrt10) * fkr.LT() - (kSqrt3/kSqrt5) * fkr.L2R();
     break;

  case kF15_1680 :
     ampl = -(3./kSqrt5) * fkr.LT() + (kSqrt3/kSqrt10) * fkr.L2R();
     break;

  case kP31_1910 :
     ampl = -(1./kSqrt15) * fkr.L2R();
     break;

  case kP33_1920 :
     ampl = (1./kSqrt15) * fkr.L2R();
     break;

  case kF35_1905 :
     ampl = (1./kSqrt35) * fkr.L2R();
     break;

  case kF37_1950 :
     ampl = -(kSqrt6/kSqrt35) * fkr.L2R();
     break;

  case kP11_1710 :
     ampl = (kSqrt3/kSqrt8) * fkr.L2R();
     break;

  case kF17_1970 :
     ampl = 0;
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
double RSHelicityAmplModelEMp::AmpPlus1(
                       const Interaction * interaction, const FKR & fkr) const
{
  double ampl = 0;

  Resonance_t res = utils::res::FromInteraction(interaction);

  switch(res) {

  case kP33_1232 :
     ampl  = kSqrt2 * fkr.R();
     break;

  case kS11_1535 :
     ampl  = -1*kSqrt3 * fkr.T() - (kSqrt3/kSqrt2) * fkr.LR();
     break;

  case kD13_1520 :
     ampl  = (kSqrt3/kSqrt2) * fkr.T() - kSqrt3 * fkr.LR();
     break;

  case kS11_1650 :
     ampl  = 0;
     break;

  case kD13_1700 :
     ampl  = 0;
     break;

  case kD15_1675 :
     ampl  = 0;
     break;

  case kS31_1620 :
     ampl  = -1*kSqrt3* fkr.T() + (1./kSqrt6)* fkr.LR();
     break;

  case kD33_1700 :
     ampl  = (kSqrt3/kSqrt2) * fkr.T() + (1./kSqrt3) * fkr.LR();
     break;

  case kP11_1440 :
     ampl  = 0.5 * kSqrt3 * fkr.L2R();
     break;

  case kP33_1600 :
     ampl  = -(1./kSqrt6) * fkr.L2R();
     break;

  case kP13_1720 :
     ampl  = (kSqrt27/kSqrt10) * fkr.LT() + (kSqrt3/kSqrt5) * fkr.L2R();
     break;

  case kF15_1680 :
     ampl  = -(3./kSqrt5) * fkr.LT() + (3./kSqrt10) * fkr.L2R();
     break;

  case kP31_1910 :
     ampl  = -(1./kSqrt15) * fkr.L2R();
     break;

  case kP33_1920 :
     ampl  = -(1./kSqrt15) * fkr.L2R();
     break;

  case kF35_1905 :
     ampl  = (1./kSqrt35) * fkr.L2R();
     break;

  case kF37_1950 :
     ampl  = (kSqrt6/kSqrt35) * fkr.L2R();
     break;

  case kP11_1710 :
     ampl  = -(kSqrt3/kSqrt8) * fkr.L2R();
     break;

  case kF17_1970 :
     ampl  = 0;
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
double RSHelicityAmplModelEMp::AmpMinus3(
                       const Interaction * interaction, const FKR & fkr) const
{
  double ampl = 0;

  Resonance_t res = utils::res::FromInteraction(interaction);

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
     ampl = 0;
     break;

  case kD15_1675 :
     ampl = 0;
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
     ampl = (3./kSqrt10) * fkr.LT();
     break;

  case kF15_1680 :
     ampl = -1*(kSqrt18/kSqrt5) * fkr.LT();
     break;

  case kP31_1910 :
     ampl = 0;
     break;

  case kP33_1920 :
     ampl = -(1./kSqrt5 ) * fkr.L2R();
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
     ampl = 0;
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
double RSHelicityAmplModelEMp::AmpPlus3(
                       const Interaction * interaction, const FKR & fkr) const
{
  double ampl = 0;

  Resonance_t res = utils::res::FromInteraction(interaction);

  switch(res) {

  case kP33_1232 :
     ampl  = kSqrt6 * fkr.R();
     break;

  case kS11_1535 :
     ampl  = 0;
     break;

  case kD13_1520 :
     ampl  = (3./kSqrt2) * fkr.T();
     break;

  case kS11_1650 :
     ampl  = 0;
     break;

  case kD13_1700 :
     ampl  = 0;
     break;

  case kD15_1675 :
     ampl  = 0;
     break;

  case kS31_1620 :
     ampl  = 0;
     break;

  case kD33_1700 :
     ampl  = (3./kSqrt2) * fkr.T();
     break;

  case kP11_1440 :
     ampl  = 0;
     break;

  case kP33_1600 :
     ampl  = -(1./kSqrt2) * fkr.L2R();
     break;

  case kP13_1720 :
     ampl  = -(3./kSqrt10) * fkr.LT();
     break;

  case kF15_1680 :
     ampl  = -1*(kSqrt18/kSqrt5) * fkr.LT();
     break;

  case kP31_1910 :
     ampl  = 0;
     break;

  case kP33_1920 :
     ampl  = (1./kSqrt5 ) * fkr.L2R();
     break;

  case kF35_1905 :
     ampl  = (kSqrt18/kSqrt35) * fkr.L2R();
     break;

  case kF37_1950 :
     ampl  = (kSqrt2/kSqrt7) * fkr.L2R();
     break;

  case kP11_1710 :
     ampl  = 0;
     break;

  case kF17_1970 :
     ampl  = 0;
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
double RSHelicityAmplModelEMp::Amp0Minus(
                       const Interaction * interaction, const FKR & fkr) const
{
  double ampl = 0;

  Resonance_t res = utils::res::FromInteraction(interaction);

  switch(res) {

  case kP33_1232 :
     ampl = 0;
     break;

  case kS11_1535 :
     ampl = -1*(kSqrt3/kSqrt2) * fkr.LS();
     break;

  case kD13_1520 :
     ampl = -1*kSqrt3 * fkr.LS();
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
     ampl = -(kSqrt3/kSqrt2)* fkr.LS();
     break;

  case kD33_1700 :
     ampl = -kSqrt3 * fkr.LS();
     break;

  case kP11_1440 :
     ampl = 0.5*kSqrt3 * fkr.L2S();
     break;

  case kP33_1600 :
     ampl = 0;
     break;

  case kP13_1720 :
     ampl = (kSqrt3/kSqrt5) * fkr.L2S();
     break;

  case kF15_1680 :
     ampl = (3./kSqrt10)* fkr.L2S();
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
     ampl = (kSqrt3/kSqrt8) * fkr.L2S();
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
double RSHelicityAmplModelEMp::Amp0Plus(
                       const Interaction * interaction, const FKR & fkr) const
{
  double ampl = 0;

  Resonance_t res = utils::res::FromInteraction(interaction);

  switch(res) {

  case kP33_1232 :
     ampl  = 0;
     break;

  case kS11_1535 :
     ampl  = (kSqrt3/kSqrt2) * fkr.LS();
     break;

  case kD13_1520 :
     ampl  = -kSqrt3 * fkr.LS();
     break;

  case kS11_1650 :
     ampl  = 0;
     break;

  case kD13_1700 :
     ampl  = 0;
     break;

  case kD15_1675 :
     ampl  = 0;
     break;

  case kS31_1620 :
     ampl  = (kSqrt3/kSqrt2)* fkr.LS();
     break;

  case kD33_1700 :
     ampl  = -kSqrt3 * fkr.LS();
     break;

  case kP11_1440 :
     ampl  = -0.5*(kSqrt3) * fkr.L2S();
     break;

  case kP33_1600 :
     ampl  = 0;
     break;

  case kP13_1720 :
     ampl  = -1*(kSqrt3/kSqrt5) * fkr.L2S();
     break;

  case kF15_1680 :
     ampl  = (3./kSqrt10) * fkr.L2S();
     break;

  case kP31_1910 :
     ampl  = 0;
     break;

  case kP33_1920 :
     ampl  = 0;
     break;

  case kF35_1905 :
     ampl  = 0;
     break;

  case kF37_1950 :
     ampl  = 0;
     break;

  case kP11_1710 :
     ampl  = (kSqrt3/kSqrt8) * fkr.L2S();
     break;

  case kF17_1970 :
     ampl  = 0;
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


