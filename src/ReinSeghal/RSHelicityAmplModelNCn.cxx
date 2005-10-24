//____________________________________________________________________________
/*!

\class    genie::RSHelicityAmplModelNCn

\brief    The Helicity Amplitudes, for all baryon resonances, for NC neutrino
          interactions on free neutrons, as computed in the Rein-Seghal's paper.

          Concrete implementation of the SPPHelicityAmplModelI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004

*/
//____________________________________________________________________________

#include "BaryonResonance/BaryonResUtils.h"
#include "Conventions/Constants.h"
#include "ReinSeghal/RSHelicityAmplModelNCn.h"
#include "Messenger/Messenger.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
RSHelicityAmplModelNCn::RSHelicityAmplModelNCn() :
RSHelicityAmplModelI("genie::RSHelicityAmplModelNCn")
{

}
//____________________________________________________________________________
RSHelicityAmplModelNCn::RSHelicityAmplModelNCn(string config) :
RSHelicityAmplModelI("genie::RSHelicityAmplModelNCn", config)
{

}
//____________________________________________________________________________
RSHelicityAmplModelNCn::~RSHelicityAmplModelNCn()
{

}
//____________________________________________________________________________
double RSHelicityAmplModelNCn::AmpMinus1(
                       const Interaction * interaction, const FKR & fkr) const
{
  double ampl = 0;
  double xi   = kSin8w_2;

  Resonance_t res = utils::res::FromInteraction(interaction);

  switch(res) {

  case kP33_1232 :
     ampl = -1 * kSqrt2 * fkr.Rminus_2wR();
     break;

  case kS11_1535 :
     ampl = -1*kSqrt3* fkr.Tminus_2wTv() -
                           (kSqrt2/kSqrt3) * fkr.Lamda() * fkr.Rminus_wR();
     break;

  case kD13_1520 :
     ampl = -(kSqrt3/kSqrt2) * fkr.Tminus_2wTv() +
                            (kSqrt4/kSqrt3) * fkr.Lamda() * fkr.Rminus_wR();
     break;

  case kS11_1650 :
     ampl = -(0.5/kSqrt6) * fkr.Lamda() * fkr.Rminus_4wR();
     break;

  case kD13_1700 :
     ampl = -0.5*(1./kSqrt30) * fkr.Lamda() * fkr.Rminus_4wR();
     break;

  case kD15_1675 :
     ampl = 0.5*(kSqrt3/kSqrt10) * fkr.Lamda() * fkr.Rminus_4wR();
     break;

  case kS31_1620 :
     ampl = kSqrt3* fkr.Tminus_2wTv() -
                                  (1./kSqrt6)* fkr.Lamda()*fkr.Rminus_2wR();
     break;

  case kD33_1700 :
     ampl = (kSqrt3/kSqrt2)* fkr.Tminus_2wTv() +
                              (1./kSqrt3) * fkr.Lamda() * fkr.Rminus_2wR();
     break;

  case kP11_1440 :
     ampl = (5./12.)*kSqrt3*fkr.L2() *(fkr.Rminus()+ 2*xi*(4./5.)*fkr.R());
     break;

  case kP33_1600 :
     ampl = 0;
     break;

  case kP13_1720 :
     ampl = 0.5*(kSqrt27/kSqrt10) * fkr.LTminus() +
                            0.5*(kSqrt5/kSqrt3) * fkr.L2() *
                                       (fkr.Rminus() + (8./5.)*xi*fkr.R());
     break;

  case kF15_1680 :
     ampl = 0.5*(3./kSqrt5) * fkr.LTminus() -
                   0.5*(kSqrt5/kSqrt2) * fkr.L2() *
                                     (fkr.Rminus() + 2*xi*(4./5.)*fkr.R());
     break;

  case kP31_1910 :
     ampl = -(1./kSqrt15) * fkr.L2() * fkr.Rminus_2wR();
     break;

  case kP33_1920 :
     ampl = (1./kSqrt15) * fkr.L2() * fkr.Rminus_2wR();
     break;

  case kF35_1905 :
     ampl = (1./kSqrt35) * fkr.L2() * fkr.Rminus_2wR();
     break;

  case kF37_1950 :
     ampl = -(kSqrt6/kSqrt35) * fkr.L2() * fkr.Rminus_2wR();
     break;

  case kP11_1710 :
     ampl = -(1./kSqrt6) * fkr.L2() * fkr.Rminus_wR();
     break;

  case kF17_1970 :
     ampl = -2.5 * (fkr.Rminus() + 1.6 * xi * fkr.R());
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
double RSHelicityAmplModelNCn::AmpPlus1(
                       const Interaction * interaction, const FKR & fkr) const
{
  double ampl = 0;
  double xi   = kSin8w_2;

  Resonance_t res = utils::res::FromInteraction(interaction);

  switch(res) {

  case kP33_1232 :
     ampl  = kSqrt2 * fkr.Rplus_2wR();
     break;

  case kS11_1535 :
     ampl  = kSqrt3 * fkr.Tplus_2wTv() +
                            (kSqrt2/kSqrt3) * fkr.Lamda() * fkr.Rplus_wR();
     break;

  case kD13_1520 :
     ampl  = -(kSqrt3/kSqrt2) * fkr.Tplus_2wTv() +
                              (kSqrt4/kSqrt3) * fkr.Lamda()* fkr.Rplus_wR();
     break;

  case kS11_1650 :
     ampl  = (0.5/kSqrt6) * fkr.Lamda() * fkr.Rplus_4wR();
     break;

  case kD13_1700 :
     ampl  = -0.5*(1./kSqrt30) * fkr.Lamda() * fkr.Rplus_4wR();
     break;

  case kD15_1675 :
     ampl  = -0.5*(kSqrt3/kSqrt10) * fkr.Lamda() * fkr.Rplus_4wR();
     break;

  case kS31_1620 :
     ampl  = -kSqrt3 * fkr.Tplus_2wTv()  +
                                (1./kSqrt6) * fkr.Lamda() * fkr.Rplus_2wR();
     break;

  case kD33_1700 :
     ampl  = (kSqrt3/kSqrt2) * fkr.Tplus_2wTv() +
                               (1./kSqrt3) * fkr.Lamda() * fkr.Rplus_2wR();
     break;

  case kP11_1440 :
     ampl  = (5./12.) * kSqrt3 * fkr.L2() *
                                      (fkr.Rplus()  + 2*xi*(4./5.)*fkr.R());
     break;

  case kP33_1600 :
     ampl  = 0;
     break;

  case kP13_1720 :
     ampl  = -0.5*(kSqrt27/kSqrt10) * fkr.LTplus()
                              -0.5*( kSqrt5/kSqrt3 ) * fkr.L2() *
                                       (fkr.Rplus()  + (8./5.)*xi*fkr.R());
     break;

  case kF15_1680 :
     ampl  = 0.5*(3./kSqrt5) * fkr.LTplus()  -
                              0.5*( kSqrt5/kSqrt2) * fkr.L2() *
                                     (fkr.Rplus()  + 2*xi*(4./5.)*fkr.R());
     break;

  case kP31_1910 :
     ampl  = -(1./kSqrt15) * fkr.L2() * fkr.Rplus_2wR();
     break;

  case kP33_1920 :
     ampl  = -(1./kSqrt15) * fkr.L2() * fkr.Rplus_2wR();
     break;

  case kF35_1905 :
     ampl  = (1./kSqrt35) * fkr.L2() * fkr.Rplus_2wR();
     break;

  case kF37_1950 :
     ampl  = (kSqrt6/kSqrt35) * fkr.L2() * fkr.Rplus_2wR();
     break;

  case kP11_1710 :
     ampl  = -(1./kSqrt6) * fkr.L2() * fkr.Rplus_wR();
     break;

  case kF17_1970 :
     ampl  = -2.5 * (fkr.Rplus() + 2.4*xi*fkr.R());
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
double RSHelicityAmplModelNCn::AmpMinus3(
                       const Interaction * interaction, const FKR & fkr) const
{
  double ampl = 0;

  Resonance_t res = utils::res::FromInteraction(interaction);

  switch(res) {

  case kP33_1232 :
     ampl = -1 * kSqrt6 * fkr.Rminus_2wR();
     break;

  case kS11_1535 :
     ampl = 0;
     break;

  case kD13_1520 :
     ampl = -(3./kSqrt2) * fkr.Tminus_2wTv();
     break;

  case kS11_1650 :
     ampl = 0;
     break;

  case kD13_1700 :
     ampl = -0.5*(3./kSqrt10) * fkr.Lamda() * fkr.Rminus_4wR();
     break;

  case kD15_1675 :
     ampl = 0.5*(kSqrt3/kSqrt5) * fkr.Lamda() * fkr.Rminus_4wR();
     break;

  case kS31_1620 :
     ampl = 0;
     break;

  case kD33_1700 :
     ampl = (3./kSqrt2) * fkr.Tminus_2wTv();
     break;

  case kP11_1440 :
     ampl = 0;
     break;

  case kP33_1600 :
     ampl = 0;
     break;

  case kP13_1720 :
     ampl = -0.5*(3./kSqrt10) * fkr.LTminus();
     break;

  case kF15_1680 :
     ampl = 0.5*(kSqrt18/kSqrt5) * fkr.LTminus();
     break;

  case kP31_1910 :
     ampl = 0;
     break;

  case kP33_1920 :
     ampl = -(1./kSqrt5 ) * fkr.L2() * fkr.Rminus_2wR();
     break;

  case kF35_1905 :
     ampl = (kSqrt18/kSqrt35) * fkr.L2() * fkr.Rminus_2wR();
     break;

  case kF37_1950 :
     ampl = -(kSqrt2/kSqrt7) * fkr.L2() * fkr.Rminus_2wR();
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
double RSHelicityAmplModelNCn::AmpPlus3(
                       const Interaction * interaction, const FKR & fkr) const
{
  double ampl = 0;

  Resonance_t res = utils::res::FromInteraction(interaction);

  switch(res) {

  case kP33_1232 :
     ampl  = kSqrt6 * fkr.Rplus_2wR();
     break;

  case kS11_1535 :
     ampl  = 0;
     break;

  case kD13_1520 :
     ampl  = -(3./kSqrt2) * fkr.Tplus_2wTv();
     break;

  case kS11_1650 :
     ampl  = 0;
     break;

  case kD13_1700 :
     ampl  = -0.5*(3./kSqrt10) * fkr.Lamda() * fkr.Rplus_4wR();
     break;

  case kD15_1675 :
     ampl  = -0.5*(kSqrt3/kSqrt5) * fkr.Lamda() * fkr.Rplus_4wR();
     break;

  case kS31_1620 :
     ampl  = 0;
     break;

  case kD33_1700 :
     ampl  = (3./kSqrt2) * fkr.Tplus_2wTv();
     break;

  case kP11_1440 :
     ampl  = 0;
     break;

  case kP33_1600 :
     ampl  = 0;
     break;

  case kP13_1720 :
     ampl  = 0.5*(3./kSqrt10) * fkr.LTplus();
     break;

  case kF15_1680 :
     ampl  = 0.5*(kSqrt18/kSqrt5) * fkr.LTplus();
     break;

  case kP31_1910 :
     ampl  = 0;
     break;

  case kP33_1920 :
     ampl  = (1./kSqrt5 ) * fkr.L2() * fkr.Rplus_2wR();
     break;

  case kF35_1905 :
     ampl  = (kSqrt18/kSqrt35) * fkr.L2() * fkr.Rplus_2wR();
     break;

  case kF37_1950 :
     ampl  = (kSqrt2/kSqrt7) * fkr.L2() * fkr.Rplus_2wR();
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
double RSHelicityAmplModelNCn::Amp0Minus(
                       const Interaction * interaction, const FKR & fkr) const
{
  double ampl = 0;
  double xi   = kSin8w_2;

  Resonance_t res = utils::res::FromInteraction(interaction);

  switch(res) {

  case kP33_1232 :
     ampl = 2 * kSqrt2 * fkr.C();
     break;

  case kS11_1535 :
     ampl = (kSqrt3/kSqrt2) * fkr.LS() * (1-2*xi) -
                                            (kSqrt2/kSqrt3) * fkr.LC_3B();
     break;

  case kD13_1520 :
     ampl = (kSqrt3) * fkr.LS() * (1-2*xi) - (2./kSqrt3) * fkr.LC();
     break;

  case kS11_1650 :
     ampl = (1./kSqrt6) * fkr.LC_3B();
     break;

  case kD13_1700 :
     ampl = -0.5*(kSqrt2/kSqrt15) * fkr.LC();
     break;

  case kD15_1675 :
     ampl = -(kSqrt3/kSqrt10) * fkr.LC();
     break;

  case kS31_1620 :
     ampl = -(kSqrt3/kSqrt2)* fkr.LS()*(1-2*xi) - (1./kSqrt6)* fkr.LC_3B();
     break;

  case kD33_1700 :
     ampl = -kSqrt3 * fkr.LS() * (1-2*xi) - (1./kSqrt3) * fkr.LC();
     break;

  case kP11_1440 :
     ampl = 0.5*(kSqrt3/2.) * fkr.L2S() -
                                (5./12.)*kSqrt3 * fkr.Lamda() * fkr.LC_2B();
     break;

  case kP33_1600 :
     ampl = 0;
     break;

  case kP13_1720 :
     ampl = -0.5*(kSqrt3/kSqrt5) * fkr.L2S()  +
                                  0.5*(kSqrt5/kSqrt3) * fkr.Lamda()* fkr.LC_5B();
     break;

  case kF15_1680 :
     ampl = -0.5*(3./kSqrt10)* fkr.L2S() + 0.5*(kSqrt5/kSqrt2) * fkr.L2C();
     break;

  case kP31_1910 :
     ampl = -(2./kSqrt15) * fkr.Lamda() * fkr.LC_5B();
     break;

  case kP33_1920 :
     ampl = -(2./kSqrt15) * fkr.Lamda() * fkr.LC_5B();
     break;

  case kF35_1905 :
     ampl = (2./kSqrt35) * fkr.L2C();
     break;

  case kF37_1950 :
     ampl = 2*(kSqrt6/kSqrt35) * fkr.L2C();
     break;

  case kP11_1710 :
     ampl =  -(kSqrt3/kSqrt8) * fkr.L2S() * (1-2*xi) +
                                    (1./kSqrt6) * fkr.Lamda() * fkr.LC_2B();
     break;

  case kF17_1970 :
     ampl = -1.5 * fkr.S() + 2.5 * fkr.C();
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
double RSHelicityAmplModelNCn::Amp0Plus(
                       const Interaction * interaction, const FKR & fkr) const
{
  double ampl = 0;
  double xi   = kSin8w_2;

  Resonance_t res = utils::res::FromInteraction(interaction);

  switch(res) {

  case kP33_1232 :
     ampl  = 2 * kSqrt2 * fkr.C();
     break;

  case kS11_1535 :
     ampl  = -1*(kSqrt3/kSqrt2) * fkr.LS() * (1-2*xi) -
                                               (kSqrt2/kSqrt3) * fkr.LC_3B();
     break;

  case kD13_1520 :
     ampl  = (kSqrt3) * fkr.LS() * (1-2*xi) + (2./kSqrt3) * fkr.LC();
     break;

  case kS11_1650 :
     ampl  = (1./kSqrt6) * fkr.LC_3B();
     break;

  case kD13_1700 :
     ampl  = 0.5*(kSqrt2/kSqrt15)* fkr.LC();
     break;

  case kD15_1675 :
     ampl  = -(kSqrt3/kSqrt10) * fkr.LC();
     break;

  case kS31_1620 :
     ampl  = (kSqrt3/kSqrt2) * fkr.LS()* (1-2*xi) - (1./kSqrt6) * fkr.LC_2B();
     break;

  case kD33_1700 :
     ampl  = -kSqrt3 * fkr.LS() * (1-2*xi) + (1./kSqrt3) * fkr.LC();
     break;

  case kP11_1440 :
     ampl  = 0.5*(kSqrt3/2.) * fkr.L2S() +
                                   (5./12.)*kSqrt3 * fkr.Lamda() * fkr.LC_2B();
     break;

  case kP33_1600 :
     ampl  = 0;
     break;

  case kP13_1720 :
     ampl  = 0.5*(kSqrt3/kSqrt5) * fkr.L2S()  +
                                    0.5*(kSqrt5/kSqrt3) * fkr.Lamda() * fkr.LC_5B();
     break;

  case kF15_1680 :
     ampl  = -0.5*(3./kSqrt10) * fkr.L2S() - 0.5*(kSqrt5/kSqrt2) * fkr.L2C();
     break;

  case kP31_1910 :
     ampl  =  (2./kSqrt15) * fkr.Lamda() * fkr.LC_5B();
     break;

  case kP33_1920 :
     ampl  = -(2./kSqrt15) * fkr.Lamda() * fkr.LC_5B();
     break;

  case kF35_1905 :
     ampl  = -(2./kSqrt35) * fkr.L2C();
     break;

  case kF37_1950 :
     ampl  = 2*(kSqrt6/kSqrt35) * fkr.L2C();
     break;

  case kP11_1710 :
     ampl  = -(kSqrt3/kSqrt8) * fkr.L2S() * (1-2*xi) -
                                     (1./kSqrt6) * fkr.Lamda() * fkr.LC_2B();
     break;

  case kF17_1970 :
     ampl  = -1.5*fkr.S() - 2.5 * fkr.C();
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
