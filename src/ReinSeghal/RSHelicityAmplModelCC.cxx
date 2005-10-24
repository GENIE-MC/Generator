//____________________________________________________________________________
/*!

\class    genie::RSHelicityAmplModelCC

\brief    The Helicity Amplitudes, for all baryon resonances, for CC neutrino
          interactions on free nucleons, as computed in the Rein-Seghal's paper.

          Concrete implementation of the SPPHelicityAmplModelI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004

*/
//____________________________________________________________________________

#include "BaryonResonance/BaryonResUtils.h"
#include "Conventions/Constants.h"
#include "ReinSeghal/RSHelicityAmplModelCC.h"
#include "Messenger/Messenger.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
RSHelicityAmplModelCC::RSHelicityAmplModelCC() :
RSHelicityAmplModelI("genie::RSHelicityAmplModelCC")
{

}
//____________________________________________________________________________
RSHelicityAmplModelCC::RSHelicityAmplModelCC(string config) :
RSHelicityAmplModelI("genie::RSHelicityAmplModelCC", config)
{

}
//____________________________________________________________________________
RSHelicityAmplModelCC::~RSHelicityAmplModelCC()
{

}
//____________________________________________________________________________
double RSHelicityAmplModelCC::AmpMinus1(
                       const Interaction * interaction, const FKR & fkr) const
{
   double ampl = 0;

   Resonance_t res = utils::res::FromInteraction(interaction);

   switch(res) {

   case (kP33_1232) :

     ampl = kSqrt2 * fkr.Rminus();
     break;

   case (kS11_1535) :
     ampl = (2.*kSqrt3) * fkr.Tminus() + (4./kSqrt6) * fkr.LRminus();
     break;

   case (kD13_1520) :
     ampl =  kSqrt6 * fkr.Tminus() - (4./kSqrt3) * fkr.LRminus();
     break;

   case (kS11_1650) :
     ampl =  (1./kSqrt6) * fkr.LRminus();
     break;

   case (kD13_1700) :
     ampl =  (1./kSqrt30) * fkr.LRminus();
     break;

   case (kD15_1675) :
     ampl =  -(kSqrt3/kSqrt10) * fkr.LRminus();
     break;

   case (kS31_1620) :
     ampl =  -kSqrt3 * fkr.Tminus() + (1./kSqrt6) * fkr.LRminus();
     break;

   case (kD33_1700) :
     ampl =  -(kSqrt3/kSqrt2) * fkr.Tminus() - (1./kSqrt3) * fkr.LRminus();
     break;

   case (kP11_1440) :
     ampl = -(5.*kSqrt3/6.) * fkr.L2Rminus();
     break;

   case (kP33_1600) :
     ampl =  -(1./kSqrt6) * fkr.L2Rminus();
     break;

   case (kP13_1720) :
     ampl =  -(3.*kSqrt3/kSqrt10) * fkr.Lamda() * fkr.Tminus() - (kSqrt5/kSqrt3) * fkr.L2Rminus();
     break;

   case (kF15_1680) :
     ampl =  -(3./kSqrt5) * fkr.Lamda() * fkr.Tminus() + (kSqrt5/kSqrt2) * fkr.L2Rminus();
     break;

   case (kP31_1910) :
     ampl =  (1./kSqrt15) * fkr.L2Rminus();
     break;

   case (kP33_1920) :
     ampl =  -(1./kSqrt15) * fkr.L2Rminus();
     break;

   case (kF35_1905) :
     ampl =  -(1./kSqrt35) * fkr.L2Rminus();
     break;

   case (kF37_1950) :
     ampl = (kSqrt6/kSqrt35)  * fkr.L2Rminus();
     break;

   case (kP11_1710) :
     ampl = (kSqrt2/kSqrt3) * fkr.L2Rminus();
     break;

   case (kF17_1970) :
     ampl =  -(kSqrt3/kSqrt35) * fkr.L2Rminus();
     break;

   default:
     LOG("RSHAmpl", pWARN) << "*** UNRECOGNIZED RESONANCE!";
     ampl = 0;
     break;
   }

   LOG("RSHAmpl", pDEBUG)
            << "A(-1) for RES: " << utils::res::AsString(res) << " = " << ampl;

   return ampl;
}
//____________________________________________________________________________
double RSHelicityAmplModelCC::AmpPlus1(
                       const Interaction * interaction, const FKR & fkr) const
{
   double ampl = 0;

   Resonance_t res = utils::res::FromInteraction(interaction);

   switch(res) {

   case (kP33_1232) :
     ampl = -kSqrt2 * fkr.Rplus();
     break;

   case (kS11_1535) :
     ampl = (-2.*kSqrt3) * fkr.Tplus() - (4./kSqrt6) * fkr.LRplus();
     break;

   case (kD13_1520) :
     ampl = (kSqrt6) * fkr.Tplus() - (4./kSqrt3) * fkr.Rplus();
     break;

   case (kS11_1650) :
     ampl = -(1./kSqrt6) * fkr.LRplus();
     break;

   case (kD13_1700) :
     ampl = (1./kSqrt30) * fkr.LRplus();
     break;

   case (kD15_1675) :
     ampl = (kSqrt3/kSqrt10) * fkr.LRplus();
     break;

   case (kS31_1620) :
     ampl = (kSqrt3) * fkr.Tplus() + (1./kSqrt6) * fkr.LRplus();
     break;

   case (kD33_1700) :
     ampl = -(kSqrt3/kSqrt2) * fkr.Tplus() - (1/kSqrt3) * fkr.LRplus();
     break;

   case (kP11_1440) :
     ampl = -(5.*kSqrt3/6.) * fkr.L2Rplus();
     break;

   case (kP33_1600) :
     ampl = (1./kSqrt6) * fkr.L2Rplus();
     break;

   case (kP13_1720) :
     ampl = (3.*kSqrt3/kSqrt10 ) * fkr.Lamda() *  fkr.Tminus() + (kSqrt5/kSqrt3) * fkr.L2Rminus();
     break;

   case (kF15_1680) :
     ampl = -(3./kSqrt5) * fkr.Lamda() * fkr.Tplus() + (kSqrt5/kSqrt2) * fkr.L2Rplus();
     break;

   case (kP31_1910) :
     ampl = (1./kSqrt15) * fkr.L2Rplus();
     break;

   case (kP33_1920) :
     ampl = (1./kSqrt15) * fkr.L2Rplus();
     break;

   case (kF35_1905) :
     ampl = -(1./kSqrt35) * fkr.L2Rplus();
     break;

   case (kF37_1950) :
     ampl = -(kSqrt6/kSqrt35) * fkr.L2Rplus();
     break;

   case (kP11_1710) :
     ampl =  (kSqrt2/kSqrt3) * fkr.L2Rplus();
     break;

   case (kF17_1970) :
     ampl =  (kSqrt3/kSqrt35) * fkr.L2Rplus();
     break;

   default:
     LOG("RSHAmpl", pWARN) << "*** UNRECOGNIZED RESONANCE!";
     ampl = 0;
     break;
   }
   LOG("RSHAmpl", pDEBUG)
            << "A(+1) for RES: " << utils::res::AsString(res) << " = " << ampl;

   return ampl;
}
//____________________________________________________________________________
double RSHelicityAmplModelCC::AmpMinus3(
                       const Interaction * interaction, const FKR & fkr) const
{
   double ampl = 0;

   Resonance_t res = utils::res::FromInteraction(interaction);

   switch(res) {

   case (kP33_1232) :
     ampl = kSqrt6 * fkr.Rminus();
     break;

   case (kS11_1535) :
     ampl = 0;
     break;

   case (kD13_1520) :
     ampl =  (2.*kSqrt9/kSqrt2) * fkr.Tminus();
     break;

   case (kS11_1650) :
     ampl = 0;
     break;

   case (kD13_1700) :
     ampl = (3./kSqrt10) * fkr.LRminus();
     break;

   case (kD15_1675) :
     ampl = -(kSqrt3/kSqrt5) * fkr.LRminus();
     break;

   case (kS31_1620) :
     ampl = 0;
     break;

   case (kD33_1700) :
     ampl = -(3./kSqrt2) * fkr.Tminus();
     break;

   case (kP11_1440) :
     ampl = 0;
     break;

   case (kP33_1600) :
     ampl = -(1./kSqrt2) * fkr.Lamda() * fkr.Lamda() * fkr.Rminus();
     break;

   case (kP13_1720) :
     ampl = (3./kSqrt10) * fkr.Lamda() * fkr.Tminus();
     break;

   case (kF15_1680) :
     ampl = -(3.*kSqrt2/kSqrt5) * fkr.Lamda() * fkr.Tminus();
     break;

   case (kP31_1910) :
     ampl =   0;
     break;

   case (kP33_1920) :
     ampl = (1./kSqrt5) * fkr.L2Rminus();
     break;

   case (kF35_1905) :
     ampl = -(3.*kSqrt2/kSqrt35) * fkr.L2Rminus();
     break;

   case (kF37_1950) :
     ampl = (kSqrt2/kSqrt7) * fkr.L2Rminus();
     break;

   case (kP11_1710) :
     ampl = 0;
     break;

   case (kF17_1970) :
     ampl = -(1./kSqrt7) * fkr.L2Rminus();
     break;

   default:
     LOG("RSHAmpl", pWARN) << "*** UNRECOGNIZED RESONANCE!";
     ampl = 0;
     break;
   }

   LOG("RSHAmpl", pDEBUG)
            << "A(-3) for RES: " << utils::res::AsString(res) << " = " << ampl;

  return ampl;
}
//____________________________________________________________________________
double RSHelicityAmplModelCC::AmpPlus3(
                       const Interaction * interaction, const FKR & fkr) const
{
   double ampl = 0;

   Resonance_t res = utils::res::FromInteraction(interaction);

   switch(res) {

   case (kP33_1232) :
     ampl = -kSqrt6 * fkr.Rplus();
     break;

   case (kS11_1535) :
     ampl = 0;
     break;

   case (kD13_1520) :
     ampl = (2.*kSqrt9/kSqrt2) * fkr.Tplus();
     break;

   case (kS11_1650) :
     ampl = 0;
     break;

   case (kD13_1700) :
     ampl = (3./kSqrt10) * fkr.Lamda() * fkr.Rplus();
     break;

   case (kD15_1675) :
     ampl = (kSqrt3/kSqrt5) * fkr.Lamda() * fkr.Rplus();
     break;

   case (kS31_1620) :
     ampl = 0;
     break;

   case (kD33_1700) :
     ampl = -(3./kSqrt2) * fkr.Tplus();
     break;

   case (kP11_1440) :
     ampl = 0;
     break;

   case (kP33_1600) :
     ampl = (1./kSqrt2) * fkr.L2Rplus();
     break;

   case (kP13_1720) :
     ampl = -(3./kSqrt10) * fkr.Lamda() * fkr.Tplus();
     break;

   case (kF15_1680) :
     ampl = -(3.*kSqrt2/kSqrt5) * fkr.Lamda() * fkr.Tplus();
     break;

   case (kP31_1910) :
     ampl = 0;
     break;

   case (kP33_1920) :
     ampl =  -(1./kSqrt5) * fkr.L2Rplus();
     break;

   case (kF35_1905) :
     ampl = -(3.*kSqrt2/kSqrt35) * fkr.L2Rplus();
     break;

   case (kF37_1950) :
     ampl = -(kSqrt2/kSqrt7) * fkr.L2Rplus();
     break;

   case (kP11_1710) :
     ampl = 0;
     break;

   case (kF17_1970) :
     ampl = (1./kSqrt7) * fkr.L2Rplus();
     break;

   default:
     LOG("RSHAmpl", pWARN) << "*** UNRECOGNIZED RESONANCE!";
     ampl = 0;
     break;
   }

   LOG("RSHAmpl", pDEBUG)
            << "A(+3) for RES: " << utils::res::AsString(res) << " = " << ampl;

   return ampl;
}
//____________________________________________________________________________
double RSHelicityAmplModelCC::Amp0Minus(
                       const Interaction * interaction, const FKR & fkr) const
{
   double ampl = 0;

   Resonance_t res = utils::res::FromInteraction(interaction);

   switch(res) {

   case (kP33_1232) :
     ampl = - 2*kSqrt2 * fkr.C();
     break;

   case (kS11_1535) :
     ampl = (-kSqrt6) * fkr.LS() + (2.*kSqrt2/kSqrt6) * fkr.LC_3B();
     break;

   case (kD13_1520) :
     ampl = (-2.*kSqrt3) * fkr.LS() + (4./kSqrt3) * fkr.LC();
     break;

   case (kS11_1650) :
     ampl = -(kSqrt2/kSqrt3) * fkr.LC_3B();
     break;

   case (kD13_1700) :
     ampl = (kSqrt2/kSqrt15) * fkr.LC();
     break;

   case (kD15_1675) :
     ampl = (kSqrt6/kSqrt5) * fkr.LC();
     break;

   case (kS31_1620) :
     ampl = (kSqrt3/kSqrt2) * fkr.LS() + (1./kSqrt6) * fkr.LC_3B();
     break;

   case (kD33_1700) :
     ampl = (kSqrt3) * fkr.LS() + (1./kSqrt3) * fkr.LC();
     break;

   case (kP11_1440) :
     ampl = -(kSqrt3/2.) * fkr.L2S() + (5.*kSqrt3/6.) * fkr.Lamda() * fkr.LC_2B();
     break;

   case (kP33_1600) :
     ampl = (kSqrt2/kSqrt3) * fkr.Lamda() * fkr.LC_2B();
     break;

   case (kP13_1720) :
     ampl = (kSqrt3/kSqrt5) * fkr.L2S() - (kSqrt5/kSqrt3) * fkr.Lamda() * fkr.LC_5B();
     break;

   case (kF15_1680) :
     ampl = (3./kSqrt10) * fkr.L2S() - (kSqrt5/kSqrt2) * fkr.L2C();
     break;

   case (kP31_1910) :
     ampl = (2./kSqrt15) * fkr.Lamda() * fkr.LC_5B();
     break;

   case (kP33_1920) :
     ampl = (2./kSqrt15) * fkr.Lamda() * fkr.LC_5B();
     break;

   case (kF35_1905) :
     ampl = -(2./kSqrt35) * fkr.L2C();
     break;

   case (kF37_1950) :
     ampl = -(2.*kSqrt6/kSqrt35) * fkr.L2C();
     break;

   case (kP11_1710) :
     ampl = (kSqrt3/kSqrt2) * fkr.L2S() -
                               (kSqrt2/kSqrt3) * fkr.Lamda() * fkr.LC_2B();
     break;

   case (kF17_1970) :
     ampl = (kSqrt6/kSqrt35) * fkr.L2C();
     break;

   default:
     LOG("RSHAmpl", pWARN) << "*** UNRECOGNIZED RESONANCE!";
     ampl =   0;
     break;
   }
   LOG("RSHAmpl", pDEBUG)
            << "A(0-) for RES: " << utils::res::AsString(res) << " = " << ampl;

   return ampl;
}
//____________________________________________________________________________
double RSHelicityAmplModelCC::Amp0Plus(
                       const Interaction * interaction, const FKR & fkr) const
{
   double ampl = 0;

   Resonance_t res = utils::res::FromInteraction(interaction);

   switch(res) {

   case (kP33_1232) :
     ampl = -(2*kSqrt2) * fkr.C();
     break;

   case (kS11_1535) :
     ampl = (kSqrt6) * fkr.LS() + (2.*kSqrt2/kSqrt6) * fkr.LC_3B();
     break;

   case (kD13_1520) :
     ampl = (-2.*kSqrt3) * fkr.LS() - (4./kSqrt3) * fkr.LC();
     break;

   case (kS11_1650) :
     ampl = -(kSqrt2/kSqrt3) * fkr.LC_3B();
     break;

   case (kD13_1700) :
     ampl = -(kSqrt2/kSqrt15) * fkr.LC();
     break;

   case (kD15_1675) :
     ampl = (kSqrt6/kSqrt5) * fkr.LC();
     break;

   case (kS31_1620) :
     ampl = -(kSqrt3/kSqrt2) * fkr.LS() + (1./kSqrt6) * fkr.LC_3B();
     break;

   case (kD33_1700) :
     ampl = (kSqrt3) * fkr.LS() - (1./kSqrt3) * fkr.LC();
     break;

   case (kP11_1440) :
     ampl = -(kSqrt3/2.) * fkr.L2S() - (5.*kSqrt3/6.) * fkr.Lamda() * fkr.LC_2B();
     break;

   case (kP33_1600) :
     ampl = (kSqrt2/kSqrt3) * fkr.Lamda() * fkr.LC_2B();
     break;

   case (kP13_1720) :
     ampl = -(kSqrt3/kSqrt5) * fkr.L2S() - (kSqrt5/kSqrt3) * fkr.Lamda() * fkr.LC_5B();
     break;

   case (kF15_1680) :
     ampl = (3./kSqrt10) * fkr.L2S() + (kSqrt5/kSqrt2) * fkr.L2C();
     break;

   case (kP31_1910) :
     ampl = -(2./kSqrt15) * fkr.Lamda() * fkr.LC_5B();
     break;

   case (kP33_1920) :
     ampl = (2./kSqrt15) * fkr.Lamda() * fkr.LC_5B();
     break;

   case (kF35_1905) :
     ampl = (2./kSqrt35) * fkr.L2C();
     break;

   case (kF37_1950) :
     ampl = -(2.*kSqrt6/kSqrt35) * fkr.L2C();
     break;

   case (kP11_1710) :
     ampl = (kSqrt3/kSqrt2) * fkr.L2S() + (kSqrt2/kSqrt3) * fkr.Lamda() * fkr.LC_2B();
     break;

   case (kF17_1970) :
     ampl = (kSqrt6/kSqrt35) * fkr.L2C();
     break;

   default:
     LOG("RSHAmpl", pWARN) << "*** UNRECOGNIZED RESONANCE!";
     ampl = 0;
     break;
   }

   LOG("RSHAmpl", pDEBUG)
            << "A(0+) for RES: " << utils::res::AsString(res) << " = " << ampl;

   return ampl;
}
//____________________________________________________________________________

