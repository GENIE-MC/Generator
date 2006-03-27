//____________________________________________________________________________
/*!

\class    genie::RSHelicityAmplModelCC

\brief    The Helicity Amplitudes, for all baryon resonances, for CC neutrino
          interactions on free nucleons, as computed in the Rein-Seghal's paper.

          Concrete implementation of the RSHelicityAmplModelI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004

*/
//____________________________________________________________________________

#include "BaryonResonance/BaryonResUtils.h"
#include "Conventions/Constants.h"
#include "ReinSeghal/RSHelicityAmplModelCC.h"
#include "ReinSeghal/RSHelicityAmpl.h"
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
RSHelicityAmpl * RSHelicityAmplModelCC::Compute(
                                       Resonance_t res, const FKR & fkr) const
{
  RSHelicityAmpl * hampl = new RSHelicityAmpl;

  switch(res) {

   case (kP33_1232) :
     hampl->fMinus1 =    kSqrt2 * fkr.Rminus();
     hampl->fPlus1  =   -kSqrt2 * fkr.Rplus();
     hampl->fMinus3 =    kSqrt6 * fkr.Rminus();
     hampl->fPlus3  =   -kSqrt6 * fkr.Rplus();
     hampl->f0Minus = -2*kSqrt2 * fkr.C();
     hampl->f0Plus  =    hampl->f0Minus;
     break;

   case (kS11_1535) :
     hampl->fMinus1 =  2.*kSqrt3 * fkr.Tminus() + (4./kSqrt6) * fkr.LRminus();
     hampl->fPlus1  = -2.*kSqrt3 * fkr.Tplus()  - (4./kSqrt6) * fkr.LRplus();
     hampl->fMinus3 =  0;
     hampl->fPlus3  =  0;
     hampl->f0Minus = -kSqrt6 * fkr.LS() + (2.*kSqrt2/kSqrt6) * fkr.LC_3B();
     hampl->f0Plus  =  kSqrt6 * fkr.LS() + (2.*kSqrt2/kSqrt6) * fkr.LC_3B();
     break;

   case (kD13_1520) :
     hampl->fMinus1 =  kSqrt6 * fkr.Tminus() - (4./kSqrt3) * fkr.LRminus();
     hampl->fPlus1  =  kSqrt6 * fkr.Tplus()  - (4./kSqrt3) * fkr.Rplus();
     hampl->fMinus3 = (2.*kSqrt9/kSqrt2) * fkr.Tminus();
     hampl->fPlus3  = (2.*kSqrt9/kSqrt2) * fkr.Tplus();
     hampl->f0Minus = -2.*kSqrt3 * fkr.LS() + (4./kSqrt3) * fkr.LC();
     hampl->f0Plus  = -2.*kSqrt3 * fkr.LS() - (4./kSqrt3) * fkr.LC();
     break;

   case (kS11_1650) :
     hampl->fMinus1 =  (1./kSqrt6) * fkr.LRminus();
     hampl->fPlus1  = -(1./kSqrt6) * fkr.LRplus();
     hampl->fMinus3 =  0;
     hampl->fPlus3  =  0;
     hampl->f0Minus = -(kSqrt2/kSqrt3) * fkr.LC_3B();
     hampl->f0Plus  =  hampl->f0Minus;
     break;

   case (kD13_1700) :
     hampl->fMinus1 =  (1./kSqrt30) * fkr.LRminus();
     hampl->fPlus1  =  (1./kSqrt30) * fkr.LRplus();
     hampl->fMinus3 =  (3./kSqrt10) * fkr.LRminus();
     hampl->fPlus3  =  (3./kSqrt10) * fkr.Lamda() * fkr.Rplus();
     hampl->f0Minus =  (kSqrt2/kSqrt15) * fkr.LC();
     hampl->f0Plus  = -1. * hampl->f0Minus;
     break;

   case (kD15_1675) :
     hampl->fMinus1 = -(kSqrt3/kSqrt10) * fkr.LRminus();
     hampl->fPlus1  =  (kSqrt3/kSqrt10) * fkr.LRplus();
     hampl->fMinus3 = -(kSqrt3/kSqrt5)  * fkr.LRminus();
     hampl->fPlus3  =  (kSqrt3/kSqrt5)  * fkr.Lamda() * fkr.Rplus();
     hampl->f0Minus =  (kSqrt6/kSqrt5)  * fkr.LC();
     hampl->f0Plus  =  hampl->f0Minus;
     break;

   case (kS31_1620) :
     hampl->fMinus1 =  -kSqrt3 * fkr.Tminus() + (1./kSqrt6) * fkr.LRminus();
     hampl->fPlus1  =   kSqrt3 * fkr.Tplus()  + (1./kSqrt6) * fkr.LRplus();
     hampl->fMinus3 =  0;
     hampl->fPlus3  =  0;
     hampl->f0Minus =  (kSqrt3/kSqrt2) * fkr.LS() + (1./kSqrt6) * fkr.LC_3B();
     hampl->f0Plus  = -(kSqrt3/kSqrt2) * fkr.LS() + (1./kSqrt6) * fkr.LC_3B();
     break;

   case (kD33_1700) :
     hampl->fMinus1 = -(kSqrt3/kSqrt2) * fkr.Tminus() -
                       (1./kSqrt3) * fkr.LRminus();
     hampl->fPlus1  = -(kSqrt3/kSqrt2) * fkr.Tplus()  -
                       (1./kSqrt3)  * fkr.LRplus();
     hampl->fMinus3 = -(3./kSqrt2) * fkr.Tminus();
     hampl->fPlus3  = -(3./kSqrt2) * fkr.Tplus();
     hampl->f0Minus = (kSqrt3) * fkr.LS() + (1./kSqrt3) * fkr.LC();
     hampl->f0Plus  = (kSqrt3) * fkr.LS() - (1./kSqrt3) * fkr.LC();
     break;

   case (kP11_1440) :
     hampl->fMinus1 = -(5.*kSqrt3/6.) * fkr.L2Rminus();
     hampl->fPlus1  = -(5.*kSqrt3/6.) * fkr.L2Rplus();
     hampl->fMinus3 =  0;
     hampl->fPlus3  =  0;
     hampl->f0Minus = -(kSqrt3/2.) * fkr.L2S() + (5.*kSqrt3/6.) * fkr.Lamda() * fkr.LC_2B();
     hampl->f0Plus  = -(kSqrt3/2.) * fkr.L2S() - (5.*kSqrt3/6.) * fkr.Lamda() * fkr.LC_2B();
     break;

   case (kP33_1600) :
     hampl->fMinus1 = -(1./kSqrt6) * fkr.L2Rminus();
     hampl->fPlus1  =  (1./kSqrt6) * fkr.L2Rplus();
     hampl->fMinus3 = -(1./kSqrt2) * fkr.Lamda() * fkr.Lamda() * fkr.Rminus();
     hampl->fPlus3  =  (1./kSqrt2) * fkr.L2Rplus();
     hampl->f0Minus =  (kSqrt2/kSqrt3) * fkr.Lamda() * fkr.LC_2B();
     hampl->f0Plus  =  hampl->f0Minus;
     break;

   case (kP13_1720) :
     hampl->fMinus1 =  -(3.*kSqrt3/kSqrt10)  * fkr.Lamda() * fkr.Tminus() -
                        (kSqrt5/kSqrt3) * fkr.L2Rminus();
     hampl->fPlus1  =   (3.*kSqrt3/kSqrt10 ) * fkr.Lamda() * fkr.Tminus() +
                        (kSqrt5/kSqrt3) * fkr.L2Rminus();
     hampl->fMinus3 =   (3./kSqrt10) * fkr.Lamda() * fkr.Tminus();
     hampl->fPlus3  =  -(3./kSqrt10) * fkr.Lamda() * fkr.Tplus();
     hampl->f0Minus =   (kSqrt3/kSqrt5) * fkr.L2S() -
                        (kSqrt5/kSqrt3) * fkr.Lamda() * fkr.LC_5B();
     hampl->f0Plus  =  -(kSqrt3/kSqrt5) * fkr.L2S() -
                        (kSqrt5/kSqrt3) * fkr.Lamda() * fkr.LC_5B();
     break;

   case (kF15_1680) :
     hampl->fMinus1 = -(3./kSqrt5) * fkr.Lamda() * fkr.Tminus() +
                       (kSqrt5/kSqrt2) * fkr.L2Rminus();
     hampl->fPlus1  = -(3./kSqrt5) * fkr.Lamda() * fkr.Tplus()  +
                       (kSqrt5/kSqrt2) * fkr.L2Rplus();
     hampl->fMinus3 = -(3.*kSqrt2/kSqrt5) * fkr.Lamda() * fkr.Tminus();
     hampl->fPlus3  = -(3.*kSqrt2/kSqrt5) * fkr.Lamda() * fkr.Tplus();
     hampl->f0Minus =  (3./kSqrt10) * fkr.L2S() - (kSqrt5/kSqrt2) * fkr.L2C();
     hampl->f0Plus  =  (3./kSqrt10) * fkr.L2S() + (kSqrt5/kSqrt2) * fkr.L2C();
     break;

   case (kP31_1910) :
     hampl->fMinus1 =  (1./kSqrt15) * fkr.L2Rminus();
     hampl->fPlus1  =  (1./kSqrt15) * fkr.L2Rplus();
     hampl->fMinus3 =  0;
     hampl->fPlus3  =  0;
     hampl->f0Minus =  (2./kSqrt15) * fkr.Lamda() * fkr.LC_5B();
     hampl->f0Plus  = -1. * hampl->f0Minus;
     break;

   case (kP33_1920) :
     hampl->fMinus1 = -(1./kSqrt15) * fkr.L2Rminus();
     hampl->fPlus1  =  (1./kSqrt15) * fkr.L2Rplus();
     hampl->fMinus3 =  (1./kSqrt5)  * fkr.L2Rminus();
     hampl->fPlus3  = -(1./kSqrt5)  * fkr.L2Rplus();
     hampl->f0Minus =  (2./kSqrt15) * fkr.Lamda() * fkr.LC_5B();
     hampl->f0Plus  =  hampl->f0Minus;
     break;

   case (kF35_1905) :
     hampl->fMinus1 =  -(1./kSqrt35) * fkr.L2Rminus();
     hampl->fPlus1  =  -(1./kSqrt35) * fkr.L2Rplus();
     hampl->fMinus3 =  -(3.*kSqrt2/kSqrt35) * fkr.L2Rminus();
     hampl->fPlus3  =  -(3.*kSqrt2/kSqrt35) * fkr.L2Rplus();
     hampl->f0Minus =  -(2./kSqrt35) * fkr.L2C();
     hampl->f0Plus  =  -1. * hampl->f0Minus;
     break;

   case (kF37_1950) :
     hampl->fMinus1 =  (kSqrt6/kSqrt35)    * fkr.L2Rminus();
     hampl->fPlus1  = -(kSqrt6/kSqrt35)    * fkr.L2Rplus();
     hampl->fMinus3 =  (kSqrt2/kSqrt7)     * fkr.L2Rminus();
     hampl->fPlus3  = -(kSqrt2/kSqrt7)     * fkr.L2Rplus();
     hampl->f0Minus = -(2.*kSqrt6/kSqrt35) * fkr.L2C();
     hampl->f0Plus  = hampl->f0Minus;
     break;

   case (kP11_1710) :
     hampl->fMinus1 = (kSqrt2/kSqrt3) * fkr.L2Rminus();
     hampl->fPlus1  = (kSqrt2/kSqrt3) * fkr.L2Rplus();
     hampl->fMinus3 = 0;
     hampl->fPlus3  = 0;
     hampl->f0Minus = (kSqrt3/kSqrt2) * fkr.L2S() -
                               (kSqrt2/kSqrt3) * fkr.Lamda() * fkr.LC_2B();
     hampl->f0Plus  = (kSqrt3/kSqrt2) * fkr.L2S() +
                               (kSqrt2/kSqrt3) * fkr.Lamda() * fkr.LC_2B();
     break;

   case (kF17_1970) :
     hampl->fMinus1 =  -(kSqrt3/kSqrt35) * fkr.L2Rminus();
     hampl->fPlus1  =   (kSqrt3/kSqrt35) * fkr.L2Rplus();
     hampl->fMinus3 =  -(1./kSqrt7)      * fkr.L2Rminus();
     hampl->fPlus3  =   (1./kSqrt7)      * fkr.L2Rplus();
     hampl->f0Minus =   (kSqrt6/kSqrt35) * fkr.L2C();
     hampl->f0Plus  =   hampl->f0Minus;
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


