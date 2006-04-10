//____________________________________________________________________________
/*!

\class    genie::RSHelicityAmplModelNCn

\brief    The Helicity Amplitudes, for all baryon resonances, for NC neutrino
          interactions on free neutrons, as computed in the Rein-Seghal's paper.

          Concrete implementation of the RSHelicityAmplModelI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004

*/
//____________________________________________________________________________

#include "BaryonResonance/BaryonResUtils.h"
#include "Conventions/Constants.h"
#include "ReinSeghal/RSHelicityAmplModelNCn.h"
#include "ReinSeghal/RSHelicityAmpl.h"
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
RSHelicityAmpl * RSHelicityAmplModelNCn::Compute(
                                       Resonance_t res, const FKR & fkr) const
{
  RSHelicityAmpl * hampl = new RSHelicityAmpl;
  double xi = kSin8w2;

  switch(res) {

   case (kP33_1232) :
     hampl->fMinus1 =  -kSqrt2 * fkr.Rminus_2wR();
     hampl->fPlus1  =   kSqrt2 * fkr.Rplus_2wR();
     hampl->fMinus3 =  -kSqrt6 * fkr.Rminus_2wR();
     hampl->fPlus3  =   kSqrt6 * fkr.Rplus_2wR();
     hampl->f0Minus = 2*kSqrt2 * fkr.C();
     hampl->f0Plus  =   hampl->f0Minus;
     break;

   case (kS11_1535) :
     hampl->fMinus1 = -1*kSqrt3* fkr.Tminus_2wTv() -
                             (kSqrt2/kSqrt3) * fkr.Lamda() * fkr.Rminus_wR();
     hampl->fPlus1  = kSqrt3 * fkr.Tplus_2wTv() +
                              (kSqrt2/kSqrt3) * fkr.Lamda() * fkr.Rplus_wR();
     hampl->fMinus3 = 0.;
     hampl->fPlus3  = 0.;
     hampl->f0Minus =  (kSqrt3/kSqrt2) * fkr.LS() * (1-2*xi) -
                                               (kSqrt2/kSqrt3) * fkr.LC_3B();
     hampl->f0Plus  = -(kSqrt3/kSqrt2) * fkr.LS() * (1-2*xi) -
                                               (kSqrt2/kSqrt3) * fkr.LC_3B();
     break;

   case (kD13_1520) :
     hampl->fMinus1 = -(kSqrt3/kSqrt2) * fkr.Tminus_2wTv() +
                       (kSqrt4/kSqrt3) * fkr.Lamda() * fkr.Rminus_wR();
     hampl->fPlus1  = -(kSqrt3/kSqrt2) * fkr.Tplus_2wTv() +
                       (kSqrt4/kSqrt3) * fkr.Lamda()* fkr.Rplus_wR();
     hampl->fMinus3 = -(3./kSqrt2) * fkr.Tminus_2wTv();
     hampl->fPlus3  = -(3./kSqrt2) * fkr.Tplus_2wTv();
     hampl->f0Minus =  kSqrt3 * fkr.LS() * (1-2*xi) - (2./kSqrt3) * fkr.LC();
     hampl->f0Plus  =  kSqrt3 * fkr.LS() * (1-2*xi) + (2./kSqrt3) * fkr.LC();
     break;

   case (kS11_1650) :
     hampl->fMinus1 = -(0.5/kSqrt6) * fkr.Lamda() * fkr.Rminus_4wR();
     hampl->fPlus1  =  (0.5/kSqrt6) * fkr.Lamda() * fkr.Rplus_4wR();
     hampl->fMinus3 = 0.;
     hampl->fPlus3  = 0.;
     hampl->f0Minus = (1./kSqrt6) * fkr.LC_3B();
     hampl->f0Plus  = hampl->f0Minus;
     break;

   case (kD13_1700) :
     hampl->fMinus1 = -0.5*(1./kSqrt30) * fkr.Lamda() * fkr.Rminus_4wR();
     hampl->fPlus1  = -0.5*(1./kSqrt30) * fkr.Lamda() * fkr.Rplus_4wR();
     hampl->fMinus3 = -0.5*(3./kSqrt10) * fkr.Lamda() * fkr.Rminus_4wR();
     hampl->fPlus3  = -0.5*(3./kSqrt10) * fkr.Lamda() * fkr.Rplus_4wR();
     hampl->f0Minus = -0.5*(kSqrt2/kSqrt15) * fkr.LC();
     hampl->f0Plus  = -1. * hampl->f0Minus;
     break;

   case (kD15_1675) :
     hampl->fMinus1 =  0.5*(kSqrt3/kSqrt10) * fkr.Lamda() * fkr.Rminus_4wR();
     hampl->fPlus1  = -0.5*(kSqrt3/kSqrt10) * fkr.Lamda() * fkr.Rplus_4wR();
     hampl->fMinus3 =  0.5*(kSqrt3/kSqrt5) * fkr.Lamda() * fkr.Rminus_4wR();
     hampl->fPlus3  = -0.5*(kSqrt3/kSqrt5) * fkr.Lamda() * fkr.Rplus_4wR();
     hampl->f0Minus = -(kSqrt3/kSqrt10) * fkr.LC();
     hampl->f0Plus  = hampl->f0Minus;
     break;

   case (kS31_1620) :
     hampl->fMinus1 = kSqrt3* fkr.Tminus_2wTv() -
                             (1./kSqrt6)* fkr.Lamda()*fkr.Rminus_2wR();
     hampl->fPlus1  = -kSqrt3 * fkr.Tplus_2wTv()  +
                             (1./kSqrt6)* fkr.Lamda()*fkr.Rplus_2wR();
     hampl->fMinus3 = 0.;
     hampl->fPlus3  = 0.;
     hampl->f0Minus = -(kSqrt3/kSqrt2)* fkr.LS()*(1-2*xi) - (1./kSqrt6)* fkr.LC_3B();
     hampl->f0Plus  =  (kSqrt3/kSqrt2)* fkr.LS()*(1-2*xi) - (1./kSqrt6)* fkr.LC_2B();
     break;

   case (kD33_1700) :
     hampl->fMinus1 = (kSqrt3/kSqrt2)* fkr.Tminus_2wTv() +
                             (1./kSqrt3) * fkr.Lamda() * fkr.Rminus_2wR();
     hampl->fPlus1  = (kSqrt3/kSqrt2) * fkr.Tplus_2wTv() +
                             (1./kSqrt3) * fkr.Lamda() * fkr.Rplus_2wR();
     hampl->fMinus3 = (3./kSqrt2) * fkr.Tminus_2wTv();
     hampl->fPlus3  = (3./kSqrt2) * fkr.Tplus_2wTv();
     hampl->f0Minus = -kSqrt3 * fkr.LS() * (1-2*xi) - (1./kSqrt3) * fkr.LC();
     hampl->f0Plus  = -kSqrt3 * fkr.LS() * (1-2*xi) + (1./kSqrt3) * fkr.LC();
     break;

   case (kP11_1440) :
     hampl->fMinus1 = (5./12.)*kSqrt3*fkr.L2() *(fkr.Rminus()+ 2*xi*(4./5.)*fkr.R());
     hampl->fPlus1  = (5./12.)*kSqrt3*fkr.L2() *(fkr.Rplus() + 2*xi*(4./5.)*fkr.R());
     hampl->fMinus3 = 0.;
     hampl->fPlus3  = 0.;
     hampl->f0Minus = 0.5*(kSqrt3/2.) * fkr.L2S() -
                                (5./12.)*kSqrt3 * fkr.Lamda() * fkr.LC_2B();
     hampl->f0Plus  = 0.5*(kSqrt3/2.) * fkr.L2S() +
                                (5./12.)*kSqrt3 * fkr.Lamda() * fkr.LC_2B();
     break;

   case (kP33_1600) :
     hampl->fMinus1 = 0.;
     hampl->fPlus1  = 0.;
     hampl->fMinus3 = 0.;
     hampl->fPlus3  = 0.;
     hampl->f0Minus = 0.;
     hampl->f0Plus  = 0.;
     break;

   case (kP13_1720) :
     hampl->fMinus1 =  0.5*(kSqrt27/kSqrt10) * fkr.LTminus() +
                       0.5*(kSqrt5/kSqrt3) * fkr.L2() * (fkr.Rminus() + (8./5.)*xi*fkr.R());
     hampl->fPlus1  = -0.5*(kSqrt27/kSqrt10) * fkr.LTplus() -
                       0.5*(kSqrt5/kSqrt3) * fkr.L2() * (fkr.Rplus()  + (8./5.)*xi*fkr.R());
     hampl->fMinus3 = -0.5*(3./kSqrt10) * fkr.LTminus();
     hampl->fPlus3  =  0.5*(3./kSqrt10) * fkr.LTplus();
     hampl->f0Minus = -0.5*(kSqrt3/kSqrt5) * fkr.L2S()  +
                                0.5*(kSqrt5/kSqrt3) * fkr.Lamda()* fkr.LC_5B();
     hampl->f0Plus  =  0.5*(kSqrt3/kSqrt5) * fkr.L2S()  +
                                0.5*(kSqrt5/kSqrt3) * fkr.Lamda()* fkr.LC_5B();
     break;

   case (kF15_1680) :
     hampl->fMinus1 =  0.5*(3./kSqrt5) * fkr.LTminus() -
                       0.5*(kSqrt5/kSqrt2) * fkr.L2() * (fkr.Rminus() + 2*xi*(4./5.)*fkr.R());
     hampl->fPlus1  =  0.5*(3./kSqrt5) * fkr.LTplus()  -
                       0.5*(kSqrt5/kSqrt2) * fkr.L2() * (fkr.Rplus()  + 2*xi*(4./5.)*fkr.R());
     hampl->fMinus3 =  0.5*(kSqrt18/kSqrt5) * fkr.LTminus();
     hampl->fPlus3  =  0.5*(kSqrt18/kSqrt5) * fkr.LTplus();
     hampl->f0Minus = -0.5*(3./kSqrt10)* fkr.L2S() + 0.5*(kSqrt5/kSqrt2) * fkr.L2C();
     hampl->f0Plus  = -0.5*(3./kSqrt10)* fkr.L2S() - 0.5*(kSqrt5/kSqrt2) * fkr.L2C();
     break;

   case (kP31_1910) :
     hampl->fMinus1 = -(1./kSqrt15) * fkr.L2() * fkr.Rminus_2wR();
     hampl->fPlus1  = -(1./kSqrt15) * fkr.L2() * fkr.Rplus_2wR();
     hampl->fMinus3 = 0.;
     hampl->fPlus3  = 0.;
     hampl->f0Minus = -(2./kSqrt15) * fkr.Lamda() * fkr.LC_5B();
     hampl->f0Plus  = -1. * hampl->f0Minus;
     break;

   case (kP33_1920) :
     hampl->fMinus1 =  (1./kSqrt15) * fkr.L2() * fkr.Rminus_2wR();
     hampl->fPlus1  = -(1./kSqrt15) * fkr.L2() * fkr.Rplus_2wR();
     hampl->fMinus3 = -(1./kSqrt5 ) * fkr.L2() * fkr.Rminus_2wR();
     hampl->fPlus3  =  (1./kSqrt5 ) * fkr.L2() * fkr.Rplus_2wR();
     hampl->f0Minus = -(2./kSqrt15) * fkr.Lamda() * fkr.LC_5B();
     hampl->f0Plus  = hampl->f0Minus;
     break;

   case (kF35_1905) :
     hampl->fMinus1 = (1./kSqrt35) * fkr.L2() * fkr.Rminus_2wR();
     hampl->fPlus1  = (1./kSqrt35) * fkr.L2() * fkr.Rplus_2wR();
     hampl->fMinus3 = (kSqrt18/kSqrt35) * fkr.L2() * fkr.Rminus_2wR();
     hampl->fPlus3  = (kSqrt18/kSqrt35) * fkr.L2() * fkr.Rplus_2wR();
     hampl->f0Minus = (2./kSqrt35) * fkr.L2C();
     hampl->f0Plus  = -1. * hampl->f0Minus;
     break;

   case (kF37_1950) :
     hampl->fMinus1 =  -(kSqrt6/kSqrt35) * fkr.L2() * fkr.Rminus_2wR();
     hampl->fPlus1  =   (kSqrt6/kSqrt35) * fkr.L2() * fkr.Rplus_2wR();
     hampl->fMinus3 =  -(kSqrt2/kSqrt7) * fkr.L2() * fkr.Rminus_2wR();
     hampl->fPlus3  =   (kSqrt2/kSqrt7) * fkr.L2() * fkr.Rplus_2wR();
     hampl->f0Minus = 2*(kSqrt6/kSqrt35) * fkr.L2C();
     hampl->f0Plus  =  hampl->f0Minus;
     break;

   case (kP11_1710) :
     hampl->fMinus1 = -(1./kSqrt6) * fkr.L2() * fkr.Rminus_wR();
     hampl->fPlus1  = -(1./kSqrt6) * fkr.L2() * fkr.Rplus_wR();
     hampl->fMinus3 = 0.;
     hampl->fPlus3  = 0.;
     hampl->f0Minus = -(kSqrt3/kSqrt8) * fkr.L2S() * (1-2*xi) +
                               (1./kSqrt6) * fkr.Lamda() * fkr.LC_2B();
     hampl->f0Plus  = -(kSqrt3/kSqrt8) * fkr.L2S() * (1-2*xi) -
                               (1./kSqrt6) * fkr.Lamda() * fkr.LC_2B();
     break;

   case (kF17_1970) :
     hampl->fMinus1 = -2.5 * (fkr.Rminus() + 1.6 * xi * fkr.R());
     hampl->fPlus1  = -2.5 * (fkr.Rplus() + 2.4*xi*fkr.R());
     hampl->fMinus3 = 0.;
     hampl->fPlus3  = 0.;
     hampl->f0Minus = -1.5 * fkr.S() + 2.5 * fkr.C();
     hampl->f0Plus  = -1.5 * fkr.S() - 2.5 * fkr.C();
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

