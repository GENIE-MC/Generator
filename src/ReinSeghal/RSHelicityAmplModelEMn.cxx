//____________________________________________________________________________
/*!

\class    genie::RSHelicityAmplModelEMn

\brief    The Helicity Amplitudes, for all baryon resonances, for Electro-
          Magnetic (EM) interactions on free neutrons, as computed in the
          Rein-Seghal's paper.

          Concrete implementation of the RSHelicityAmplModelI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  March 30, 2005

*/
//____________________________________________________________________________

#include "BaryonResonance/BaryonResUtils.h"
#include "Conventions/Constants.h"
#include "ReinSeghal/RSHelicityAmplModelEMn.h"
#include "ReinSeghal/RSHelicityAmpl.h"
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
RSHelicityAmpl * RSHelicityAmplModelEMn::Compute(
                                       Resonance_t res, const FKR & fkr) const
{
  RSHelicityAmpl * hampl = new RSHelicityAmpl;

  switch(res) {

   case (kP33_1232) :
     hampl->fMinus1 = -kSqrt2 * fkr.R();
     hampl->fPlus1  =  kSqrt2 * fkr.R();
     hampl->fMinus3 = -kSqrt6 * fkr.R();
     hampl->fPlus3  =  kSqrt6 * fkr.R();
     hampl->f0Minus =  0.;
     hampl->f0Plus  =  0.;
     break;

   case (kS11_1535) :
     hampl->fMinus1 = -kSqrt3 * fkr.T() - (1./kSqrt6) * fkr.LR();
     hampl->fPlus1  =  kSqrt3 * fkr.T() + (1./kSqrt6) * fkr.LR();
     hampl->fMinus3 =   0.;
     hampl->fPlus3  =   0.;
     hampl->f0Minus =  (kSqrt3/kSqrt2) * fkr.LS();
     hampl->f0Plus  = -(kSqrt3/kSqrt2) * fkr.LS();
     break;

   case (kD13_1520) :
     hampl->fMinus1 = -(kSqrt3/kSqrt2) * fkr.T() + (1./kSqrt3) * fkr.LR();
     hampl->fPlus1  = -(kSqrt3/kSqrt2) * fkr.T() + (1./kSqrt3) * fkr.LR();
     hampl->fMinus3 =  (3./kSqrt2) * fkr.T();
     hampl->fPlus3  = -(3./kSqrt2) * fkr.T();
     hampl->f0Minus =  kSqrt3 * fkr.LS();
     hampl->f0Plus  =  kSqrt3 * fkr.LS();
     break;

   case (kS11_1650) :
     hampl->fMinus1 = -(1./kSqrt6) * fkr.LR();
     hampl->fPlus1  =  (1./kSqrt6) * fkr.LR();
     hampl->fMinus3 =  0.;
     hampl->fPlus3  =  0.;
     hampl->f0Minus =  0.;
     hampl->f0Plus  =  0.;
     break;

   case (kD13_1700) :
     hampl->fMinus1 = -(1./kSqrt30) * fkr.LR();
     hampl->fPlus1  = -(1./kSqrt30) * fkr.LR();
     hampl->fMinus3 = -(3./kSqrt10) * fkr.LR();
     hampl->fPlus3  = -(3./kSqrt10) * fkr.LR();
     hampl->f0Minus =  0.;
     hampl->f0Plus  =  0.;
     break;

   case (kD15_1675) :
     hampl->fMinus1 =  (kSqrt3/kSqrt10) * fkr.LR();
     hampl->fPlus1  = -(kSqrt3/kSqrt10) * fkr.LR();
     hampl->fMinus3 =  (kSqrt3/kSqrt5)  * fkr.LR();
     hampl->fPlus3  = -(kSqrt3/kSqrt5)  * fkr.LR();
     hampl->f0Minus =  0.;
     hampl->f0Plus  =  0.;
     break;

   case (kS31_1620) :
     hampl->fMinus1 =  kSqrt3 * fkr.T() - (1./kSqrt6) * fkr.LR();
     hampl->fPlus1  = -kSqrt3 * fkr.T() + (1./kSqrt6) * fkr.LR();
     hampl->fMinus3 =  0.;
     hampl->fPlus3  =  0.;
     hampl->f0Minus = -(kSqrt3/kSqrt2) * fkr.LS();
     hampl->f0Plus  =  (kSqrt3/kSqrt2) * fkr.LS();
     break;

   case (kD33_1700) :
     hampl->fMinus1 = (kSqrt3/kSqrt2) * fkr.T() + (1./kSqrt3) * fkr.LR();
     hampl->fPlus1  = (kSqrt3/kSqrt2) * fkr.T() + (1./kSqrt3) * fkr.LR();
     hampl->fMinus3 =  (3./kSqrt2) * fkr.T();
     hampl->fPlus3  =  (3./kSqrt2) * fkr.T();
     hampl->f0Minus = -kSqrt3 * fkr.LS();
     hampl->f0Plus  = -kSqrt3 * fkr.L2S();
     break;

   case (kP11_1440) :
     hampl->fMinus1 = (1./kSqrt3) * fkr.L2R();
     hampl->fPlus1  = (1./kSqrt3)*fkr.L2R();
     hampl->fMinus3 =  0.;
     hampl->fPlus3  =  0.;
     hampl->f0Minus =  0.;
     hampl->f0Plus  =  0.;
     break;

   case (kP33_1600) :
     hampl->fMinus1 =  (1./kSqrt6) * fkr.L2R();
     hampl->fPlus1  = -(1./kSqrt6) * fkr.L2R();
     hampl->fMinus3 =  (1./kSqrt2) * fkr.L2R();
     hampl->fPlus3  = -(1./kSqrt2) * fkr.L2R();
     hampl->f0Minus =  0;
     hampl->f0Plus  =  0;
     break;

   case (kP13_1720) :
     hampl->fMinus1 =  (2./kSqrt15) * fkr.L2R();
     hampl->fPlus1  = -(2./kSqrt15) * fkr.L2R();
     hampl->fMinus3 =  0.;
     hampl->fPlus3  =  0.;
     hampl->f0Minus =  0.;
     hampl->f0Plus  =  0.;
     break;

   case (kF15_1680) :
     hampl->fMinus1 = -(kSqrt2/kSqrt5) * fkr.L2R();
     hampl->fPlus1  = -(kSqrt2/kSqrt5) * fkr.L2R();
     hampl->fMinus3 =  0.;
     hampl->fPlus3  =  0.;
     hampl->f0Minus =  0.;
     hampl->f0Plus  =  0.;
     break;

   case (kP31_1910) :
     hampl->fMinus1 = -(1./kSqrt15) * fkr.L2R();
     hampl->fPlus1  = -(1./kSqrt15) * fkr.L2R();
     hampl->fMinus3 =  0.;
     hampl->fPlus3  =  0.;
     hampl->f0Minus =  0.;
     hampl->f0Plus  =  0.;
     break;

   case (kP33_1920) :
     hampl->fMinus1 =  (1./kSqrt15) * fkr.L2R();
     hampl->fPlus1  = -(1./kSqrt15) * fkr.L2R();
     hampl->fMinus3 = -(1./kSqrt5)  * fkr.L2R();
     hampl->fPlus3  =  (1./kSqrt5)  * fkr.L2R();
     hampl->f0Minus =  0.;
     hampl->f0Plus  =  0.;
     break;

   case (kF35_1905) :
     hampl->fMinus1 = (1./kSqrt35) * fkr.L2R();
     hampl->fPlus1  = (1./kSqrt35) * fkr.L2R();
     hampl->fMinus3 = (kSqrt18/kSqrt35) * fkr.L2R();
     hampl->fPlus3  = (kSqrt18/kSqrt35) * fkr.L2R();
     hampl->f0Minus =  0.;
     hampl->f0Plus  =  0.;
     break;

   case (kF37_1950) :
     hampl->fMinus1 = -(kSqrt6/kSqrt35) * fkr.L2R();
     hampl->fPlus1  =  (kSqrt6/kSqrt35) * fkr.L2() * fkr.Rplus_2wR();
     hampl->fMinus3 = -(kSqrt2/kSqrt7)  * fkr.L2R();
     hampl->fPlus3  =  (kSqrt2/kSqrt7)  * fkr.L2R();
     hampl->f0Minus =  0.;
     hampl->f0Plus  =  0.;
     break;

   case (kP11_1710) :
     hampl->fMinus1 = -(1./kSqrt24) * fkr.L2R();
     hampl->fPlus1  = -(1./kSqrt24) * fkr.L2R();
     hampl->fMinus3 =  0.;
     hampl->fPlus3  =  0.;
     hampl->f0Minus = -(kSqrt3/kSqrt8) * fkr.L2S();
     hampl->f0Plus  = -(kSqrt3/kSqrt8) * fkr.L2S();
     break;

   case (kF17_1970) :
     hampl->fMinus1 =  (kSqrt3/kSqrt35) * fkr.L2R();
     hampl->fPlus1  = -(kSqrt3/kSqrt35) * fkr.L2R();
     hampl->fMinus3 =  (1./kSqrt7) * fkr.L2R();
     hampl->fPlus3  =  0.;
     hampl->f0Minus =  0.;
     hampl->f0Plus  =  0.;
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

