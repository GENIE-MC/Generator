//____________________________________________________________________________
/*
  Copyright (c) 2023-2025, The GENIE Collaboration
  For the full text of the license visit http://copyright.genie-mc.org

  Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
  University of Liverpool & STFC Rutherford Appleton Laboratory
*/
//____________________________________________________________________________

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/Resonance/XSection/MAIDHelicityAmplModelEMp.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
MAIDHelicityAmplModelEMp::MAIDHelicityAmplModelEMp() :
  MAIDHelicityAmplModelI("genie::MAIDHelicityAmplModelEMp")
{

}
//____________________________________________________________________________
MAIDHelicityAmplModelEMp::MAIDHelicityAmplModelEMp(string config) :
  MAIDHelicityAmplModelI("genie::MAIDHelicityAmplModelEMp", config)
{

}
//____________________________________________________________________________
MAIDHelicityAmplModelEMp::~MAIDHelicityAmplModelEMp()
{

}
//____________________________________________________________________________
const MAIDHelicityAmpl & MAIDHelicityAmplModelEMp::Compute( const Interaction * interaction ) const {

  const InitialState & init_state = interaction -> InitState();
  const Target & target = init_state.Tgt();
  const Resonance_t res = interaction->ExclTag().Resonance();
  
  //Get kinematical parameters
  const Kinematics & kinematics = interaction -> Kine();
  double W  = kinematics.W();
  double q2 = kinematics.q2();

  //Get Baryon resonance parameters
  double MR  = utils::res::Mass(res);

  //Compute auxiliary & kinematical factors
  double Mnuc   = target.HitNucMass();
  double W2     = TMath::Power(W,2);
  double Mnuc2  = TMath::Power(Mnuc,2);
  double kR = (MR*MR - Mnuc2)/(2.*MR);  //LAR choice
                                        //k(W=MR,Q2=0) Q2=0 real photon
  double kgcm0 = (W2 - Mnuc2)/(2.*MR); // kW at equation 5
  double egcm = (W2+q2-Mnuc2)/(2.*MR); //photon energy at the W center of mass frame
  double qcm = TMath::Sqrt(egcm*egcm-q2); //photon momentum k in Equation 3
  double tau = -q2/(4.*Mnuc2); //tau is 

  switch(res) {
  case (kP33_1232) :
    {
      double Fq = 1./TMath::Power(1-q2/0.71,2)*qcm/kgcm0;
      double AM = 300.*(1. - 0.01*q2)*TMath::Exp(0.23*q2)*Fq;
      double AE = -6.37 * (1. + 0.021*q2)*TMath::Exp(0.16*q2)*Fq;
      double AS = -12.40 * (1. - 0.12*q2) / (1. + 4.9*tau)*qcm/kR*TMath::Exp(0.23*q2)*Fq;
      fAmpl.fA12= -(3.*AE+AM)/2./1000.;
      fAmpl.fA32= TMath::Sqrt(3.)/2.*(AE-AM)/1000.;
      fAmpl.fS12 = TMath::Sqrt(2.)*AS/1000.;
      break ; 
    }
  case (kS11_1535) :
    {
      fAmpl.fA12 = 66.*(1.-1.61*q2)*exp(+0.70*q2)/1000.;
      fAmpl.fA32 = 0 ;
      fAmpl.fS12 = -2.*(1.-23.9*q2)*exp(+0.81*q2)/1000.;
      break ;
    }
  case (kD13_1520) :
    {
      fAmpl.fA12 = -27.*(1.-7.77*q2)*exp(1.09*q2)/1000.;
      fAmpl.fA32 = 161.*(1.-0.69*q2)*exp(2.1*q2)/1000.;
      fAmpl.fS12 = -63.6*(1.-4.19*q2)*exp(3.40*q2)/1000.;
      break ;
    }
  case (kS11_1650) :
    {
      fAmpl.fA12 = 33.*(1.-1.45*q2)*exp(+0.62*q2)/1000.;
      fAmpl.fA32 = 0;
      fAmpl.fS12 = -3.5*(1.-2.88*q2)*exp(+0.76*q2)/1000.;
      break ; 
    }
  case (kD13_1700) :
    {
      // Not implemented 
      break ; 
    }
  case (kD15_1675) :
    {
      fAmpl.fA12 = 15.*(1.-0.1*q2)*exp(+2.00*q2)/1000.;
      fAmpl.fA32 = 22.*(1.-0.1*q2)*exp(+2.00*q2)/1000.;
      fAmpl.fS12 = 0.;
      break ;
    }
  case (kS31_1620) :
    {
      fAmpl.fA12 =  66.*(1.-1.86*q2)*exp(+2.5*q2)/1000.;
      fAmpl.fA32 = 0;
      fAmpl.fS12 = 16.2*(1.-2.83*q2)*exp(+2.0*q2)/1000.;
      break ; 
    }
  case (kD33_1700) :
    {
      fAmpl.fA12 = 226.*(1.-1.91*q2)*exp(+1.77*q2)/1000.;
      fAmpl.fA32 =  210.*(1.-1.97*q2)*exp(+2.2*q2)/1000.;
      fAmpl.fS12 = 0.;
      // ??
      fAmpl.fA12=0;
      fAmpl.fA32=0; //set to 0 to get agreement with (e,e')
      break ; 
    }
  case (kP11_1440) :
    {
      fAmpl.fA12 = -61.4*(1.+1.22*q2-0.55*q2*q2*q2*q2)*exp(+1.51*q2)/1000.;
      fAmpl.fS12 = 4.2*(1.-40*q2+1.5*q2*q2*q2*q2)*exp(+1.75*q2)/1000.;
      break ; 
    }
  case (kP33_1600) :
    {
      // Not implemented 
      break ; 
    }
  case (kP13_1720) :
    {
      fAmpl.fA12 = 73.*(1.-1.89*q2)*exp(+1.55*q2)/1000.;
      fAmpl.fA32 = -11.*(1.-16.*q2)*exp(+1.55*q2)/1000.;
      fAmpl.fS12 = -53.*(1.-2.46*q2)*exp(+1.55*q2)/1000.;
      break ; 
    }
  case (kF15_1680) :
    {
      fAmpl.fA12 =  -25.*(1.-3.98*q2)*exp(+1.20*q2)/1000.;
      fAmpl.fA32 =  134.*(1.-1.*q2)*exp(+2.22*q2)/1000.;
      fAmpl.fS12 =  -44.*(1.-3.14*q2)*exp(+1.68*q2)/1000.;
      break ; 
    }
  case (kP31_1910) :
    {
      // Not implemented
      break ;
    }
  case (kP33_1920) :
    {
      // Not implemented 
      break ; 
    }
  case (kF35_1905) :
    {
      // Not implemented
      break ; 
    }
  case (kF37_1950) :
    {
      // Not implemented
      break ; 
    }
  case (kP11_1710) :
    {
      // Not implemented
      break ; 
    }
  case (kF17_1970) :
    {
      // Not implemented
      break ; 
    }
  default:
    {
      LOG("MAIDHAmpl", pWARN) << "*** UNRECOGNIZED RESONANCE!";
      fAmpl.fA12 = 0.;
      fAmpl.fA32  = 0.;
      fAmpl.fS12 = 0.;
      break;
    }

      
  }//switch

  return fAmpl;
}
//____________________________________________________________________________
