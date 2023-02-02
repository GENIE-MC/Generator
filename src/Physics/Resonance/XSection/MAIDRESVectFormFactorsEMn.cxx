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
#include "Physics/Resonance/XSection/MAIDRESVectFormFactorsEMn.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
MAIDRESVectFormFactorsEMn::MAIDRESVectFormFactorsEMn() :
  RESVectFormFactorsI("genie::MAIDRESVectFormFactorsEMn")
{

}
//____________________________________________________________________________
MAIDRESVectFormFactorsEMn::MAIDRESVectFormFactorsEMn(string config) :
  RESVectFormFactorsI("genie::MAIDRESVectFormFactorsEMn", config)
{

}
//____________________________________________________________________________
MAIDRESVectFormFactorsEMn::~MAIDRESVectFormFactorsEMn()
{

}
//____________________________________________________________________________
RESVectFFAmplitude MAIDRESVectFormFactorsEMn::Compute( const Interaction interaction ) const {
  RESVectFFAmplitude ampl ; 

  const InitialState & init_state = interaction.InitState();
  const ProcessInfo &  proc_info  = interaction.ProcInfo();
  const Target & target = init_state.Tgt();
  const Resonance_t res = interaction.ExclTag().Resonance();
  
  //Get kinematical parameters
  const Kinematics & kinematics = interaction.Kine();
  double W  = kinematics.W();
  double q2 = kinematics.q2();

  //Compute auxiliary & kinematical factors
  double Mnuc   = target.HitNucMass();
  double W2     = TMath::Power(W,    2);
  double Mnuc2  = TMath::Power(Mnuc, 2);

  double MR  = utils::res::Mass(res);
  double kR = (MR*MR - Mnuc2)/(2.*MR);  //LAR choice
                                        //k(W=MR,Q2=0) Q2=0 real photon
  double kgcm0 = (W2 - Mnuc2)/(2.*MR); // kW at equation 5
  double egcm = (W2+q2-Mnuc2)/(2.*MR); //photon energy at the W center of mass frame
  double qcm = TMath::Sqrt(egcm*egcm-q2); //photon momentum k in Equation 3

  double tau = -q2/(4.*Mnuc2);  

  switch(res) {

  case (kP33_1232) :
    {
      double Fq = 1./TMath::Power(1-q2/0.71,2)*qcm/kgcm0;
      double AM = 300.*(1. - 0.01*q2)*TMath::Exp(0.23*q2)*Fq;
      double AE = -6.37 * (1. + 0.021*q2)*TMath::Exp(0.16*q2)*Fq;
      double AS = -12.40 * (1. - 0.12*q2) / (1. + 4.9*tau)*qcm/kR*TMath::Exp(0.23*q2)*Fq;
      ampl.SetAmplA12( (3.*AE+AM)/2./1000. );
      ampl.SetAmplA32( TMath::Sqrt(3.)/2.*(AE-AM)/1000. );
      ampl.SetAmplS12( TMath::Sqrt(2.)*AS/1000. );
      break;
    }
  case (kS11_1535) :
    {
      ampl.SetAmplA12( -51.*(1.-4.75*q2)*exp(1.69*q2)/1000. );
      ampl.SetAmplS12( 28.5*(1.-0.36*q2)*exp(1.55*q2)/1000. );
      break;
    }
  case (kD13_1520) :
    {
      ampl.SetAmplA12( -77.*(1.+0.53*q2)*exp(1.55*q2)/1000. );
      ampl.SetAmplA32( -154.*(1.-0.58*q2)*exp(1.75*q2)/1000. );
      ampl.SetAmplS12( 13.6*(1.-15.7*q2)*exp(1.57*q2)/1000. );
      break;
    }
  case (kS11_1650) :
    {
      ampl.SetAmplA12( 9.*(1.-0.13*q2)*exp(1.55*q2)/1000. );
      ampl.SetAmplS12( 10.1*(1.+0.50*q2)*exp(1.55*q2)/1000. );
      break;
    }
  case (kD15_1675) :
    {
      ampl.SetAmplA12( -62.*(1.-0.01*q2)*exp(+2.00*q2)/1000. );
      ampl.SetAmplA32( -84.*(1.-0.01*q2)*exp(+2.00*q2)/1000. );
      ampl.SetAmplS12( 0. );
      break;
    }
  case (kS31_1620) :
    {
      ampl.SetAmplA12( 66.*(1.-1.86*q2)*exp(+2.5*q2)/1000. );
      ampl.SetAmplS12( 16.2*(1.-2.83*q2)*exp(+2.0*q2)/1000. );
      break;
    }
  case (kD33_1700) :
    {
      ampl.SetAmplA12( 226.*(1.-1.91*q2)*exp(+1.77*q2)/1000. );
      ampl.SetAmplA32( 210.*(1.-1.97*q2)*exp(+2.2*q2)/1000. );
      ampl.SetAmplS12( 0. );
      break;
    }
  case (kP11_1440) :
    {
      ampl.SetAmplA12( 54.1*(1.-0.95*q2)*exp(+1.77*q2)/1000. );
      ampl.SetAmplS12( -41.5*(1.-2.98*q2)*exp(+1.55*q2)/1000. );
      ampl.SetAmplA32( 0. ) ;// ??
      break;
    }
  case (kP13_1720) :
    {
      ampl.SetAmplA12( -3*(1.-12.7*q2)*exp(+1.55*q2)/1000. );
      ampl.SetAmplA32( -31*(1.-4.99*q2)*exp(+1.55*q2)/1000. );
      ampl.SetAmplS12( 0. );
      break;
    }
  case (kF15_1680) :
    {
      ampl.SetAmplA12( 28*(1.-0.*q2)*exp(+1.2*q2)/1000. );
      ampl.SetAmplA32( -38*(1.-4.09*q2)*exp(+1.75*q2)/1000. );
      ampl.SetAmplS12( 0. );
      break;
    }
  case (kD13_1700) :
    {
      // Not implemented
      break;
    }
  case (kP33_1600) :
    {
      // Not implemented 
      break;
    }
  case (kP31_1910) :
    {
      // Not implemented
      break;
    }
  case (kP33_1920) :
    {
      // Not implemented
      break;
    }
  case (kF35_1905) :
    {
      // Not implemented
      break;
    }
  case (kF37_1950) :
    {
      // Not implemented
      break;
    }
  case (kP11_1710) :
    {
      // Not implemented
      break;
    }
  case (kF17_1970) :
    {
      // Not implemented
      break;
    }
  default:
    {
      LOG("MAIDRESVectFormFactorsEMn", pWARN) << "*** UNRECOGNIZED RESONANCE!";
      break;
    }

  }//switch

  return ampl;
}
//____________________________________________________________________________
void MAIDRESVectFormFactorsEMn::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void MAIDRESVectFormFactorsEMn::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfig();
}
//____________________________________________________________________________
void MAIDRESVectFormFactorsEMn::LoadConfig(void)
{
  // LOAD PARAMETERS HERE 
}
