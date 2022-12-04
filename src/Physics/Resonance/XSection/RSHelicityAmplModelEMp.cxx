//____________________________________________________________________________
/*
 Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Oct 05, 2009 - CA
   Compute() now returns a `const RSHelicityAmpl &' and avoids creating a new
   RSHelicityAmpl at each call.                      

*/
//____________________________________________________________________________

#include "BaryonResonance/BaryonResUtils.h"
#include "Conventions/Constants.h"
#include "ReinSehgal/RSHelicityAmplModelEMp.h"
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
const RSHelicityAmpl & 
    RSHelicityAmplModelEMp::Compute(
          Resonance_t res, const FKR & fkr, const Interaction * interaction) const
{
  // need kinematical factors to compute MAID vector form factors
  const InitialState & init_state = interaction -> InitState();
  const ProcessInfo &  proc_info  = interaction -> ProcInfo();
  const Target & target = init_state.Tgt();
 
  bool is_EM     = proc_info.IsEM();
 
  //Get kinematical parameters
  const Kinematics & kinematics = interaction -> Kine();
  double W  = kinematics.W();
  double q2 = kinematics.q2();

  //Get Baryon resonance parameters
  int    IR  = utils::res::ResonanceIndex    (res);
  int    LR  = utils::res::OrbitalAngularMom (res);
  double MR  = utils::res::Mass              (res);
  double WR  = utils::res::Width             (res);
  double NR  = utils::res::BWNorm            (res);

  //Compute auxiliary & kinematical factors
  double E      = init_state.ProbeE(kRfHitNucRest);
  double Mnuc   = target.HitNucMass();
  double W2     = TMath::Power(W,    2);
  double Mnuc2  = TMath::Power(Mnuc, 2);
  //qvec is k(W,Q2) at W=MR 
  double qvec   = TMath::Sqrt(TMath::Power((MR*MR-Mnuc2+q2),2)/(4*MR*MR) - q2); 

  double kR = (MR*MR - Mnuc2)/(2.*MR);  //LAR choice
                                        //k(W=MR,Q2=0) Q2=0 real photon
  double kgcm0 = (W2 - Mnuc2)/(2.*MR); // kW at equation 5
  double egcm = (W2+q2-Mnuc2)/(2.*MR); //photon energy at the W center of mass frame
  double qcm = TMath::Sqrt(egcm*egcm-q2); //photon momentum k in Equation 3
                                          //qvec equals to qcm at W=MR
  double tau = -q2/(4.*Mnuc2); //tau is 
  LOG("ReinSeghalRes", pWARN) << "q2= " << q2 << " IR= "<< IR <<" LR= "<< LR <<" WR= "<< WR <<" NR= "<< NR <<" E= "<< E <<" W2= "<< W2 <<" is_EM "<<is_EM;


  //---------------------------------------



  switch(res) {

   case (kP33_1232) :
   {
     fAmpl.fPlus1  =  kSqrt2 * fkr.R;
     fAmpl.fPlus3  =  kSqrt6 * fkr.R;
     fAmpl.fMinus1 = -1 * fAmpl.fPlus1;
     fAmpl.fMinus3 = -1 * fAmpl.fPlus3;
     fAmpl.f0Minus =  0.;
     fAmpl.f0Plus  =  0.;
     //recalculae the parameters 
     //here MAID multipoles are used
     double Fq = 1./TMath::Power(1-q2/0.71,2)*qcm/kgcm0;
     double AM = 300.*(1. - 0.01*q2)*TMath::Exp(0.23*q2)*Fq;
     double AE = -6.37 * (1. + 0.021*q2)*TMath::Exp(0.16*q2)*Fq;
     double AS = -12.40 * (1. - 0.12*q2) / (1. + 4.9*tau)*qcm/kR*TMath::Exp(0.23*q2)*Fq;
     double A12= -(3.*AE+AM)/2./1000.;
     double A32= TMath::Sqrt(3.)/2.*(AE-AM)/1000.;
     double S12 = TMath::Sqrt(2.)*AS/1000.;

     LOG("ReinSeghalRes", pWARN) << "q2= " << q2 << " A12, A32, S12 =" <<A12<<" "<<A32<<"  "<<S12;

      
     fAmpl.fMinus1 =  TMath::Sqrt(115./130.)*TMath::Sqrt((Mnuc/MR)*kgcm0/(2.*kPi*kAem))*A12;
     fAmpl.fPlus1  = -fAmpl.fMinus1;
     fAmpl.fMinus3 =  TMath::Sqrt(115./130.)*TMath::Sqrt((Mnuc/MR)*kgcm0/(2.*kPi*kAem))*A32;
     fAmpl.fPlus3  = -fAmpl.fMinus3;
     fAmpl.f0Minus =  TMath::Sqrt(115./130.)*TMath::Sqrt((Mnuc/MR)*kgcm0/(2.*kPi*kAem))*(-q2/(qvec*qvec))*S12;
     fAmpl.f0Plus  = fAmpl.f0Minus;

     LOG("ReinSeghalRes", pWARN) << "MAID fminus1= " << fAmpl.fMinus1 << " fminus3 =" << fAmpl.fMinus3 ;


     break;
   }
   case (kS11_1535) :
   {
     fAmpl.fMinus1 =  kSqrt3 * fkr.T + kSqrt3_2 * fkr.Lamda * fkr.R;
     fAmpl.f0Minus = -kSqrt3_2 * fkr.Lamda * fkr.S;
     fAmpl.fPlus1  = -1. * fAmpl.fMinus1;
     fAmpl.f0Plus  = -1. * fAmpl.f0Minus;
     fAmpl.fMinus3 =  0.;
     fAmpl.fPlus3  =  0.;
     //change to MAID parameterizations of vector amplitude
     //numbers are from Equation 47 and Table 13 and table 7
     double A12 =  66.*(1.-1.61*q2)*exp(+0.70*q2)/1000.;
     double S12 = -2.*(1.-23.9*q2)*exp(+0.81*q2)/1000.;
     //euqation 24 
     fAmpl.fMinus1 =  TMath::Sqrt((Mnuc/MR)*kgcm0/(2.*kPi*kAem))*A12;
     fAmpl.f0Minus =  TMath::Sqrt((Mnuc/MR)*kgcm0/(2.*kPi*kAem))*(-q2/(qvec*qvec))*S12;
     fAmpl.fPlus1  = -1. * fAmpl.fMinus1;
     fAmpl.f0Plus  = -1. * fAmpl.f0Minus;
     fAmpl.fMinus3 =  0.;
     fAmpl.fPlus3  =  0.;
      



     //-------------------------------------------
     break;
   }
   case (kD13_1520) :
   {
     fAmpl.fMinus1 =  kSqrt3_2 * fkr.T - kSqrt3 * fkr.Lamda * fkr.R;
     fAmpl.fMinus3 =  k3_Sqrt2 * fkr.T;
     fAmpl.f0Minus = -kSqrt3 * fkr.Lamda * fkr.S;
     fAmpl.fPlus1  =  fAmpl.fMinus1;
     fAmpl.fPlus3  =  fAmpl.fMinus3;
     fAmpl.f0Plus  =  fAmpl.f0Minus;
     //change to MAID parameterization of vector amplitude
     double A12 = -27.*(1.-7.77*q2)*exp(1.09*q2)/1000.;
     double A32 =  161.*(1.-0.69*q2)*exp(2.1*q2)/1000.;
     double S12 = -63.6*(1.-4.19*q2)*exp(3.40*q2)/1000.;
     fAmpl.fMinus1 =  TMath::Sqrt((Mnuc/MR)*kgcm0/(2.*kPi*kAem))*A12;
     fAmpl.fPlus1  = -fAmpl.fMinus1;
     fAmpl.fMinus3 =  TMath::Sqrt((Mnuc/MR)*kgcm0/(2.*kPi*kAem))*A32;
     fAmpl.fPlus3  = -fAmpl.fMinus3;
     fAmpl.f0Minus =  TMath::Sqrt((Mnuc/MR)*kgcm0/(2.*kPi*kAem))*(-q2/(qvec*qvec))*S12;
     fAmpl.f0Plus  = fAmpl.f0Minus;



     //-------------------------------------------
     break;
   }
   case (kS11_1650) :
   {
     fAmpl.fMinus1 = 0.;
     fAmpl.fPlus1  = 0.;
     fAmpl.fMinus3 = 0.;
     fAmpl.fPlus3  = 0.;
     fAmpl.f0Minus = 0.;
     fAmpl.f0Plus  = 0.;
     //change to MAID parameterization of vector amplitude
     double A12 =  33.*(1.-1.45*q2)*exp(+0.62*q2)/1000.;
     double S12 = -3.5*(1.-2.88*q2)*exp(+0.76*q2)/1000.;
     fAmpl.fMinus1 =  TMath::Sqrt((Mnuc/MR)*kgcm0/(2.*kPi*kAem))*A12;
     fAmpl.f0Minus =  TMath::Sqrt((Mnuc/MR)*kgcm0/(2.*kPi*kAem))*(-q2/(qvec*qvec))*S12;
     fAmpl.fPlus1  = -1. * fAmpl.fMinus1;
     fAmpl.f0Plus  = -1. * fAmpl.f0Minus;
     fAmpl.fMinus3 =  0.;
     fAmpl.fPlus3  =  0.;



     //----------------------------------------------
     break;
   }
   case (kD13_1700) :
   {
     fAmpl.fMinus1 = 0.;
     fAmpl.fPlus1  = 0.;
     fAmpl.fMinus3 = 0.;
     fAmpl.fPlus3  = 0.;
     fAmpl.f0Minus = 0.;
     fAmpl.f0Plus  = 0.;
     // change to MAID parameterization of vector amplitude
     

     //---------------------------------------------------- 
     break;
   }
   case (kD15_1675) :
   {
     fAmpl.fMinus1 = 0.;
     fAmpl.fPlus1  = 0.;
     fAmpl.fMinus3 = 0.;
     fAmpl.fPlus3  = 0.;
     fAmpl.f0Minus = 0.;
     fAmpl.f0Plus  = 0.;
     //change to MAID parameterization of vector amplitude
     double A12 =  15.*(1.-0.1*q2)*exp(+2.00*q2)/1000.;
     double A32 =  22.*(1.-0.1*q2)*exp(+2.00*q2)/1000.;
     double S12 = 0.;
     fAmpl.fMinus1 =  TMath::Sqrt((Mnuc/MR)*kgcm0/(2.*kPi*kAem))*A12;
     fAmpl.fPlus1  = -fAmpl.fMinus1;
     fAmpl.fMinus3 =  TMath::Sqrt((Mnuc/MR)*kR/(2.*kPi*kAem))*A32;
     fAmpl.fPlus3  = -fAmpl.fMinus3;
     fAmpl.f0Minus =  TMath::Sqrt((Mnuc/MR)*kR/(2.*kPi*kAem))*(-q2/(qvec*qvec))*S12;
     fAmpl.f0Plus  = fAmpl.f0Minus;


     
     break;
   }
   case (kS31_1620) :
   {
     fAmpl.fMinus1 =  kSqrt3 * fkr.T - k1_Sqrt6 * fkr.Lamda * fkr.R;
     fAmpl.f0Minus = -kSqrt3_2 * fkr.Lamda * fkr.S;
     fAmpl.fPlus1  = -1. * fAmpl.fMinus1;
     fAmpl.f0Plus  = -1. * fAmpl.f0Minus;
     fAmpl.fMinus3 = 0.;
     fAmpl.fPlus3  = 0.;

     //change to MAID parameterization of vector amplitude
     double A12 =  66.*(1.-1.86*q2)*exp(+2.5*q2)/1000.;
     double S12 = 16.2*(1.-2.83*q2)*exp(+2.0*q2)/1000.;
     fAmpl.fMinus1 =  TMath::Sqrt((Mnuc/MR)*kgcm0/(2.*kPi*kAem))*A12;
     fAmpl.f0Minus =  TMath::Sqrt((Mnuc/MR)*kgcm0/(2.*kPi*kAem))*(-q2/(qvec*qvec))*S12;
     double A32=0;
     fAmpl.fPlus1  = -1. * fAmpl.fMinus1;
     fAmpl.f0Plus  = -1. * fAmpl.f0Minus;
     fAmpl.fMinus3 =  TMath::Sqrt((Mnuc/MR)*kgcm0/(2.*kPi*kAem))*A32;
     fAmpl.fPlus3  = -1. * fAmpl.fMinus3;



     break;
   }
   case (kD33_1700) :
   {
     fAmpl.fMinus1 =  kSqrt3_2 * fkr.T + k1_Sqrt3 * fkr.Lamda * fkr.R;
     fAmpl.fMinus3 =  k3_Sqrt2 * fkr.T;
     fAmpl.f0Minus = -kSqrt3 * fkr.Lamda * fkr.S;
     fAmpl.fPlus1  =  fAmpl.fMinus1;
     fAmpl.fPlus3  =  fAmpl.fMinus3;
     fAmpl.f0Plus  =  fAmpl.f0Minus;
     //change to MAID parameterization of vector amplitude
     double A12 = 226.*(1.-1.91*q2)*exp(+1.77*q2)/1000.;
     double A32 =  210.*(1.-1.97*q2)*exp(+2.2*q2)/1000.;
     double S12 = 0.;
     A12=0;
     A32=0; //set to 0 to get agreement with (e,e')
     fAmpl.fMinus1 =  TMath::Sqrt((Mnuc/MR)*kgcm0/(2.*kPi*kAem))*A12;
     fAmpl.fPlus1  = -fAmpl.fMinus1;
     fAmpl.fMinus3 =  TMath::Sqrt((Mnuc/MR)*kgcm0/(2.*kPi*kAem))*A32;
     fAmpl.fPlus3  = -fAmpl.fMinus3;
     fAmpl.f0Minus =  TMath::Sqrt((Mnuc/MR)*kgcm0/(2.*kPi*kAem))*(-q2/(qvec*qvec))*S12;
     fAmpl.f0Plus  = -fAmpl.f0Minus;

     //------------------------------------------------------------
     break;
   }
   case (kP11_1440) :
   {
     double L2  = TMath::Power(fkr.Lamda, 2);

     fAmpl.fMinus1 = -0.5*kSqrt3 * L2 * fkr.R;
     fAmpl.fPlus1  =  fAmpl.fMinus1;
     fAmpl.fMinus3 =  0.;
     fAmpl.fPlus3  =  0.;
     fAmpl.f0Minus = -0.5*kSqrt3 * L2 * fkr.S;
     fAmpl.f0Plus  =  fAmpl.f0Minus;
     //change to MAID parameterization of vector amplitude
     double A12 =  -61.4*(1.+1.22*q2-0.55*q2*q2*q2*q2)*exp(+1.51*q2)/1000.;
     double S12 = 4.2*(1.-40*q2+1.5*q2*q2*q2*q2)*exp(+1.75*q2)/1000.;
     fAmpl.fMinus1 =  TMath::Sqrt((Mnuc/MR)*kgcm0/(2.*kPi*kAem))*A12;
     fAmpl.f0Minus =  -TMath::Sqrt((Mnuc/MR)*kgcm0/(2.*kPi*kAem))*(-q2/(qvec*qvec))*S12;
     fAmpl.fPlus1  = -1. * fAmpl.fMinus1;
     fAmpl.f0Plus  = -1. * fAmpl.f0Minus;
     fAmpl.fMinus3 =  0.;
     fAmpl.fPlus3  =  0.;

 


     //----------------------------------------------------------
     break;
   }
   case (kP33_1600) :
   {
     double L2R  = TMath::Power(fkr.Lamda, 2) * fkr.R;

     fAmpl.fMinus1 = k1_Sqrt6 * L2R;
     fAmpl.fMinus3 = k1_Sqrt2 * L2R;
     fAmpl.fPlus1  = -1. * fAmpl.fMinus1;
     fAmpl.fPlus3  = -1. * fAmpl.fMinus3;
     fAmpl.f0Minus = 0.;
     fAmpl.f0Plus  = 0.;
     //MAID
     fAmpl.fMinus3 =  0.;
     fAmpl.fPlus3  =  0.;
     fAmpl.fMinus1 =  0.;
     fAmpl.fPlus1  =  0.;
     fAmpl.f0Minus =  0.;
     fAmpl.f0Plus  =  0.;


     break;
   }
   case (kP13_1720) :
   {
     double L2  = TMath::Power(fkr.Lamda, 2);
     double LT  = fkr.Lamda * fkr.T;

     fAmpl.fMinus1 = -kSqrt27_10 * LT - kSqrt3_5 * L2 * fkr.R;
     fAmpl.fMinus3 =  k3_Sqrt10 * LT;
     fAmpl.f0Minus =  kSqrt3_5  * L2 * fkr.S;
     fAmpl.fPlus1  = -1. * fAmpl.fMinus1;
     fAmpl.fPlus3  = -1. * fAmpl.fMinus3;
     fAmpl.f0Plus  = -1. * fAmpl.f0Minus;
      //change to MAID parameterization of vector amplitude
     double A12 =  73.*(1.-1.89*q2)*exp(+1.55*q2)/1000.;
     double A32 = -11.*(1.-16.*q2)*exp(+1.55*q2)/1000.;
     double S12 = -53.*(1.-2.46*q2)*exp(+1.55*q2)/1000.;
     
     fAmpl.fMinus1 =  TMath::Sqrt((Mnuc/MR)*kgcm0/(2.*kPi*kAem))*A12;
     fAmpl.f0Minus =  TMath::Sqrt((Mnuc/MR)*kgcm0/(2.*kPi*kAem))*(-q2/(qvec*qvec))*S12;
     fAmpl.fPlus1  = -1. * fAmpl.fMinus1;
     fAmpl.f0Plus  = 1. * fAmpl.f0Minus;
     fAmpl.fMinus3 =  TMath::Sqrt((Mnuc/MR)*kgcm0/(2.*kPi*kAem))*A32;
     fAmpl.fPlus3  = -1. * fAmpl.fMinus3;


     //-----------------------------------------------------------------
     break;
   }
   case (kF15_1680) :
   {
     double L2  = TMath::Power(fkr.Lamda, 2);
     double LT  = fkr.Lamda * fkr.T;

     fAmpl.fMinus1 =  -k3_Sqrt5  * LT + k3_Sqrt10 * L2 * fkr.R;
     fAmpl.fMinus3 =  -kSqrt18_5 * LT;
     fAmpl.f0Minus =   k3_Sqrt10 * L2 * fkr.S;
     fAmpl.fPlus1  =  fAmpl.fMinus1;
     fAmpl.fPlus3  =  fAmpl.fMinus3;
     fAmpl.f0Plus  =  fAmpl.f0Minus;
     //change to MAID parameterization of vector amplitude
     double A12 =  -25.*(1.-3.98*q2)*exp(+1.20*q2)/1000.;
     double A32 =  134.*(1.-1.*q2)*exp(+2.22*q2)/1000.;
     double S12 =  -44.*(1.-3.14*q2)*exp(+1.68*q2)/1000.;
     
     fAmpl.fMinus1 =  -TMath::Sqrt((Mnuc/MR)*kgcm0/(2.*kPi*kAem))*A12;
     fAmpl.f0Minus =  -TMath::Sqrt((Mnuc/MR)*kgcm0/(2.*kPi*kAem))*(-q2/(qvec*qvec))*S12;
     fAmpl.fPlus1  = 1. * fAmpl.fMinus1;
     fAmpl.f0Plus  = 1. * fAmpl.f0Minus;
     fAmpl.fMinus3 =  -TMath::Sqrt((Mnuc/MR)*kgcm0/(2.*kPi*kAem))*A32;
     fAmpl.fPlus3  = 1. * fAmpl.fMinus3;

 



     break;
   }
   case (kP31_1910) :
   {
     fAmpl.fMinus1 = -k1_Sqrt15 * TMath::Power(fkr.Lamda, 2) * fkr.R;
     fAmpl.fPlus1  = fAmpl.fMinus1;
     fAmpl.fMinus3 = 0.;
     fAmpl.fPlus3  = 0.;
     fAmpl.f0Minus = 0.;
     fAmpl.f0Plus  = 0.;
     //MAID
     fAmpl.fMinus3 =  0.;
     fAmpl.fPlus3  =  0.;
     fAmpl.fMinus1 =  0.;
     fAmpl.fPlus1  =  0.;
     fAmpl.f0Minus =  0.;
     fAmpl.f0Plus  =  0.;

     break;
   }
   case (kP33_1920) :
   {
     double L2R  = TMath::Power(fkr.Lamda, 2) * fkr.R;

     fAmpl.fMinus1 =  k1_Sqrt15 * L2R;
     fAmpl.fMinus3 = -k1_Sqrt5  * L2R;
     fAmpl.fPlus1  = -1.* fAmpl.fMinus1;
     fAmpl.fPlus3  = -1.* fAmpl.fMinus3;
     fAmpl.f0Minus =  0.;
     fAmpl.f0Plus  =  0.;
     break;
   }
   case (kF35_1905) :
   {
     double L2R  = TMath::Power(fkr.Lamda, 2) * fkr.R;

     fAmpl.fMinus1 = k1_Sqrt35  * L2R;
     fAmpl.fMinus3 = kSqrt18_35 * L2R;
     fAmpl.fPlus1  = fAmpl.fMinus1;
     fAmpl.fPlus3  = fAmpl.fMinus3;
     fAmpl.f0Minus = 0.;
     fAmpl.f0Plus  = 0.;
     //MAID
     fAmpl.fMinus3 =  0.;
     fAmpl.fPlus3  =  0.;
     fAmpl.fMinus1 =  0.;
     fAmpl.fPlus1  =  0.;
     fAmpl.f0Minus =  0.;
     fAmpl.f0Plus  =  0.;

     break;
   }
   case (kF37_1950) :
   {
     double L2R  = TMath::Power(fkr.Lamda, 2) * fkr.R;

     fAmpl.fMinus1 = -kSqrt6_35 * L2R;
     fAmpl.fMinus3 = -kSqrt2_7  * L2R;
     fAmpl.fPlus1  = -1. * fAmpl.fMinus1;
     fAmpl.fPlus3  = -1. * fAmpl.fMinus3;
     fAmpl.f0Minus = 0.;
     fAmpl.f0Plus  = 0.;
     //MAID
     fAmpl.fMinus3 =  0.;
     fAmpl.fPlus3  =  0.;
     fAmpl.fMinus1 =  0.;
     fAmpl.fPlus1  =  0.;
     fAmpl.f0Minus =  0.;
     fAmpl.f0Plus  =  0.;


     break;
   }
   case (kP11_1710) :
   {
     double L2  = TMath::Power(fkr.Lamda, 2);

     fAmpl.fMinus1 = kSqrt3_8 * L2 * fkr.R;
     fAmpl.f0Minus = kSqrt3_8 * L2 * fkr.S;
     fAmpl.fPlus1  = fAmpl.fMinus1;
     fAmpl.f0Plus  = fAmpl.f0Minus;
     fAmpl.fMinus3 = 0.;
     fAmpl.fPlus3  = 0.;
     //MAID
     fAmpl.fMinus3 =  0.;
     fAmpl.fPlus3  =  0.;
     fAmpl.fMinus1 =  0.;
     fAmpl.fPlus1  =  0.;
     fAmpl.f0Minus =  0.;
     fAmpl.f0Plus  =  0.;


     break;
   }
   case (kF17_1970) :
   {
     fAmpl.fMinus1 = 0.;
     fAmpl.fPlus1  = 0.;
     fAmpl.fMinus3 = 0.;
     fAmpl.fPlus3  = 0.;
     fAmpl.f0Minus = 0.;
     fAmpl.f0Plus  = 0.;
     //MAID
     fAmpl.fMinus3 =  0.;
     fAmpl.fPlus3  =  0.;
     fAmpl.fMinus1 =  0.;
     fAmpl.fPlus1  =  0.;
     fAmpl.f0Minus =  0.;
     fAmpl.f0Plus  =  0.;


     break;
   }
   default:
   {
     LOG("RSHAmpl", pWARN) << "*** UNRECOGNIZED RESONANCE!";
     fAmpl.fMinus1 = 0.;
     fAmpl.fPlus1  = 0.;
     fAmpl.fMinus3 = 0.;
     fAmpl.fPlus3  = 0.;
     fAmpl.f0Minus = 0.;
     fAmpl.f0Plus  = 0.;
     break;
   }

  }//switch

  return fAmpl;
}
//____________________________________________________________________________


