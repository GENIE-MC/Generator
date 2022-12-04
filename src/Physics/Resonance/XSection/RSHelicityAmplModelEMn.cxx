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

#include "Algorithm/AlgConfigPool.h"
#include "BaryonResonance/BaryonResUtils.h"
#include "Conventions/Constants.h"
#include "ReinSehgal/RSHelicityAmplModelEMn.h"
#include "ReinSehgal/RSHelicityAmpl.h"
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
const RSHelicityAmpl & 
   RSHelicityAmplModelEMn::Compute(
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
     
     //here MAID multipoles are used
     double Fq = 1./TMath::Power(1-q2/0.71,2)*qcm/kgcm0;
     double AM = 300.*(1. - 0.01*q2)*TMath::Exp(0.23*q2)*Fq;
     double AE = -6.37 * (1. + 0.021*q2)*TMath::Exp(0.16*q2)*Fq;
     double AS = -12.40 * (1. - 0.12*q2) / (1. + 4.9*tau)*qcm/kR*TMath::Exp(0.23*q2)*Fq;
     double A12= -(3.*AE+AM)/2./1000.;
     double A32= TMath::Sqrt(3.)/2.*(AE-AM)/1000.;
     double S12 = TMath::Sqrt(2.)*AS/1000.;
     LOG("ReinSeghalRes", pWARN) << "q2= " << q2 << " A12, A32, S12 =" <<A12<<" "<<A32<<"  "<<S12;

     fAmpl.fMinus1 =  TMath::Sqrt(115/130)*TMath::Sqrt((Mnuc/MR)*kgcm0/(2.*kPi*kAem))*A12;
     fAmpl.fPlus1  = -fAmpl.fMinus1;
     fAmpl.fMinus3 =  TMath::Sqrt(115/130)*(Mnuc/MR)*TMath::Sqrt((Mnuc/MR)*kgcm0/(2.*kPi*kAem))*A32;
     fAmpl.fPlus3  = -fAmpl.fMinus3;
     fAmpl.f0Minus =  TMath::Sqrt(115/130)*(Mnuc/MR)*TMath::Sqrt((Mnuc/MR)*kgcm0/(2.*kPi*kAem))*(-q2/(qvec*qvec))*S12;
     fAmpl.f0Plus  = fAmpl.f0Minus;
     LOG("ReinSeghalRes", pWARN) << "MAID fminus1= " << fAmpl.fMinus1 << " fminus3 =" << fAmpl.fMinus3 ;    

     break;
   }
   case (kS11_1535) :
   {
     fAmpl.fPlus1  =  kSqrt3   * fkr.T + k1_Sqrt6 * fkr.Lamda * fkr.R;
     fAmpl.f0Minus =  kSqrt3_2 * fkr.Lamda * fkr.S;
     fAmpl.fMinus1 = -1 * fAmpl.fPlus1;
     fAmpl.f0Plus  = -1 * fAmpl.f0Minus;
     fAmpl.fMinus3 =  0.;
     fAmpl.fPlus3  =  0.;
     
     // change to MAID parameterization of vector amplitude
     double A12= -51.*(1.-4.75*q2)*exp(1.69*q2)/1000.;
     double S12=  28.5*(1.-0.36*q2)*exp(1.55*q2)/1000.;
     //double A32=  0.;
     fAmpl.fMinus1 =  TMath::Sqrt(kgcm0/(2.*kPi*kAem))*A12;
     fAmpl.f0Minus =  TMath::Sqrt(kgcm0/(2.*kPi*kAem))*(-q2/(qvec*qvec))*S12;
     fAmpl.fPlus1  = -1. * fAmpl.fMinus1;
     fAmpl.f0Plus  = -1. * fAmpl.f0Minus;
     fAmpl.fMinus3 =  0.;
     fAmpl.fPlus3  =  0.;
     fAmpl.fMinus3 =  0.;
     fAmpl.fPlus3  =  0.;
     


     //-------------------------------------------------
     break;
   }
   case (kD13_1520) :
   {
     fAmpl.fMinus1 = -kSqrt3_2 * fkr.T + k1_Sqrt3 * fkr.Lamda * fkr.R;
     fAmpl.fMinus3 = -k3_Sqrt2 * fkr.T;
     fAmpl.f0Minus =  kSqrt3 * fkr.Lamda * fkr.S;
     fAmpl.fPlus1  =  fAmpl.fMinus1;
     fAmpl.fPlus3  =  fAmpl.fMinus3;
     fAmpl.f0Plus  =  fAmpl.f0Minus;
     //change to MAID parameterizations of vector amplitude
     double A12= -77.*(1.+0.53*q2)*exp(1.55*q2)/1000.;
     double A32= -154.*(1.-0.58*q2)*exp(1.75*q2)/1000.;
     double S12=  13.6*(1.-15.7*q2)*exp(1.57*q2)/1000.;
     fAmpl.fMinus1 =  TMath::Sqrt(kgcm0/(2.*kPi*kAem))*A12;
     fAmpl.fPlus1  = -fAmpl.fMinus1;
     fAmpl.fMinus3 =  TMath::Sqrt(kgcm0/(2.*kPi*kAem))*A32;
     fAmpl.fPlus3  = -fAmpl.fMinus3;
     fAmpl.f0Minus =  TMath::Sqrt(kgcm0/(2.*kPi*kAem))*(-q2/(qvec*qvec))*S12;
     fAmpl.f0Plus  = fAmpl.f0Minus;

     
     //----------------------------------------------------
     break;
   }
   case (kS11_1650) :
   {
     fAmpl.fPlus1  =  k1_Sqrt6 * fkr.Lamda * fkr.R;
     fAmpl.fMinus1 = -1 * fAmpl.fPlus1;
     fAmpl.fMinus3 =  0.;
     fAmpl.fPlus3  =  0.;
     fAmpl.f0Minus =  0.;
     fAmpl.f0Plus  =  0.;
     //change to MAID parameterizations of vector amplitude
     double A12=  9.*(1.-0.13*q2)*exp(1.55*q2)/1000.;
     double S12=  10.1*(1.+0.50*q2)*exp(1.55*q2)/1000.;
     //double A32=  0.;
     fAmpl.fMinus1 =  TMath::Sqrt(kgcm0/(2.*kPi*kAem))*A12;
     fAmpl.f0Minus =  TMath::Sqrt(kgcm0/(2.*kPi*kAem))*(-q2/(qvec*qvec))*S12;
     fAmpl.fPlus1  = -1. * fAmpl.fMinus1;
     fAmpl.f0Plus  = -1. * fAmpl.f0Minus;
     fAmpl.fMinus3 =  0.;
     fAmpl.fPlus3  =  0.;



     //---------------------------------------------------
     break;
   }
   case (kD13_1700) :
   {
     double LR = fkr.Lamda * fkr.R;

     fAmpl.fMinus1 = -(1./kSqrt30) * LR;
     fAmpl.fMinus3 = -(3./kSqrt10) * LR;
     fAmpl.fPlus1  =  fAmpl.fMinus1;
     fAmpl.fPlus3  =  fAmpl.fMinus3;
     fAmpl.f0Minus =  0.;
     fAmpl.f0Plus  =  0.;
     break;
   }
   case (kD15_1675) :
   {
     double LR = fkr.Lamda * fkr.R;

     fAmpl.fMinus1 = kSqrt3_10 * LR;
     fAmpl.fMinus3 = kSqrt3_5  * LR;
     fAmpl.fPlus1  = -1 * fAmpl.fMinus1;
     fAmpl.fPlus3  = -1 * fAmpl.fMinus3;
     fAmpl.f0Minus =  0.;
     fAmpl.f0Plus  =  0.;
     //change to MAID parameterizations of vector amplitude
     double A12= -62.*(1.-0.01*q2)*exp(+2.00*q2)/1000.;
     double A32= -84.*(1.-0.01*q2)*exp(+2.00*q2)/1000.;
     double S12= 0.;
     fAmpl.fMinus1 =  TMath::Sqrt(kgcm0/(2.*kPi*kAem))*A12;
     fAmpl.fPlus1  = -fAmpl.fMinus1;
     fAmpl.fMinus3 =  TMath::Sqrt(kgcm0/(2.*kPi*kAem))*A32;
     fAmpl.fPlus3  = -fAmpl.fMinus3;
     fAmpl.f0Minus =  TMath::Sqrt(kgcm0/(2.*kPi*kAem))*(-q2/(qvec*qvec))*S12;
     fAmpl.f0Plus  = fAmpl.f0Minus;

 

     //-----------------------------------------------------------------
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
     fAmpl.fMinus1 =  TMath::Sqrt((Mnuc/MR)*kgcm0/(2.*kPi*kAem))*A12;
     fAmpl.fPlus1  = -fAmpl.fMinus1;
     fAmpl.fMinus3 =  TMath::Sqrt((Mnuc/MR)*kgcm0/(2.*kPi*kAem))*A32;
     fAmpl.fPlus3  = -fAmpl.fMinus3;
     fAmpl.f0Minus =  TMath::Sqrt((Mnuc/MR)*kgcm0/(2.*kPi*kAem))*(-q2/(qvec*qvec))*S12;
     fAmpl.f0Plus  = -fAmpl.f0Minus;

     //------------------------------------------------------------
     //



     break;
   }
   case (kP11_1440) :
   {
     fAmpl.fMinus1 = k1_Sqrt3 * TMath::Power(fkr.Lamda, 2) * fkr.R;
     fAmpl.fPlus1  = fAmpl.fMinus1;
     fAmpl.fMinus3 =  0.;
     fAmpl.fPlus3  =  0.;
     fAmpl.f0Minus =  0.;
     fAmpl.f0Plus  =  0.;
     
     //change to MAID parameterization of vector amplitude
     double A12 =  54.1*(1.-0.95*q2)*exp(+1.77*q2)/1000.;
     double S12 = -41.5*(1.-2.98*q2)*exp(+1.55*q2)/1000.;
     double A32=0;
     fAmpl.fMinus1 =  TMath::Sqrt((Mnuc/MR)*kgcm0/(2.*kPi*kAem))*A12;
     fAmpl.f0Minus =  -TMath::Sqrt((Mnuc/MR)*kgcm0/(2.*kPi*kAem))*(-q2/(qvec*qvec))*S12;
     fAmpl.fPlus1  = 1. * fAmpl.fMinus1;
     fAmpl.f0Plus  = 1. * fAmpl.f0Minus;
     fAmpl.fMinus3 =  0.;
     fAmpl.fPlus3  =  0.;
     
      

     //-------------------------------------------------------------
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
     break;
   }
   case (kP13_1720) :
   {
     fAmpl.fMinus1 = k2_Sqrt15 * TMath::Power(fkr.Lamda, 2) * fkr.R;
     fAmpl.fPlus1  = -1 * fAmpl.fMinus1;
     fAmpl.fMinus3 =  0.;
     fAmpl.fPlus3  =  0.;
     fAmpl.f0Minus =  0.;
     fAmpl.f0Plus  =  0.;
     //change to MAID parameterization of vector amplitude
     double A12 =  -3*(1.-12.7*q2)*exp(+1.55*q2)/1000.;
     double A32 = -31*(1.-4.99*q2)*exp(+1.55*q2)/1000.;
     double S12 = 0.;
     fAmpl.fMinus1 =  -TMath::Sqrt((Mnuc/MR)*kgcm0/(2.*kPi*kAem))*A12;
     fAmpl.f0Minus =  -TMath::Sqrt((Mnuc/MR)*kgcm0/(2.*kPi*kAem))*(-q2/(qvec*qvec))*S12;
     fAmpl.fPlus1  = -1. * fAmpl.fMinus1;
     fAmpl.f0Plus  = 1. * fAmpl.f0Minus;
     fAmpl.fMinus3 =  TMath::Sqrt((Mnuc/MR)*kgcm0/(2.*kPi*kAem))*A32;
     fAmpl.fPlus3  = -1. * fAmpl.fMinus3;

     //----------------------------------------------------------------      
     break;
   }
   case (kF15_1680) :
   {
     fAmpl.fMinus1 =  -kSqrt2_5 * TMath::Power(fkr.Lamda, 2) * fkr.R;
     fAmpl.fPlus1  =  fAmpl.fMinus1;
     fAmpl.fMinus3 =  0.;
     fAmpl.fPlus3  =  0.;
     fAmpl.f0Minus =  0.;
     fAmpl.f0Plus  =  0.;
     //change to MAID parameterization of vector amplitude
     double A12 =  28*(1.-0.*q2)*exp(+1.2*q2)/1000.;
     double A32 = -38*(1.-4.09*q2)*exp(+1.75*q2)/1000.;
     double S12 = 0.;
     fAmpl.fMinus1 =  -TMath::Sqrt((Mnuc/MR)*kgcm0/(2.*kPi*kAem))*A12;
     fAmpl.f0Minus =  -TMath::Sqrt((Mnuc/MR)*kgcm0/(2.*kPi*kAem))*(-q2/(qvec*qvec))*S12;
     fAmpl.fPlus1  = 1. * fAmpl.fMinus1;
     fAmpl.f0Plus  = 1. * fAmpl.f0Minus;
     fAmpl.fMinus3 =  -TMath::Sqrt((Mnuc/MR)*kgcm0/(2.*kPi*kAem))*A32;
     fAmpl.fPlus3  = 1. * fAmpl.fMinus3;



     //--------------------------------------------------------
     break;
   }
   case (kP31_1910) :
   {
     fAmpl.fMinus1 =  -k1_Sqrt15 * TMath::Power(fkr.Lamda, 2) * fkr.R;
     fAmpl.fPlus1  =  fAmpl.fMinus1;
     fAmpl.fMinus3 =  0.;
     fAmpl.fPlus3  =  0.;
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
     //MAID
     fAmpl.fMinus3 =  0.;
     fAmpl.fPlus3  =  0.;
     fAmpl.fMinus1 =  0.;
     fAmpl.fPlus1  =  0.;
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
     double L2 = TMath::Power(fkr.Lamda, 2);

     fAmpl.fMinus1 = -k1_Sqrt24 * L2 * fkr.R;
     fAmpl.f0Minus = -kSqrt3_8  * L2 * fkr.S;
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
     double L2R = TMath::Power(fkr.Lamda, 2) * fkr.R;

     fAmpl.fMinus1 = kSqrt3_35 * L2R;
     fAmpl.fPlus1  = -1 * fAmpl.fMinus1;
     fAmpl.fMinus3 = k1_Sqrt7  * L2R;
     fAmpl.fPlus3  = -1 * fAmpl.fMinus3;
     fAmpl.f0Minus =  0.;
     fAmpl.f0Plus  =  0.;
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
