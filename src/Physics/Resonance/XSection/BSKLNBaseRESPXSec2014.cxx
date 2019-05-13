//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

          Afroditi Papadopoulou <apapadop \at mit.edu>
          Massachusetts Institute of Technology

          Adi Ashkenazi <adishka \at gmail.com>
          Massachusetts Institute of Technology

 @ July 4, 2018 - Afroditi Papadopoulou
   For electromagnetic (EM) interactions, the weak g2 was still used for the
   calculation of the helicity amplitude. Fixed by replacing with the correct EM g2

*/
//____________________________________________________________________________

#include <TMath.h>
#include <TSystem.h>

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/Spline.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/Utils/Range1.h"
#include "Framework/Utils/BWFunc.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Physics/Resonance/XSection/BSKLNBaseRESPXSec2014.h"
#include "Physics/Resonance/XSection/RSHelicityAmplModelI.h"
#include "Physics/Resonance/XSection/RSHelicityAmpl.h"
#include "Physics/NuclearState/FermiMomentumTablePool.h"
#include "Physics/NuclearState/FermiMomentumTable.h"
#include "Physics/NuclearState/NuclearUtils.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
BSKLNBaseRESPXSec2014::BSKLNBaseRESPXSec2014(string name) :
XSecAlgorithmI(name)
{

}
//____________________________________________________________________________
BSKLNBaseRESPXSec2014::BSKLNBaseRESPXSec2014(string name, string config) :
XSecAlgorithmI(name, config)
{

}
//____________________________________________________________________________
BSKLNBaseRESPXSec2014::~BSKLNBaseRESPXSec2014()
{

}
//____________________________________________________________________________
double BSKLNBaseRESPXSec2014::XSec(
    const Interaction * interaction, KinePhaseSpace_t kps) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  const InitialState & init_state = interaction -> InitState();
  const ProcessInfo &  proc_info  = interaction -> ProcInfo();
  const Target & target = init_state.Tgt();

  // Get kinematical parameters
  const Kinematics & kinematics = interaction -> Kine();
  double W  = kinematics.W();
  double q2 = kinematics.q2();
  double costh = kinematics.FSLeptonP4().CosTheta();

  // Under the DIS/RES joining scheme, xsec(RES)=0 for W>=Wcut
  if(fUsingDisResJoin) {
    if(W>=fWcut) {
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
      LOG("BSKLNBaseRESPXSec2014", pDEBUG)
        << "RES/DIS Join Scheme: XSec[RES, W=" << W
        << " >= Wcut=" << fWcut << "] = 0";
#endif
      return 0;
    }
  }

  // Get the input baryon resonance
  Resonance_t resonance = interaction->ExclTag().Resonance();
  string      resname   = utils::res::AsString(resonance);
  bool        is_delta  = utils::res::IsDelta (resonance);

  // Get the neutrino, hit nucleon & weak current
  int  nucpdgc   = target.HitNucPdg();
  int  probepdgc = init_state.ProbePdg();
  bool is_nu     = pdg::IsNeutrino         (probepdgc);
  bool is_nubar  = pdg::IsAntiNeutrino     (probepdgc);
  bool is_lplus  = pdg::IsPosChargedLepton (probepdgc);
  bool is_lminus = pdg::IsNegChargedLepton (probepdgc);
  bool is_p      = pdg::IsProton  (nucpdgc);
  bool is_n      = pdg::IsNeutron (nucpdgc);
  bool is_CC     = proc_info.IsWeakCC();
  bool is_NC     = proc_info.IsWeakNC();
  bool is_EM     = proc_info.IsEM();

  //  bool new_GV = fGA; //JN
  //  bool new_GA = fGV; //JN


  if(is_CC && !is_delta) {
    if((is_nu && is_p) || (is_nubar && is_n)) return 0;
  }

  // Get baryon resonance parameters
  int    IR  = utils::res::ResonanceIndex    (resonance);
  int    LR  = utils::res::OrbitalAngularMom (resonance);
  double MR  = utils::res::Mass              (resonance);
  double WR  = utils::res::Width             (resonance);
   double NR  = fNormBW?utils::res::BWNorm    (resonance,fN0ResMaxNWidths,fN2ResMaxNWidths,fGnResMaxNWidths):1;

  // Following NeuGEN, avoid problems with underlying unphysical
  // model assumptions by restricting the allowed W phase space
  // around the resonance peak
 if (fNormBW) {
        if      (W > MR + fN0ResMaxNWidths * WR && IR==0) return 0.;
        else if (W > MR + fN2ResMaxNWidths * WR && IR==2) return 0.;
        else if (W > MR + fGnResMaxNWidths * WR)          return 0.;
  }

  // Compute auxiliary & kinematical factors
  double E      = init_state.ProbeE(kRfHitNucRest);
  double Mnuc   = target.HitNucMass();
  double W2     = TMath::Power(W,    2);
  double Mnuc2  = TMath::Power(Mnuc, 2);
  double k      = 0.5 * (W2 - Mnuc2)/Mnuc;
  double v      = k - 0.5 * q2/Mnuc;
  double v2     = TMath::Power(v, 2);
  double Q2     = v2 - q2;
  double Q      = TMath::Sqrt(Q2);
  double Eprime = E - v;
  double U      = 0.5 * (E + Eprime + Q) / E;
  double V      = 0.5 * (E + Eprime - Q) / E;
  double U2     = TMath::Power(U, 2);
  double V2     = TMath::Power(V, 2);
  double UV     = U*V;


  //JN parameter from the KUZMIN et al.

  //  bool is_RS  = true;
  bool is_KLN = false;
  if(fKLN && is_CC) is_KLN=true;

  bool is_BRS = false;
  if(fBRS && is_CC) is_BRS=true;

  double ml    = interaction->FSPrimLepton()->Mass();
  double Pl    = TMath::Sqrt(Eprime*Eprime - ml*ml);

  double vstar = (Mnuc*v + q2)/W;  //missing W
  double Qstar = TMath::Sqrt(-q2 + vstar*vstar);
  double sqrtq2 = TMath::Sqrt(-q2);
  double a = 1. + 0.5*(W2-q2+Mnuc2)/Mnuc/W;

  double KNL_Alambda_plus  = 0;
  double KNL_Alambda_minus = 0;
  double KNL_j0_plus  = 0;
  double KNL_j0_minus = 0;
  double KNL_jx_plus  = 0;
  double KNL_jx_minus = 0;
  double KNL_jy_plus  = 0;
  double KNL_jy_minus = 0;
  double KNL_jz_plus  = 0;
  double KNL_jz_minus = 0;
  double KNL_Qstar_plus =0;
  double KNL_Qstar_minus =0;

  double KNL_K = Q/E/TMath::Sqrt(2*(-q2));

  double KNL_cL_plus  = 0;
  double KNL_cL_minus = 0;

  double KNL_cR_plus  = 0;
  double KNL_cR_minus = 0;

  double KNL_cS_plus  = 0;
  double KNL_cS_minus = 0;

  double KNL_vstar_plus  = 0;
  double KNL_vstar_minus = 0;

  if(is_CC && (is_KLN || is_BRS)){

    LOG("BSKLNBaseRESPXSec2014",pINFO) "costh1="<<costh;
    costh = (q2 - ml*ml + 2.*E*Eprime)/2./E/Pl;
    //ml=0;
    LOG("BSKLNBaseRESPXSec2014",pINFO) "q2="<<q2<< "m2="<<ml*ml<<" 2.*E*Eprime="<<2.*E*Eprime<<" nom="<< (q2 - ml*ml + 2.*E*Eprime)<<" den="<<2.*E*Pl;
    LOG("BSKLNBaseRESPXSec2014",pINFO) "costh2="<<costh;

    KNL_Alambda_plus  = TMath::Sqrt(E*(Eprime - Pl));
    KNL_Alambda_minus = TMath::Sqrt(E*(Eprime + Pl));
    LOG("BSKLNBaseRESPXSec2014",pINFO)
       << "\n+++++++++++++++++++++++ \n"
       << "E="<<E << " K= "<<KNL_K << "\n"
       << "El="<<Eprime<<" Pl="<<Pl<<" ml="<<ml << "\n"
       << "W="<<W<<" Q="<<Q<<" q2="<<q2 << "\n"
       << "A-="<<KNL_Alambda_minus<<" A+="<<KNL_Alambda_plus << "\n"
       << "xxxxxxxxxxxxxxxxxxxxxxx";

    KNL_j0_plus  = KNL_Alambda_plus /W * TMath::Sqrt(1 - costh) * (Mnuc - Eprime - Pl);
    KNL_j0_minus = KNL_Alambda_minus/W * TMath::Sqrt(1 + costh) * (Mnuc - Eprime + Pl);

    KNL_jx_plus  = KNL_Alambda_plus/ Q * TMath::Sqrt(1 + costh) * (Pl - E);
    KNL_jx_minus = KNL_Alambda_minus/Q * TMath::Sqrt(1 - costh) * (Pl + E);

    KNL_jy_plus  =  KNL_Alambda_plus  * TMath::Sqrt(1 + costh);
    KNL_jy_minus = -KNL_Alambda_minus * TMath::Sqrt(1 - costh);

    KNL_jz_plus  = KNL_Alambda_plus /W/Q *  TMath::Sqrt(1 - costh) * ( (E + Pl)*(Mnuc -Eprime) + Pl*(  E + 2*E*costh -Pl) );
    KNL_jz_minus = KNL_Alambda_minus/W/Q *  TMath::Sqrt(1 + costh) * ( (E - Pl)*(Mnuc -Eprime) + Pl*( -E + 2*E*costh -Pl) );

    if (is_nu || is_lminus) {
      KNL_Qstar_plus  = sqrtq2 * KNL_j0_plus  / TMath::Sqrt(TMath::Abs(KNL_j0_plus*KNL_j0_plus - KNL_jz_plus*KNL_jz_plus) );
      KNL_Qstar_minus = sqrtq2 * KNL_j0_minus / TMath::Sqrt(TMath::Abs(KNL_j0_minus*KNL_j0_minus - KNL_jz_minus*KNL_jz_minus) );
    }

    else if (is_nubar || is_lplus){
      KNL_Qstar_plus  = sqrtq2 * KNL_j0_minus / TMath::Sqrt(TMath::Abs(KNL_j0_minus*KNL_j0_minus - KNL_jz_minus*KNL_jz_minus) );
      KNL_Qstar_minus = sqrtq2 * KNL_j0_plus / TMath::Sqrt(TMath::Abs(KNL_j0_plus*KNL_j0_plus - KNL_jz_plus*KNL_jz_plus) );
    }

    if (is_nu || is_lminus) {
      KNL_vstar_plus  = sqrtq2 * KNL_jz_plus  / TMath::Sqrt(TMath::Abs(KNL_j0_plus*KNL_j0_plus - KNL_jz_plus*KNL_jz_plus) );
      KNL_vstar_minus = sqrtq2 * KNL_jz_minus / TMath::Sqrt(TMath::Abs(KNL_j0_minus*KNL_j0_minus - KNL_jz_minus*KNL_jz_minus) );
    }
    else if (is_nubar || is_lplus) {
      KNL_vstar_minus  = sqrtq2 * KNL_jz_plus / TMath::Sqrt(TMath::Abs(KNL_j0_plus*KNL_j0_plus - KNL_jz_plus*KNL_jz_plus) );
      KNL_vstar_plus   = sqrtq2 * KNL_jz_minus / TMath::Sqrt(TMath::Abs(KNL_j0_minus*KNL_j0_minus - KNL_jz_minus*KNL_jz_minus) );
    }

    if(is_nu || is_lminus){
      KNL_cL_plus  = TMath::Sqrt(0.5)* KNL_K * (KNL_jx_plus  - KNL_jy_plus);
      KNL_cL_minus = TMath::Sqrt(0.5)* KNL_K * (KNL_jx_minus - KNL_jy_minus);

      KNL_cR_plus  = TMath::Sqrt(0.5)* KNL_K * (KNL_jx_plus  + KNL_jy_plus);
      KNL_cR_minus = TMath::Sqrt(0.5)* KNL_K * (KNL_jx_minus + KNL_jy_minus);

      KNL_cS_plus   = KNL_K *  TMath::Sqrt(TMath::Abs(KNL_j0_plus *KNL_j0_plus  - KNL_jz_plus *KNL_jz_plus ) );
      KNL_cS_minus  = KNL_K *  TMath::Sqrt(TMath::Abs(KNL_j0_minus*KNL_j0_minus - KNL_jz_minus*KNL_jz_minus) );
    }

    if (is_nubar || is_lplus) {
      KNL_cL_plus  =  1 * TMath::Sqrt(0.5)* KNL_K * (KNL_jx_minus + KNL_jy_minus);
      KNL_cL_minus = -1 * TMath::Sqrt(0.5)* KNL_K * (KNL_jx_plus  + KNL_jy_plus);

      KNL_cR_plus  =  1 * TMath::Sqrt(0.5)* KNL_K * (KNL_jx_minus - KNL_jy_minus);
      KNL_cR_minus = -1 * TMath::Sqrt(0.5)* KNL_K * (KNL_jx_plus  - KNL_jy_plus);

      KNL_cS_plus  = -1 * KNL_K *  TMath::Sqrt(TMath::Abs(KNL_j0_minus*KNL_j0_minus - KNL_jz_minus*KNL_jz_minus) );
      KNL_cS_minus =  1 * KNL_K *  TMath::Sqrt(TMath::Abs(KNL_j0_plus*KNL_j0_plus - KNL_jz_plus*KNL_jz_plus) );
    }
  }

     LOG("BSKLNBaseRESPXSec2014",pINFO) <<"j0-="<<KNL_j0_minus<<" j0+="<<KNL_j0_plus;
     LOG("BSKLNBaseRESPXSec2014",pINFO) <<"jx-="<<KNL_jx_minus<<" jx+="<<KNL_jx_plus;
     LOG("BSKLNBaseRESPXSec2014",pINFO) <<"jy-="<<KNL_jy_minus<<" jy+="<<KNL_jy_plus;
     LOG("BSKLNBaseRESPXSec2014",pINFO) <<"jz-="<<KNL_jz_minus<<" jz+="<<KNL_jz_plus;

     LOG("BSKLNBaseRESPXSec2014",pINFO) "sqrt2="<<sqrtq2<<" jz+=:"<<KNL_jz_plus<<" j0+="<<KNL_j0_plus<<" denom="<<TMath::Sqrt(TMath::Abs(KNL_j0_plus*KNL_j0_plus - KNL_jz_plus*KNL_jz_plus) );

     LOG("BSKLNBaseRESPXSec2014",pINFO) <<"vstar-="<<KNL_vstar_minus<<" vstar+="<<KNL_vstar_plus;
     LOG("BSKLNBaseRESPXSec2014",pINFO) <<"Qstar-="<<KNL_Qstar_minus<<" Qstar+="<<KNL_Qstar_plus;

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("BSKLNBaseRESPXSec2014", pDEBUG)
    << "Kinematical params V = " << V << ", U = " << U;
#endif

  // Calculate the Feynman-Kislinger-Ravndall parameters

  double Go  = TMath::Power(1 - 0.25 * q2/Mnuc2, 0.5-IR);
  double GV  = Go * TMath::Power( 1./(1-q2/fMv2), 2);
  double GA  = Go * TMath::Power( 1./(1-q2/fMa2), 2);

  if(fGV){

    LOG("BSKLNBaseRESPXSec2014",pDEBUG) <<"Using new GV";
    double CV0 =  1./(1-q2/fMv2/4.);
    double CV3 =  2.13 * CV0 * TMath::Power( 1-q2/fMv2,-2);
    double CV4 = -1.51 * CV0 * TMath::Power( 1-q2/fMv2,-2);
    double CV5 =  0.48 * CV0 * TMath::Power( 1-q2/fMv2/0.766, -2);

    double GV3 =  0.5 / TMath::Sqrt(3) * ( CV3 * (W + Mnuc)/Mnuc
                  + CV4 * (W2 + q2 -Mnuc2)/2./Mnuc2
                  + CV5 * (W2 - q2 -Mnuc2)/2./Mnuc2 );

    double GV1 = - 0.5 / TMath::Sqrt(3) * ( CV3 * (Mnuc2 -q2 +Mnuc*W)/W/Mnuc
                 + CV4 * (W2 +q2 - Mnuc2)/2./Mnuc2
                 + CV5 * (W2 -q2 - Mnuc2)/2./Mnuc2 );

    GV = 0.5 * TMath::Power( 1 - q2/(Mnuc + W)/(Mnuc + W), 0.5-IR)
         * TMath::Sqrt( 3 * GV3*GV3 + GV1*GV1);
  }

  if(fGA){
    LOG("BSKLNBaseRESPXSec2014",pDEBUG) << "Using new GA";

    double CA5_0 = 1.2;
    double CA5 = CA5_0 *  TMath::Power( 1./(1-q2/fMa2), 2);
    //  GA = 0.5 * TMath::Sqrt(3.) * TMath::Power( 1 - q2/(Mnuc + W)/(Mnuc + W), 0.5-IR) * (1- (W2 +q2 -Mnuc2)/8./Mnuc2) * CA5/fZeta;
    GA = 0.5 * TMath::Sqrt(3.) * TMath::Power( 1 - q2/(Mnuc + W)/(Mnuc + W), 0.5-IR) * (1- (W2 +q2 -Mnuc2)/8./Mnuc2) * CA5;

    LOG("BSKLNBaseRESPXSec2014",pINFO) <<"GA= " <<GA << "  C5A= " <<CA5;
  }
  //JN end of new form factors code

  if(is_EM) {
    GA = 0.; // zero the axial term for EM scattering
  }

  double d      = TMath::Power(W+Mnuc,2.) - q2;
  double sq2omg = TMath::Sqrt(2./fOmega);
  double nomg   = IR * fOmega;
  double mq_w   = Mnuc*Q/W;

  fFKR.Lamda  = sq2omg * mq_w;
  fFKR.Tv     = GV / (3.*W*sq2omg);
  fFKR.Rv     = kSqrt2 * mq_w*(W+Mnuc)*GV / d;
  fFKR.S      = (-q2/Q2) * (3*W*Mnuc + q2 - Mnuc2) * GV / (6*Mnuc2);
  fFKR.Ta     = (2./3.) * (fZeta/sq2omg) * mq_w * GA / d;
  fFKR.Ra     = (kSqrt2/6.) * fZeta * (GA/W) * (W+Mnuc + 2*nomg*W/d );
  fFKR.B      = fZeta/(3.*W*sq2omg) * (1 + (W2-Mnuc2+q2)/ d) * GA;
  fFKR.C      = fZeta/(6.*Q) * (W2 - Mnuc2 + nomg*(W2-Mnuc2+q2)/d) * (GA/Mnuc);
  fFKR.R      = fFKR.Rv;
  fFKR.Rplus  = - (fFKR.Rv + fFKR.Ra);
  fFKR.Rminus = - (fFKR.Rv - fFKR.Ra);
  fFKR.T      = fFKR.Tv;
  fFKR.Tplus  = - (fFKR.Tv + fFKR.Ta);
  fFKR.Tminus = - (fFKR.Tv - fFKR.Ta);

  //JN KNL
  double KNL_S_plus = 0;
  double KNL_S_minus = 0;
  double KNL_B_plus = 0;
  double KNL_B_minus = 0;
  double KNL_C_plus = 0;
  double KNL_C_minus = 0;

  if(is_CC && is_KLN){
    KNL_S_plus  = (KNL_vstar_plus*vstar  - KNL_Qstar_plus *Qstar )* (Mnuc2 -q2 - 3*W*Mnuc ) * GV / (6*Mnuc2)/Q2; //possibly missing minus sign ()
    KNL_S_minus = (KNL_vstar_minus*vstar - KNL_Qstar_minus*Qstar )* (Mnuc2 -q2 - 3*W*Mnuc ) * GV / (6*Mnuc2)/Q2;

    LOG("BSKLNBaseRESPXSec2014",pINFO) <<"KNL S= " <<KNL_S_plus<<"\t"<<KNL_S_minus<<"\t"<<fFKR.S;

    KNL_B_plus  = fZeta/(3.*W*sq2omg)/Qstar * (KNL_Qstar_plus  + KNL_vstar_plus *Qstar/a/Mnuc ) * GA;
    KNL_B_minus = fZeta/(3.*W*sq2omg)/Qstar * (KNL_Qstar_minus + KNL_vstar_minus*Qstar/a/Mnuc ) * GA;
    LOG("BSKLNBaseRESPXSec2014",pINFO) <<"KNL B= " <<KNL_B_plus<<"\t"<<KNL_B_minus<<"\t"<<fFKR.B;

    KNL_C_plus = ( (KNL_Qstar_plus*Qstar - KNL_vstar_plus*vstar ) * ( 1./3. + vstar/a/Mnuc)
        + KNL_vstar_plus*(2./3.*W +q2/a/Mnuc + nomg/3./a/Mnuc) )* fZeta * (GA/2./W/Qstar);

    KNL_C_minus = ( (KNL_Qstar_minus*Qstar - KNL_vstar_minus*vstar ) * ( 1./3. + vstar/a/Mnuc)
        + KNL_vstar_minus*(2./3.*W +q2/a/Mnuc + nomg/3./a/Mnuc) )* fZeta * (GA/2./W/Qstar);

    LOG("BSKLNBaseRESPXSec2014",pINFO)  <<"KNL C= "<<KNL_C_plus<<"\t"<<KNL_C_minus<<"\t"<<fFKR.C;
  }
  double BRS_S_plus = 0;
  double BRS_S_minus = 0;
  double BRS_B_plus = 0;
  double BRS_B_minus = 0;
  double BRS_C_plus = 0;
  double BRS_C_minus = 0;


  if(is_CC && is_BRS){

    KNL_S_plus  = (KNL_vstar_plus*vstar  - KNL_Qstar_plus *Qstar )* (Mnuc2 -q2 - 3*W*Mnuc ) * GV / (6*Mnuc2)/Q2;
    KNL_S_minus = (KNL_vstar_minus*vstar - KNL_Qstar_minus*Qstar )* (Mnuc2 -q2 - 3*W*Mnuc ) * GV / (6*Mnuc2)/Q2;


    KNL_B_plus  = fZeta/(3.*W*sq2omg)/Qstar * (KNL_Qstar_plus  + KNL_vstar_plus *Qstar/a/Mnuc ) * GA;
    KNL_B_minus = fZeta/(3.*W*sq2omg)/Qstar * (KNL_Qstar_minus + KNL_vstar_minus*Qstar/a/Mnuc ) * GA;


    KNL_C_plus = ( (KNL_Qstar_plus*Qstar - KNL_vstar_plus*vstar ) * ( 1./3. + vstar/a/Mnuc)
        + KNL_vstar_plus*(2./3.*W +q2/a/Mnuc + nomg/3./a/Mnuc) )* fZeta * (GA/2./W/Qstar);

    KNL_C_minus = ( (KNL_Qstar_minus*Qstar - KNL_vstar_minus*vstar ) * ( 1./3. + vstar/a/Mnuc)
        + KNL_vstar_minus*(2./3.*W +q2/a/Mnuc + nomg/3./a/Mnuc) )* fZeta * (GA/2./W/Qstar);

    BRS_S_plus = KNL_S_plus;
    BRS_S_minus = KNL_S_minus;
    LOG("BSKLNBaseRESPXSec2014",pINFO) <<"BRS S= " <<KNL_S_plus<<"\t"<<KNL_S_minus<<"\t"<<fFKR.S;

    BRS_B_plus = KNL_B_plus + fZeta*GA/2./W/Qstar*( KNL_Qstar_plus*vstar - KNL_vstar_plus*Qstar)
      *( 2./3 /sq2omg *(vstar + Qstar*Qstar/Mnuc/a))/(kPionMass2 -q2);

    BRS_B_minus = KNL_B_minus + fZeta*GA/2./W/Qstar*( KNL_Qstar_minus*vstar - KNL_vstar_minus*Qstar)
      *( 2./3 /sq2omg *(vstar + Qstar*Qstar/Mnuc/a))/(kPionMass2 -q2);
    LOG("BSKLNBaseRESPXSec2014",pINFO) <<"BRS B= " <<KNL_B_plus<<"\t"<<KNL_B_minus<<"\t"<<fFKR.B;

    BRS_C_plus = KNL_C_plus  + fZeta*GA/2./W/Qstar*( KNL_Qstar_plus*vstar - KNL_vstar_plus*Qstar)
      * Qstar*(2./3.*W +q2/Mnuc/a +nomg/3./a/Mnuc)/(kPionMass2 -q2);

    BRS_C_minus = KNL_C_minus  + fZeta*GA/2./W/Qstar*( KNL_Qstar_minus*vstar - KNL_vstar_minus*Qstar)
      * Qstar*(2./3.*W +q2/Mnuc/a +nomg/3./a/Mnuc)/(kPionMass2 -q2);
    LOG("BSKLNBaseRESPXSec2014",pINFO) <<"BRS C= " <<KNL_C_plus<<"\t"<<KNL_C_minus<<"\t"<<fFKR.C;
  }

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("FKR", pDEBUG)
    << "FKR params for RES = " << resname << " : " << fFKR;
#endif

  // Calculate the Rein-Sehgal Helicity Amplitudes
  double sigL_minus = 0;
  double sigR_minus = 0;
  double sigS_minus = 0;

  double sigL_plus = 0;
  double sigR_plus = 0;
  double sigS_plus = 0;

  const RSHelicityAmplModelI * hamplmod = 0;
  const RSHelicityAmplModelI * hamplmod_KNL_minus = 0;
  const RSHelicityAmplModelI * hamplmod_KNL_plus = 0;
  const RSHelicityAmplModelI * hamplmod_BRS_minus = 0;
  const RSHelicityAmplModelI * hamplmod_BRS_plus = 0;

  // These lines were ~ 100 lines below, which means that, for EM interactions, the coefficients below were still calculated using the weak coupling constant - Afro
  double g2 = kGF2;

  // For EM interaction replace  G_{Fermi} with :
  // a_{em} * pi / ( sqrt(2) * sin^2(theta_weinberg) * Mass_{W}^2 }
  // See C.Quigg, Gauge Theories of the Strong, Weak and E/M Interactions,
  // ISBN 0-8053-6021-2, p.112 (6.3.57)
  // Also, take int account that the photon propagator is 1/p^2 but the
  // W propagator is 1/(p^2-Mass_{W}^2), so weight the EM case with
  // Mass_{W}^4 / q^4
  // So, overall:
  // G_{Fermi}^2 --> a_{em}^2 * pi^2 / (2 * sin^4(theta_weinberg) * q^{4})
  //

  if(is_EM) {
    double q4 = q2*q2;
    g2 = kAem2 * kPi2 / (2.0 * fSin48w * q4);
  }

  if(is_CC) g2 = kGF2*fVud2;

  double sig0 = 0.125*(g2/kPi)*(-q2/Q2)*(W/Mnuc);
  double scLR = W/Mnuc;
  double scS  = (Mnuc/W)*(-Q2/q2);

  double sigL =0;
  double sigR =0;
  double sigS =0;

  double sigRSL =0;
  double sigRSR =0;
  double sigRSS =0;

  if(is_CC && !(is_KLN || is_BRS) ) {

    hamplmod = fHAmplModelCC;
  }
  else
    if(is_NC) {
      if (is_p) { hamplmod = fHAmplModelNCp;}
      else      { hamplmod = fHAmplModelNCn;}
    }
    else
      if(is_EM) {
        if (is_p) { hamplmod = fHAmplModelEMp;}
        else      { hamplmod = fHAmplModelEMn;}
      }
      else
        if(is_CC && is_KLN ){
          fFKR.S = KNL_S_minus;        //2 times fFKR.S?
          fFKR.B = KNL_B_minus;
          fFKR.C = KNL_C_minus;

          hamplmod_KNL_minus = fHAmplModelCC;

          assert(hamplmod_KNL_minus);

          const RSHelicityAmpl & hampl_KNL_minus = hamplmod_KNL_minus->Compute(resonance, fFKR);

          sigL_minus = (hampl_KNL_minus.Amp2Plus3 () + hampl_KNL_minus.Amp2Plus1 ());
          sigR_minus = (hampl_KNL_minus.Amp2Minus3() + hampl_KNL_minus.Amp2Minus1());
          sigS_minus = (hampl_KNL_minus.Amp20Plus () + hampl_KNL_minus.Amp20Minus());


          fFKR.S = KNL_S_plus;
          fFKR.B = KNL_B_plus;
          fFKR.C = KNL_C_plus;
          hamplmod_KNL_plus = fHAmplModelCC;
          assert(hamplmod_KNL_plus);

          const RSHelicityAmpl & hampl_KNL_plus = hamplmod_KNL_plus->Compute(resonance, fFKR);

          sigL_plus = (hampl_KNL_plus.Amp2Plus3 () + hampl_KNL_plus.Amp2Plus1 ());
          sigR_plus = (hampl_KNL_plus.Amp2Minus3() + hampl_KNL_plus.Amp2Minus1());
          sigS_plus = (hampl_KNL_plus.Amp20Plus () + hampl_KNL_plus.Amp20Minus());

        }
        else
          if(is_CC && is_BRS ){
            fFKR.S = BRS_S_minus;
            fFKR.B = BRS_B_minus;
            fFKR.C = BRS_C_minus;

            hamplmod_BRS_minus = fHAmplModelCC;
            assert(hamplmod_BRS_minus);

            const RSHelicityAmpl & hampl_BRS_minus = hamplmod_BRS_minus->Compute(resonance, fFKR);

            sigL_minus = (hampl_BRS_minus.Amp2Plus3 () + hampl_BRS_minus.Amp2Plus1 ());
            sigR_minus = (hampl_BRS_minus.Amp2Minus3() + hampl_BRS_minus.Amp2Minus1());
            sigS_minus = (hampl_BRS_minus.Amp20Plus () + hampl_BRS_minus.Amp20Minus());

            fFKR.S = BRS_S_plus;
            fFKR.B = BRS_B_plus;
            fFKR.C = BRS_C_plus;
            hamplmod_BRS_plus = fHAmplModelCC;
            assert(hamplmod_BRS_plus);

            const RSHelicityAmpl & hampl_BRS_plus = hamplmod_BRS_plus->Compute(resonance, fFKR);

            sigL_plus = (hampl_BRS_plus.Amp2Plus3 () + hampl_BRS_plus.Amp2Plus1 ());
            sigR_plus = (hampl_BRS_plus.Amp2Minus3() + hampl_BRS_plus.Amp2Minus1());
            sigS_plus = (hampl_BRS_plus.Amp20Plus () + hampl_BRS_plus.Amp20Minus());
          }

  // Compute the cross section
  if(is_KLN || is_BRS) {

     sigL_minus *= scLR;
     sigR_minus *= scLR;
     sigS_minus *= scS;
     sigL_plus  *= scLR;
     sigR_plus  *= scLR;
     sigS_plus  *= scS;

     LOG("BSKLNBaseRESPXSec2014", pINFO)
         << "sL,R,S minus = " << sigL_minus << "," << sigR_minus << "," << sigS_minus;
     LOG("BSKLNBaseRESPXSec2014", pINFO)
         << "sL,R,S plus = " << sigL_plus << "," << sigR_plus << "," << sigS_plus;
  }
  else {
     assert(hamplmod);

     const RSHelicityAmpl & hampl = hamplmod->Compute(resonance, fFKR);

     sigL = scLR* (hampl.Amp2Plus3 () + hampl.Amp2Plus1 ());
     sigR = scLR* (hampl.Amp2Minus3() + hampl.Amp2Minus1());
     sigS = scS * (hampl.Amp20Plus () + hampl.Amp20Minus());
  }

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("BSKLNBaseRESPXSec2014", pDEBUG) << "sig_{0} = " << sig0;
  LOG("BSKLNBaseRESPXSec2014", pDEBUG) << "sig_{L} = " << sigL;
  LOG("BSKLNBaseRESPXSec2014", pDEBUG) << "sig_{R} = " << sigR;
  LOG("BSKLNBaseRESPXSec2014", pDEBUG) << "sig_{S} = " << sigS;
#endif

  double xsec = 0.0;

  if(is_KLN || is_BRS) {
      xsec =   TMath::Power(KNL_cL_minus,2)*sigL_minus + TMath::Power(KNL_cL_plus,2)*sigL_plus
             + TMath::Power(KNL_cR_minus,2)*sigR_minus + TMath::Power(KNL_cR_plus,2)*sigR_plus
             + TMath::Power(KNL_cS_minus,2)*sigS_minus + TMath::Power(KNL_cS_plus,2)*sigS_plus;
      xsec *=sig0;

      LOG("BSKLNBaseRESPXSec2014",pINFO) << "A-="<<KNL_Alambda_minus<<" A+="<<KNL_Alambda_plus;
      // protect against sigRSR=sigRSL=sigRSS=0
      LOG("BSKLNBaseRESPXSec2014",pINFO) <<q2<<"\t"<<xsec<<"\t"<<sig0*(V2*sigR + U2*sigL + 2*UV*sigS)<<"\t"<<xsec/TMath::Max(sig0*(V2*sigRSR + U2*sigRSL + 2*UV*sigRSS),1.0e-100);
      LOG("BSKLNBaseRESPXSec2014",pINFO) <<"fFKR.B="<<fFKR.B<<" fFKR.C="<<fFKR.C<<" fFKR.S="<<fFKR.S;
      LOG("BSKLNBaseRESPXSec2014",pINFO) <<"CL-="<<TMath::Power(KNL_cL_minus,2)<<" CL+="<<TMath::Power(KNL_cL_plus,2)<<" U2="<<U2;
      LOG("BSKLNBaseRESPXSec2014",pINFO) <<"SL-="<<sigL_minus<<" SL+="<<sigL_plus<<" SL="<<sigRSL;

      LOG("BSKLNBaseRESPXSec2014",pINFO) <<"CR-="<<TMath::Power(KNL_cR_minus,2)<<" CR+="<<TMath::Power(KNL_cR_plus,2)<<" V2="<<V2;
      LOG("BSKLNBaseRESPXSec2014",pINFO) <<"SR-="<<sigR_minus<<" SR+="<<sigR_plus<<" sR="<<sigRSR;

      LOG("BSKLNBaseRESPXSec2014",pINFO) <<"CS-="<<TMath::Power(KNL_cS_minus,2)<<" CS+="<<TMath::Power(KNL_cS_plus,2)<<" UV="<<UV;
      LOG("BSKLNBaseRESPXSec2014",pINFO) <<"SS-="<<sigL_minus<<" SS+="<<sigS_plus<<" sS="<<sigRSS;
  }
  else {
     if (is_nu || is_lminus) {
         xsec = sig0*(V2*sigR + U2*sigL + 2*UV*sigS);
     }
     else
     if (is_nubar || is_lplus) {
         xsec = sig0*(U2*sigR + V2*sigL + 2*UV*sigS);
     }
     xsec = TMath::Max(0.,xsec);
  }
  double mult = 1.0;
  if ( is_CC && is_delta ) {
     if ( (is_nu && is_p) || (is_nubar && is_n) ) mult=3.0;
  }
  xsec *= mult;

  // Check whether the cross section is to be weighted with a Breit-Wigner distribution
  // (default: true)
  double bw = 1.0;
  if ( fWghtBW ) {
     bw = utils::bwfunc::BreitWignerL(W,LR,MR,WR,NR);
  }
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("BSKLNBaseRESPXSec2014", pDEBUG)
      << "BreitWigner(RES=" << resname << ", W=" << W << ") = " << bw;
#endif
  xsec *= bw;

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("BSKLNBaseRESPXSec2014", pINFO)
      << "\n d2xsec/dQ2dW"  << "[" << interaction->AsString()
      << "](W=" << W << ", q2=" << q2 << ", E=" << E << ") = " << xsec;
#endif

  // The algorithm computes d^2xsec/dWdQ2
  // Check whether variable tranformation is needed
  if ( kps != kPSWQ2fE ) {
     double J = utils::kinematics::Jacobian(interaction,kPSWQ2fE,kps);
     xsec *= J;
  }

  // Apply given scaling factor
  if      (is_CC) { xsec *= fXSecScaleCC; }
  else if (is_NC) { xsec *= fXSecScaleNC; }

  // If requested return the free nucleon xsec even for input nuclear tgt
  if ( interaction->TestBit(kIAssumeFreeNucleon) ) return xsec;

  int Z = target.Z();
  int A = target.A();
  int N = A-Z;

  // Take into account the number of scattering centers in the target
  int NNucl = (is_p) ? Z : N;
  xsec*=NNucl; // nuclear xsec (no nuclear suppression factor)

  if ( fUsePauliBlocking && A!=1 )
  {
    // Calculation of Pauli blocking according references:
    //
    //     [1] S.L. Adler,  S. Nussinov,  and  E.A.  Paschos,  "Nuclear
    //         charge exchange corrections to leptonic pion  production
    //         in  the (3,3) resonance  region,"  Phys. Rev. D 9 (1974)
    //         2125-2143 [Erratum Phys. Rev. D 10 (1974) 1669].
    //     [2] J.Y. Yu, "Neutrino interactions and  nuclear  effects in
    //         oscillation experiments and the  nonperturbative disper-
    //         sive  sector in strong (quasi-)abelian  fields,"  Ph. D.
    //         Thesis, Dortmund U., Dortmund, 2002 (unpublished).
    //     [3] E.A. Paschos, J.Y. Yu,  and  M. Sakuda,  "Neutrino  pro-
    //         duction  of  resonances,"  Phys. Rev. D 69 (2004) 014013
    //         [arXiv: hep-ph/0308130].

    double P_Fermi = 0.0;

    // Maximum value of Fermi momentum of target nucleon (GeV)
    if ( A<6 || ! fUseRFGParametrization )
    {
        // look up the Fermi momentum for this target
        FermiMomentumTablePool * kftp = FermiMomentumTablePool::Instance();
        const FermiMomentumTable * kft = kftp->GetTable(fKFTable);
        P_Fermi = kft->FindClosestKF(pdg::IonPdgCode(A, Z), nucpdgc);
     }
     else {
        // define the Fermi momentum for this target
        P_Fermi = utils::nuclear::FermiMomentumForIsoscalarNucleonParametrization(target);
        // correct the Fermi momentum for the struck nucleon
        if(is_p) { P_Fermi *= TMath::Power( 2.*Z/A, 1./3); }
        else     { P_Fermi *= TMath::Power( 2.*N/A, 1./3); }
     }

     double FactorPauli_RES = 1.0;

     double k0 = 0., q = 0., q0 = 0.;

     if (P_Fermi > 0.)
     {
        k0 = (W2-Mnuc2-Q2)/(2*W);
        k = TMath::Sqrt(k0*k0+Q2);  // previous value of k is overridden
        q0 = (W2-Mnuc2+kPionMass2)/(2*W);
        q = TMath::Sqrt(q0*q0-kPionMass2);
     }

     if ( 2*P_Fermi < k-q )
        FactorPauli_RES = 1.0;
     if ( 2*P_Fermi >= k+q )
        FactorPauli_RES = ((3*k*k+q*q)/(2*P_Fermi)-(5*TMath::Power(k,4)+TMath::Power(q,4)+10*k*k*q*q)/(40*TMath::Power(P_Fermi,3)))/(2*k);
     if ( 2*P_Fermi >= k-q && 2*P_Fermi <= k+q )
        FactorPauli_RES = ((q+k)*(q+k)-4*P_Fermi*P_Fermi/5-TMath::Power(k-q, 3)/(2*P_Fermi)+TMath::Power(k-q, 5)/(40*TMath::Power(P_Fermi, 3)))/(4*q*k);

     xsec *= FactorPauli_RES;
  }
  return xsec;
}
//____________________________________________________________________________
double BSKLNBaseRESPXSec2014::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool BSKLNBaseRESPXSec2014::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const InitialState & init_state = interaction->InitState();
  const ProcessInfo &  proc_info  = interaction->ProcInfo();
  const XclsTag &      xcls       = interaction->ExclTag();

  if(!proc_info.IsResonant()) return false;
  if(!xcls.KnownResonance())  return false;

  int  hitnuc = init_state.Tgt().HitNucPdg();
  bool is_pn = (pdg::IsProton(hitnuc) || pdg::IsNeutron(hitnuc));

  if (!is_pn) return false;

  int  probe   = init_state.ProbePdg();
  bool is_weak = proc_info.IsWeak();
  bool is_em   = proc_info.IsEM();
  bool nu_weak = (pdg::IsNeutralLepton(probe) && is_weak);
  bool l_em    = (pdg::IsChargedLepton(probe) && is_em  );

  if (!nu_weak && !l_em) return false;

  return true;
}
//____________________________________________________________________________
void BSKLNBaseRESPXSec2014::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void BSKLNBaseRESPXSec2014::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void BSKLNBaseRESPXSec2014::LoadConfig(void)
{
  // Cross section scaling factors
  this->GetParam( "RES-CC-XSecScale", fXSecScaleCC ) ;
  this->GetParam( "RES-NC-XSecScale", fXSecScaleNC ) ;

  // Load all configuration data or set defaults

  this->GetParam( "RES-Zeta"   , fZeta  ) ;
  this->GetParam( "RES-Omega"  , fOmega ) ;
  this->GetParam( "minibooneGA", fGA    ) ;
  this->GetParam( "minibooneGV", fGV    ) ;

  double ma, mv ;
  this->GetParam( "RES-Ma", ma ) ;
  this->GetParam( "RES-Mv", mv ) ;
  fMa2 = TMath::Power(ma,2);
  fMv2 = TMath::Power(mv,2);

  this->GetParamDef( "BreitWignerWeight", fWghtBW, true ) ;
  this->GetParamDef( "BreitWignerNorm",   fNormBW, true);
  double thw ;
  this->GetParam( "WeinbergAngle", thw ) ;
  fSin48w = TMath::Power( TMath::Sin(thw), 4 );
  double Vud;
  this->GetParam("CKM-Vud", Vud );
  fVud2 = TMath::Power( Vud, 2 );
  this->GetParam("FermiMomentumTable", fKFTable);
  this->GetParam("RFG-UseParametrization", fUseRFGParametrization);
  this->GetParam("UsePauliBlockingForRES", fUsePauliBlocking);

  // Load all the sub-algorithms needed

  fHAmplModelCC     = 0;
  fHAmplModelNCp    = 0;
  fHAmplModelNCn    = 0;
  fHAmplModelEMp    = 0;
  fHAmplModelEMn    = 0;

  AlgFactory * algf = AlgFactory::Instance();

  fHAmplModelCC  = dynamic_cast<const RSHelicityAmplModelI *> (
      algf->GetAlgorithm("genie::RSHelicityAmplModelCC","Default"));
  fHAmplModelNCp = dynamic_cast<const RSHelicityAmplModelI *> (
      algf->GetAlgorithm("genie::RSHelicityAmplModelNCp","Default"));
  fHAmplModelNCn = dynamic_cast<const RSHelicityAmplModelI *> (
      algf->GetAlgorithm("genie::RSHelicityAmplModelNCn","Default"));
  fHAmplModelEMp = dynamic_cast<const RSHelicityAmplModelI *> (
      algf->GetAlgorithm("genie::RSHelicityAmplModelEMp","Default"));
  fHAmplModelEMn = dynamic_cast<const RSHelicityAmplModelI *> (
      algf->GetAlgorithm("genie::RSHelicityAmplModelEMn","Default"));

  assert( fHAmplModelCC  );
  assert( fHAmplModelNCp );
  assert( fHAmplModelNCn );
  assert( fHAmplModelEMp );
  assert( fHAmplModelEMn );

  // Use algorithm within a DIS/RES join scheme. If yes get Wcut
  this->GetParam( "UseDRJoinScheme", fUsingDisResJoin ) ;
  fWcut = 999999;
  if(fUsingDisResJoin) {
    this->GetParam( "Wcut", fWcut ) ;
  }

  // NeuGEN limits in the allowed resonance phase space:
  // W < min{ Wmin(physical), (res mass) + x * (res width) }
  // It limits the integration area around the peak and avoids the
  // problem with huge xsec increase at low Q2 and high W.
  // In correspondence with Hugh, Rein said that the underlying problem
  // are unphysical assumptions in the model.
  this->GetParamDef( "MaxNWidthForN2Res", fN2ResMaxNWidths, 2.0 ) ;
  this->GetParamDef( "MaxNWidthForN0Res", fN0ResMaxNWidths, 6.0 ) ;
  this->GetParamDef( "MaxNWidthForGNRes", fGnResMaxNWidths, 4.0 ) ;

  // Load the differential cross section integrator
  fXSecIntegrator =
    dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);
}
//____________________________________________________________________________
