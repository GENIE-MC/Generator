//____________________________________________________________________________
/*
 Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

 @ Oct 05, 2009 - CA
   Modified code to handle charged lepton scattering too.
   Also, the helicity amplitude code now returns a `const RSHelicityAmpl &'.
 @ July 23, 2010 - CA
   BaryonResParams, and BreitWignerI, BaryonResDataSetI implementations are
   now redundant. Get resonance parameters from BaryonResUtils and use the
   Breit-Weigner functions from utils::bwfunc.
 @ December 15, 2014 - JN, SD
   Add new version due to Jarek Nowak.  Based on parameters set in
   UserPhysicsOptions.xml.  Default is Berger-Sehgal (Phys. Rev. D76, 
   113004 (2007)) with new GA and new GV.  Additional model due to Kuzmin, 
   Lyubushkin, and Naumov (Mod. Phys. Lett. A19 (2004) 2815 and Phys. Part. 
   Nucl. 35 (2004) S133) also available.  Each of these models
   includes effect of lepton mass.  BS adds a new pole diagram.  
   New GA and GV form factors are based on studies with MiniBooNE data.
 @ Dec 22, 2014 - GP
   Incorporating changes from J. Nowak into a new class (was 
   ReinSehgalRESPXSec, now BergerSehgalRESPXSec2014).
*/
//____________________________________________________________________________

#include <TMath.h>
#include <TSystem.h>

#include "Algorithm/AlgFactory.h"
#include "Algorithm/AlgConfigPool.h"
#include "Base/XSecIntegratorI.h"
#include "BaryonResonance/BaryonResUtils.h"
#include "Conventions/GBuild.h"
#include "Conventions/Constants.h"
#include "Conventions/RefFrame.h"
#include "Conventions/KineVar.h"
#include "Conventions/Units.h"
#include "Messenger/Messenger.h"
#include "Numerical/Spline.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "ReinSehgal/BergerSehgalRESPXSec2014.h"
#include "ReinSehgal/RSHelicityAmplModelI.h"
#include "ReinSehgal/RSHelicityAmpl.h"
#include "Utils/KineUtils.h"
#include "Utils/MathUtils.h"
#include "Utils/Range1.h"
#include "Utils/BWFunc.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
BergerSehgalRESPXSec2014::BergerSehgalRESPXSec2014() :
  XSecAlgorithmI("genie::BergerSehgalRESPXSec2014")
{
  fNuTauRdSpl    = 0;
  fNuTauBarRdSpl = 0;
}
//____________________________________________________________________________
BergerSehgalRESPXSec2014::BergerSehgalRESPXSec2014(string config) :
  XSecAlgorithmI("genie::BergerSehgalRESPXSec2014", config)
{
  fNuTauRdSpl    = 0;
  fNuTauBarRdSpl = 0;
}
//____________________________________________________________________________
BergerSehgalRESPXSec2014::~BergerSehgalRESPXSec2014()
{
  if(fNuTauRdSpl)    delete fNuTauRdSpl;
  if(fNuTauBarRdSpl) delete fNuTauBarRdSpl;
}
//____________________________________________________________________________
double BergerSehgalRESPXSec2014::XSec(
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
      LOG("BergerSehgalRESPXSec2014", pDEBUG)
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
  double NR  = utils::res::BWNorm            (resonance);

  // Following NeuGEN, avoid problems with underlying unphysical
  // model assumptions by restricting the allowed W phase space
  // around the resonance peak
  if      (W > MR + fN0ResMaxNWidths * WR && IR==0) return 0.;
  else if (W > MR + fN2ResMaxNWidths * WR && IR==2) return 0.;
  else if (W > MR + fGnResMaxNWidths * WR)          return 0.;

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

  double KNL_vstar_plus = 0 ;
  double KNL_vstar_minus = 0 ;



  if(is_CC && (is_KLN || is_BRS)){

    LOG("BergerSehgalRESPXSec2014",pINFO) "costh1="<<costh;    
    costh = (q2 - ml*ml + 2.*E*Eprime)/2./E/Pl;
    //ml=0;
    LOG("BergerSehgalRESPXSec2014",pINFO) "q2="<<q2<< "m2="<<ml*ml<<" 2.*E*Eprime="<<2.*E*Eprime<<" nom="<< (q2 - ml*ml + 2.*E*Eprime)<<" den="<<2.*E*Pl;
    LOG("BergerSehgalRESPXSec2014",pINFO) "costh2="<<costh;
    Pl    = TMath::Sqrt(Eprime*Eprime - ml*ml);


    if(costh <= -1. + 1e-7) {
      LOG("BergerSehgalRESPXSec2014", pDEBUG)
        << "Changing costh = " << costh << " to -1";
      costh = -1 + 1e-6;
    }
    if(costh >= 1. - 1e-7){
      LOG("BergerSehgalRESPXSec2014", pDEBUG)
        << "Changing costh = " << costh << " to 1";
      costh = 1 - 1e-6;
    }
    vstar = (Mnuc*v + q2)/W;//missing W
    Qstar = TMath::Sqrt(-q2 + vstar*vstar);


    KNL_Alambda_plus  = TMath::Sqrt(E*(Eprime - Pl));
    KNL_Alambda_minus = TMath::Sqrt(E*(Eprime + Pl));
    LOG("BergerSehgalRESPXSec2014",pINFO) <<"+++++++++++++++++++++++";
    LOG("BergerSehgalRESPXSec2014",pINFO) <<"E="<<E << " K= "<<KNL_K;
    LOG("BergerSehgalRESPXSec2014",pINFO) <<"El="<<Eprime<<" Pl="<<Pl<<" ml="<<ml;
    LOG("BergerSehgalRESPXSec2014",pINFO) <<"W="<<W<<" Q="<<Q<<" q2="<<q2;
    LOG("BergerSehgalRESPXSec2014",pINFO) <<"A-="<<KNL_Alambda_minus<<" A+="<<KNL_Alambda_plus;
    LOG("BergerSehgalRESPXSec2014",pINFO) <<"xxxxxxxxxxxxxxxxxxxxxxx";

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
      KNL_cL_plus  =  1 * TMath::Sqrt(0.5)* KNL_K * (KNL_jx_minus - KNL_jy_minus);
      KNL_cL_minus = -1 * TMath::Sqrt(0.5)* KNL_K * (KNL_jx_plus  - KNL_jy_plus);

      KNL_cR_plus  =  1 * TMath::Sqrt(0.5)* KNL_K * (KNL_jx_minus + KNL_jy_minus);
      KNL_cR_minus = -1 * TMath::Sqrt(0.5)* KNL_K * (KNL_jx_plus  + KNL_jy_plus);

      KNL_cS_plus  = -1 * KNL_K *  TMath::Sqrt(TMath::Abs(KNL_j0_minus*KNL_j0_minus - KNL_jz_minus*KNL_jz_minus) );
      KNL_cS_minus =  1 * KNL_K *  TMath::Sqrt(TMath::Abs(KNL_j0_plus*KNL_j0_plus - KNL_jz_plus*KNL_jz_plus) );
    }
  }

     LOG("BergerSehgalRESPXSec2014",pINFO) <<"j0-="<<KNL_j0_minus<<" j0+="<<KNL_j0_plus;
     LOG("BergerSehgalRESPXSec2014",pINFO) <<"jx-="<<KNL_jx_minus<<" jx+="<<KNL_jx_plus;
     LOG("BergerSehgalRESPXSec2014",pINFO) <<"jy-="<<KNL_jy_minus<<" jy+="<<KNL_jy_plus;
     LOG("BergerSehgalRESPXSec2014",pINFO) <<"jz-="<<KNL_jz_minus<<" jz+="<<KNL_jz_plus;

     LOG("BergerSehgalRESPXSec2014",pINFO) "sqrt2="<<sqrtq2<<" jz+=:"<<KNL_jz_plus<<" j0+="<<KNL_j0_plus<<" denom="<<TMath::Sqrt(TMath::Abs(KNL_j0_plus*KNL_j0_plus - KNL_jz_plus*KNL_jz_plus) );

     LOG("BergerSehgalRESPXSec2014",pINFO) <<"vstar-="<<KNL_vstar_minus<<" vstar+="<<KNL_vstar_plus;
     LOG("BergerSehgalRESPXSec2014",pINFO) <<"Qstar-="<<KNL_Qstar_minus<<" Qstar+="<<KNL_Qstar_plus;

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("BergerSehgalRESPXSec2014", pDEBUG) 
    << "Kinematical params V = " << V << ", U = " << U;
#endif

  // Calculate the Feynman-Kislinger-Ravndall parameters

  double Go  = TMath::Power(1 - 0.25 * q2/Mnuc2, 0.5-IR);
  double GV  = Go * TMath::Power( 1./(1-q2/fMv2), 2);
  double GA  = Go * TMath::Power( 1./(1-q2/fMa2), 2);


  //JN New form factors code
  //  new_GV = false; //JN
  //  new_GA = true; //JN

  if(fGV){

        LOG("BergerSehgal2014",pDEBUG) <<"Using new GV";
    double CV0 =  1./(1-q2/fMv2/4.);
    double CV3 =  2.13 * CV0 * TMath::Power( 1-q2/fMv2,-2);
    double CV4 = -1.51 * CV0 * TMath::Power( 1-q2/fMv2,-2);
    double CV5 =  0.48 * CV0 * TMath::Power( 1-q2/fMv2/0.766, -2);


    double GV3 =   0.5 / TMath::Sqrt(3) * ( CV3 * (W + Mnuc)/Mnuc
        + CV4 * (W2 + q2 -Mnuc2)/2./Mnuc2 
        + CV5 * (W2 - q2 -Mnuc2)/2./Mnuc2 );

    double GV1 = - 0.5 / TMath::Sqrt(3) * ( CV3 * (Mnuc2 -q2 +Mnuc*W)/W/Mnuc
        + CV4 * (W2 +q2 - Mnuc2)/2./Mnuc2 
        + CV5 * (W2 -q2 - Mnuc2)/2./Mnuc2 );

    GV = 0.5 * TMath::Power( 1 - q2/(Mnuc + W)/(Mnuc + W), 0.5-IR)
      * TMath::Sqrt( 3 * GV3*GV3 + GV1*GV1);


  }

  if(fGA){

    LOG("BergerSehgalRESPXSec2014",pDEBUG) <<"Using new GA";

    double CA5_0 = 1.2;
    double CA5 = CA5_0 *  TMath::Power( 1./(1-q2/fMa2), 2);
    //  GA = 0.5 * TMath::Sqrt(3.) * TMath::Power( 1 - q2/(Mnuc + W)/(Mnuc + W), 0.5-IR) * (1- (W2 +q2 -Mnuc2)/8./Mnuc2) * CA5/fZeta;
    GA = 0.5 * TMath::Sqrt(3.) * TMath::Power( 1 - q2/(Mnuc + W)/(Mnuc + W), 0.5-IR) * (1- (W2 +q2 -Mnuc2)/8./Mnuc2) * CA5;

    //    LOG("BergerSehgal2014",pINFO) << 0.5 * TMath::Sqrt(3.) * TMath::Power( 1 - q2/(Mnuc + W)/(Mnuc + W), 0.5-IR)* (1- (W2 +q2 -Mnuc2)/8./Mnuc2);
    LOG("BergerSehgalRESPXSec2014",pINFO) <<"GA= " <<GA << "  C5A= " <<CA5;
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

    LOG("BergerSehgalRESPXSec2014",pINFO) <<"KNL S= " <<KNL_S_plus<<"\t"<<KNL_S_minus<<"\t"<<fFKR.S;

    KNL_B_plus  = fZeta/(3.*W*sq2omg)/Qstar * (KNL_Qstar_plus  + KNL_vstar_plus *Qstar/a/Mnuc ) * GA;
    KNL_B_minus = fZeta/(3.*W*sq2omg)/Qstar * (KNL_Qstar_minus + KNL_vstar_minus*Qstar/a/Mnuc ) * GA;
    LOG("BergerSehgalRESPXSec2014",pINFO) <<"KNL B= " <<KNL_B_plus<<"\t"<<KNL_B_minus<<"\t"<<fFKR.B;

    KNL_C_plus = ( (KNL_Qstar_plus*Qstar - KNL_vstar_plus*vstar ) * ( 1./3. + vstar/a/Mnuc)
        + KNL_vstar_plus*(2./3.*W +q2/a/Mnuc + nomg/3./a/Mnuc) )* fZeta * (GA/2./W/Qstar);

    KNL_C_minus = ( (KNL_Qstar_minus*Qstar - KNL_vstar_minus*vstar ) * ( 1./3. + vstar/a/Mnuc)
        + KNL_vstar_minus*(2./3.*W +q2/a/Mnuc + nomg/3./a/Mnuc) )* fZeta * (GA/2./W/Qstar);

    LOG("BergerSehgalRESPXSec2014",pINFO)  <<"KNL C= "<<KNL_C_plus<<"\t"<<KNL_C_minus<<"\t"<<fFKR.C;
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
    LOG("BergerSehgalRESPXSec2014",pINFO) <<"BRS S= " <<KNL_S_plus<<"\t"<<KNL_S_minus<<"\t"<<fFKR.S;

    BRS_B_plus = KNL_B_plus + fZeta*GA/2./W/Qstar*( KNL_Qstar_plus*vstar - KNL_vstar_plus*Qstar)
      *( 2./3 /sq2omg *(vstar + Qstar*Qstar/Mnuc/a))/(kPionMass2 -q2);

    BRS_B_minus = KNL_B_minus + fZeta*GA/2./W/Qstar*( KNL_Qstar_minus*vstar - KNL_vstar_minus*Qstar)
      *( 2./3 /sq2omg *(vstar + Qstar*Qstar/Mnuc/a))/(kPionMass2 -q2);
    LOG("BergerSehgalRESPXSec2014",pINFO) <<"BRS B= " <<KNL_B_plus<<"\t"<<KNL_B_minus<<"\t"<<fFKR.B;

    BRS_C_plus = KNL_C_plus  + fZeta*GA/2./W/Qstar*( KNL_Qstar_plus*vstar - KNL_vstar_plus*Qstar)
      * Qstar*(2./3.*W +q2/Mnuc/a +nomg/3./a/Mnuc)/(kPionMass2 -q2);

    BRS_C_minus = KNL_C_minus  + fZeta*GA/2./W/Qstar*( KNL_Qstar_minus*vstar - KNL_vstar_minus*Qstar)
      * Qstar*(2./3.*W +q2/Mnuc/a +nomg/3./a/Mnuc)/(kPionMass2 -q2);
    LOG("BergerSehgalRESPXSec2014",pINFO) <<"BRS C= " <<KNL_C_plus<<"\t"<<KNL_C_minus<<"\t"<<fFKR.C;

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

  double g2 = kGF2;

  double sig0 = 0.125*(g2/kPi)*(-q2/Q2)*(W/Mnuc);
  double scLR = W/Mnuc;
  double scS  = (Mnuc/W)*(-Q2/q2);

  double sigL =0;
  double sigR =0;
  double sigS =0;

  double sigRSL =0;
  double sigRSR =0;
  double sigRSS =0;

  /*
     hamplmod = fHAmplModelCC;
     const RSHelicityAmpl & hampl = hamplmod->Compute(resonance, fFKR);

     sigRSL = scLR* (hampl.Amp2Plus3 () + hampl.Amp2Plus1 ());
     sigRSR = scLR* (hampl.Amp2Minus3() + hampl.Amp2Minus1());
     sigRSS = scS * (hampl.Amp20Plus () + hampl.Amp20Minus());

*/
  //<<<<<<<<<
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
          fFKR.S = KNL_C_minus;

          hamplmod_KNL_minus = fHAmplModelCC;

          assert(hamplmod_KNL_minus);

          const RSHelicityAmpl & hampl_KNL_minus = hamplmod_KNL_minus->Compute(resonance, fFKR);

          sigL_minus = (hampl_KNL_minus.Amp2Plus3 () + hampl_KNL_minus.Amp2Plus1 ());
          sigR_minus = (hampl_KNL_minus.Amp2Minus3() + hampl_KNL_minus.Amp2Minus1());
          sigS_minus = (hampl_KNL_minus.Amp20Plus () + hampl_KNL_minus.Amp20Minus());


          fFKR.S = KNL_S_plus;
          fFKR.B = KNL_B_plus;
          fFKR.S = KNL_C_plus;
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
            fFKR.S = BRS_C_minus;

            hamplmod_BRS_minus = fHAmplModelCC;
            assert(hamplmod_BRS_minus);

            const RSHelicityAmpl & hampl_BRS_minus = hamplmod_BRS_minus->Compute(resonance, fFKR);

            sigL_minus = (hampl_BRS_minus.Amp2Plus3 () + hampl_BRS_minus.Amp2Plus1 ());
            sigR_minus = (hampl_BRS_minus.Amp2Minus3() + hampl_BRS_minus.Amp2Minus1());
            sigS_minus = (hampl_BRS_minus.Amp20Plus () + hampl_BRS_minus.Amp20Minus());

            fFKR.S = BRS_S_plus;
            fFKR.B = BRS_B_plus;
            fFKR.S = BRS_C_plus;
            hamplmod_BRS_plus = fHAmplModelCC;
            assert(hamplmod_BRS_plus);

            const RSHelicityAmpl & hampl_BRS_plus = hamplmod_BRS_plus->Compute(resonance, fFKR);

            sigL_plus = (hampl_BRS_plus.Amp2Plus3 () + hampl_BRS_plus.Amp2Plus1 ());
            sigR_plus = (hampl_BRS_plus.Amp2Minus3() + hampl_BRS_plus.Amp2Minus1());
            sigS_plus = (hampl_BRS_plus.Amp20Plus () + hampl_BRS_plus.Amp20Minus());
          }
  /*
     if(is_CC) { 
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
     assert(hamplmod);

     const RSHelicityAmpl & hampl = hamplmod->Compute(resonance, fFKR); 
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("RSHAmpl", pDEBUG)
    << "Helicity Amplitudes for RES = " << resname << " : " << hampl;
#endif
     */

  //JNtest double g2 = kGF2;
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

// Compute the cross section

/*
   double sig0 = 0.125*(g2/kPi)*(-q2/Q2)*(W/Mnuc);
   double scLR = W/Mnuc;
   double scS  = (Mnuc/W)*(-Q2/q2);

   JN test */
//  double sigL = scLR* (hampl.Amp2Plus3 () + hampl.Amp2Plus1 ());
//  double sigR = scLR* (hampl.Amp2Minus3() + hampl.Amp2Minus1());
//  double sigS = scS * (hampl.Amp20Plus () + hampl.Amp20Minus());

/*
   double sigL =0;
   double sigR =0;
   double sigS =0;
   JN tests  */

if( is_KLN || is_BRS){

  sigL_minus *= scLR;
  sigR_minus *= scLR;
  sigS_minus *= scS;

  sigL_plus *= scLR;
  sigR_plus *= scLR;
  sigS_plus *= scS;
  LOG("BergerSehgalRESPXSec2014",pINFO) <<"sL,R,S minus="<<sigL_minus<<","<<sigR_minus<<","<<sigS_minus;
  LOG("BergerSehgalRESPXSec2014",pINFO) <<"sL,R,S plus ="<<sigL_plus <<","<<sigR_plus <<","<<sigS_plus;
}
else {
  assert(hamplmod);

  const RSHelicityAmpl & hampl = hamplmod->Compute(resonance, fFKR);

  sigL = scLR* (hampl.Amp2Plus3 () + hampl.Amp2Plus1 ());
  sigR = scLR* (hampl.Amp2Minus3() + hampl.Amp2Minus1());
  sigS = scS * (hampl.Amp20Plus () + hampl.Amp20Minus());
}

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
LOG("BergerSehgalRESPXSec2014", pDEBUG) << "sig_{0} = " << sig0;
LOG("BergerSehgalRESPXSec2014", pDEBUG) << "sig_{L} = " << sigL;
LOG("BergerSehgalRESPXSec2014", pDEBUG) << "sig_{R} = " << sigR;
LOG("BergerSehgalRESPXSec2014", pDEBUG) << "sig_{S} = " << sigS;
#endif

double xsec = 0.0;

if(is_KLN || is_BRS){
  xsec =  TMath::Power(KNL_cL_minus,2)*sigL_minus + TMath::Power(KNL_cL_plus,2)*sigL_plus
    + TMath::Power(KNL_cR_minus,2)*sigR_minus + TMath::Power(KNL_cR_plus,2)*sigR_plus
    + TMath::Power(KNL_cS_minus,2)*sigS_minus + TMath::Power(KNL_cS_plus,2)*sigS_plus;
  xsec *=sig0;

     LOG("BergerSehgalRESPXSec2014",pINFO) <<"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n";

     LOG("BergerSehgalRESPXSec2014",pINFO) <<"A-="<<KNL_Alambda_minus<<" A+="<<KNL_Alambda_plus;
     LOG("BergerSehgalRESPXSec2014",pINFO) <<q2<<"\t"<<xsec<<"\t"<<sig0*(V2*sigR + U2*sigL + 2*UV*sigS)<<"\t"<<xsec/(sig0*(V2*sigRSR + U2*sigRSL + 2*UV*sigRSS));
     LOG("BergerSehgalRESPXSec2014",pINFO) <<"fFKR.B="<<fFKR.B<<" fFKR.C="<<fFKR.C<<" fFKR.S="<<fFKR.S;
     LOG("BergerSehgalRESPXSec2014",pINFO) <<"CL-="<<TMath::Power(KNL_cL_minus,2)<<" CL+="<<TMath::Power(KNL_cL_plus,2)<<" U2="<<U2;
     LOG("BergerSehgalRESPXSec2014",pINFO) <<"SL-="<<sigL_minus<<" SL+="<<sigL_plus<<" SL="<<sigRSL;

     LOG("BergerSehgalRESPXSec2014",pINFO) <<"CR-="<<TMath::Power(KNL_cR_minus,2)<<" CR+="<<TMath::Power(KNL_cR_plus,2)<<" V2="<<V2;
     LOG("BergerSehgalRESPXSec2014",pINFO) <<"SR-="<<sigR_minus<<" SR+="<<sigR_plus<<" sR="<<sigRSR;

     LOG("BergerSehgalRESPXSec2014",pINFO) <<"CS-="<<TMath::Power(KNL_cS_minus,2)<<" CS+="<<TMath::Power(KNL_cS_plus,2)<<" UV="<<UV;
     LOG("BergerSehgalRESPXSec2014",pINFO) <<"SS-="<<sigL_minus<<" SS+="<<sigS_plus<<" sS="<<sigRSS;
     LOG("BergerSehgalRESPXSec2014",pINFO) <<">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";

}
else{
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
if(is_CC && is_delta) {
  if((is_nu && is_p) || (is_nubar && is_n)) mult=3.0;
}
xsec *= mult;

// Check whether the cross section is to be weighted with a
// Breit-Wigner distribution (default: true)
double bw = 1.0;
if(fWghtBW) {
  bw = utils::bwfunc::BreitWignerL(W,LR,MR,WR,NR); 
} 
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
LOG("BergerSehgalRESPXSec2014", pDEBUG) << "BreitWigner(RES=" 
  << resname << ", W=" << W << ") = " << bw;
#endif
xsec *= bw; 

// Apply NeuGEN nutau cross section reduction factors
double rf = 1.0;
Spline * spl = 0;
if (is_CC && fUsingNuTauScaling) {
  if (pdg::IsNuTau(probepdgc)) {
    spl = fNuTauRdSpl;
  }
  else if (pdg::IsAntiNuTau(probepdgc)) {
    spl = fNuTauBarRdSpl;
  }
  if(spl) {
    if(E <spl->XMax()) rf = spl->Evaluate(E);
  }
}
xsec *= rf;

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
LOG("BergerSehgalRESPXSec2014", pINFO) 
  << "\n d2xsec/dQ2dW"  << "[" << interaction->AsString()
  << "](W=" << W << ", q2=" << q2 << ", E=" << E << ") = " << xsec;
#endif

// The algorithm computes d^2xsec/dWdQ2
// Check whether variable tranformation is needed
if(kps!=kPSWQ2fE) {
  double J = utils::kinematics::Jacobian(interaction,kPSWQ2fE,kps);
  xsec *= J;
}

// If requested return the free nucleon xsec even for input nuclear tgt
if( interaction->TestBit(kIAssumeFreeNucleon) ) return xsec;

// Take into account the number of scattering centers in the target
int NNucl = (is_p) ? target.Z() : target.N();

xsec*=NNucl; // nuclear xsec (no nuclear suppression factor)

return xsec;
}
//____________________________________________________________________________
double BergerSehgalRESPXSec2014::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool BergerSehgalRESPXSec2014::ValidProcess(const Interaction * interaction) const
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
void BergerSehgalRESPXSec2014::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void BergerSehgalRESPXSec2014::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void BergerSehgalRESPXSec2014::LoadConfig(void)
{
  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();

  // Load all configuration data or set defaults

  LOG("RSHAmpl", pWARN)
    << "get to beg = ";
  fZeta  = fConfig->GetDoubleDef( "Zeta",  gc->GetDouble("RS-Zeta")  );
  LOG("RSHAmpl", pWARN)
    << "load fZeta";
  fOmega = fConfig->GetDoubleDef( "Omega", gc->GetDouble("RS-Omega") );
  LOG("RSHAmpl", pWARN)
    << "load Omega";
  fKLN = fConfig->GetBoolDef("is_KLN", gc->GetBool("is_KLN"));
  LOG("RSHAmpl", pWARN)
    << "load is_KLN";
  fBRS = fConfig->GetBoolDef("is_BRS", gc->GetBool("is_BRS"));
  fGA  = fConfig->GetBoolDef("minibooneGA", gc->GetBool("minibooneGA"));
  fGV  = fConfig->GetBoolDef("minibooneGV", gc->GetBool("minibooneGV"));

  double ma  = fConfig->GetDoubleDef( "Ma", gc->GetDouble("RES-Ma") );
  double mv  = fConfig->GetDoubleDef( "Mv", gc->GetDouble("RES-Mv") );

  fMa2 = TMath::Power(ma,2);
  fMv2 = TMath::Power(mv,2);

  fWghtBW = fConfig->GetBoolDef("BreitWignerWeight", true);

  double thw = fConfig->GetDoubleDef(
      "weinberg-angle", gc->GetDouble("WeinbergAngle"));

  fSin48w = TMath::Power( TMath::Sin(thw), 4 );

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
  fUsingDisResJoin = fConfig->GetBoolDef(
      "UseDRJoinScheme", gc->GetBool("UseDRJoinScheme"));
  fWcut = 999999;
  if(fUsingDisResJoin) {
    fWcut = fConfig->GetDoubleDef("Wcut",gc->GetDouble("Wcut"));
  }

  // NeuGEN limits in the allowed resonance phase space:
  // W < min{ Wmin(physical), (res mass) + x * (res width) }
  // It limits the integration area around the peak and avoids the
  // problem with huge xsec increase at low Q2 and high W.
  // In correspondence with Hugh, Rein said that the underlying problem
  // are unphysical assumptions in the model. 
  fN2ResMaxNWidths = fConfig->GetDoubleDef("MaxNWidthForN2Res", 2.0);
  fN0ResMaxNWidths = fConfig->GetDoubleDef("MaxNWidthForN0Res", 6.0);
  fGnResMaxNWidths = fConfig->GetDoubleDef("MaxNWidthForGNRes", 4.0);

  // NeuGEN reduction factors for nu_tau: a gross estimate of the effect of
  // neglected form factors in the R/S model
  fUsingNuTauScaling = fConfig->GetBoolDef("UseNuTauScalingFactors", true);
  if(fUsingNuTauScaling) {
    if(fNuTauRdSpl)    delete fNuTauRdSpl;
    if(fNuTauBarRdSpl) delete fNuTauBarRdSpl;

    assert(gSystem->Getenv("GENIE"));
    string base = gSystem->Getenv("GENIE");

    string filename = base + "/data/evgen/rein_sehgal/res/nutau_xsec_scaling_factors.dat";
    LOG("BergerSehgalRESPXSec2014", pINFO) 
      << "Loading nu_tau xsec reduction spline from: " << filename;
    fNuTauRdSpl = new Spline(filename);

    filename = base + "/data/evgen/rein_sehgal/res/nutaubar_xsec_scaling_factors.dat";
    LOG("BergerSehgalRESPXSec2014", pINFO) 
      << "Loading bar{nu_tau} xsec reduction spline from: " << filename;
    fNuTauBarRdSpl = new Spline(filename);
  }

  // load the differential cross section integrator
  fXSecIntegrator =
    dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);
}
//____________________________________________________________________________
