//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Authors: Igor Kakorin <kakorin@jinr.ru>, Joint Institute for Nuclear Research
          Vadim Naumov <vnaumov@theor.jinr.ru >,  Joint Institute for Nuclear Research \n
          adapted from code provided by
          Minoo Kabirnezhad <minoo.kabirnezhad@physics.ox.ac.uk>
          University of Oxford, Department of Physics
          based on code of Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool


  For the class documentation see the corresponding header file.
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
#include <Framework/Conventions/KinePhaseSpace.h>
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Interaction/SppChannel.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Utils/Range1.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Physics/Resonance/XSection/MKSPPPXSec2020.h"
#include "Physics/NuclearState/FermiMomentumTablePool.h"
#include "Physics/NuclearState/FermiMomentumTable.h"
#include "Physics/NuclearState/NuclearUtils.h"

#include <algorithm>


using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
MKSPPPXSec2020::MKSPPPXSec2020() :
XSecAlgorithmI("genie::MKSPPPXSec2020")
{

}
//____________________________________________________________________________
MKSPPPXSec2020::MKSPPPXSec2020(string config) :
XSecAlgorithmI("genie::MKSPPPXSec2020", config)
{

}
//____________________________________________________________________________
MKSPPPXSec2020::~MKSPPPXSec2020()
{

}
//____________________________________________________________________________
double MKSPPPXSec2020::XSec(const Interaction * interaction, KinePhaseSpace_t kps) const
{
  if(! this -> ValidProcess    (interaction) ) return 0;
  if(! this -> ValidKinematics (interaction) ) return 0;

  const InitialState & init_state = interaction -> InitState();
  const ProcessInfo &  proc_info  = interaction -> ProcInfo();
  const Target & target = init_state.Tgt();

  double xsec = 0;
  //-- Get 1pi exclusive channel
  SppChannel_t spp_channel = SppChannel::FromInteraction(interaction);

  int  nucpdgc   = target.HitNucPdg();
  int  probepdgc = init_state.ProbePdg();
  bool is_nu     = pdg::IsNeutrino      (probepdgc);
  bool is_nubar  = pdg::IsAntiNeutrino  (probepdgc);
  bool is_CC     = proc_info.IsWeakCC();
  bool is_NC     = proc_info.IsWeakNC();
  bool is_p      = pdg::IsProton  (nucpdgc);


  // Get kinematical parameters
  const Kinematics & kinematics = interaction -> Kine();
  double E   = init_state.ProbeE(kRfHitNucRest);
  double ml  = interaction->FSPrimLepton()->Mass();
  double ml2 = ml*ml;
  double Q2  = kinematics.Q2();
  double W   = kinematics.W();
  double W2  = W*W;
  double Wt2 = W*2;


  // dimension of kine phase space
  std::string s = KinePhaseSpace::AsString(kps);
  int kpsdim = s!="<|E>"?1 + std::count(s.begin(), s.begin()+s.find('}'), ','):0;
  if (kpsdim < 3 || kpsdim > 4) return 0;
  // Pion angles should be given in Adler frame
  double CosTheta = kinematics.GetKV(kKVctp);
  double Phi = 0;
  if (kpsdim == 4)
    Phi = kinematics.GetKV(kKVphip);

  double SinTheta   = TMath::Sqrt(1 - CosTheta*CosTheta);
  double SinHalfTheta = TMath::Sqrt((1 - CosTheta)/2);
  double CosHalfTheta = TMath::Sqrt((1 + CosTheta)/2);

  PDGLibrary * pdglib = PDGLibrary::Instance();

  // imply isospin symmetry
  double m_pi  = (pdglib->Find(kPdgPiP)->Mass() + pdglib->Find(kPdgPi0)->Mass() + pdglib->Find(kPdgPiM)->Mass())/3;
  double m_pi2 = m_pi*m_pi;

  double M = (pdglib->Find(kPdgProton)->Mass() + pdglib->Find(kPdgNeutron)->Mass())/2;
  double M2  = M*M;
  double Mt2 = M*2;

  double q_0               = (W2 - M2 + m_pi2)/Wt2;
  double q_02              = q_0*q_0;
  double k_0               = (W2 - M2 - Q2)/Wt2;
  double abs_mom_q         = TMath::Sqrt(TMath::Max(0., q_02 - m_pi2));
  double abs_mom_k         = TMath::Sqrt(k_0*k_0 + Q2);
  //double E_2L              = (M2 - W2 - Q2 + Mt2*E)/Mt2;
  double abs_mom_k_L       = W*abs_mom_k/M;
  double abs_mom_k_L2      = abs_mom_k_L*abs_mom_k_L;
  double k_2               = (M2 + Mt2*E - W2 - ml2)/Wt2;
  double k_1               =  k_2 + k_0;
  double p_10              = (W2 + M2 + Q2)/Wt2;
  double qk                = q_0*k_0 - abs_mom_k*abs_mom_q*CosTheta;
  //double k_2L              = TMath::Sqrt(E_2L*E_2L - ml2);               //magnitude of lepton momentum in lab frame
  
  // There is no check of positivity under root expression in the original code
  double k_2_iso           = TMath::Sqrt(TMath::Max(0., k_2*k_2 - ml2));   //magnitude of lepton momentum in isobaric frame
  double cos_theta         = k_2_iso>0?(2*k_1*k_2 - Q2 - ml2)/2/k_1/k_2_iso:0;
  // There are no checks presented below in the original code
  if (cos_theta > 1)
    cos_theta = 1;
  if (cos_theta < -1)
    cos_theta = -1;

  // Eq. 7 of ref. 1
  double A_plus            = TMath::Sqrt( k_1*(k_2 - k_2_iso) );
  double A_minus           = TMath::Sqrt( k_1*(k_2 + k_2_iso) );
  //Eq. 6 of ref. 1
  double eps_1_plus        =   2*A_plus *(k_1 - k_2_iso)/abs_mom_k*TMath::Sqrt(1 + cos_theta);
  double eps_1_minus       =  -2*A_minus*(k_1 + k_2_iso)/abs_mom_k*TMath::Sqrt(1 - cos_theta);
  double eps_2_plus        =   2*A_plus *TMath::Sqrt(1 + cos_theta);
  double eps_2_minus       =   2*A_minus*TMath::Sqrt(1 - cos_theta);
  //Eq. 9 of ref. 1
  double eps_zero_L        =  -2*A_minus*TMath::Sqrt(1 + cos_theta);       // L->lambda = -1
  double eps_zero_R        =   2*A_plus *TMath::Sqrt(1 - cos_theta);       // R->lambda = +1
  double eps_z_L           =  -2*A_minus*(k_1 - k_2_iso)/abs_mom_k*TMath::Sqrt(1 + cos_theta);
  double eps_z_R           =   2*A_plus *(k_1 + k_2_iso)/abs_mom_k*TMath::Sqrt(1 - cos_theta);
  
  if (is_nubar)
  {
     if (fUseAuthorCode)
     {
         //This ``recipe'' of transition from neutrino to antineutrino case is promoted by Minoo Kabirnezhad.
         //However, it is not correct in our opinion. All details can be found in Ref. [12], see
         //section "Problem with transition from neutrino to antineutrino case".
         Phi = -Phi;
         eps_1_plus        =  -2*A_minus *(k_1 + k_2_iso)/abs_mom_k*TMath::Sqrt(1 - cos_theta);
         eps_1_minus       =  -2*A_plus*(k_1 - k_2_iso)/abs_mom_k*TMath::Sqrt(1 + cos_theta);
         eps_2_plus        =  -2*A_minus *TMath::Sqrt(1 - cos_theta);
         eps_2_minus       =   2*A_plus*TMath::Sqrt(1 + cos_theta);
         eps_zero_L        =   2*A_plus*TMath::Sqrt(1 - cos_theta);       // L->lambda = -1
         eps_zero_R        =   2*A_minus *TMath::Sqrt(1 + cos_theta);       // R->lambda = +1
         eps_z_L           =   2*A_plus*(k_1 + k_2_iso)/abs_mom_k*TMath::Sqrt(1 - cos_theta);
         eps_z_R           =   2*A_minus *(k_1 - k_2_iso)/abs_mom_k*TMath::Sqrt(1 + cos_theta);
     }
     else
     {
         eps_1_plus        =   2*A_minus *(k_1 + k_2_iso)/abs_mom_k*TMath::Sqrt(1 - cos_theta);
         eps_1_minus       =   2*A_plus*(k_1 - k_2_iso)/abs_mom_k*TMath::Sqrt(1 + cos_theta);
         eps_2_plus        =   2*A_minus *TMath::Sqrt(1 - cos_theta);
         eps_2_minus       =  -2*A_plus*TMath::Sqrt(1 + cos_theta);
         eps_zero_L        =   2*A_plus*TMath::Sqrt(1 - cos_theta);       // L->lambda = -1
         eps_zero_R        =   2*A_minus *TMath::Sqrt(1 + cos_theta);     // R->lambda = +1
         eps_z_L           =   2*A_plus*(k_1 + k_2_iso)/abs_mom_k*TMath::Sqrt(1 - cos_theta);
         eps_z_R           =   2*A_minus *(k_1 - k_2_iso)/abs_mom_k*TMath::Sqrt(1 + cos_theta);
     }
  }
  //Eq. 10 of ref. 1
  double C_L_plus          =  k1_Sqrt2*(eps_1_plus  - eps_2_plus);
  double C_L_minus         =  k1_Sqrt2*(eps_1_minus - eps_2_minus);
  double C_R_plus          = -k1_Sqrt2*(eps_1_plus  + eps_2_plus);
  double C_R_minus         = -k1_Sqrt2*(eps_1_minus + eps_2_minus);
  double C_S_plus          =  TMath::Sqrt(TMath::Abs(eps_zero_R*eps_zero_R -  eps_z_R*eps_z_R));
  double C_S_minus         =  TMath::Sqrt(TMath::Abs(eps_zero_L*eps_zero_L -  eps_z_L*eps_z_L));

  // if C_S_plus or C_S_minus are equal to zero, then the singularity appers
  // in the expressions for Hbkg and Hres at BosonPolarization=PLUS0 or BosonPolarization=MINUS0
  // which further eliminated by corresponding factors in the sum for calculating of the cross section.
  // Therefore we can just set them equal to 1 in this case. In our opinion it is simplier to eliminate
  // singularity in such way, because it allows to keep intact original formulas from paper.
  C_S_plus  = C_S_plus==0?1:C_S_plus;
  C_S_minus = C_S_minus==0?1:C_S_minus;


  /*    Bkg Contribution   */
  //Auxiliary functions
  double W_plus            = W + M;
  double W_minus           = W - M;
  double W_plus2           = W_plus*W_plus;
  double W_minus2          = W_minus*W_minus;
  double s_minus_M2        = W_plus*W_minus;
  double u_minus_M2        = m_pi2 - 2*(q_0*p_10 + abs_mom_q*abs_mom_k*CosTheta);
  double t_minus_mpi2      = -(Q2 + 2*qk);
  double mpi2_minus_k2     = Q2 + m_pi2;
  double Q                 = TMath::Sqrt(Q2);

  // Helicity amplitudes for bkg
  HelicityBkgAmp<std::complex<double>> Hbkg;

  //Vector_helicity amplituds for bkg
  // vector cut
  double FV_cut = 0;
  //virtual form symmetry_factor to kill the bkg smothly
  if (W < fBkgVWmin)
    FV_cut = 1;
  else if (W >= fBkgVWmin && W < fBkgVWmax)
    FV_cut = fBkgV0 + W*(fBkgV1 + W*(fBkgV2 + W*fBkgV3));

  fFormFactors.Calculate(interaction);
  double g_A     = -FA0;
  double FF_pi   = fFormFactors.F1V();
  double FF_1    = 0.5*FF_pi;
  double FF_2    = 0.5*fFormFactors.xiF2V();
  


  // Eq. 38 and Table 5 of ref. 1
  double V_1 = 0, V_2 = 0, V_3 = 0, V_4 = 0, V_5 = 0;

  // [p,n,pi+,pi0,pi-]
  switch (spp_channel) {

    //CC
    case (kSpp_vp_cc_10100) :
    case (kSpp_vbn_cc_01001):
        V_1  = g_A*kSqrt2*FV_cut/f_pi*(Mt2*FF_1/u_minus_M2 + FF_2/M);
        V_2  = g_A*kSqrt2*FV_cut/f_pi*Mt2*FF_1/u_minus_M2/qk;
        V_3  = g_A*kSqrt2*FV_cut/f_pi*2*FF_2/u_minus_M2;
        V_4  = -V_3;
        V_5  = g_A*kSqrt2*FV_cut/f_pi*M*FF_pi/t_minus_mpi2/qk;
        break;

    case (kSpp_vn_cc_01100) :
    case (kSpp_vbp_cc_10001):
        V_1 =  g_A*kSqrt2*FV_cut/f_pi*(Mt2*FF_1/s_minus_M2 + FF_2/M);
        V_2 =  g_A*kSqrt2*FV_cut/f_pi*Mt2*FF_1/s_minus_M2/qk;
        V_3 = -g_A*kSqrt2*FV_cut/f_pi*2*FF_2/s_minus_M2;
        V_4 =  V_3;
        V_5 = -g_A*kSqrt2*FV_cut/f_pi*M*FF_pi/t_minus_mpi2/qk;
        break;

    case (kSpp_vn_cc_10010) :
    case (kSpp_vbp_cc_01010):
        V_1 =  g_A*FV_cut/f_pi*Mt2*FF_1*(1/s_minus_M2 - 1/u_minus_M2);
        V_2 =  g_A*FV_cut/f_pi*Mt2*FF_1/qk*(1/s_minus_M2 - 1/u_minus_M2);
        V_3 = -g_A*FV_cut/f_pi*2*FF_2*(1/s_minus_M2 + 1/u_minus_M2);
        V_4 = -g_A*FV_cut/f_pi*2*FF_2*(1/s_minus_M2 - 1/u_minus_M2);;
        V_5 = -g_A*FV_cut/f_pi*Mt2*FF_pi/t_minus_mpi2/qk;
        break;

    //NC
    case (kSpp_vp_nc_10010) :
    case (kSpp_vbp_nc_10010):
    case (kSpp_vn_nc_01010) :
    case (kSpp_vbn_nc_01010):
        V_1 =  g_A*FV_cut*M/f_pi*(FF_1*(1/u_minus_M2 + 1/s_minus_M2) + FF_2/M2);
        V_2 =  g_A*FV_cut*M/f_pi*FF_1*(1/u_minus_M2 + 1/s_minus_M2)/qk;
        V_3 =  g_A*FV_cut/2/f_pi*FF_2*(1/u_minus_M2 - 1/s_minus_M2);
        V_4 = -g_A*FV_cut/2/f_pi*FF_2*(1/u_minus_M2 + 1/s_minus_M2);
        break;

    case (kSpp_vp_nc_01100) :
    case (kSpp_vbp_nc_01100):
        V_1 = -k1_Sqrt2*g_A*FV_cut*Mt2/f_pi*FF_1*(1/u_minus_M2 - 1/s_minus_M2);
        V_2 = -k1_Sqrt2*g_A*FV_cut*Mt2/f_pi*FF_1*(1/u_minus_M2 - 1/s_minus_M2)/qk;
        V_3 = -k1_Sqrt2*g_A*FV_cut/f_pi*FF_2*(1/u_minus_M2 + 1/s_minus_M2);
        V_4 =  k1_Sqrt2*g_A*FV_cut/f_pi*FF_2*(1/u_minus_M2 - 1/s_minus_M2);
        V_5 = -k1_Sqrt2*g_A*FV_cut*Mt2/f_pi*FF_pi/t_minus_mpi2/qk;
        break;

    case (kSpp_vn_nc_10001) :
    case (kSpp_vbn_nc_10001):
        V_1 =  k1_Sqrt2*g_A*FV_cut*Mt2/f_pi*FF_1*(1/u_minus_M2 - 1/s_minus_M2);
        V_2 =  k1_Sqrt2*g_A*FV_cut*Mt2/f_pi*FF_1*(1/u_minus_M2 - 1/s_minus_M2)/qk;
        V_3 =  k1_Sqrt2*g_A*FV_cut/f_pi*FF_2*(1/u_minus_M2 + 1/s_minus_M2);
        V_4 = -k1_Sqrt2*g_A*FV_cut/f_pi*FF_2*(1/u_minus_M2 - 1/s_minus_M2);
        V_5 =  k1_Sqrt2*g_A*FV_cut*Mt2/f_pi*FF_pi/t_minus_mpi2/qk;
        break;
    default:
                 // should not be here - meaningless to return anything
                 gAbortingInErr = true;
                 LOG("MKSPPPXSec2020", pFATAL) << "Unknown resonance production channel: " << spp_channel;
                 exit(1);


  }
  double V_25 = (W_plus*W_minus)*V_2 - Q2*V_5;

  
  double O_1_plus          = TMath::Sqrt((W_plus2 + Q2)*( W_plus2 - m_pi2 ))/Wt2;
  // There is no check of positivity under root expression in the original code
  double O_1_minus         = TMath::Sqrt(TMath::Max(0., (W_minus2 + Q2)*(W_minus2 - m_pi2)))/Wt2;
  double O_2_plus          = TMath::Sqrt((W_plus2 + Q2)/(W_plus2 - m_pi2));
  // There is no check of positivity under root expression in the original code
  double O_2_minus         = W_minus2>m_pi2?TMath::Sqrt((W_minus2 + Q2)/(W_minus2 - m_pi2)):0;
  
  double K_1_V = W_minus*O_1_plus;
  double K_2_V = W_plus*O_1_minus;
  // There is no check of positivity under root expression in the original code
  double K_3_V = W_minus2>m_pi2?abs_mom_q*abs_mom_q*W_plus*O_2_minus:0;
  double K_4_V = abs_mom_q*abs_mom_q*W_minus*O_2_plus;
  double K_5_V = 1/O_2_plus;
  // the of K_6_V = 1/O_2_minus differ from original code
  double K_6_V = TMath::Sqrt(TMath::Max(0., (W_minus2 - m_pi2))/(W_minus2 + Q2));

  double F_1 =   V_1 + qk*(V_3 - V_4)/W_minus + W_minus*V_4;
  double F_2 =  -V_1 + qk*(V_3 - V_4)/W_plus  + W_plus *V_4;
  double F_3 =  (V_3 - V_4) + V_25/W_plus;
  double F_4 =  (V_3 - V_4) - V_25/W_minus;
  double F_5 = (W_plus2  + Q2)*(V_1 + W_minus*V_4)/Wt2 + (W_plus*q_0  - qk)*(V_3 - V_4) + q_0*V_25 - k_0*qk*V_5 -
                                                                         qk*(W_plus2 + Q2 + 2*W*W_minus)*V_2/Wt2;
  double F_6 =-(W_minus2 + Q2)*(V_1 - W_plus *V_4)/Wt2 + (W_minus*q_0 - qk)*(V_3 - V_4) - q_0*V_25 + k_0*qk*V_5 +
                                                                         qk*(W_plus2 + Q2 + 2*W*W_minus)*V_2/Wt2;

  double sF_1 = K_1_V*F_1/Mt2;
  double sF_2 = K_2_V*F_2/Mt2;
  double sF_3 = K_3_V*F_3/Mt2;
  double sF_4 = K_4_V*F_4/Mt2;
  double sF_5 = K_5_V*F_5/Mt2;
  double sF_6 = K_6_V*F_6/Mt2;

  Hbkg(Current::VECTOR,          BosonPolarization::LEFT, NucleonPolarization::PLUS,  NucleonPolarization::PLUS)  =
      k1_Sqrt2*SinTheta*CosHalfTheta*(sF_3 + sF_4)*std::complex<double>
     (TMath::Cos(Phi*PhaseFactor(BosonPolarization::LEFT, NucleonPolarization::PLUS,  NucleonPolarization::PLUS)),
      TMath::Sin(Phi*PhaseFactor(BosonPolarization::LEFT, NucleonPolarization::PLUS,  NucleonPolarization::PLUS)));

  Hbkg(Current::VECTOR,          BosonPolarization::LEFT,  NucleonPolarization::MINUS, NucleonPolarization::PLUS)  =
     -k1_Sqrt2*SinTheta*SinHalfTheta*(sF_3 - sF_4)*std::complex<double>
     (TMath::Cos(Phi*PhaseFactor(BosonPolarization::LEFT,  NucleonPolarization::MINUS, NucleonPolarization::PLUS)),
      TMath::Sin(Phi*PhaseFactor(BosonPolarization::LEFT,  NucleonPolarization::MINUS, NucleonPolarization::PLUS)));

  Hbkg(Current::VECTOR,          BosonPolarization::LEFT,  NucleonPolarization::PLUS,  NucleonPolarization::MINUS) =
      kSqrt2*(CosHalfTheta*(sF_1 - sF_2) - 0.5*SinTheta*SinHalfTheta*(sF_3 - sF_4))*std::complex<double>
     (TMath::Cos(Phi*PhaseFactor(BosonPolarization::LEFT,  NucleonPolarization::PLUS,  NucleonPolarization::MINUS)),
      TMath::Sin(Phi*PhaseFactor(BosonPolarization::LEFT,  NucleonPolarization::PLUS,  NucleonPolarization::MINUS)));

  Hbkg(Current::VECTOR,          BosonPolarization::LEFT,  NucleonPolarization::MINUS, NucleonPolarization::MINUS) =
     -kSqrt2*(SinHalfTheta*(sF_1 + sF_2) + 0.5*SinTheta*CosHalfTheta*(sF_3 + sF_4))*std::complex<double>
     (TMath::Cos(Phi*PhaseFactor(BosonPolarization::LEFT,  NucleonPolarization::MINUS, NucleonPolarization::MINUS)),
      TMath::Sin(Phi*PhaseFactor(BosonPolarization::LEFT,  NucleonPolarization::MINUS, NucleonPolarization::MINUS)));


  Hbkg(Current::VECTOR,          BosonPolarization::RIGHT, NucleonPolarization::PLUS,  NucleonPolarization::PLUS)  =
     -kSqrt2*(SinHalfTheta*(sF_1 + sF_2) + 0.5*SinTheta*CosHalfTheta*(sF_3 + sF_4))*std::complex<double>
     (TMath::Cos(Phi*PhaseFactor(BosonPolarization::RIGHT, NucleonPolarization::PLUS,  NucleonPolarization::PLUS)),
      TMath::Sin(Phi*PhaseFactor(BosonPolarization::RIGHT, NucleonPolarization::PLUS,  NucleonPolarization::PLUS)));

  Hbkg(Current::VECTOR,          BosonPolarization::RIGHT, NucleonPolarization::MINUS, NucleonPolarization::PLUS)  =
     -kSqrt2*(CosHalfTheta*(sF_1 - sF_2) - 0.5*SinTheta*SinHalfTheta*(sF_3 - sF_4))*std::complex<double>
     (TMath::Cos(Phi*PhaseFactor(BosonPolarization::RIGHT, NucleonPolarization::MINUS, NucleonPolarization::PLUS)),
      TMath::Sin(Phi*PhaseFactor(BosonPolarization::RIGHT, NucleonPolarization::MINUS, NucleonPolarization::PLUS)));

  Hbkg(Current::VECTOR,          BosonPolarization::RIGHT, NucleonPolarization::PLUS,  NucleonPolarization::MINUS) =
      k1_Sqrt2*SinTheta*SinHalfTheta*(sF_3 - sF_4)*std::complex<double>
     (TMath::Cos(Phi*PhaseFactor(BosonPolarization::RIGHT, NucleonPolarization::PLUS,  NucleonPolarization::MINUS)),
      TMath::Sin(Phi*PhaseFactor(BosonPolarization::RIGHT, NucleonPolarization::PLUS,  NucleonPolarization::MINUS)));

  Hbkg(Current::VECTOR,          BosonPolarization::RIGHT, NucleonPolarization::MINUS, NucleonPolarization::MINUS) =
      k1_Sqrt2*SinTheta*CosHalfTheta*(sF_3 + sF_4)*std::complex<double>
     (TMath::Cos(Phi*PhaseFactor(BosonPolarization::RIGHT, NucleonPolarization::MINUS, NucleonPolarization::MINUS)),
      TMath::Sin(Phi*PhaseFactor(BosonPolarization::RIGHT, NucleonPolarization::MINUS, NucleonPolarization::MINUS)));

      
  if (is_CC || (is_NC && is_nu) )
  {
     Hbkg(Current::VECTOR,          BosonPolarization::MINUS0, NucleonPolarization::PLUS,  NucleonPolarization::PLUS)  =
         CosHalfTheta*(k_0*eps_z_L - abs_mom_k*eps_zero_L)*(sF_5 + sF_6)/C_S_minus*std::complex<double>
        (TMath::Cos(Phi*PhaseFactor(BosonPolarization::MINUS0, NucleonPolarization::PLUS,  NucleonPolarization::PLUS)),
         TMath::Sin(Phi*PhaseFactor(BosonPolarization::MINUS0, NucleonPolarization::PLUS,  NucleonPolarization::PLUS)));
     
     Hbkg(Current::VECTOR,          BosonPolarization::MINUS0, NucleonPolarization::MINUS, NucleonPolarization::PLUS)  =
        -SinHalfTheta*(k_0*eps_z_L - abs_mom_k*eps_zero_L)*(sF_5 - sF_6)/C_S_minus*std::complex<double>
        (TMath::Cos(Phi*PhaseFactor(BosonPolarization::MINUS0, NucleonPolarization::MINUS, NucleonPolarization::PLUS)),
         TMath::Sin(Phi*PhaseFactor(BosonPolarization::MINUS0, NucleonPolarization::MINUS, NucleonPolarization::PLUS)));
     
     Hbkg(Current::VECTOR,          BosonPolarization::MINUS0, NucleonPolarization::PLUS,  NucleonPolarization::MINUS) =
        -SinHalfTheta*(k_0*eps_z_L - abs_mom_k*eps_zero_L)*(sF_5 - sF_6)/C_S_minus*std::complex<double>
        (TMath::Cos(Phi*PhaseFactor(BosonPolarization::MINUS0, NucleonPolarization::PLUS,  NucleonPolarization::MINUS)),
         TMath::Sin(Phi*PhaseFactor(BosonPolarization::MINUS0, NucleonPolarization::PLUS,  NucleonPolarization::MINUS)));
     
     Hbkg(Current::VECTOR,          BosonPolarization::MINUS0, NucleonPolarization::MINUS, NucleonPolarization::MINUS) =
        -CosHalfTheta*(k_0*eps_z_L - abs_mom_k*eps_zero_L)*(sF_5 + sF_6)/C_S_minus*std::complex<double>
        (TMath::Cos(Phi*PhaseFactor(BosonPolarization::MINUS0, NucleonPolarization::MINUS, NucleonPolarization::MINUS)),
         TMath::Sin(Phi*PhaseFactor(BosonPolarization::MINUS0, NucleonPolarization::MINUS, NucleonPolarization::MINUS)));
  }
  if (is_CC || (is_NC && is_nubar) )
  {
     Hbkg(Current::VECTOR,          BosonPolarization::PLUS0,  NucleonPolarization::PLUS,  NucleonPolarization::PLUS)  =
         CosHalfTheta*(k_0*eps_z_R - abs_mom_k*eps_zero_R)*(sF_5 + sF_6)/C_S_plus*std::complex<double>
        (TMath::Cos(Phi*PhaseFactor(BosonPolarization::PLUS0,  NucleonPolarization::PLUS,  NucleonPolarization::PLUS)),
         TMath::Sin(Phi*PhaseFactor(BosonPolarization::PLUS0,  NucleonPolarization::PLUS,  NucleonPolarization::PLUS)));
  
     Hbkg(Current::VECTOR,          BosonPolarization::PLUS0,  NucleonPolarization::MINUS, NucleonPolarization::PLUS)  =
        -SinHalfTheta*(k_0*eps_z_R - abs_mom_k*eps_zero_R)*(sF_5 - sF_6)/C_S_plus*std::complex<double>
        (TMath::Cos(Phi*PhaseFactor(BosonPolarization::PLUS0,  NucleonPolarization::MINUS, NucleonPolarization::PLUS)),
         TMath::Sin(Phi*PhaseFactor(BosonPolarization::PLUS0,  NucleonPolarization::MINUS, NucleonPolarization::PLUS)));
  
     Hbkg(Current::VECTOR,          BosonPolarization::PLUS0,  NucleonPolarization::PLUS,  NucleonPolarization::MINUS) =
        -SinHalfTheta*(k_0*eps_z_R - abs_mom_k*eps_zero_R)*(sF_5 - sF_6)/C_S_plus*std::complex<double>
        (TMath::Cos(Phi*PhaseFactor(BosonPolarization::PLUS0,  NucleonPolarization::PLUS,  NucleonPolarization::MINUS)),
         TMath::Sin(Phi*PhaseFactor(BosonPolarization::PLUS0,  NucleonPolarization::PLUS,  NucleonPolarization::MINUS)));
  
     Hbkg(Current::VECTOR,          BosonPolarization::PLUS0,  NucleonPolarization::MINUS, NucleonPolarization::MINUS) =
        -CosHalfTheta*(k_0*eps_z_R - abs_mom_k*eps_zero_R)*(sF_5 + sF_6)/C_S_plus*std::complex<double>
        (TMath::Cos(Phi*PhaseFactor(BosonPolarization::PLUS0,  NucleonPolarization::MINUS, NucleonPolarization::MINUS)),
         TMath::Sin(Phi*PhaseFactor(BosonPolarization::PLUS0,  NucleonPolarization::MINUS, NucleonPolarization::MINUS)));
  }

  if (is_NC)
  {
    // isoscalar electromagnetic amplitudes eq.68 of ref. 8
    fEMFormFactors.Calculate(interaction);
    double FF_em_1 = 0.5*fEMFormFactors.F1V();
    double FF_em_2 = 0.5*fEMFormFactors.xiF2V();

    double Vem_1 = 0, Vem_2 = 0, Vem_3 = 0, Vem_4 = 0, Vem_5 = 0;

    // [p,n,pi+,pi0,pi-]
    switch (spp_channel) {
       //NC
       case (kSpp_vp_nc_10010) :
       case (kSpp_vbp_nc_10010):
           Vem_1 = -k1_Sqrt2*g_A*FV_cut*M/f_pi*(FF_em_1*(1/u_minus_M2 + 1/s_minus_M2) + FF_em_2/M2);
           Vem_2 = -k1_Sqrt2*g_A*FV_cut*M/f_pi*FF_em_1*(1/u_minus_M2  + 1/s_minus_M2)/qk;
           Vem_3 = -k1_Sqrt2*g_A*FV_cut/2/f_pi*FF_em_2*(1/u_minus_M2 - 1/s_minus_M2);
           Vem_4 =  k1_Sqrt2*g_A*FV_cut/2/f_pi*FF_em_2*(1/u_minus_M2 + 1/s_minus_M2);
           break;
       case (kSpp_vn_nc_01010) :
       case (kSpp_vbn_nc_01010):
           Vem_1 = k1_Sqrt2*g_A*FV_cut*M/f_pi*(FF_em_1*(1/u_minus_M2 + 1/s_minus_M2) + FF_em_2/M2);
           Vem_2 = k1_Sqrt2*g_A*FV_cut*M/f_pi*FF_em_1*(1/u_minus_M2  + 1/s_minus_M2)/qk;
           Vem_3 = k1_Sqrt2*g_A*FV_cut/2/f_pi*FF_em_2*(1/u_minus_M2 - 1/s_minus_M2);
           Vem_4 =-k1_Sqrt2*g_A*FV_cut/2/f_pi*FF_em_2*(1/u_minus_M2 + 1/s_minus_M2);
           break;

       case (kSpp_vp_nc_01100) :
       case (kSpp_vbp_nc_01100):
       case (kSpp_vn_nc_10001) :
       case (kSpp_vbn_nc_10001):
           Vem_1 = g_A*FV_cut*M/f_pi*(FF_em_1*(1/u_minus_M2 + 1/s_minus_M2) + FF_em_2/M2/2);
           Vem_2 = g_A*FV_cut*M/f_pi*FF_em_1*(1/u_minus_M2  + 1/s_minus_M2)/qk;
           Vem_3 = g_A*FV_cut/2/f_pi*FF_em_2*(1/u_minus_M2 - 1/s_minus_M2);
           Vem_4 =-g_A*FV_cut/2/f_pi*FF_em_2*(1/u_minus_M2 + 1/s_minus_M2);
           break;
       default:
             // should not be here - meaningless to return anything
             gAbortingInErr = true;
             LOG("MKSPPPXSec2020", pFATAL) << "CC channel in NC mode " << spp_channel;
             exit(1);
    }
    double Vem_25 = (W_plus*W_minus)*Vem_2 + Vem_5;

    double Fem_1 = Vem_1 + qk*(Vem_3 - Vem_4)/W_minus + W_minus*Vem_4;
    double Fem_2 = -Vem_1 + qk*(Vem_3 - Vem_4)/W_plus  + W_plus*Vem_4;
    double Fem_3 = (Vem_3 - Vem_4) + Vem_25/W_plus;
    double Fem_4 = (Vem_3 - Vem_4) - Vem_25/W_minus;
    double Fem_5 = (W_plus2  + Q2)*(Vem_1 + W_minus*Vem_4)/Wt2 + (W_plus*q_0  - qk)*(Vem_3 - Vem_4) + q_0*Vem_25 -
                                                                       qk*(W_plus2 + Q2 + 2*W*W_minus)*Vem_2/Wt2;
    double Fem_6 =-(W_minus2 + Q2)*(Vem_1 - W_plus*Vem_4)/Wt2 + (W_minus*q_0 - qk)*(Vem_3 - Vem_4) - q_0*Vem_25 +
                                                                       qk*(W_plus2 + Q2 + 2*W*W_minus)*Vem_2/Wt2;

    double sFem_1 = K_1_V*Fem_1/Mt2;
    double sFem_2 = K_2_V*Fem_2/Mt2;
    double sFem_3 = K_3_V*Fem_3/Mt2;
    double sFem_4 = K_4_V*Fem_4/Mt2;
    double sFem_5 = K_5_V*Fem_5/Mt2;
    double sFem_6 = K_6_V*Fem_6/Mt2;

    HelicityBkgAmp<std::complex<double>> HbkgEM;

    HbkgEM(Current::VECTOR,        BosonPolarization::LEFT,  NucleonPolarization::PLUS,  NucleonPolarization::PLUS)  =
        k1_Sqrt2*SinTheta*CosHalfTheta*(sFem_3 + sFem_4)*std::complex<double>
       (TMath::Cos(Phi*PhaseFactor(BosonPolarization::LEFT, NucleonPolarization::PLUS,  NucleonPolarization::PLUS)),
        TMath::Sin(Phi*PhaseFactor(BosonPolarization::LEFT, NucleonPolarization::PLUS,  NucleonPolarization::PLUS)));

    HbkgEM(Current::VECTOR,        BosonPolarization::LEFT,  NucleonPolarization::MINUS, NucleonPolarization::PLUS)  =
       -k1_Sqrt2*SinTheta*SinHalfTheta*(sFem_3 - sFem_4)*std::complex<double>
       (TMath::Cos(Phi*PhaseFactor(BosonPolarization::LEFT,  NucleonPolarization::MINUS, NucleonPolarization::PLUS)),
        TMath::Sin(Phi*PhaseFactor(BosonPolarization::LEFT,  NucleonPolarization::MINUS, NucleonPolarization::PLUS)));

    HbkgEM(Current::VECTOR,        BosonPolarization::LEFT,  NucleonPolarization::PLUS,  NucleonPolarization::MINUS) =
        kSqrt2*(CosHalfTheta*(sFem_1 - sFem_2) - 0.5*SinTheta*SinHalfTheta*(sFem_3 - sFem_4))*std::complex<double>
       (TMath::Cos(Phi*PhaseFactor(BosonPolarization::LEFT,  NucleonPolarization::PLUS,  NucleonPolarization::MINUS)),
        TMath::Sin(Phi*PhaseFactor(BosonPolarization::LEFT,  NucleonPolarization::PLUS,  NucleonPolarization::MINUS)));

    HbkgEM(Current::VECTOR,        BosonPolarization::LEFT,  NucleonPolarization::MINUS, NucleonPolarization::MINUS) =
       -kSqrt2*(SinHalfTheta*(sFem_1 + sFem_2) + 0.5*SinTheta*CosHalfTheta*(sFem_3 + sFem_4))*std::complex<double>
       (TMath::Cos(Phi*PhaseFactor(BosonPolarization::LEFT,  NucleonPolarization::MINUS, NucleonPolarization::MINUS)),
        TMath::Sin(Phi*PhaseFactor(BosonPolarization::LEFT,  NucleonPolarization::MINUS, NucleonPolarization::MINUS)));


    HbkgEM(Current::VECTOR,        BosonPolarization::RIGHT, NucleonPolarization::PLUS,  NucleonPolarization::PLUS)  =
       -kSqrt2*(SinHalfTheta*(sFem_1 + sFem_2) + 0.5*SinTheta*CosHalfTheta*(sFem_3 + sFem_4))*std::complex<double>
       (TMath::Cos(Phi*PhaseFactor(BosonPolarization::RIGHT, NucleonPolarization::PLUS,  NucleonPolarization::PLUS)),
        TMath::Sin(Phi*PhaseFactor(BosonPolarization::RIGHT, NucleonPolarization::PLUS,  NucleonPolarization::PLUS)));

    HbkgEM(Current::VECTOR,        BosonPolarization::RIGHT, NucleonPolarization::MINUS, NucleonPolarization::PLUS)  =
       -kSqrt2*(CosHalfTheta*(sFem_1 - sFem_2) - 0.5*SinTheta*SinHalfTheta*(sFem_3 - sFem_4))*std::complex<double>
       (TMath::Cos(Phi*PhaseFactor(BosonPolarization::RIGHT, NucleonPolarization::MINUS, NucleonPolarization::PLUS)),
        TMath::Sin(Phi*PhaseFactor(BosonPolarization::RIGHT, NucleonPolarization::MINUS, NucleonPolarization::PLUS)));

    HbkgEM(Current::VECTOR,        BosonPolarization::RIGHT, NucleonPolarization::PLUS,  NucleonPolarization::MINUS) =
        k1_Sqrt2*SinTheta*SinHalfTheta*(sFem_3 - sFem_4)*std::complex<double>
       (TMath::Cos(Phi*PhaseFactor(BosonPolarization::RIGHT, NucleonPolarization::PLUS,  NucleonPolarization::MINUS)),
        TMath::Sin(Phi*PhaseFactor(BosonPolarization::RIGHT, NucleonPolarization::PLUS,  NucleonPolarization::MINUS)));

    HbkgEM(Current::VECTOR,        BosonPolarization::RIGHT, NucleonPolarization::MINUS, NucleonPolarization::MINUS) =
        k1_Sqrt2*SinTheta*CosHalfTheta*(sFem_3 + sFem_4)*std::complex<double>
       (TMath::Cos(Phi*PhaseFactor(BosonPolarization::RIGHT, NucleonPolarization::MINUS, NucleonPolarization::MINUS)),
        TMath::Sin(Phi*PhaseFactor(BosonPolarization::RIGHT, NucleonPolarization::MINUS, NucleonPolarization::MINUS)));


    if (is_nu)
    {
      HbkgEM(Current::VECTOR,       BosonPolarization::MINUS0, NucleonPolarization::PLUS,  NucleonPolarization::PLUS)  =
          CosHalfTheta*(k_0*eps_z_L - abs_mom_k*eps_zero_L)*(sFem_5 + sFem_6)/C_S_minus*std::complex<double>
         (TMath::Cos(Phi*PhaseFactor(BosonPolarization::MINUS0, NucleonPolarization::PLUS,  NucleonPolarization::PLUS)),
          TMath::Sin(Phi*PhaseFactor(BosonPolarization::MINUS0, NucleonPolarization::PLUS,  NucleonPolarization::PLUS)));
  
      HbkgEM(Current::VECTOR,        BosonPolarization::MINUS0, NucleonPolarization::MINUS, NucleonPolarization::PLUS)  =
         -SinHalfTheta*(k_0*eps_z_L - abs_mom_k*eps_zero_L)*(sFem_5 - sFem_6)/C_S_minus*std::complex<double>
         (TMath::Cos(Phi*PhaseFactor(BosonPolarization::MINUS0, NucleonPolarization::MINUS, NucleonPolarization::PLUS)),
          TMath::Sin(Phi*PhaseFactor(BosonPolarization::MINUS0, NucleonPolarization::MINUS, NucleonPolarization::PLUS)));
  
      HbkgEM(Current::VECTOR,        BosonPolarization::MINUS0, NucleonPolarization::PLUS,  NucleonPolarization::MINUS) =
         -SinHalfTheta*(k_0*eps_z_L - abs_mom_k*eps_zero_L)*(sFem_5 - sFem_6)/C_S_minus*std::complex<double>
         (TMath::Cos(Phi*PhaseFactor(BosonPolarization::MINUS0, NucleonPolarization::PLUS,  NucleonPolarization::MINUS)),
          TMath::Sin(Phi*PhaseFactor(BosonPolarization::MINUS0, NucleonPolarization::PLUS,  NucleonPolarization::MINUS)));
  
      HbkgEM(Current::VECTOR,        BosonPolarization::MINUS0, NucleonPolarization::MINUS, NucleonPolarization::MINUS) =
         -CosHalfTheta*(k_0*eps_z_L - abs_mom_k*eps_zero_L)*(sFem_5 + sFem_6)/C_S_minus*std::complex<double>
         (TMath::Cos(Phi*PhaseFactor(BosonPolarization::MINUS0, NucleonPolarization::MINUS, NucleonPolarization::MINUS)),
          TMath::Sin(Phi*PhaseFactor(BosonPolarization::MINUS0, NucleonPolarization::MINUS, NucleonPolarization::MINUS)));
    }
    else
    {
      HbkgEM(Current::VECTOR,       BosonPolarization::PLUS0, NucleonPolarization::PLUS,  NucleonPolarization::PLUS)  =
          CosHalfTheta*(k_0*eps_z_R - abs_mom_k*eps_zero_R)*(sFem_5 + sFem_6)/C_S_plus*std::complex<double>
         (TMath::Cos(Phi*PhaseFactor(BosonPolarization::MINUS0, NucleonPolarization::PLUS,  NucleonPolarization::PLUS)),
          TMath::Sin(Phi*PhaseFactor(BosonPolarization::MINUS0, NucleonPolarization::PLUS,  NucleonPolarization::PLUS)));
  
      HbkgEM(Current::VECTOR,        BosonPolarization::PLUS0, NucleonPolarization::MINUS, NucleonPolarization::PLUS)  =
         -SinHalfTheta*(k_0*eps_z_R - abs_mom_k*eps_zero_R)*(sFem_5 - sFem_6)/C_S_plus*std::complex<double>
         (TMath::Cos(Phi*PhaseFactor(BosonPolarization::MINUS0, NucleonPolarization::MINUS, NucleonPolarization::PLUS)),
          TMath::Sin(Phi*PhaseFactor(BosonPolarization::MINUS0, NucleonPolarization::MINUS, NucleonPolarization::PLUS)));
  
      HbkgEM(Current::VECTOR,        BosonPolarization::PLUS0, NucleonPolarization::PLUS,  NucleonPolarization::MINUS) =
         -SinHalfTheta*(k_0*eps_z_R - abs_mom_k*eps_zero_R)*(sFem_5 - sFem_6)/C_S_plus*std::complex<double>
         (TMath::Cos(Phi*PhaseFactor(BosonPolarization::MINUS0, NucleonPolarization::PLUS,  NucleonPolarization::MINUS)),
          TMath::Sin(Phi*PhaseFactor(BosonPolarization::MINUS0, NucleonPolarization::PLUS,  NucleonPolarization::MINUS)));
  
      HbkgEM(Current::VECTOR,        BosonPolarization::PLUS0, NucleonPolarization::MINUS, NucleonPolarization::MINUS) =
         -CosHalfTheta*(k_0*eps_z_R - abs_mom_k*eps_zero_R)*(sFem_5 + sFem_6)/C_S_plus*std::complex<double>
         (TMath::Cos(Phi*PhaseFactor(BosonPolarization::MINUS0, NucleonPolarization::MINUS, NucleonPolarization::MINUS)),
          TMath::Sin(Phi*PhaseFactor(BosonPolarization::MINUS0, NucleonPolarization::MINUS, NucleonPolarization::MINUS)));
    }

    Hbkg = (1 - 2*fSin2Wein)*Hbkg - 4*fSin2Wein*HbkgEM;
  }

  //Axial helicity amp for bkg
  //Axial cut is same as FV_cut in the latest version
  double FA_cut = FV_cut;
  
  // FF_A here is equal to -g_A*FF_A of original code
  double FF_A    = fFormFactors.FA();
  double F_rho   = Frho0/(1 - (t_minus_mpi2 + m_pi2 )/fRho770Mass/fRho770Mass);
  
    // Eq. 38 and Table 5 of ref. 1
  double A_1 = 0, A_3 = 0, A_4 = 0, A_7 = 0, A_8 = 0;

  // [p,n,pi+,pi0,pi-]
  switch (spp_channel) {

    //CC
    case (kSpp_vp_cc_10100) :
    case (kSpp_vbn_cc_01001):
        A_1 = -kSqrt2*FA_cut*M*FF_A/f_pi/u_minus_M2;
        A_3 = -A_1;
        A_4 =  k1_Sqrt2*FA_cut/f_pi*(FF_A + F_rho)/M;
        A_7 =  kSqrt2*FA_cut/f_pi*M*FF_A/mpi2_minus_k2;
        A_8 = -k1_Sqrt2*FA_cut/f_pi*(FF_A*(1 +  4*M2/u_minus_M2) + F_rho)/mpi2_minus_k2;
        break;

    case (kSpp_vn_cc_01100) :
    case (kSpp_vbp_cc_10001):
        A_1 =  kSqrt2*FA_cut/f_pi*M*FF_A/s_minus_M2;
        A_3 =  A_1;
        A_4 = -k1_Sqrt2*FA_cut/f_pi/M*(FF_A + F_rho);
        A_7 =  kSqrt2*FA_cut/f_pi*M*FF_A/mpi2_minus_k2;
        A_8 =  k1_Sqrt2/f_pi*FA_cut*(FF_A/mpi2_minus_k2*(1 + 4*M2/s_minus_M2) + F_rho/mpi2_minus_k2);
        break;

    case (kSpp_vn_cc_10010) :
    case (kSpp_vbp_cc_01010):
        A_1 =  FA_cut/f_pi*M*FF_A*(1/u_minus_M2 + 1/s_minus_M2);
        A_3 = -FA_cut/f_pi*M*FF_A*(1/u_minus_M2 - 1/s_minus_M2) ;
        A_4 = -FA_cut/f_pi/M*(FF_A + F_rho);
        A_8 = FA_cut/f_pi*(F_rho + FF_A*(1 + 2*M2*(1/u_minus_M2 + 1/s_minus_M2)))/mpi2_minus_k2;
        break;

    //NC
    case (kSpp_vp_nc_10010) :
    case (kSpp_vbp_nc_10010):
    case (kSpp_vn_nc_01010) :
    case (kSpp_vbn_nc_01010):
        A_1 = -FA_cut/f_pi*M*FF_A*(1/u_minus_M2 - 1/s_minus_M2)/2;
        A_3 =  FA_cut/f_pi*M*FF_A*(1/u_minus_M2 + 1/s_minus_M2)/2;
        A_7 =  FA_cut/f_pi*M*FF_A/mpi2_minus_k2;
        break;

    case (kSpp_vp_nc_01100) :
    case (kSpp_vbp_nc_01100):
        A_1 =  k1_Sqrt2*FA_cut/f_pi*M*FF_A*(1/u_minus_M2 + 1/s_minus_M2);
        A_3 = -k1_Sqrt2*FA_cut/f_pi*M*FF_A*(1/u_minus_M2 - 1/s_minus_M2);
        A_4 = -k1_Sqrt2*FA_cut/f_pi/M*(FF_A + F_rho);
        A_8 =  k1_Sqrt2*FA_cut/f_pi/mpi2_minus_k2*(F_rho + FF_A*(1 + 2*M2*(1/u_minus_M2 + 1/s_minus_M2)));
        break;

    case (kSpp_vn_nc_10001) :
    case (kSpp_vbn_nc_10001):
        A_1 = -k1_Sqrt2*FA_cut/f_pi*M*FF_A*(1/u_minus_M2 + 1/s_minus_M2);
        A_3 =  k1_Sqrt2*FA_cut/f_pi*M*FF_A*(1/u_minus_M2 - 1/s_minus_M2);
        A_4 =  k1_Sqrt2*FA_cut/f_pi/M*(FF_A + F_rho);
        A_8 =  k1_Sqrt2*FA_cut/f_pi/mpi2_minus_k2*(FF_A*(1 +  2*M2*(1/u_minus_M2  + 1/s_minus_M2)) + F_rho);
        break;
    default:
        // should not be here - meaningless to return anything
        gAbortingInErr = true;
        LOG("MKSPPPXSec2020", pFATAL) << "Unknown resonance production channel: " << spp_channel;
        exit(1);
  }

  double K_1_A = abs_mom_q*O_2_plus;
  // There is no check of positivity under root expression in the original code
  double K_2_A = W_minus2>m_pi2?abs_mom_q*O_2_minus:1;
  double K_3_A = abs_mom_q*O_1_minus;
  double K_4_A = abs_mom_q* O_1_plus;
  double K_5_A = O_1_minus;
  double K_6_A = O_1_plus;
  double abs_mom_k2 = abs_mom_k*abs_mom_k;
  double Delta = k_0*(q_0*k_0 - qk)/abs_mom_k2;

  double G_1 =  W_plus* A_1 - M*A_4;
  double G_2 = -W_minus*A_1 - M*A_4;
  double G_3 =  A_1 - A_3;
  double G_4 = -G_3;
  double G_5 = (Delta + (W_plus2  - m_pi2)/Wt2 + Wt2*k_0*W_plus/(W_minus2 + Q2))*A_1  + (q_0 - Delta)*A_3 - M*W_minus*A_4/(p_10 - M);
  double G_6 =-(Delta + (W_minus2 - m_pi2)/Wt2 + Wt2*k_0*W_minus/(W_plus2 + Q2))*A_1  - (q_0 - Delta)*A_3 - M*W_plus*A_4 /(p_10 + M);
  double G_7=   (W_plus2  - m_pi2)*A_1/Wt2 + q_0*A_3 - M*A_4 + k_0*A_7  + W_plus *k_0*A_8;
  double G_8 = -(W_minus2 - m_pi2)*A_1/Wt2 - q_0*A_3 - M*A_4 - k_0*A_7  + W_minus*k_0*A_8;

  //Isobaric amplitud (Scribal G)
  double sG_1 = K_1_A*G_1/Mt2;
  double sG_2 = K_2_A*G_2/Mt2;
  double sG_3 = K_3_A*G_3/Mt2;
  double sG_4 = K_4_A*G_4/Mt2;
  double sG_5 = K_5_A*G_5/Mt2;
  double sG_6 = K_6_A*G_6/Mt2;
  double sG_7 = K_5_A*G_7/Mt2;
  double sG_8 = K_6_A*G_8/Mt2;

  // scalar part eq. 71 of ref.8
  if (is_NC)
  {
     double g_s =  0.15;
     if (spp_channel == kSpp_vp_nc_01100 || spp_channel == kSpp_vbp_nc_01100 || spp_channel == kSpp_vn_nc_10001 || spp_channel == kSpp_vbn_nc_10001)
       g_s *= kSqrt2;

     double Gs_A = g_s/g_A*FF_A;
     double As_1 = k1_Sqrt2*g_A*FA_cut/f_pi*M*Gs_A*(1/u_minus_M2 - 1/s_minus_M2);
     double As_3 =-k1_Sqrt2*g_A*FA_cut/f_pi*M*Gs_A*(1/u_minus_M2 + 1/s_minus_M2);
     double As_4 = 0;
     double As_7 = 0;
     double As_8 = 0;
     if (spp_channel == kSpp_vn_nc_01010 || spp_channel == kSpp_vbn_nc_01010)
     {
       As_1 *= -1;
       As_3 *= -1;
     }

     double Gs_1 = W_plus* As_1 - M*As_4;
     double Gs_2 = -W_minus*As_1 - M*As_4;
     double Gs_3 = As_1 - As_3;
     double Gs_4 = -Gs_3;
     double Gs_5 = (Delta + (W_plus2 - m_pi2)/Wt2  + 2*W*k_0*W_plus /(W_minus2 + Q2))*As_1 + (q_0 - Delta)*As_3 - M*W_minus*As_4/(p_10-M);
     double Gs_6 =-(Delta + (W_minus2 - m_pi2)/Wt2 + 2*W*k_0*W_minus/(W_plus2 + Q2)) *As_1 - (q_0 - Delta)*As_3 - M*W_plus*As_4/(p_10 + M);
     double Gs_7 = (W_plus2  - m_pi2)*As_1/Wt2 + q_0*As_3 - M*As_4 + k_0*As_7 + W_plus* k_0*As_8;
     double Gs_8 =-(W_minus2 - m_pi2)*As_1/Wt2 - q_0*As_3 - M*As_4 - k_0*As_7 + W_minus* k_0*As_8;

     sG_1 -= K_1_A*Gs_1/Mt2;
     sG_2 -= K_2_A*Gs_2/Mt2;
     sG_3 -= K_3_A*Gs_3/Mt2;
     sG_4 -= K_4_A*Gs_4/Mt2;
     sG_5 -= K_5_A*Gs_5/Mt2;
     sG_6 -= K_6_A*Gs_6/Mt2;
     sG_7 -= K_5_A*Gs_7/Mt2;
     sG_8 -= K_6_A*Gs_8/Mt2;
  }

  // helicity amplitudes
  Hbkg(Current::AXIAL,           BosonPolarization::LEFT,  NucleonPolarization::PLUS,  NucleonPolarization::PLUS)  =
      k1_Sqrt2*SinTheta*CosHalfTheta*(sG_3 + sG_4)*std::complex<double>
     (TMath::Cos(Phi*PhaseFactor(BosonPolarization::LEFT, NucleonPolarization::PLUS,  NucleonPolarization::PLUS)),
      TMath::Sin(Phi*PhaseFactor(BosonPolarization::LEFT, NucleonPolarization::PLUS,  NucleonPolarization::PLUS)));

  Hbkg(Current::AXIAL,           BosonPolarization::LEFT,  NucleonPolarization::MINUS, NucleonPolarization::PLUS)  =
      k1_Sqrt2*SinTheta*SinHalfTheta*(sG_3 - sG_4)*std::complex<double>
     (TMath::Cos(Phi*PhaseFactor(BosonPolarization::LEFT,  NucleonPolarization::MINUS, NucleonPolarization::PLUS)),
      TMath::Sin(Phi*PhaseFactor(BosonPolarization::LEFT,  NucleonPolarization::MINUS, NucleonPolarization::PLUS)));

  Hbkg(Current::AXIAL,           BosonPolarization::LEFT,  NucleonPolarization::PLUS,  NucleonPolarization::MINUS) =
      kSqrt2*(CosHalfTheta*(sG_1 - sG_2) - 0.5*SinTheta*SinHalfTheta*(sG_3 - sG_4))*std::complex<double>
     (TMath::Cos(Phi*PhaseFactor(BosonPolarization::LEFT,  NucleonPolarization::PLUS,  NucleonPolarization::MINUS)),
      TMath::Sin(Phi*PhaseFactor(BosonPolarization::LEFT,  NucleonPolarization::PLUS,  NucleonPolarization::MINUS)));

  Hbkg(Current::AXIAL,           BosonPolarization::LEFT,  NucleonPolarization::MINUS, NucleonPolarization::MINUS) =
      kSqrt2*(SinHalfTheta*(sG_1 + sG_2) + 0.5*SinTheta*CosHalfTheta*(sG_3 + sG_4))*std::complex<double>
     (TMath::Cos(Phi*PhaseFactor(BosonPolarization::LEFT,  NucleonPolarization::MINUS, NucleonPolarization::MINUS)),
      TMath::Sin(Phi*PhaseFactor(BosonPolarization::LEFT,  NucleonPolarization::MINUS, NucleonPolarization::MINUS)));


   Hbkg(Current::AXIAL,          BosonPolarization::RIGHT, NucleonPolarization::PLUS,  NucleonPolarization::PLUS)  =
     -kSqrt2*(SinHalfTheta*(sG_1 + sG_2) + 0.5*SinTheta*CosHalfTheta*(sG_3 + sG_4))*std::complex<double>
     (TMath::Cos(Phi*PhaseFactor(BosonPolarization::RIGHT, NucleonPolarization::PLUS,  NucleonPolarization::PLUS)),
      TMath::Sin(Phi*PhaseFactor(BosonPolarization::RIGHT, NucleonPolarization::PLUS,  NucleonPolarization::PLUS)));

  Hbkg(Current::AXIAL,           BosonPolarization::RIGHT, NucleonPolarization::MINUS, NucleonPolarization::PLUS)  =
      kSqrt2*(CosHalfTheta*(sG_1 - sG_2) - 0.5*SinTheta*SinHalfTheta*(sG_3 - sG_4))*std::complex<double>
     (TMath::Cos(Phi*PhaseFactor(BosonPolarization::RIGHT, NucleonPolarization::MINUS, NucleonPolarization::PLUS)),
      TMath::Sin(Phi*PhaseFactor(BosonPolarization::RIGHT, NucleonPolarization::MINUS, NucleonPolarization::PLUS)));

  Hbkg(Current::AXIAL,           BosonPolarization::RIGHT, NucleonPolarization::PLUS,  NucleonPolarization::MINUS) =
      k1_Sqrt2*SinTheta*SinHalfTheta*(sG_3 - sG_4)*std::complex<double>
     (TMath::Cos(Phi*PhaseFactor(BosonPolarization::RIGHT, NucleonPolarization::PLUS,  NucleonPolarization::MINUS)),
      TMath::Sin(Phi*PhaseFactor(BosonPolarization::RIGHT, NucleonPolarization::PLUS,  NucleonPolarization::MINUS)));

  Hbkg(Current::AXIAL,           BosonPolarization::RIGHT, NucleonPolarization::MINUS, NucleonPolarization::MINUS) =
     -k1_Sqrt2*SinTheta*CosHalfTheta*(sG_3 + sG_4)*std::complex<double>
     (TMath::Cos(Phi*PhaseFactor(BosonPolarization::RIGHT, NucleonPolarization::MINUS, NucleonPolarization::MINUS)),
      TMath::Sin(Phi*PhaseFactor(BosonPolarization::RIGHT, NucleonPolarization::MINUS, NucleonPolarization::MINUS)));

  if (is_CC || (is_NC && is_nu) )
  {
    Hbkg(Current::AXIAL,           BosonPolarization::MINUS0, NucleonPolarization::PLUS,  NucleonPolarization::PLUS)  =
        CosHalfTheta*(abs_mom_k*eps_z_L*(sG_5 + sG_6) + (k_0*eps_zero_L - abs_mom_k*eps_z_L)*(sG_7 + sG_8))/k_0/C_S_minus*std::complex<double>
       (TMath::Cos(Phi*PhaseFactor(BosonPolarization::MINUS0, NucleonPolarization::PLUS,  NucleonPolarization::PLUS)),
        TMath::Sin(Phi*PhaseFactor(BosonPolarization::MINUS0, NucleonPolarization::PLUS,  NucleonPolarization::PLUS)));
    
    Hbkg(Current::AXIAL,           BosonPolarization::MINUS0, NucleonPolarization::MINUS, NucleonPolarization::PLUS)  =
        SinHalfTheta*(abs_mom_k*eps_z_L*(sG_5 - sG_6) + (k_0*eps_zero_L - abs_mom_k*eps_z_L)*(sG_7 - sG_8))/k_0/C_S_minus*std::complex<double>
       (TMath::Cos(Phi*PhaseFactor(BosonPolarization::MINUS0, NucleonPolarization::MINUS, NucleonPolarization::PLUS)),
        TMath::Sin(Phi*PhaseFactor(BosonPolarization::MINUS0, NucleonPolarization::MINUS, NucleonPolarization::PLUS)));
    
    Hbkg(Current::AXIAL,           BosonPolarization::MINUS0, NucleonPolarization::PLUS,  NucleonPolarization::MINUS) =
       -SinHalfTheta*(abs_mom_k*eps_z_L*(sG_5 - sG_6) + (k_0*eps_zero_L - abs_mom_k*eps_z_L)*(sG_7 - sG_8))/k_0/C_S_minus*std::complex<double>
       (TMath::Cos(Phi*PhaseFactor(BosonPolarization::MINUS0, NucleonPolarization::PLUS,  NucleonPolarization::MINUS)),
        TMath::Sin(Phi*PhaseFactor(BosonPolarization::MINUS0, NucleonPolarization::PLUS,  NucleonPolarization::MINUS)));
    
    Hbkg(Current::AXIAL,           BosonPolarization::MINUS0, NucleonPolarization::MINUS, NucleonPolarization::MINUS) =
        CosHalfTheta*(abs_mom_k*eps_z_L*(sG_5 + sG_6) + (k_0*eps_zero_L - abs_mom_k*eps_z_L)*(sG_7 + sG_8))/k_0/C_S_minus*std::complex<double>
       (TMath::Cos(Phi*PhaseFactor(BosonPolarization::MINUS0, NucleonPolarization::MINUS, NucleonPolarization::MINUS)),
        TMath::Sin(Phi*PhaseFactor(BosonPolarization::MINUS0, NucleonPolarization::MINUS, NucleonPolarization::MINUS)));
  }
  if (is_CC || (is_NC && is_nubar) )
  {
    Hbkg(Current::AXIAL,           BosonPolarization::PLUS0,  NucleonPolarization::PLUS,  NucleonPolarization::PLUS)  =
        CosHalfTheta*(abs_mom_k*eps_z_R*(sG_5 + sG_6) + (k_0*eps_zero_R - abs_mom_k*eps_z_R)*(sG_7 + sG_8))/k_0/C_S_plus*std::complex<double>
       (TMath::Cos(Phi*PhaseFactor(BosonPolarization::PLUS0,  NucleonPolarization::PLUS,  NucleonPolarization::PLUS)),
        TMath::Sin(Phi*PhaseFactor(BosonPolarization::PLUS0,  NucleonPolarization::PLUS,  NucleonPolarization::PLUS)));

    Hbkg(Current::AXIAL,           BosonPolarization::PLUS0,  NucleonPolarization::MINUS, NucleonPolarization::PLUS)  =
        SinHalfTheta*(abs_mom_k*eps_z_R*(sG_5 - sG_6) + (k_0*eps_zero_R - abs_mom_k*eps_z_R)*(sG_7 - sG_8))/k_0/C_S_plus*std::complex<double>
       (TMath::Cos(Phi*PhaseFactor(BosonPolarization::PLUS0,  NucleonPolarization::MINUS, NucleonPolarization::PLUS)),
        TMath::Sin(Phi*PhaseFactor(BosonPolarization::PLUS0,  NucleonPolarization::MINUS, NucleonPolarization::PLUS)));

    Hbkg(Current::AXIAL,           BosonPolarization::PLUS0,  NucleonPolarization::PLUS,  NucleonPolarization::MINUS) =
       -SinHalfTheta*(abs_mom_k*eps_z_R*(sG_5 - sG_6) + (k_0*eps_zero_R - abs_mom_k*eps_z_R)*(sG_7 - sG_8))/k_0/C_S_plus*std::complex<double>
       (TMath::Cos(Phi*PhaseFactor(BosonPolarization::PLUS0,  NucleonPolarization::PLUS,  NucleonPolarization::MINUS)),
        TMath::Sin(Phi*PhaseFactor(BosonPolarization::PLUS0,  NucleonPolarization::PLUS,  NucleonPolarization::MINUS)));

    Hbkg(Current::AXIAL,           BosonPolarization::PLUS0,  NucleonPolarization::MINUS, NucleonPolarization::MINUS) =
        CosHalfTheta*(abs_mom_k*eps_z_R*(sG_5 + sG_6) + (k_0*eps_zero_R - abs_mom_k*eps_z_R)*(sG_7 + sG_8))/k_0/C_S_plus*std::complex<double>
       (TMath::Cos(Phi*PhaseFactor(BosonPolarization::PLUS0,  NucleonPolarization::MINUS, NucleonPolarization::MINUS)),
        TMath::Sin(Phi*PhaseFactor(BosonPolarization::PLUS0,  NucleonPolarization::MINUS, NucleonPolarization::MINUS)));
  }

  /*   Resonace contribution   */
  HelicityAmpVminusARes<std::complex<double>>  Hres;
  //f_BW_A is same as f_BW_V in the latest version, so rename f_BW_V = f_BW_A -> f_BW
  std::complex<double>  kappa_f_BW;

  for (auto res : fResList)
  {

    if (utils::res::Isospin(res) == 1 && SppChannel::FinStateIsospin(spp_channel) == 3)  // skip resonances with I=1/2 if isospin of final state is 3/2
      continue;

    int    NR         = utils::res::ResonanceIndex    (res);
    int    LR         = utils::res::OrbitalAngularMom (res);
    int    JR         = (utils::res::AngularMom       (res) + 1)/2;
    double MR         = utils::res::Mass              (res);
    double WR         = utils::res::Width             (res);
    double Cjsgn_plus = utils::res::Cjsgn_plus        (res);
    double Dsgn       = utils::res::Dsgn              (res);
    double BR         = SppChannel::BranchingRatio    (res);


    double d = W_plus2 + Q2;
    double sq2omg = TMath::Sqrt(2/fOmega);
    double nomg = NR*fOmega;

    //Graczyk and Sobczyk vector form-factors
    double CV_factor = 1/(1 + Q2/fMv2/4);
    // Eq. 29 of ref. 6
    double CV3 =  fCv3*CV_factor/(1 + Q2/fMv2)/(1 + Q2/fMv2);
    // Eq. 30 of ref. 6
    double CV4 = -1. * fCv4 / fCv3 * CV3;
    // Eq. 31 of ref. 6
    double CV5 =  fCv51*CV_factor/(1 + Q2/fMv2/fCv52)/(1 + Q2/fMv2/fCv52);

    // Eq. 38 of ref. 6
    double GV3 =  0.5*k1_Sqrt3*(CV4*(W2 - M2 - Q2)/2/M2 + CV5*(W2 - M2 + Q2)/2/M2 + CV3*W_plus/M);
    // Eq. 39 of ref. 6
    double GV1 = -0.5*k1_Sqrt3*(CV4*(W2 - M2 - Q2)/2/M2 + CV5*(W2 - M2 + Q2)/2/M2 - CV3*(W_plus*M + Q2)/W/M);
    // Eq. 36 of ref. 6, which is implied to use for EM-production
    double GV  =  0.5*TMath::Sqrt(1 + Q2/W_plus2)/TMath::Power(1 + Q2/4/M2, 0.5*NR)*TMath::Sqrt(3*GV3*GV3 + GV1*GV1);
    // Eq. 37 of ref. 6, which is implied to use for neutrino-production
    // double GV  =  0.5*TMath::Sqrt(1 + Q2/W_plus2)/TMath::Power(1 + Q2/4/M2, NR)*TMath::Sqrt(3*GV3*GV3 + GV1*GV1);


    //Graczyk and Sobczyk axial form-symmetry_factor
    // Eq. 52 of ref. 6
    double CA5 = fCA50/(1 + Q2/fMa2)/(1 + Q2/fMa2);

    // The form is the same like in Eq. 54 of ref. 6, but differ from it by index, which in ref. 6 is equal to NR.
    double GA = 0.5*kSqrt3*TMath::Sqrt(1 + Q2/W_plus2)*(1 - (W2 - Q2 -M2)/8/M2)*CA5/TMath::Power(1+ Q2/4/M2, 0.5*NR);

    double qMR_0            = (MR*MR - M2 + m_pi2)/(2*MR);
    double abs_mom_qMR      = TMath::Sqrt( qMR_0*qMR_0 - m_pi2);
    double Gamma            = WR*TMath::Power((abs_mom_q/abs_mom_qMR), 2*LR + 1);

    // denominator of Breit-Wigner function
    std::complex<double> denom(W - MR, Gamma/2);
    // Breit-Wigner amplitude multiplied by kappa*sqrt(BR) to avoid singularity at abs_mom_q=0, where BR = chi_E (see eq. 25 and 27 of ref. 1)
    //   double kappa            = kPi*W*TMath::Sqrt(2/JR/abs_mom_q)/M;
    //   f_BW                    = TMath::Sqrt(BR*Gamma/2/kPi)/denom;
    kappa_f_BW = W*TMath::Sqrt(kPi*BR*WR/JR/abs_mom_qMR)*TMath::Power((abs_mom_q/abs_mom_qMR), LR)/denom/M;

    fFKR.Lamda  = sq2omg*abs_mom_k;
    fFKR.Tv     = GV/3/W/sq2omg;
    fFKR.Ta     = 2./3/sq2omg*abs_mom_k*GA/d;
    fFKR.Rv     = kSqrt2*abs_mom_k*W_plus*GV/d;
    fFKR.Ra     = kSqrt2/6*(W_plus + 2*nomg*W/d)*GA/W;
    fFKR.R      = fFKR.Rv;
    fFKR.T      = fFKR.Tv;
    fFKR.Rplus  = - (fFKR.Rv + fFKR.Ra);
    fFKR.Rminus = - (fFKR.Rv - fFKR.Ra);
    fFKR.Tplus  = - (fFKR.Tv + fFKR.Ta);
    fFKR.Tminus = - (fFKR.Tv - fFKR.Ta);

    double a_aux = 1 + ((W2 + Q2 + M2)/(Mt2*W));
    // if Q2=Q=sqrt(Q2)=0 then the singularity appears in k_sqrtQ2
    // which is eliminated when in Hres at BosonPolarization=PLUS0 or BosonPolarization=MINUS0
    // when k_sqrtQ2 is multiplied by C, B and S. Therefore in this case we put Q equal to 1.
    // In our opinion it is simplier to eliminate singularity in such way, because it allows to
    // keep intact original formulas from paper.
    if (Q == 0) Q = 1;

    double C_plus = Q/C_S_plus*((eps_zero_R*abs_mom_k - eps_z_R*k_0)*(1./3 + k_0/a_aux/M) +
                   (Wt2/3 - Q2/a_aux/M + nomg/a_aux/M/3)*(eps_z_R + (eps_zero_R*k_0 - eps_z_R*abs_mom_k)*abs_mom_k/mpi2_minus_k2))*GA/Wt2/abs_mom_k;
    double C_minus = Q/C_S_minus*((eps_zero_L*abs_mom_k - eps_z_L*k_0)*(1./3 + k_0/a_aux/M) +
                   (Wt2/3 - Q2/a_aux/M + nomg/a_aux/M/3)*(eps_z_L + (eps_zero_L*k_0 - eps_z_L*abs_mom_k)*abs_mom_k/mpi2_minus_k2))*GA/Wt2/abs_mom_k;

    double B_plus =  Q/C_S_plus *(eps_zero_R + eps_z_R*abs_mom_k/a_aux/M +
                    (eps_zero_R*k_0 - eps_z_R*abs_mom_k)*(k_0 + abs_mom_k*abs_mom_k/M/a_aux)/mpi2_minus_k2)*GA/W/3/sq2omg/abs_mom_k;
    double B_minus = Q/C_S_minus*(eps_zero_L + eps_z_L*abs_mom_k/a_aux/M +
                    (eps_zero_L*k_0 - eps_z_L*abs_mom_k)*(k_0 + abs_mom_k*abs_mom_k/M/a_aux)/mpi2_minus_k2)*GA/W/3/sq2omg/abs_mom_k;

    double S_plus  = Q/C_S_plus *(eps_z_R*k_0 - eps_zero_R*abs_mom_k)*(1 + Q2/M2 - 3*W/M)*GV/abs_mom_k_L2/6;
    double S_minus = Q/C_S_minus*(eps_z_L*k_0 - eps_zero_L*abs_mom_k)*(1 + Q2/M2 - 3*W/M)*GV/abs_mom_k_L2/6;
    double k_sqrtQ2 = abs_mom_k/Q;

    const RSHelicityAmplModelI * hamplmod = 0;

    if (is_CC)
      hamplmod = fHAmplModelCC;
    else if(is_NC)
    {
      if (is_p)
        hamplmod = fHAmplModelNCp;
      else
        hamplmod = fHAmplModelNCn;
    }

    fFKR.S = S_minus;
    fFKR.B = B_minus;
    fFKR.C = C_minus;
    const RSHelicityAmpl & hampl = hamplmod->Compute(res, fFKR);
    double fp3 = hampl.AmpPlus3();
    double fp1 = hampl.AmpPlus1();
    double fm3 = hampl.AmpMinus3();
    double fm1 = hampl.AmpMinus1();
    double fm0m = 0, fm0p = 0, fp0m = 0, fp0p = 0;
    if (is_CC || (is_NC && is_nu) )
    {
      fm0m = hampl.Amp0Minus();
      fm0p = hampl.Amp0Plus();
    }
    if (is_CC || (is_NC && is_nubar) )
    {
      fFKR.S = S_plus;
      fFKR.B = B_plus;
      fFKR.C = C_plus;
      const RSHelicityAmpl & hampl_plus = hamplmod->Compute(res, fFKR);
      fp0m = hampl_plus.Amp0Minus();
      fp0p = hampl_plus.Amp0Plus();
    }

    double JRtSqrt2 = kSqrt2*JR;

    // The signs of Hres are not correct. Now we consider these signs are exactly equal to ones in the latest version
    // of code provided by Minoo Kabirnezhad, however pay attention to Cjsgn_plus
    // which differ from what stated in refs. 1 and 2.
    // More details about other problems can be found in Ref. 12, see section "Ambiguity in calculation of signs of the amplitudes"
    // (most important conclusion in subsection "Paradox and probable explanation").

    // Helicity amplitudes V-A, eq. 23 - 25 and Table 3 of ref. 1
    // Cjsgn_plus are C^j_{lS2z} for S2z=1/2 and given in Table 9 of ref. 4
    // Cjsgn_minus are C^j_{lS2z} for S2z=-1/2 and all equal to 1 (see Table 7 of ref. 4)

    // The sign of the following amplitude is opposite to one from original code, because of it will be change
    // when it will be multiplied by Wigner functions d^j_{\lambda\mu}
    Hres(res,                      BosonPolarization::LEFT, NucleonPolarization::PLUS,  NucleonPolarization::PLUS) =
        JRtSqrt2*Dsgn*Cjsgn_plus*kappa_f_BW*fp3*std::complex<double>
       (TMath::Cos(Phi*PhaseFactor(BosonPolarization::LEFT, NucleonPolarization::PLUS,  NucleonPolarization::PLUS)),
        TMath::Sin(Phi*PhaseFactor(BosonPolarization::LEFT, NucleonPolarization::PLUS,  NucleonPolarization::PLUS)));

    Hres(res,                      BosonPolarization::LEFT, NucleonPolarization::MINUS,  NucleonPolarization::PLUS) =
       -JRtSqrt2*Dsgn*Cjsgn_plus*kappa_f_BW*fp3*std::complex<double>
       (TMath::Cos(Phi*PhaseFactor(BosonPolarization::LEFT, NucleonPolarization::MINUS,  NucleonPolarization::PLUS)),
        TMath::Sin(Phi*PhaseFactor(BosonPolarization::LEFT, NucleonPolarization::MINUS,  NucleonPolarization::PLUS)));

    Hres(res,                      BosonPolarization::LEFT, NucleonPolarization::PLUS,  NucleonPolarization::MINUS) =
        JRtSqrt2*Dsgn*Cjsgn_plus*kappa_f_BW*fp1*std::complex<double>
       (TMath::Cos(Phi*PhaseFactor(BosonPolarization::LEFT, NucleonPolarization::PLUS,  NucleonPolarization::MINUS)),
        TMath::Sin(Phi*PhaseFactor(BosonPolarization::LEFT, NucleonPolarization::PLUS,  NucleonPolarization::MINUS)));


    // The sign of the following amplitude is opposite to one from original code, because of it will be change
    // when it will be multiplied by Wigner functions d^j_{\lambda\mu}
    Hres(res,                      BosonPolarization::LEFT, NucleonPolarization::MINUS,  NucleonPolarization::MINUS) =
       -JRtSqrt2*Dsgn*Cjsgn_plus*kappa_f_BW*fp1*std::complex<double>
       (TMath::Cos(Phi*PhaseFactor(BosonPolarization::LEFT, NucleonPolarization::MINUS,  NucleonPolarization::MINUS)),
        TMath::Sin(Phi*PhaseFactor(BosonPolarization::LEFT, NucleonPolarization::MINUS,  NucleonPolarization::MINUS)));



    Hres(res,                      BosonPolarization::RIGHT, NucleonPolarization::PLUS,  NucleonPolarization::PLUS) =
       -JRtSqrt2*Dsgn*Cjsgn_plus*kappa_f_BW*fm1*std::complex<double>
       (TMath::Cos(Phi*PhaseFactor(BosonPolarization::RIGHT, NucleonPolarization::PLUS,  NucleonPolarization::PLUS)),
        TMath::Sin(Phi*PhaseFactor(BosonPolarization::RIGHT, NucleonPolarization::PLUS,  NucleonPolarization::PLUS)));

    Hres(res,                      BosonPolarization::RIGHT, NucleonPolarization::MINUS,  NucleonPolarization::PLUS) =
        JRtSqrt2*Dsgn*Cjsgn_plus*kappa_f_BW*fm1*std::complex<double>
       (TMath::Cos(Phi*PhaseFactor(BosonPolarization::RIGHT, NucleonPolarization::MINUS,  NucleonPolarization::PLUS)),
        TMath::Sin(Phi*PhaseFactor(BosonPolarization::RIGHT, NucleonPolarization::MINUS,  NucleonPolarization::PLUS)));

    Hres(res,                      BosonPolarization::RIGHT, NucleonPolarization::PLUS,  NucleonPolarization::MINUS) =
       -JRtSqrt2*Dsgn*Cjsgn_plus*kappa_f_BW*fm3*std::complex<double>
       (TMath::Cos(Phi*PhaseFactor(BosonPolarization::RIGHT, NucleonPolarization::PLUS,  NucleonPolarization::MINUS)),
        TMath::Sin(Phi*PhaseFactor(BosonPolarization::RIGHT, NucleonPolarization::PLUS,  NucleonPolarization::MINUS)));

    Hres(res,                      BosonPolarization::RIGHT, NucleonPolarization::MINUS,  NucleonPolarization::MINUS) =
        JRtSqrt2*Dsgn*Cjsgn_plus*kappa_f_BW*fm3*std::complex<double>
       (TMath::Cos(Phi*PhaseFactor(BosonPolarization::RIGHT, NucleonPolarization::MINUS,  NucleonPolarization::MINUS)),
        TMath::Sin(Phi*PhaseFactor(BosonPolarization::RIGHT, NucleonPolarization::MINUS,  NucleonPolarization::MINUS)));



    Hres(res,                      BosonPolarization::MINUS0, NucleonPolarization::PLUS,  NucleonPolarization::PLUS) =
       -k_sqrtQ2*JRtSqrt2*Dsgn*kappa_f_BW*fm0m*std::complex<double>
       (TMath::Cos(Phi*PhaseFactor(BosonPolarization::MINUS0, NucleonPolarization::PLUS,  NucleonPolarization::PLUS)),
        TMath::Sin(Phi*PhaseFactor(BosonPolarization::MINUS0, NucleonPolarization::PLUS,  NucleonPolarization::PLUS)));

    // The sign of the following amplitude is opposite to one from original code, because of it will be change
    // when it will be multiplied by Wigner functions d^j_{\lambda\mu}
    Hres(res,                      BosonPolarization::MINUS0, NucleonPolarization::MINUS,  NucleonPolarization::PLUS) =
        k_sqrtQ2*JRtSqrt2*Dsgn*kappa_f_BW*fm0m*std::complex<double>
       (TMath::Cos(Phi*PhaseFactor(BosonPolarization::MINUS0, NucleonPolarization::MINUS,  NucleonPolarization::PLUS)),
        TMath::Sin(Phi*PhaseFactor(BosonPolarization::MINUS0, NucleonPolarization::MINUS,  NucleonPolarization::PLUS)));

    Hres(res,                      BosonPolarization::MINUS0, NucleonPolarization::PLUS,  NucleonPolarization::MINUS) =
        k_sqrtQ2*JRtSqrt2*Dsgn*kappa_f_BW*fm0p*std::complex<double>
       (TMath::Cos(Phi*PhaseFactor(BosonPolarization::MINUS0, NucleonPolarization::PLUS,  NucleonPolarization::MINUS)),
        TMath::Sin(Phi*PhaseFactor(BosonPolarization::MINUS0, NucleonPolarization::PLUS,  NucleonPolarization::MINUS)));

    Hres(res,                      BosonPolarization::MINUS0, NucleonPolarization::MINUS,  NucleonPolarization::MINUS) =
       -k_sqrtQ2*JRtSqrt2*Dsgn*kappa_f_BW*fm0p*std::complex<double>
       (TMath::Cos(Phi*PhaseFactor(BosonPolarization::MINUS0, NucleonPolarization::MINUS,  NucleonPolarization::MINUS)),
        TMath::Sin(Phi*PhaseFactor(BosonPolarization::MINUS0, NucleonPolarization::MINUS,  NucleonPolarization::MINUS)));



    Hres(res,                      BosonPolarization::PLUS0, NucleonPolarization::PLUS,  NucleonPolarization::PLUS) =
       -k_sqrtQ2*JRtSqrt2*Dsgn*kappa_f_BW*fp0m*std::complex<double>
       (TMath::Cos(Phi*PhaseFactor(BosonPolarization::PLUS0, NucleonPolarization::PLUS,  NucleonPolarization::PLUS)),
        TMath::Sin(Phi*PhaseFactor(BosonPolarization::PLUS0, NucleonPolarization::PLUS,  NucleonPolarization::PLUS)));

    // The sign of the following amplitude is opposite to one from original code, because of it will be change
    // when it will be multiplied by Wigner functions d^j_{\lambda\mu}
    Hres(res,                      BosonPolarization::PLUS0, NucleonPolarization::MINUS,  NucleonPolarization::PLUS) =
        k_sqrtQ2*JRtSqrt2*Dsgn*kappa_f_BW*fp0m*std::complex<double>
       (TMath::Cos(Phi*PhaseFactor(BosonPolarization::PLUS0, NucleonPolarization::MINUS,  NucleonPolarization::PLUS)),
        TMath::Sin(Phi*PhaseFactor(BosonPolarization::PLUS0, NucleonPolarization::MINUS,  NucleonPolarization::PLUS)));

    Hres(res,                      BosonPolarization::PLUS0, NucleonPolarization::PLUS,  NucleonPolarization::MINUS) =
        k_sqrtQ2*JRtSqrt2*Dsgn*kappa_f_BW*fp0p*std::complex<double>
       (TMath::Cos(Phi*PhaseFactor(BosonPolarization::PLUS0, NucleonPolarization::PLUS,  NucleonPolarization::MINUS)),
        TMath::Sin(Phi*PhaseFactor(BosonPolarization::PLUS0, NucleonPolarization::PLUS,  NucleonPolarization::MINUS)));

    Hres(res,                      BosonPolarization::PLUS0, NucleonPolarization::MINUS,  NucleonPolarization::MINUS) =
       -k_sqrtQ2*JRtSqrt2*Dsgn*kappa_f_BW*fp0p*std::complex<double>
       (TMath::Cos(Phi*PhaseFactor(BosonPolarization::PLUS0, NucleonPolarization::MINUS,  NucleonPolarization::MINUS)),
        TMath::Sin(Phi*PhaseFactor(BosonPolarization::PLUS0, NucleonPolarization::MINUS,  NucleonPolarization::MINUS)));
  } //end resonances loop

  // Sum of all helicities
  SumHelicityAmpVminusARes<std::complex<double>> sum3;
  SumHelicityAmpVminusARes<std::complex<double>> sum1;

  double pt_1 = 1;
  double pt_2 = 3*CosTheta;
  double pt_3 = 15./2*CosTheta*CosTheta - 3./2;
  double pt_4 = 35./2*CosTheta*CosTheta*CosTheta - 15./2*CosTheta;

  // Wigner functions in terms of pion polar angle (eq. B4 of ref. 1)
  double d[4][2][2];
  d[0][0][1] = CosHalfTheta;
  d[0][0][0] =-SinHalfTheta;
  d[0][1][1] = 0;
  d[0][1][0] = 0;

  d[1][0][1] = CosHalfTheta*(pt_2 - pt_1)/2;
  d[1][0][0] =-SinHalfTheta*(pt_2 + pt_1)/2;
  d[1][1][1] =-SinHalfTheta*(k1_Sqrt3*pt_2 + kSqrt3*pt_1)/2;
  d[1][1][0] = CosHalfTheta*(-k1_Sqrt3*pt_2 + kSqrt3*pt_1)/2;

  d[2][0][1] = CosHalfTheta*(pt_3 - pt_2)/3;
  d[2][0][0] =-SinHalfTheta*(pt_3 + pt_2)/3;
  d[2][1][1] =-SinHalfTheta*(k1_Sqrt2*pt_3 + kSqrt2*pt_2)/3;
  d[2][1][0] = CosHalfTheta*(-k1_Sqrt2*pt_3 + kSqrt2*pt_2)/3;

  d[3][0][1] = CosHalfTheta*(pt_4 - pt_3)/4;
  d[3][0][0] =-SinHalfTheta*(pt_4 + pt_3)/4;
  d[3][1][1] =-SinHalfTheta*(kSqrt3_5*pt_4 + kSqrt5_3*pt_3)/4;
  d[3][1][0] = CosHalfTheta*(-kSqrt3_5*pt_4 + kSqrt5_3*pt_3)/4;


  for (BosonPolarization lk : BosonPolarizationIterator() )
  {
    int lambda_k = Lambda(lk);
    for (NucleonPolarization l2 : NucleonPolarizationIterator() )
    {
       int lambda_2 = Lambda(l2);
       for (NucleonPolarization l1 : NucleonPolarizationIterator() )
       {
         int lambda_1 = Lambda(l1);
         int lambda = lambda_k - lambda_1;
         int mu = -lambda_2;
         // Symmetry properties of Wigner d-functions, see eq. A1 from ref. 7
         int symmetry_factor = 1;
         int itemp;
         if (lambda < 0)
         {
           if (TMath::Abs(lambda)>TMath::Abs(mu)) //swap 1 + swap 2
           {
             symmetry_factor = TMath::Power(-1, (lambda - mu)/2);
             lambda*=-1;
             mu*=-1;
           }
           else
           {
             if (mu<0) //swap 1
             {
               itemp = lambda;
               lambda = -mu;
               mu = -itemp;
             }
             else     //swap 2
             {
               symmetry_factor = TMath::Power(-1, (lambda - mu)/2);
               itemp = lambda;
               lambda = mu;
               mu = itemp;
             }
           }
         }

         for (auto res : fResList)
         {
            // Get baryon resonance parameters
            int    IR   = utils::res::Isospin           (res);
            int    JR   = (utils::res::AngularMom       (res) + 1)/2;
            if (SppChannel::FinStateIsospin(spp_channel) == 3 && IR == 1)  // skip resonances with I=1/2 if isospin of final state is 3/2
              continue;
            // Eq. 24 of ref. 1
            if (IR == 3)
               sum3(lk, l2, l1) += symmetry_factor*d[JR - 1][(lambda - 1)/2][(mu + 1)/2]*Hres(res, lk, l2, l1);


            if (IR == 1)
               sum1(lk, l2, l1) += symmetry_factor*d[JR - 1][(lambda - 1)/2][(mu + 1)/2]*Hres(res, lk, l2, l1);
         }
       }
    }
  }

  // Isospin Clebsch–Gordan coefficients to sum amplitudes for I=1/2 and I=3/2, see eq.25 and Table 2 of ref. 1
  double C1         = SppChannel::Isospin1Coefficients(spp_channel);
  double C3         = SppChannel::Isospin3Coefficients(spp_channel);


  double g = fFermiConstant ;
  if(is_CC) g *= fVud;

  double Lcoeff= abs_mom_k_L/2/kSqrt2/E;
  // Eq. 3.59 of ref. 2, which is multiplied by Q to avoid singularity at Q2=0
  double xsec0 = TMath::Power(g/8/kPi/kPi, 2)*abs_mom_q*Lcoeff*Lcoeff/abs_mom_k_L2/2;
  // We divide xsec0 by Q2 due to redifinition of Lcoeff


  if (kpsdim == 4)
  {
    for (NucleonPolarization l2 : NucleonPolarizationIterator() )
    {
       for (NucleonPolarization l1 : NucleonPolarizationIterator() )
       {
         // Eq. 18 ref. 1
         xsec +=
         TMath::Power(
            std::abs(
                C_L_minus*(
                   Hbkg(Current::VECTOR, BosonPolarization::LEFT,   l2, l1) - Hbkg(Current::AXIAL, BosonPolarization::LEFT,   l2, l1)  +
                    C1*sum1(             BosonPolarization::LEFT,   l2, l1) + C3*sum3(             BosonPolarization::LEFT,   l2, l1)
                ) +
                C_R_minus*(
                   Hbkg(Current::VECTOR, BosonPolarization::RIGHT,  l2, l1) - Hbkg(Current::AXIAL, BosonPolarization::RIGHT,  l2, l1) +
                 C1*sum1(                BosonPolarization::RIGHT,  l2, l1) + C3*sum3(             BosonPolarization::RIGHT,  l2, l1)
                ) +
                C_S_minus*(
                   Hbkg(Current::VECTOR, BosonPolarization::MINUS0, l2, l1) - Hbkg(Current::AXIAL, BosonPolarization::MINUS0, l2, l1) +
                   C1*sum1(              BosonPolarization::MINUS0, l2, l1) + C3*sum3(             BosonPolarization::MINUS0, l2, l1)
                )
              )
        , 2)
                                                                      +
         TMath::Power(
            std::abs(
               C_L_plus*(
                   Hbkg(Current::VECTOR, BosonPolarization::LEFT,   l2, l1) - Hbkg(Current::AXIAL,  BosonPolarization::LEFT,  l2, l1)   +
                   C1*sum1(              BosonPolarization::LEFT,   l2, l1) + C3*sum3(              BosonPolarization::LEFT,  l2, l1)
               ) +
               C_R_plus*(
                   Hbkg(Current::VECTOR, BosonPolarization::RIGHT,  l2, l1) - Hbkg(Current::AXIAL,  BosonPolarization::RIGHT, l2, l1)  +
                   C1*sum1(              BosonPolarization::RIGHT,  l2, l1) + C3*sum3(              BosonPolarization::RIGHT, l2, l1)
               ) +
               C_S_plus*(
                   Hbkg(Current::VECTOR, BosonPolarization::PLUS0,  l2, l1) - Hbkg(Current::AXIAL,  BosonPolarization::PLUS0, l2, l1)  +
                   C1*sum1(              BosonPolarization::PLUS0,  l2, l1) + C3*sum3(              BosonPolarization::PLUS0, l2, l1)
               )
             )
        , 2);
       }
    }
  }
  else if (kpsdim == 3)
  {
    for (NucleonPolarization l2 : NucleonPolarizationIterator() )
    {
       for (NucleonPolarization l1 : NucleonPolarizationIterator() )
       {
          xsec +=
                   (C_L_minus*C_L_minus + C_L_plus*C_L_plus)*(
                   TMath::Power((Hbkg(Current::VECTOR, BosonPolarization::LEFT,   l2, l1) - Hbkg(Current::AXIAL,  BosonPolarization::LEFT,   l2, l1)).real(),  2) +
                   TMath::Power(std::abs(C1*sum1(      BosonPolarization::LEFT,   l2, l1) + C3*sum3(              BosonPolarization::LEFT,   l2, l1)), 2) +
                   2*(          Hbkg(Current::VECTOR,  BosonPolarization::LEFT,   l2, l1) - Hbkg(Current::AXIAL,  BosonPolarization::LEFT,   l2, l1)).real()*
                                        (C1*sum1(      BosonPolarization::LEFT,   l2, l1) + C3*sum3(              BosonPolarization::LEFT,   l2, l1)).real()
                                                             )+
                   (C_R_minus*C_R_minus + C_R_plus*C_R_plus)*(
                   TMath::Power((Hbkg(Current::VECTOR, BosonPolarization::RIGHT,  l2, l1) - Hbkg(Current::AXIAL,  BosonPolarization::RIGHT,  l2, l1)).real(),  2) +
                   TMath::Power(std::abs(C1*sum1(      BosonPolarization::RIGHT,  l2, l1) + C3*sum3(              BosonPolarization::RIGHT,  l2, l1)), 2) +
                   2*(          Hbkg(Current::VECTOR,  BosonPolarization::RIGHT,  l2, l1) - Hbkg(Current::AXIAL,  BosonPolarization::RIGHT,  l2, l1)).real()*
                                        (C1*sum1(      BosonPolarization::RIGHT,  l2, l1) + C3*sum3(              BosonPolarization::RIGHT,  l2, l1)).real()
                                                             )+
                   (C_S_minus*C_S_minus)*(
                   TMath::Power((Hbkg(Current::VECTOR, BosonPolarization::MINUS0, l2, l1) - Hbkg(Current::AXIAL,  BosonPolarization::MINUS0, l2, l1)).real(),  2) +
                   TMath::Power(std::abs(C1*sum1(      BosonPolarization::MINUS0, l2, l1) + C3*sum3(              BosonPolarization::MINUS0, l2, l1)), 2) +
                   2*(          Hbkg(Current::VECTOR,  BosonPolarization::MINUS0, l2, l1) - Hbkg(Current::AXIAL,  BosonPolarization::MINUS0, l2, l1)).real()*
                                        (C1*sum1(      BosonPolarization::MINUS0, l2, l1) + C3*sum3(              BosonPolarization::MINUS0, l2, l1)).real()
                                                             )+
                   (C_S_plus*C_S_plus)*(
                   TMath::Power((Hbkg(Current::VECTOR, BosonPolarization::PLUS0,  l2, l1) - Hbkg(Current::AXIAL,  BosonPolarization::PLUS0,  l2, l1)).real(),  2) +
                   TMath::Power(std::abs(C1*sum1(      BosonPolarization::PLUS0,  l2, l1) + C3*sum3(              BosonPolarization::PLUS0,  l2, l1)), 2) +
                   2*(          Hbkg(Current::VECTOR,  BosonPolarization::PLUS0,  l2, l1) - Hbkg(Current::AXIAL,  BosonPolarization::PLUS0,  l2, l1)).real()*
                                        (C1*sum1(      BosonPolarization::PLUS0,  l2, l1) + C3*sum3(              BosonPolarization::PLUS0,  l2, l1)).real()
                                                             );
       }
   }
   xsec*= 2*kPi;
 }

  xsec*=xsec0;


  // The algorithm computes d^4xsec/dWdQ2dCosTheta_pidPhi_pi or d^3xsec/dWdQ2dCosTheta_pi
  // Check whether variable tranformation is needed
  if ( kps != kPSWQ2ctpphipfE && kps != kPSWQ2ctpfE )
  {
     double J = 1;
     if (kpsdim == 3)
       J = utils::kinematics::Jacobian(interaction, kPSWQ2ctpfE, kps);
     else if (kpsdim == 4)
       J = utils::kinematics::Jacobian(interaction, kPSWQ2ctpphipfE, kps);
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
  int NNucl = (SppChannel::InitStateNucleon(spp_channel) == kPdgProton) ? Z : N;
  xsec*=NNucl; // nuclear xsec (no nuclear suppression symmetry_factor)

  if ( fUsePauliBlocking && A!=1 && kps == kPSWQ2ctpfE )
  {
    // Calculation of Pauli blocking according to refs. 9-11
    double P_Fermi = 0;

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

     double FactorPauli_RES = 1;

     if (P_Fermi > 0)
     {
        if ( 2*P_Fermi < abs_mom_k-abs_mom_q )
           FactorPauli_RES = 1;
        if ( 2*P_Fermi >= abs_mom_k+abs_mom_q )
           FactorPauli_RES = ((3*abs_mom_k*abs_mom_k + abs_mom_q*abs_mom_q)/(2*P_Fermi) - (5*TMath::Power(abs_mom_k,4) + TMath::Power(abs_mom_q,4) + 10*abs_mom_k*abs_mom_k*abs_mom_q*abs_mom_q)/(40*TMath::Power(P_Fermi,3)))/(2*abs_mom_k);
        if ( 2*P_Fermi >= abs_mom_k-abs_mom_q && 2*P_Fermi <= abs_mom_k+abs_mom_q )
           FactorPauli_RES = ((abs_mom_q + abs_mom_k)*(abs_mom_q + abs_mom_k) - 4*P_Fermi*P_Fermi/5 - TMath::Power(abs_mom_k - abs_mom_q, 3)/(2*P_Fermi)+TMath::Power(abs_mom_k - abs_mom_q, 5)/(40*TMath::Power(P_Fermi, 3)))/(4*abs_mom_q*abs_mom_k);
     }
     xsec *= FactorPauli_RES;
  }
  return xsec;

}

//____________________________________________________________________________
double MKSPPPXSec2020::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool MKSPPPXSec2020::ValidProcess(const Interaction * interaction) const
{

  if(interaction->TestBit(kISkipProcessChk)) return true;

  //-- Get the requested SPP channel
  SppChannel_t spp_channel = SppChannel::FromInteraction(interaction);
  if( spp_channel == kSppNull ) {
    return false;
  }

  return true;

}
//____________________________________________________________________________
inline int MKSPPPXSec2020::Lambda (NucleonPolarization l) const
{
  return 2*l - 1;
}
//____________________________________________________________________________
inline int MKSPPPXSec2020::Lambda (BosonPolarization l) const
{
  return (l*(l*(4*l-21)+29)-6)/3;
}
//____________________________________________________________________________
inline int MKSPPPXSec2020::PhaseFactor(BosonPolarization lk, NucleonPolarization l1, NucleonPolarization l2) const
{
  return lk*(lk*(4*lk - 21) + 29)/6 - l1 - l2;
}
//____________________________________________________________________________
bool MKSPPPXSec2020::ValidKinematics(const Interaction * interaction) const
{
  // call only after ValidProcess
  if ( interaction->TestBit(kISkipKinematicChk) ) return true;

  const KPhaseSpace& kps = interaction->PhaseSpace();

  // Get kinematical parameters
  const InitialState & init_state = interaction -> InitState();
  const Kinematics & kinematics = interaction -> Kine();
  double Enu = init_state.ProbeE(kRfHitNucRest);
  double W    = kinematics.W();
  double Q2   = kinematics.Q2();

  if (Enu < kps.Threshold_SPP_iso())
    return false;

  Range1D_t Wl  = kps.WLim_SPP_iso();
  if (fWmax >= Wl.min)
    Wl.max = TMath::Min(fWmax, Wl.max);
  Range1D_t Q2l = kps.Q2Lim_W_SPP_iso();

  if (W < Wl.min || W > Wl.max)
    return false;

  if (Q2 < Q2l.min || Q2 > Q2l.max)
    return false;

  return true;

}
//____________________________________________________________________________
void MKSPPPXSec2020::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void MKSPPPXSec2020::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void MKSPPPXSec2020::LoadConfig(void)
{

  fResList.Clear();
  string resonances ;
  GetParam( "ResonanceNameList", resonances ) ;
  fResList.DecodeFromNameList(resonances);

  // Cross section scaling symmetry_factors
  this->GetParam( "RES-CC-XSecScale", fXSecScaleCC );
  this->GetParam( "RES-NC-XSecScale", fXSecScaleNC );

  // Load all configuration data or set defaults
  this->GetParamDef( "RES-Omega"  , fOmega, 1.05);

  double ma, mv;
  this->GetParam( "RES-Ma"   , ma);
  this->GetParam( "RES-Mv"   , mv);
  fMa2 = ma*ma;
  fMv2 = mv*mv;
  this->GetParam( "RES-CA50" , fCA50) ;

  this->GetParam( "GVCAL-Cv3"  , fCv3)  ;
  this->GetParam( "GVCAL-Cv4"  , fCv4)  ;
  this->GetParam( "GVCAL-Cv51" , fCv51) ;
  this->GetParam( "GVCAL-Cv52" , fCv52) ;

  this->GetParam( "FermiConstant", fFermiConstant );
  
  double thw;
  this->GetParam( "WeinbergAngle", thw );
  fSin2Wein = TMath::Power( TMath::Sin(thw), 2 );

  this->GetParam("CKM-Vud", fVud );

  // Load all the sub-algorithms needed
  fHAmplModelCC     = 0;
  fHAmplModelNCp    = 0;
  fHAmplModelNCn    = 0;

  fHAmplModelCC  = dynamic_cast<const RSHelicityAmplModelI *> ( this -> SubAlg( "HelictyAmplCCAlg" ) );
  fHAmplModelNCp = dynamic_cast<const RSHelicityAmplModelI *> ( this -> SubAlg( "HelictyAmplNCpAlg" ));
  fHAmplModelNCn = dynamic_cast<const RSHelicityAmplModelI *> ( this -> SubAlg( "HelictyAmplNCnAlg" ));

  assert(fHAmplModelCC);
  assert(fHAmplModelNCp);
  assert(fHAmplModelNCn);

  // Parameters for calculation of background contribution
  fFormFactorsModel = dynamic_cast<const QELFormFactorsModelI *> (this->SubAlg("CCFormFactorsAlg"));
  assert(fFormFactorsModel);
  fFormFactors.SetModel(fFormFactorsModel); // <-- attach algorithm

  fEMFormFactorsModel = dynamic_cast<const QELFormFactorsModelI *> (this->SubAlg("EMFormFactorsAlg"));
  assert(fEMFormFactorsModel);
  fEMFormFactors.SetModel(fEMFormFactorsModel); // <-- attach algorithm

  this->GetParam("FermiMomentumTable", fKFTable);
  this->GetParam("RFG-UseParametrization", fUseRFGParametrization);
  this->GetParam("UsePauliBlockingForRES", fUsePauliBlocking);


  this->GetParamDef("MK-fPi",             f_pi,           0.093 );
  this->GetParamDef("QEL-FA0",            FA0,            -1.26 );
  this->GetParamDef("MK-Frho0",           Frho0,            1.0 );

  this->GetParamDef("MK-NonResBkg-VWmin", fBkgVWmin,       1.30 );
  this->GetParamDef("MK-NonResBkg-VWmax", fBkgVWmax,       1.60 );
  this->GetParamDef("MK-NonResBkg-V3",    fBkgV3,       8.08497 );
  this->GetParamDef("MK-NonResBkg-V2",    fBkgV2,      -41.6767 );
  this->GetParamDef("MK-NonResBkg-V1",    fBkgV1,       66.3419 );
  this->GetParamDef("MK-NonResBkg-V0",    fBkgV0,      -32.5733 );
  this->GetParamDef("Mass-rho770",        fRho770Mass,   0.7758 );
  this->GetParamDef("MK-WMax",            fWmax,           -1.0 );
  
  this->GetParamDef("UseAuthorCode", fUseAuthorCode, false );

  // Load the differential cross section integrator
  fXSecIntegrator = dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);

}
//____________________________________________________________________________
