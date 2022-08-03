//____________________________________________________________________________
/*
 Copyright (c) 2003-2022, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Authors: Igor Kakorin <kakorin@jinr.ru>, Joint Institute for Nuclear Research

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
#include "Physics/Resonance/XSection/DCCEMSPPPXSec.h"
#include "Physics/NuclearState/FermiMomentumTablePool.h"
#include "Physics/NuclearState/FermiMomentumTable.h"
#include "Physics/NuclearState/NuclearUtils.h"

#include <algorithm>


using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
DCCEMSPPPXSec::DCCEMSPPPXSec() :
XSecAlgorithmI("genie::DCCEMSPPPXSec")
{

}
//____________________________________________________________________________
DCCEMSPPPXSec::DCCEMSPPPXSec(string config) :
XSecAlgorithmI("genie::DCCEMSPPPXSec", config)
{

}
//____________________________________________________________________________
DCCEMSPPPXSec::~DCCEMSPPPXSec()
{

}
//____________________________________________________________________________
double DCCEMSPPPXSec::XSec(const Interaction * interaction, KinePhaseSpace_t kps) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  const InitialState & init_state = interaction -> InitState();
  const ProcessInfo &  proc_info  = interaction -> ProcInfo();
  const Target & target = init_state.Tgt();

  double xsec = 0;
  //-- Get 1pi exclusive channel
  SppChannel_t spp_channel = SppChannel::FromInteraction(interaction);

  int  nucpdgc   = target.HitNucPdg();
  int  probepdgc = init_state.ProbePdg();
  bool is_nubar  = pdg::IsAntiNeutrino     (probepdgc);
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
  if (kpsdim < 3 || kpsdim > 4) return 0.;
  // Pion angles should be given in Adler frame
  double CosTheta = kinematics.GetKV(kKVctp);
  double Phi = 0.;
  if (kpsdim == 4)
    Phi = kinematics.GetKV(kKVphip);

  double SinTheta   = TMath::Sqrt(1. - CosTheta*CosTheta);
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
  double abs_mom_q         = TMath::Sqrt(q_02 - m_pi2);
  double abs_mom_k         = TMath::Sqrt(k_0*k_0 + Q2);
  //double E_2L              = (M2 - W2 - Q2 + Mt2*E)/Mt2;
  double abs_mom_k_L       = W*abs_mom_k/M;
  double abs_mom_k_L2      = abs_mom_k_L*abs_mom_k_L;
  double k_2               = (M2 + Mt2*E - W2 -ml2)/Wt2;
  double k_1               =  k_2 + k_0;
  double p_10              = (W2 + M2 + Q2)/Wt2;
  double qk                = q_0*k_0 - abs_mom_k*abs_mom_q*CosTheta;
  //double k_2L              = TMath::Sqrt(E_2L*E_2L - ml2);         //magnitude of lepton momentum in lab frame
  double k_2_iso           = TMath::Sqrt(k_2*k_2 - ml2);   //magnitude of lepton momentum in isobaric frame
  double cos_theta         = (2*k_1*k_2 - Q2 - ml2)/2/k_1/k_2_iso;

  // Eq. 7 of ref. 1
  double A_plus            = TMath::Sqrt( k_1*(k_2 - k_2_iso) );
  double A_minus           = TMath::Sqrt( k_1*(k_2 + k_2_iso) );
  //Eq. 6 of ref. 1
  double eps_1_plus        =  2.*A_plus *(k_1 - k_2_iso)/abs_mom_k*TMath::Sqrt(1. + cos_theta);
  double eps_1_minus       =  -2.*A_minus*(k_1 + k_2_iso)/abs_mom_k*TMath::Sqrt(1. - cos_theta);
  double eps_2_plus        =  2.*A_plus *TMath::Sqrt(1. + cos_theta);
  double eps_2_minus       =  2.*A_minus*TMath::Sqrt(1. - cos_theta);
  //Eq. 9 of ref. 1
  double eps_zero_L        = -2.*A_minus*TMath::Sqrt(1. + cos_theta);       // L->lambda = -1
  double eps_zero_R        = 2.*A_plus *TMath::Sqrt(1. - cos_theta);       // R->lambda = +1
  double eps_z_L           = -2.*A_minus*(k_1 - k_2_iso)/abs_mom_k*TMath::Sqrt(1. + cos_theta);
  double eps_z_R           = 2.*A_plus *(k_1 + k_2_iso)/abs_mom_k*TMath::Sqrt(1. - cos_theta);
  // This ``recipe'' of transition from neutrino to antineutrino case is promoted by Minoo Kabirnezhad.
  // However, it is not correct in our opinion. All details can be found in Ref. [12], see
  // section "Problem with transition from neutrino to antineutrino case".
  if (is_nubar)
  {
     Phi = -Phi;
     //Eq. 6 of ref. 1
     eps_1_plus        =  -2.*A_minus *(k_1 + k_2_iso)/abs_mom_k*TMath::Sqrt(1. - cos_theta);
     eps_1_minus       =  -2.*A_plus*(k_1 - k_2_iso)/abs_mom_k*TMath::Sqrt(1. + cos_theta);
     eps_2_plus        =  -2.*A_minus *TMath::Sqrt(1. - cos_theta);
     eps_2_minus       =  2.*A_plus*TMath::Sqrt(1. + cos_theta);
     //Eq. 9 of ref. 1
     eps_zero_L        = 2.*A_plus*TMath::Sqrt(1. - cos_theta);       // L->lambda = -1
     eps_zero_R        = 2.*A_minus *TMath::Sqrt(1. + cos_theta);       // R->lambda = +1
     eps_z_L           = 2.*A_plus*(k_1 + k_2_iso)/abs_mom_k*TMath::Sqrt(1. - cos_theta);
     eps_z_R           = 2.*A_minus *(k_1 - k_2_iso)/abs_mom_k*TMath::Sqrt(1. + cos_theta);
  }



  /*    Bkg Contribution   */
  //Auxiliary functions
  double W_plus            = W + M;
  double W_minus           = W - M;
  double W_plus2           = W_plus*W_plus;
  double W_minus2          = W_minus*W_minus;
  double O_1_plus          = TMath::Sqrt((W_plus2 + Q2)*( W_plus2 - m_pi2 ))/Wt2;
  double O_1_minus         = TMath::Sqrt((W_minus2 + Q2)*(W_minus2 - m_pi2))/Wt2;
  double O_2_plus          = TMath::Sqrt((W_plus2 + Q2)/(W_plus2 - m_pi2));
  double O_2_minus         = TMath::Sqrt((W_minus2 + Q2)/(W_minus2 - m_pi2));
  double s_minus_M2        = W*W - M*M;
  double u_minus_M2        = m_pi2 - 2*(q_0*p_10 + abs_mom_q*abs_mom_k*CosTheta);
  double t_minus_mpi2      = -(Q2 + 2*qk);
  double mpi2_minus_k2     = Q2 + m_pi2;
  double Q                 = TMath::Sqrt(Q2);





  //f_BW_A is same as f_BW_V in the latest version, so rename f_BW_V = f_BW_A -> f_BW
  std::complex<double>  f_BW;

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
    double sq2omg = TMath::Sqrt(2./fOmega);
    double nomg = NR*fOmega;
  }



    //Graczyk and Sobczyk vector form-factors
    double CV_factor = 1./(1. + Q2/fMv2/4.);
    // Eq. 29 of ref. 6
    double CV3 =  2.13*CV_factor/(1 + Q2/fMv2)/(1 + Q2/fMv2);
    // Eq. 30 of ref. 6
    double CV4 = -1.51/2.13*CV3;
    // Eq. 31 of ref. 6
    double CV5 =  0.48*CV_factor/(1 + Q2/fMv2/0.766)/(1 + Q2/fMv2/0.766);

    // Eq. 38 of ref. 6
    double GV3 =  0.5*k1_Sqrt3*(CV4*(W2 - M2 - Q2)/2./M2 + CV5*(W2 - M2 + Q2)/2./M2 + CV3*W_plus/M);
    // Eq. 39 of ref. 6
    double GV1 = -0.5*k1_Sqrt3*(CV4*(W2 - M2 - Q2)/2./M2 + CV5*(W2 - M2 + Q2)/2./M2 - CV3*(W_plus*M + Q2)/W/M);
    // Eq. 36 of ref. 6, which is implied to use for EM-production
    double GV  =  0.5*TMath::Sqrt(1 + Q2/W_plus2)/TMath::Power(1 + Q2/4/M2, 0.5*NR)*TMath::Sqrt(3.*GV3*GV3 + GV1*GV1);
    // Eq. 37 of ref. 6, which is implied to use for neutrino-production
    // double GV  =  0.5*TMath::Sqrt(1 + Q2/W_plus2)/TMath::Power(1 + Q2/4/M2, NR)*TMath::Sqrt(3.*GV3*GV3 + GV1*GV1);



    //Graczyk and Sobczyk axial form-symmetry_factor
    // Eq. 52 of ref. 6
    double CA5 = fCA50/(1 + Q2/fMa2)/(1 + Q2/fMa2);

    // The form is the same like in Eq. 54 of ref. 6, but differ from it by index, which in ref. 6 is equal to NR.
    double GA = 0.5*kSqrt3*TMath::Sqrt(1 + Q2/W_plus2)*(1 - (W2 - Q2 -M2)/8/M2)*CA5/TMath::Power(1+ Q2/4/M2, 0.5*NR);

    double qMR_0        = (MR*MR - M2 + m_pi2)/(2.*MR);
    double abs_mom_qMR  = TMath::Sqrt( qMR_0*qMR_0 - m_pi2);
    double kappa        = kPi*W*TMath::Sqrt(2./JR/abs_mom_q)/M;
    double Gamma        = WR*TMath::Power((abs_mom_q/abs_mom_qMR), 2*LR + 1);



    // denominator of Breit-Wigner function
    std::complex<double> denom(W-MR, Gamma/2);
    //Breit-Wigner amplitude multiplied by sqrt(BR), where BR = chi_E (see eq. 25 and 27 of ref. 1)
    f_BW = TMath::Sqrt(BR*Gamma/2/kPi)/denom;




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

/*
    // Not used in the latest version
    const MKHelicityAmplModelI * hamplmod = 0;
*/
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

/*
    // Not used in the latest version
    const MKHelicityAmpl & hampl = hamplmod->Compute(res, fFKR);
*/
    fFKR.S = S_minus;
    fFKR.B = B_minus;
    fFKR.C = C_minus;
    const RSHelicityAmpl & hampl = hamplmod->Compute(res, fFKR);
    double fp3 = hampl.AmpPlus3();
    double fp1 = hampl.AmpPlus1();
    double fm3 = hampl.AmpMinus3();
    double fm1 = hampl.AmpMinus1();
    double fm0m = hampl.Amp0Minus();
    double fm0p = hampl.Amp0Plus();

    double fp0m = 0., fp0p = 0.;
    if (is_CC)
    {
       fFKR.S = S_plus;
       fFKR.B = B_plus;
       fFKR.C = C_plus;
       const RSHelicityAmpl & hampl_plus = hamplmod->Compute(res, fFKR);
       fp0m = hampl_plus.Amp0Minus();
       fp0p = hampl_plus.Amp0Plus();
    }

    double JRtSqrt2 = kSqrt2*JR;


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
    
  }
  else if (kpsdim == 3)
  {
    
   xsec*= 2*kPi;
  }

  xsec*=xsec0;


  // The algorithm computes d^4xsec/dWdQ2dCostThetadPhi or d^3xsec/dWdQ2dCostTheta
  // Check whether variable tranformation is needed
  if ( kps != kPSWQ2ctpphipfE && kps != kPSWQ2ctpfE )
  {
     double J = 1.;
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

     double k0 = 0., q = 0., q0 = 0., k = 0.;

     if (P_Fermi > 0.)
     {
        k0 = (W2 - M2 - Q2)/(2*W);
        k = TMath::Sqrt(k0*k0 + Q2);  // previous value of k is overridden
        q0 = (W2 - M2 + kPionMass2)/(2*W);
        q = TMath::Sqrt(q0*q0 - kPionMass2);
     }

     if ( 2*P_Fermi < k-q )
        FactorPauli_RES = 1.0;
     if ( 2*P_Fermi >= k+q )
        FactorPauli_RES = ((3*k*k + q*q)/(2*P_Fermi) - (5*TMath::Power(k,4) + TMath::Power(q,4) + 10*k*k*q*q)/(40*TMath::Power(P_Fermi,3)))/(2*k);
     if ( 2*P_Fermi >= k-q && 2*P_Fermi <= k+q )
        FactorPauli_RES = ((q + k)*(q + k) - 4*P_Fermi*P_Fermi/5 - TMath::Power(k - q, 3)/(2*P_Fermi)+TMath::Power(k - q, 5)/(40*TMath::Power(P_Fermi, 3)))/(4*q*k);

     xsec *= FactorPauli_RES;
  }
  return xsec;

}

//____________________________________________________________________________
double DCCEMSPPPXSec::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool DCCEMSPPPXSec::ValidProcess(const Interaction * interaction) const
{

  if(interaction->TestBit(kISkipProcessChk)) return true;

  //-- Get the requested SPP channel
  SppChannel_t spp_channel = SppChannel::FromInteraction(interaction);
  if( spp_channel == kSppNull ) {

    LOG("DCCEMSPPPXSec", pERROR)
            << "Insufficient SPP exclusive final state information!\n";
    return false;
  }

  return true;

}
//____________________________________________________________________________
bool DCCEMSPPPXSec::ValidKinematics(const Interaction * interaction) const
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

  if (Enu < kps.Threshold_RSPP())
    return false;

  Range1D_t Wl  = kps.WLim_RSPP();
  Range1D_t Q2l = kps.Q2Lim_W_RSPP();

  if (W < Wl.min || W > Wl.max)
    return false;

  if (Q2 < Q2l.min || Q2 > Q2l.max)
    return false;

  return true;

}
//____________________________________________________________________________
void DCCEMSPPPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DCCEMSPPPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DCCEMSPPPXSec::LoadConfig(void)
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

  this->GetParam( "FermiConstant", fFermiConstant );
  
  double thw;
  this->GetParam( "WeinbergAngle", thw );
  fSin2Wein = TMath::Power( TMath::Sin(thw), 2 );

  this->GetParam("CKM-Vud", fVud );

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


  this->GetParamDef("Mass-rho770 ", fRho770Mass,   0.7758 );

  // Load the differential cross section integrator
  fXSecIntegrator = dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);

}
//____________________________________________________________________________
