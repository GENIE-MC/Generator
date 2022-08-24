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
#include "Physics/Resonance/XSection/DCCSPPPXSec.h"
#include "Physics/NuclearState/FermiMomentumTablePool.h"
#include "Physics/NuclearState/FermiMomentumTable.h"
#include "Physics/NuclearState/NuclearUtils.h"

#include <algorithm>


using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
DCCSPPPXSec::DCCSPPPXSec() :
XSecAlgorithmI("genie::DCCSPPPXSec")
{

}
//____________________________________________________________________________
DCCSPPPXSec::DCCSPPPXSec(string config) :
XSecAlgorithmI("genie::DCCSPPPXSec", config)
{

}
//____________________________________________________________________________
DCCSPPPXSec::~DCCSPPPXSec()
{

}
//____________________________________________________________________________
double DCCSPPPXSec::XSec(const Interaction * interaction, KinePhaseSpace_t kps) const
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
  double Eli_L = init_state.ProbeE(kRfHitNucRest);
  double ml    = interaction->FSPrimLepton()->Mass();
  double ml2   = ml*ml;
  double Q2    = kinematics.Q2();
  double W     = kinematics.W();
  double W2    = W*W;
  
  // dimension of kine phase space
  std::string s = KinePhaseSpace::AsString(kps);
  int kpsdim = s!="<|E>"?1 + std::count(s.begin(), s.begin()+s.find('}'), ','):0;
  if (kpsdim < 3 || kpsdim > 4) return 0.;
  // Pion angles should be given in pi-N rest frame
  double cos_theta_pi = kinematics.GetKV(kKVctp);
  double phi_pi = 0.;
  if (kpsdim == 4)
    phi_pi = kinematics.GetKV(kKVphip);

  double sin_theta_pi   = TMath::Sqrt(1. - cos_theta_pi*cos_theta_pi);
//  double SinHalfTheta  = TMath::Sqrt((1 - CosTheta_pi)/2);
//  double CosHalfTheta  = TMath::Sqrt((1 + CosTheta_pi)/2);

  PDGLibrary * pdglib = PDGLibrary::Instance();
  double Mi    = pdglib->Find(SppChannel::InitStateNucleon(spp_channel))->Mass();
  double Mi2   = Mi*Mi;
  double Mf    = pdglib->Find(SppChannel::FinStateNucleon(spp_channel))->Mass();
  double Mf2   = Mf*Mf;
  double m_pi   = pdglib->Find(SppChannel::FinStatePion(spp_channel))->Mass();
  double m_pi2 = m_pi*m_pi;


  // Eq. 14 of ref. 4
  double omega             = (W2 - Mi2 - Q2)/2./W;
  double omega_L           = omega*W/Mi;
  // Eq. 15 of ref. 4
  double q                 = TMath::Sqrt(Q2 + omega*omega);
  double q_L               = TMath::Sqrt(Q2 + omega_L*omega_L);
  
  double k0                = (W2 - Mf2 + m_pi2)/2./W;
  // Eq. 16 of ref. 4
  double k                 = TMath::Sqrt(k0*k0 - m_pi2);
  // Eq. 17 of ref. 4
  double q_gamma           = (W2 - Mi2)/2./W;
  // Eq. 3 of ref. 4
  double qL_gamma          = (W2 - Mi2)/2./Mi;
  
  double Elf_L             = Eli_L - omega_L;
  double pli_L             = TMath::Sqrt(Eli_L*Eli_L - ml2);
  double plf_L             = TMath::Sqrt(Elf_L*Elf_L - ml2);
  double cos_theta_l       = (Eli_L*Elf_L - Q2/2 - ml2)/pli_L/plf_L;
  // Eq. 4 of ref. 4
  double epsilon = 0.;
  if (cos_theta_l > -1.)
  {
    double tan_half_theta_l2  = (1. - cos_theta_l)/(1. + cos_theta_l);
    epsilon                   = 1./(1. + 2.*qL*qL*tan_half_theta_l2/Q2);
  }
  // Eq. 2 of ref. 4
  double Gamma             = kAem*qL_gamma*Elf_L/2./kPi2/Q2/Eli_L/(1. - epsilon);
    
 



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
  if      (is_EM) { xsec *= fXSecScaleEM; }

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
double DCCSPPPXSec::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool DCCSPPPXSec::ValidProcess(const Interaction * interaction) const
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
bool DCCSPPPXSec::ValidKinematics(const Interaction * interaction) const
{
  // call only after ValidProcess
  if ( interaction->TestBit(kISkipKinematicChk) ) return true;

  const KPhaseSpace  & kps        = interaction->PhaseSpace();
  SppChannel_t spp_channel        = SppChannel::FromInteraction(interaction);
  
  // Get kinematical parameters
  const InitialState & init_state = interaction -> InitState();
  const Kinematics & kinematics = interaction -> Kine();
  double Enu = init_state.ProbeE(kRfHitNucRest);
  double W    = kinematics.W();
  double Q2   = kinematics.Q2();

  if (Enu < kps.Threshold())
    return false;
    
  PDGLibrary * pdglib = PDGLibrary::Instance();
  double Mi   = pdglib->Find(SppChannel::InitStateNucleon(spp_channel))->Mass();
  double Mf   = pdglib->Find(SppChannel::FinStateNucleon(spp_channel))->Mass();
  double mpi  = pdglib->Find(SppChannel::FinStatePion(spp_channel))->Mass();
  double ml   = interaction->FSPrimLepton()->Mass();
  double ml2  = ml*ml;
  
  double s = Mi*(Mi + 2*Enu);
  double sqrt_s = TMath::Sqrt(s);
 
  // kinematic W-limits
  double Wmin = Mf + mpi;
  double Wmax  = sqrt_s - ml;
  
  // kinematic Q2-limits
  double Enu_CM = (s - Mi*Mi)/2/sqrt_s;
  double El_CM  = (s + ml2 - W*W)/2/sqrt_s;
  double Pl_CM  = (El_CM - ml)<0?0:TMath::Sqrt(El_CM*El_CM - ml2);
  Q2min = (2*Enu_CM*(El_CM - Pl_CM) - ml2)*(1. + std::numeric_limits<double>::epsilon());
  Q2max = (2*Enu_CM*(El_CM + Pl_CM) - ml2)*(1. - std::numeric_limits<double>::epsilon());
  
  // model restrictions
  Wmin  = TMath::Max (Wmin,  1.08);
  Wmax  = TMath::Min (Wmax,  2.00);
  Q2min = TMath::Max (Q2min, 0.00);
  Q2max = TMath::Min (Q2max, 3.00);
  
  if (W < Wmin || W > Wmax)
    return false;

  if (Q2 < Q2min || Q2 > Q2max)
    return false;

  return true;

}
//____________________________________________________________________________
void DCCSPPPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DCCSPPPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DCCSPPPXSec::LoadConfig(void)
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
