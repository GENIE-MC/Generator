//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
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
 @ May 01, 2012 - CA
   Pick nutau/nutaubar scaling factors from new location.
 @ May 01, 2016 - Libo Jiang
   Add W dependence to Delta->N gamma

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
#include "Physics/Resonance/XSection/ReinSehgalRESPXSec.h"
#include "Physics/Resonance/XSection/RSHelicityAmplModelI.h"
#include "Physics/Resonance/XSection/RSHelicityAmpl.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Physics/NuclearState/FermiMomentumTablePool.h"
#include "Physics/NuclearState/FermiMomentumTable.h"
#include "Physics/NuclearState/NuclearUtils.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
ReinSehgalRESPXSec::ReinSehgalRESPXSec() :
XSecAlgorithmI("genie::ReinSehgalRESPXSec")
{
  fNuTauRdSpl    = 0;
  fNuTauBarRdSpl = 0;
}
//____________________________________________________________________________
ReinSehgalRESPXSec::ReinSehgalRESPXSec(string config) :
XSecAlgorithmI("genie::ReinSehgalRESPXSec", config)
{
  fNuTauRdSpl    = 0;
  fNuTauBarRdSpl = 0;
}
//____________________________________________________________________________
ReinSehgalRESPXSec::~ReinSehgalRESPXSec()
{
  if(fNuTauRdSpl)    delete fNuTauRdSpl;
  if(fNuTauBarRdSpl) delete fNuTauBarRdSpl;
}
//____________________________________________________________________________
double ReinSehgalRESPXSec::XSec(
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

  // Under the DIS/RES joining scheme, xsec(RES)=0 for W>=Wcut
  if(fUsingDisResJoin) {
    if(W>=fWcut) {
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
       LOG("ReinSehgalRes", pDEBUG)
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

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("ReinSehgalRes", pDEBUG) 
     << "Kinematical params V = " << V << ", U = " << U;
#endif

  // Calculate the Feynman-Kislinger-Ravndall parameters

  double Go  = TMath::Power(1 - 0.25 * q2/Mnuc2, 0.5-IR);
  double GV  = Go * TMath::Power( 1./(1-q2/fMv2), 2);
  double GA  = Go * TMath::Power( 1./(1-q2/fMa2), 2);

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

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("FKR", pDEBUG) 
     << "FKR params for RES = " << resname << " : " << fFKR;
#endif

  // Calculate the Rein-Sehgal Helicity Amplitudes

  const RSHelicityAmplModelI * hamplmod = 0;
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

  double g2 = kGF2;
  if(is_CC) g2 = kGF2*fVud2;
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

  double sig0 = 0.125*(g2/kPi)*(-q2/Q2)*(W/Mnuc);
  double scLR = W/Mnuc;
  double scS  = (Mnuc/W)*(-Q2/q2);
  double sigL = scLR* (hampl.Amp2Plus3 () + hampl.Amp2Plus1 ());
  double sigR = scLR* (hampl.Amp2Minus3() + hampl.Amp2Minus1());
  double sigS = scS * (hampl.Amp20Plus () + hampl.Amp20Minus());

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("ReinSehgalRes", pDEBUG) << "sig_{0} = " << sig0;
  LOG("ReinSehgalRes", pDEBUG) << "sig_{L} = " << sigL;
  LOG("ReinSehgalRes", pDEBUG) << "sig_{R} = " << sigR;
  LOG("ReinSehgalRes", pDEBUG) << "sig_{S} = " << sigS;
#endif

  double xsec = 0.0;
  if (is_nu || is_lminus) {
     xsec = sig0*(V2*sigR + U2*sigL + 2*UV*sigS);
  } 
  else 
  if (is_nubar || is_lplus) {
     xsec = sig0*(U2*sigR + V2*sigL + 2*UV*sigS);
  } 
  xsec = TMath::Max(0.,xsec);

  double mult = 1.0;
  if(is_CC && is_delta) {
    if((is_nu && is_p) || (is_nubar && is_n)) mult=3.0;
  }
  xsec *= mult;

  // Check whether the cross section is to be weighted with a
  // Breit-Wigner distribution (default: true)
  double bw = 1.0;
  if(fWghtBW) {
     //different Delta photon decay branch
     if(is_delta){
     bw = utils::bwfunc::BreitWignerLGamma(W,LR,MR,WR,NR); 
     }
     else{
     bw = utils::bwfunc::BreitWignerL(W,LR,MR,WR,NR); 
     }
  } 
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
     LOG("ReinSehgalRes", pDEBUG) 
       << "BreitWigner(RES=" << resname << ", W=" << W << ") = " << bw;
#endif
  xsec *= bw; 

  // Apply NeuGEN nutau cross section reduction factors
  double rf = 1.0;
  Spline * spl = 0;
  if (is_CC && fUsingNuTauScaling) {
    if      (pdg::IsNuTau    (probepdgc)) spl = fNuTauRdSpl;
    else if (pdg::IsAntiNuTau(probepdgc)) spl = fNuTauBarRdSpl;

    if(spl) {
      if(E <spl->XMax()) rf = spl->Evaluate(E);
    }
  }
  xsec *= rf;

  // Apply given scaling factor
  double xsec_scale = 1.;
  if      (is_CC) { xsec_scale = fXSecScaleCC; }
  else if (is_NC) { xsec_scale = fXSecScaleNC; }
  xsec *= xsec_scale;

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("ReinSehgalRes", pINFO) 
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

  
  int Z = target.Z();
  int A = target.A();
  int N = A-Z;
  
  // Take into account the number of scattering centers in the target
  int NNucl = (is_p) ? Z : N;

  xsec*=NNucl; // nuclear xsec (no nuclear suppression factor) 
  
  if (fUsePauliBlocking && A!=1)
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
     if (A<6 || !fUseRFGParametrization)
     {
         // Look up the Fermi momentum for this target
         FermiMomentumTablePool * kftp = FermiMomentumTablePool::Instance();
         const FermiMomentumTable * kft = kftp->GetTable(fKFTable);
         P_Fermi = kft->FindClosestKF(pdg::IonPdgCode(A, Z), nucpdgc);
     }
     else {
        // Define the Fermi momentum for this target
        P_Fermi = utils::nuclear::FermiMomentumForIsoscalarNucleonParametrization(target);
        // Correct the Fermi momentum for the struck nucleon
        if(is_p) { P_Fermi *= TMath::Power( 2.*Z/A, 1./3); }
        else     { P_Fermi *= TMath::Power( 2.*N/A, 1./3); }
     }
  
     double FactorPauli_RES = 1.0;
  
     double k0 = 0., q = 0., q0 = 0.;
 
     if (P_Fermi > 0.)
     {
        k0 = (W2-Mnuc2-Q2)/(2*W);
        k = TMath::Sqrt(k0*k0+Q2);                  // previous value of k is overridden
        q0 = (W2-Mnuc2+kPionMass2)/(2*W);
        q = TMath::Sqrt(q0*q0-kPionMass2);
     }
           
     if (2*P_Fermi < k-q) 
        FactorPauli_RES = 1.0;
     if (2*P_Fermi >= k+q)
        FactorPauli_RES = ((3*k*k+q*q)/(2*P_Fermi)-(5*TMath::Power(k,4)+TMath::Power(q,4)+10*k*k*q*q)/(40*TMath::Power(P_Fermi,3)))/(2*k);
     if (2*P_Fermi >= k-q && 2*P_Fermi <= k+q)
        FactorPauli_RES = ((q+k)*(q+k)-4*P_Fermi*P_Fermi/5-TMath::Power(k-q, 3)/(2*P_Fermi)+TMath::Power(k-q, 5)/(40*TMath::Power(P_Fermi, 3)))/(4*q*k);
     
     xsec *= FactorPauli_RES;
  }
  
  return xsec;
}
//____________________________________________________________________________
double ReinSehgalRESPXSec::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool ReinSehgalRESPXSec::ValidProcess(const Interaction * interaction) const
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
void ReinSehgalRESPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void ReinSehgalRESPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void ReinSehgalRESPXSec::LoadConfig(void)
{
  // Cross section scaling factors
  this->GetParam( "RES-CC-XSecScale", fXSecScaleCC ) ;
  this->GetParam( "RES-NC-XSecScale", fXSecScaleNC ) ;

  this->GetParam( "RES-Zeta", fZeta ) ;
  this->GetParam( "RES-Omega", fOmega ) ;

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

  // NeuGEN reduction factors for nu_tau: a gross estimate of the effect of
  // neglected form factors in the R/S model
  this->GetParamDef( "UseNuTauScalingFactors", fUsingNuTauScaling, true ) ;
  if(fUsingNuTauScaling) {
     if(fNuTauRdSpl)    delete fNuTauRdSpl;
     if(fNuTauBarRdSpl) delete fNuTauBarRdSpl;

     assert( std::getenv( "GENIE") );
     string base = std::getenv( "GENIE") ;

     string filename = base + "/data/evgen/rein_sehgal/res/nutau_xsec_scaling_factors.dat";
     LOG("ReinSehgalRes", pNOTICE) 
                << "Loading nu_tau xsec reduction spline from: " << filename;
     fNuTauRdSpl = new Spline(filename);

     filename = base + "/data/evgen/rein_sehgal/res/nutaubar_xsec_scaling_factors.dat";
     LOG("ReinSehgalRes", pNOTICE) 
           << "Loading bar{nu_tau} xsec reduction spline from: " << filename;
     fNuTauBarRdSpl = new Spline(filename);
  }

  // load the differential cross section integrator
  fXSecIntegrator =
      dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);
}
//____________________________________________________________________________

