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
   Modified code to handle charged lepton scattering too.
   Also, the helicity amplitude code now returns a `const RSHelicityAmpl &'.
 @ July 23, 2010 - CA
   BaryonResParams, and BreitWignerI, BaryonResDataSetI implementations are
   now redundant. Get resonance parameters from BaryonResUtils and use the
   Breit-Weigner functions from utils::bwfunc.
 @ May 01, 2012 - CA
   Pick nutau/nutaubar scaling factors from new location.

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
#include "ReinSehgal/ReinSehgalRESPXSec.h"
#include "ReinSehgal/RSHelicityAmplModelI.h"
#include "ReinSehgal/RSHelicityAmpl.h"
#include "Utils/KineUtils.h"
#include "Utils/MathUtils.h"
#include "Utils/Range1.h"
#include "Utils/BWFunc.h"

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
     bw = utils::bwfunc::BreitWignerL(W,LR,MR,WR,NR); 
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

  // Take into account the number of scattering centers in the target
  int NNucl = (is_p) ? target.Z() : target.N();

  xsec*=NNucl; // nuclear xsec (no nuclear suppression factor)

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
  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();

  // Load all configuration data or set defaults

  fZeta  = fConfig->GetDoubleDef( "Zeta",  gc->GetDouble("RS-Zeta")  );
  fOmega = fConfig->GetDoubleDef( "Omega", gc->GetDouble("RS-Omega") );

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

