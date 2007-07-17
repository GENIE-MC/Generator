//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory - May 05, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TMath.h>
#include <TSystem.h>

#include "Algorithm/AlgFactory.h"
#include "Algorithm/AlgConfigPool.h"
#include "Base/XSecIntegratorI.h"
#include "BaryonResonance/BaryonResDataSetI.h"
#include "BaryonResonance/BaryonResUtils.h"
#include "Conventions/Constants.h"
#include "Conventions/RefFrame.h"
#include "Conventions/KineVar.h"
#include "Conventions/Units.h"
#include "Messenger/Messenger.h"
#include "Numerical/Spline.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "ReinSeghal/ReinSeghalRESPXSec.h"
#include "ReinSeghal/RSHelicityAmplModelI.h"
#include "ReinSeghal/RSHelicityAmpl.h"
#include "Utils/KineUtils.h"
#include "Utils/MathUtils.h"
#include "Utils/Range1.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
ReinSeghalRESPXSec::ReinSeghalRESPXSec() :
XSecAlgorithmI("genie::ReinSeghalRESPXSec")
{
  fNuTauRdSpl    = 0;
  fNuTauBarRdSpl = 0;
}
//____________________________________________________________________________
ReinSeghalRESPXSec::ReinSeghalRESPXSec(string config) :
XSecAlgorithmI("genie::ReinSeghalRESPXSec", config)
{
  fNuTauRdSpl    = 0;
  fNuTauBarRdSpl = 0;
}
//____________________________________________________________________________
ReinSeghalRESPXSec::~ReinSeghalRESPXSec()
{
  if(fNuTauRdSpl)    delete fNuTauRdSpl;
  if(fNuTauBarRdSpl) delete fNuTauBarRdSpl;
}
//____________________________________________________________________________
double ReinSeghalRESPXSec::XSec(
                 const Interaction * interaction, KinePhaseSpace_t kps) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  const InitialState & init_state = interaction -> InitState();
  const ProcessInfo &  proc_info  = interaction -> ProcInfo();
  const Target & target = init_state.Tgt();

  //----- Get kinematical parameters
  const Kinematics & kinematics = interaction -> Kine();
  double W  = kinematics.W();
  double q2 = kinematics.q2();

  //----- Under the DIS/RES joining scheme, xsec(RES)=0 for W>=Wcut
  if(fUsingDisResJoin) {
    if(W>=fWcut) {
       LOG("ReinSeghalRes", pDEBUG)
         << "RES/DIS Join Scheme: XSec[RES, W=" << W 
                                << " >= Wcut=" << fWcut << "] = 0";
       return 0;
    }
  }

  //-- Get the input baryon resonance
  Resonance_t resonance = interaction->ExclTag().Resonance();
  string      resname   = utils::res::AsString(resonance);
  bool        is_delta  = utils::res::IsDelta (resonance);

  //-- Get the neutrino, hit nucleon & weak current
  int  nucpdgc = target.HitNucPdg();
  int  nupdgc  = init_state.ProbePdg();
  bool is_nu    = pdg::IsNeutrino     (nupdgc);
  bool is_nubar = pdg::IsAntiNeutrino (nupdgc);
  bool is_p     = pdg::IsProton       (nucpdgc);
  bool is_n     = pdg::IsNeutron      (nucpdgc);
  bool is_CC    = proc_info.IsWeakCC();

  if(is_CC && !is_delta) {
    if((is_nu && is_p) || (is_nubar && is_n)) return 0;
  }

  //-- Get baryon resonance parameters
  fBRP.RetrieveData(resonance);  
  double Mres = fBRP.Mass();
  double Gres = fBRP.Width();
  int    Nres = fBRP.ResonanceIndex();

  //-- Following NeuGEN, avoid problems with underlying unphysical
  //   model assumptions by restricting the allowed W phase space
  //   around the resonance peak
  if      (W > Mres + fN0ResMaxNWidths * Gres && Nres==0) return 0.;
  else if (W > Mres + fN2ResMaxNWidths * Gres && Nres==2) return 0.;
  else if (W > Mres + fGnResMaxNWidths * Gres)            return 0.;

  //-- Compute auxiliary & kinematical factors 
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

  LOG("ReinSeghalRes", pDEBUG) << "V = " << V << ", U = " << U;

  //-- Calculate the Feynman-Kislinger-Ravndall parameters
  LOG("ReinSeghalRes", pDEBUG) << "Computing the FKR parameters";

  double Go     = TMath::Power(1 - 0.25 * q2/Mnuc2, 0.5-Nres);
  double GV     = Go * TMath::Power( 1./(1-q2/fMv2), 2);
  double GA     = Go * TMath::Power( 1./(1-q2/fMa2), 2);

  double d      = TMath::Power(W+Mnuc,2.) - q2;
  double sq2omg = TMath::Sqrt(2./fOmega);
  double nomg   = Nres * fOmega;
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

  LOG("FKR", pDEBUG) 
           << "FKR params for RES=" << resname << " : " << fFKR;

  //-- Calculate the Rein-Seghal Helicity Amplitudes
  LOG("ReinSeghalRes", pDEBUG) << "Computing Helicity Amplitudes";

  const RSHelicityAmplModelI * hamplmod = 0;

  if(is_CC)   { hamplmod = fHAmplModelCC; }
  else {
    if (is_p) { hamplmod = fHAmplModelNCp;}
    else      { hamplmod = fHAmplModelNCn;}
  }
  assert(hamplmod);
  
  RSHelicityAmpl * hampl = hamplmod->Compute(resonance, fFKR); 
  assert(hampl);

  LOG("RSHAmpl", pDEBUG)
    << "Helicity Ampl for RES=" << resname << " : " << *hampl;

  //-- Compute the cross section

  double sig0 = 0.125*(kGF2/kPi)*(-q2/Q2)*(W/Mnuc);
  double scLR = W/Mnuc;
  double scS  = (Mnuc/W)*(-Q2/q2);
  double sigL = scLR* (hampl->Amp2Plus3 () + hampl->Amp2Plus1 ());
  double sigR = scLR* (hampl->Amp2Minus3() + hampl->Amp2Minus1());
  double sigS = scS * (hampl->Amp20Plus () + hampl->Amp20Minus());

  delete hampl;

  LOG("ReinSeghalRes", pDEBUG) << "sig_{0} = " << sig0;
  LOG("ReinSeghalRes", pDEBUG) << "sig_{L} = " << sigL;
  LOG("ReinSeghalRes", pDEBUG) << "sig_{R} = " << sigR;
  LOG("ReinSeghalRes", pDEBUG) << "sig_{S} = " << sigS;

  double xsec = 0.0;
  if (is_nu) {
     xsec = sig0*(V2*sigR + U2*sigL + 2*UV*sigS);
  } else {
     xsec = sig0*(U2*sigR + V2*sigL + 2*UV*sigS);
  }
  xsec = TMath::Max(0.,xsec);

  double mult = 1.0;
  if(is_CC && is_delta) {
    if((is_nu && is_p) || (is_nubar && is_n)) mult=3.0;
  }
  xsec *= mult;

  //-- Check whether the cross section is to be weighted with a
  //   Breit-Wigner distribution (default: true)
  double bw = 1.0;
  if(fWghtBW) {
     bw = fBreitWigner->Eval(resonance, W);
     LOG("ReinSeghalRes", pDEBUG) 
       << "BreitWigner(RES=" << resname << ", W=" << W << ") = " << bw;
  } else {
      LOG("ReinSeghalRes", pDEBUG) << "Breit-Wigner wght is turned-off";
  }
  xsec *= bw; 

  //-- Apply NeuGEN nutau cross section reduction factors
  double rf = 1.0;
  Spline * spl = 0;
  if (is_CC && fUsingNuTauScaling) {
    if      (pdg::IsNuTau(nupdgc)    ) spl = fNuTauRdSpl;
    else if (pdg::IsAntiNuTau(nupdgc)) spl = fNuTauBarRdSpl;

    if(spl) {
      if(E <spl->XMax()) rf = spl->Evaluate(E);
    }
  }
  xsec *= rf;

  LOG("ReinSeghalRes", pINFO) 
    << "\n d2xsec/dQ2dW"  << "[" << interaction->AsString()
          << "](W=" << W << ", q2=" << q2 << ", E=" << E << ") = " << xsec;

  //-- The algorithm computes d^2xsec/dWdQ2
  //   Check whether variable tranformation is needed
  if(kps!=kPSWQ2fE) {
    double J = utils::kinematics::Jacobian(interaction,kPSWQ2fE,kps);
    xsec *= J;
  }

  //-- If requested return the free nucleon xsec even for input nuclear tgt
  if( interaction->TestBit(kIAssumeFreeNucleon) ) return xsec;

  //-- number of scattering centers in the target
  int NNucl = (is_p) ? target.Z() : target.N();

  xsec*=NNucl; // nuclear xsec (no nuclear suppression factor)

  return xsec;
}
//____________________________________________________________________________
double ReinSeghalRESPXSec::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool ReinSeghalRESPXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const InitialState & init_state = interaction->InitState();
  const ProcessInfo &  proc_info  = interaction->ProcInfo();
  const XclsTag &      xcls       = interaction->ExclTag();

  if(!proc_info.IsResonant()) return false;
  if(!proc_info.IsWeak())     return false;

  int  nuc = init_state.Tgt().HitNucPdg();
  int  nu  = init_state.ProbePdg();

  if (!pdg::IsProton(nuc)  && !pdg::IsNeutron(nuc))     return false;
  if (!pdg::IsNeutrino(nu) && !pdg::IsAntiNeutrino(nu)) return false;

  if(!xcls.KnownResonance()) return false;

  return true;
}
//____________________________________________________________________________
/*
bool ReinSeghalRESPXSec::ValidKinematics(const Interaction* interaction) const
{
  if(interaction->TestBit(kISkipKinematicChk)) return true;

  const Kinematics &   kinematics = interaction -> Kine();
  const InitialState & init_state = interaction -> InitState();

  double E    = init_state.ProbeE(kRfHitNucRest);
  double W    = kinematics.W();
  double q2   = kinematics.q2();

  //-- Check energy threshold & kinematical limits in q2, W
  double EvThr = interaction->EnergyThreshold();
  if(E <= EvThr) {
    LOG("ReinSeghalRes", pINFO) << "E  = " << E << " < Ethr = " << EvThr;
    return false;
  }

  //-- Check against physical range in W and Q2
  Range1D_t rW  = utils::kinematics::KineRange(interaction, kKVW);
  Range1D_t rQ2 = utils::kinematics::KineRange(interaction, kKVQ2);

  bool in_physical_range = utils::math::IsWithinLimits(W, rW) &&
                           utils::math::IsWithinLimits(-q2, rQ2);
  if(!in_physical_range) return false;

  return true;
}*/
//____________________________________________________________________________
void ReinSeghalRESPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void ReinSeghalRESPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void ReinSeghalRESPXSec::LoadConfig(void)
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

  // Load all the sub-algorithms needed

  fBaryonResDataSet = 0;
  fHAmplModelCC     = 0;
  fHAmplModelNCp    = 0;
  fHAmplModelNCn    = 0;
  fBreitWigner      = 0;

  //-- Access a "Baryon Resonance data-set" sub-algorithm
  fBaryonResDataSet = 
         dynamic_cast<const BaryonResDataSetI *> (
                           this->SubAlg("BaryonResData"));
  assert(fBaryonResDataSet);

  fBRP.SetDataSet(fBaryonResDataSet); // <-- attach data set;

  if(fWghtBW) {
    //-- Access a "Breit-Wigner" sub-algorithm
    fBreitWigner = dynamic_cast<const BreitWignerI *> (
                             this->SubAlg("BreitWignerAlg"));
    assert(fBreitWigner);
  }

  //-- Access a "Helicity Amplitudes model" sub-algorithms

  AlgFactory * algf = AlgFactory::Instance();

  fHAmplModelCC  = dynamic_cast<const RSHelicityAmplModelI *> (
	        algf->GetAlgorithm("genie::RSHelicityAmplModelCC","Default"));
  fHAmplModelNCp = dynamic_cast<const RSHelicityAmplModelI *> (
	       algf->GetAlgorithm("genie::RSHelicityAmplModelNCp","Default"));
  fHAmplModelNCn = dynamic_cast<const RSHelicityAmplModelI *> (
	       algf->GetAlgorithm("genie::RSHelicityAmplModelNCn","Default"));

  assert( fHAmplModelCC  );
  assert( fHAmplModelNCp );
  assert( fHAmplModelNCn );

  //-- Use algorithm within a DIS/RES join scheme. If yes get Wcut
  fUsingDisResJoin = fConfig->GetBoolDef(
                           "UseDRJoinScheme", gc->GetBool("UseDRJoinScheme"));
  fWcut = 999999;
  if(fUsingDisResJoin) {
    fWcut = fConfig->GetDoubleDef("Wcut",gc->GetDouble("Wcut"));
  }

  //-- NeuGEN limits in the allowed resonance phase space:
  //           W < min{ Wmin(physical), (res mass) + x * (res width) }
  //   This limits the integration area around the peak and avoids the
  //   problem with huge xsec increase at low Q2 and high W.
  //   In correspondence with Hugh, Rein said that the underlying problem
  //   are unphysical assumptions in the model. 
  fN2ResMaxNWidths = fConfig->GetDoubleDef("MaxNWidthForN2Res", 2.0);
  fN0ResMaxNWidths = fConfig->GetDoubleDef("MaxNWidthForN0Res", 6.0);
  fGnResMaxNWidths = fConfig->GetDoubleDef("MaxNWidthForGNRes", 4.0);

  //-- NeuGEN reduction factors for nu_tau: a gross estimate of the effect of
  //   neglected form factors in the R/S model
  fUsingNuTauScaling = fConfig->GetBoolDef("UseNuTauScalingFactors", true);
  if(fUsingNuTauScaling) {
     if(fNuTauRdSpl)    delete fNuTauRdSpl;
     if(fNuTauBarRdSpl) delete fNuTauBarRdSpl;

     assert(gSystem->Getenv("GENIE"));
     string base = gSystem->Getenv("GENIE");

     string filename = base + "/data/etc/rs-res-xsec-scaling-nutau.dat";
     LOG("ReinSeghalRes", pNOTICE) 
                << "Loading nu_tau xsec reduction spline from: " << filename;
     fNuTauRdSpl = new Spline(filename);

     filename = base + "/data/etc/rs-res-xsec-scaling-nutaubar.dat";
     LOG("ReinSeghalRes", pNOTICE) 
           << "Loading bar{nu_tau} xsec reduction spline from: " << filename;
     fNuTauBarRdSpl = new Spline(filename);
  }

  //-- load the differential cross section integrator
  fXSecIntegrator =
      dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);
}
//____________________________________________________________________________

