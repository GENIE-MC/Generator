//____________________________________________________________________________
/*
 Copyright (c) 2003-2018, The GENIE Collaboration
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
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/Utils/Range1.h"
#include "Framework/Utils/BWFunc.h"
#include "Physics/BoostedDarkMatter/XSection/ReinSehgalDMRESPXSec.h"
#include "Physics/Resonance/XSection/RSHelicityAmplModelI.h"
#include "Physics/Resonance/XSection/RSHelicityAmpl.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Physics/NuclearState/FermiMomentumTablePool.h"
#include "Physics/NuclearState/FermiMomentumTable.h"
#include "Physics/NuclearState/NuclearUtils.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
ReinSehgalDMRESPXSec::ReinSehgalDMRESPXSec() :
XSecAlgorithmI("genie::ReinSehgalDMRESPXSec")
{
  fNuTauRdSpl    = 0;
  fNuTauBarRdSpl = 0;
}
//____________________________________________________________________________
ReinSehgalDMRESPXSec::ReinSehgalDMRESPXSec(string config) :
XSecAlgorithmI("genie::ReinSehgalDMRESPXSec", config)
{
  fNuTauRdSpl    = 0;
  fNuTauBarRdSpl = 0;
}
//____________________________________________________________________________
ReinSehgalDMRESPXSec::~ReinSehgalDMRESPXSec()
{
  if(fNuTauRdSpl)    delete fNuTauRdSpl;
  if(fNuTauBarRdSpl) delete fNuTauBarRdSpl;
}
//____________________________________________________________________________
double ReinSehgalDMRESPXSec::XSec(
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
  bool is_dm     = pdg::IsDarkMatter       (probepdgc);
  bool is_dmbar  = pdg::IsAntiDarkMatter   (probepdgc);
  bool is_p      = pdg::IsProton  (nucpdgc);
  bool is_n      = pdg::IsNeutron (nucpdgc);
  bool is_DM     = proc_info.IsDarkMatter();

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
  double ml     = init_state.GetProbeP4(kRfHitNucRest)->M();
  double Mnuc   = target.HitNucMass(); 
  double E2     = TMath::Power(E,    2);
  double ml2    = TMath::Power(ml,   2);
  double W2     = TMath::Power(W,    2);
  double Mnuc2  = TMath::Power(Mnuc, 2);
  double p2     = E2 - ml2;
  double p      = TMath::Sqrt(p2);
  double k      = 0.5 * (W2 - Mnuc2)/Mnuc;
  double v      = k - 0.5 * q2/Mnuc;
  double v2     = TMath::Power(v, 2);
  double Q2     = v2 - q2;
  double Q      = TMath::Sqrt(Q2);
  double Eprime = E - v;
  double U      = 0.5 * (E + Eprime + Q) / p;
  double V      = 0.5 * (E + Eprime - Q) / p;
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
  if(is_DM) { 
    if (is_p) { hamplmod = fHAmplModelDMp;}
    else      { hamplmod = fHAmplModelDMn;}
  }
  assert(hamplmod);
  
  const RSHelicityAmpl & hampl = hamplmod->Compute(resonance, fFKR); 

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("RSHAmpl", pDEBUG)
     << "Helicity Amplitudes for RES = " << resname << " : " << hampl;
#endif

  double gZp4 = TMath::Power( fgZp, 4);
  double MZ2 = TMath::Power(fMedMass,2);
  double prop = 1. / (q2 - MZ2);
  double prop2 = TMath::Power(prop,2);
  double g2 = gZp4 * prop2;

  // Compute the cross section

  double sig0 = 0.03125*(g2/kPi)*(-q2/Q2)*(W/Mnuc);
  double scLR = W/Mnuc;
  double scS  = (Mnuc/W)*(-Q2/q2);
  double sigL = scLR* (hampl.Amp2Plus3 () + hampl.Amp2Plus1 ());
  double sigR = scLR* (hampl.Amp2Minus3() + hampl.Amp2Minus1());
  double sigS = scS * (hampl.Amp20Plus () + hampl.Amp20Minus());

//#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("ReinSehgalRes", pDEBUG) << "sig_{0} = " << -sig0 * Q2 / q2;
  LOG("ReinSehgalRes", pDEBUG) << "sig_{L} = " << sigL;
  LOG("ReinSehgalRes", pDEBUG) << "sig_{R} = " << sigR;
  LOG("ReinSehgalRes", pDEBUG) << "sig_{S} = " << sigS;
//#endif

  double xsec = 0.0;
  double QL2 = TMath::Power(fQDV - fQDA,2);
  double QR2 = TMath::Power(fQDV + fQDA,2);
  double QVpA = TMath::Power(fQDV,2) + TMath::Power(fQDA,2);
  double QVmA = TMath::Power(fQDV,2) - TMath::Power(fQDA,2);
  LOG("ReinSehgalRes", pDEBUG) << "MZ = " << fMedMass;
  LOG("ReinSehgalRes", pDEBUG) << "q2 = " << q2 << ", W = " << W2 << " -> sigmaR = " << sig0*((V2*QL2 + U2*QR2 + 2 * ml2 / (E2 - ml2) * QVmA * Q2 / q2)*sigR); 
  xsec = sig0*((V2*QL2 + U2*QR2 + 2 * ml2 / (E2 - ml2) * QVmA * Q2 / q2)*sigR + (U2*QL2 + V2*QR2 + 2 * ml2 / (E2 - ml2) * QVmA * Q2 / q2)*sigL + ((4*UV + 2 * ml2 / q2) * QVpA - 2 * ml2 / (E2 - ml2) * QVmA * Q2 / q2)*sigS);
  xsec = TMath::Max(0.,xsec);

  double mult = 1.0;
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
  xsec *= rf;

  // Apply given scaling factor
  double xsec_scale = 1.;
  if (is_DM) { xsec_scale = fXSecScaleDM; }
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
double ReinSehgalDMRESPXSec::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool ReinSehgalDMRESPXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const InitialState & init_state = interaction->InitState();
  const ProcessInfo &  proc_info  = interaction->ProcInfo();
  const XclsTag &      xcls       = interaction->ExclTag();

  if(!proc_info.IsDarkMatterResonant()) return false;
  if(!xcls.KnownResonance())  return false;

  int  hitnuc = init_state.Tgt().HitNucPdg();
  bool is_pn = (pdg::IsProton(hitnuc) || pdg::IsNeutron(hitnuc));

  if (!is_pn) return false;

  int  probe   = init_state.ProbePdg();
  bool is_dm = proc_info.IsDarkMatter(); 
  bool dm_scat = (pdg::IsDarkMatter(probe) && is_dm);
  
  if (!dm_scat) {
	  LOG("ReinSehgalRes", pWARN) << "Not dark matter but " << probe << " and " << is_dm;
  }

  if (!dm_scat) return false;

  return true;
}
//____________________________________________________________________________
void ReinSehgalDMRESPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void ReinSehgalDMRESPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void ReinSehgalDMRESPXSec::LoadConfig(void)
{
  // Cross section scaling factors
  this->GetParam( "RES-NC-XSecScale", fXSecScaleDM ) ;

  this->GetParam( "RES-Zeta", fZeta ) ;
  this->GetParam( "RES-Omega", fOmega ) ;

  // velocity dependence of interaction
  this->GetParamDef("velocity-mode", fVelMode, 0 );

  // mediator coupling
  this->GetParam("ZpCoupling", fgZp ) ;

  // mediator mass
  fMedMass = PDGLibrary::Instance()->Find(kPdgMediator)->Mass();

  // dark matter charges
  double QDL, QDR;
  this->GetParam("DarkLeftCharge", QDL);
  this->GetParam("DarkRightCharge", QDR);
  fQDV = 0.5*(QDL + QDR);
  fQDA = 0.5*(-QDL + QDR);

  double ma, mv ;
  this->GetParam( "RES-Ma", ma ) ;
  this->GetParam( "RES-Mv", mv ) ;
  fMa2 = TMath::Power(ma,2);
  fMv2 = TMath::Power(mv,2);

  this->GetParamDef( "BreitWignerWeight", fWghtBW, true ) ;
  this->GetParamDef( "BreitWignerNorm",   fNormBW, true);

  this->GetParam("FermiMomentumTable", fKFTable);
  this->GetParam("RFG-UseParametrization", fUseRFGParametrization);
  this->GetParam("UsePauliBlockingForRES", fUsePauliBlocking);
  
  // Load all the sub-algorithms needed

  fHAmplModelDMp    = 0;
  fHAmplModelDMn    = 0;

  AlgFactory * algf = AlgFactory::Instance();

  fHAmplModelDMp = dynamic_cast<const RSHelicityAmplModelI *> (
      algf->GetAlgorithm("genie::RSHelicityAmplModelDMp","Default"));
  fHAmplModelDMn = dynamic_cast<const RSHelicityAmplModelI *> (
      algf->GetAlgorithm("genie::RSHelicityAmplModelDMn","Default"));

  assert( fHAmplModelDMp );
  assert( fHAmplModelDMn );

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

