//____________________________________________________________________________
/*
 Copyright (c) 2003-2023, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Zachary W. Orr <zachorrt@colostate.edu>
 Colorado State University

based on code by
 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory
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
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/Utils/Range1.h"
#include "Framework/Utils/BWFunc.h"
#include "Physics/BoostedDarkMatter/XSection/DMRESPXSec.h"
#include "Physics/BoostedDarkMatter/XSection/RSHelicityAmplModelDMI.h"
#include "Physics/BoostedDarkMatter/XSection/RSHelicityAmplDM.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Physics/NuclearState/FermiMomentumTablePool.h"
#include "Physics/NuclearState/FermiMomentumTable.h"
#include "Physics/NuclearState/NuclearUtils.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
DMRESPXSec::DMRESPXSec() :
XSecAlgorithmI("genie::DMRESPXSec")
{

}
//____________________________________________________________________________
DMRESPXSec::DMRESPXSec(string config) :
XSecAlgorithmI("genie::DMRESPXSec", config)
{

}
//____________________________________________________________________________
DMRESPXSec::~DMRESPXSec()
{

}
//____________________________________________________________________________
double DMRESPXSec::XSec(
                 const Interaction * interaction, KinePhaseSpace_t kps) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  const InitialState & init_state = interaction -> InitState();
  const ProcessInfo &  proc_info  = interaction -> ProcInfo();
  const Target & target = init_state.Tgt();

  const Kinematics & kinematics = interaction -> Kine();
  LOG("DMRESPXSec", pDEBUG) << "Using v^" << fVelMode << " dependence";

//Get Kinematic Parameters:
  double W  = kinematics.W();
  double W2     = TMath::Power(W,    2);       //W^2
  double q2 = kinematics.q2();

//Rejection Sampling from Kinematics:
// Under the DIS/RES joining scheme, xsec(RES)=0 for W>=Wcut
  if(fUsingDisResJoin) {
    if(W>=fWcut) {
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
       LOG("DMRES", pDEBUG)
         << "RES/DIS Join Scheme: XSec[RES, W=" << W
         << " >= Wcut=" << fWcut << "] = 0";
#endif
       return 0;
    }
  }
//Pass: Resonance XSec

//Probe Res Interaction:
// Get the input baryon resonance:
  Resonance_t resonance = interaction->ExclTag().Resonance();
  string      resname   = utils::res::AsString(resonance);
  bool        is_delta  = utils::res::IsDelta (resonance);

// Get the dark matter, target nucleon & DM current:
  int  nucpdgc   = target.HitNucPdg();
  int  probepdgc = init_state.ProbePdg();
  bool is_dm     = pdg::IsDarkMatter         (probepdgc);
  bool is_dmbar  = pdg::IsAntiDarkMatter    (probepdgc);


  bool is_p      = pdg::IsProton  (nucpdgc);
  bool is_n      = pdg::IsNeutron (nucpdgc);
  bool is_dmres     = proc_info.IsDarkMatterResonant();

  //Need break for is not DM
    if(!is_dmres) {
    return 0;
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

  //Auxillary Kinematic factors
  double E      = init_state.ProbeE(kRfHitNucRest);   //E1
  double E2 = TMath::Power(E,2); //E1^2
  double Mnuc   = target.HitNucMass();   //Nucleon mass
  double Mnuc2  = TMath::Power(Mnuc, 2);       //m^2

  //Scalar Charge
  double fQchiS2 = TMath::Power(fQchiS,2); //QS^2
  // L and R charges
  double fQchiL = fQchiV-fQchiA;
  double fQchiR = fQchiV+fQchiA;

  // Compute kinematic factors
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


// Masses
  double mchi   = init_state.GetProbeP4(kRfHitNucRest)->M(); //DM mass: mx
  double mchi2 = TMath::Power(mchi, 2);        // mx^2
  double mZprime2 = TMath::Power(fMedMass, 2); // mZ'^2
  double mZprime4 = TMath::Power(mZprime2, 2); //mZ'^4
//DM charge auxillary factors
  double fQchiLmR = fQchiL-fQchiR; // QL - QR
  double fQchiLR = fQchiL * fQchiR; // QL*QR
  double fQchiLmR2 = TMath::Power(fQchiLmR, 2); //(QL - QR)^2
  double fQchiL2 = TMath::Power(fQchiL, 2); //QL^2
  double fQchiR2 = TMath::Power(fQchiR, 2); //QR^2
//DM mass auxillary factors
  double mZ_q2 = TMath::Power(mZprime2 - q2, 2)/mZprime4; //(mZ'^2 - q^2)^2 / mZ'^4
  double mchiTerm = (mchi2 * Q2)/(E2 * q2); //mchi^2 * Q^2 / E^2 * q^2

 //Spin average over incident DM = 0.5
  double AL = U2*fQchiL2 + V2*fQchiR2 + 2*mchiTerm*fQchiLR;
  double AR = V2*fQchiL2 + U2*fQchiR2 + 2*mchiTerm*fQchiLR;
  double AS = 2*UV*(fQchiL2 + fQchiR2) + mchiTerm*fQchiLmR2;
  //DM term
  double AZ = 0.0;
  if (mZprime2 != 0.0){
    AZ = -mchiTerm*mZ_q2*fQchiLmR2;
  }

  //For Scalar DM: spin average = 1
 if(fVelMode == 2){
  AL = fQchiS2 * (2*UV +  2*mchiTerm);
  AR = fQchiS2 * (2*UV +  2*mchiTerm);
  AS = fQchiS2 * (TMath::Power((U+V), 2));
  AZ = 0.0;
}
  // Calculate the Feynman-Kislinger-Ravndall parameters
  double Go  = TMath::Power(1 - 0.25 * q2/Mnuc2, 0.5-IR);
  double GV  = Go * TMath::Power( 1./(1-q2/fMv2), 2);
  double GA  = Go * TMath::Power( 1./(1-q2/fMa2), 2);
  if(mZprime2 == 0.){
    GA = 0;
  }

  double d      = TMath::Power(W+Mnuc,2.) - q2;
  double d2     = W2 - Mnuc2 + q2;
  double sq2omg = TMath::Sqrt(2./fOmega);
  double nomg   = IR * fOmega;
  double mq_w   = Mnuc*Q/W;
  double BC     = (2*W*mq_w) / (W2 - Mnuc2 + q2);
  fFKRDM.Lamda  = sq2omg * mq_w;
  fFKRDM.Tv     = GV / (3.*W*sq2omg);
  fFKRDM.Rv     = kSqrt2 * mq_w*(W+Mnuc)*GV / d;
  fFKRDM.S      = (-q2/Q2) * (3*W*Mnuc + q2 - Mnuc2) * GV / (6*Mnuc2);
  fFKRDM.Ta     = (2./3.) * (fZeta/sq2omg) * mq_w * GA / d;
  fFKRDM.Ra     = (kSqrt2/6.) * fZeta * (GA/W) * (W+Mnuc + 2*nomg*W/d );
  fFKRDM.Bs     = fZeta/(3.*W*sq2omg) * (1 + (W2-Mnuc2+q2)/ d) * GA;
  fFKRDM.Cs     = fZeta/(6.*Q) * (W2 - Mnuc2 + nomg*(W2-Mnuc2+q2)/d) * (GA/Mnuc);
  fFKRDM.Bz     = ((fZeta * GA) * (W - Mnuc)) / (2 * W * mq_w * sq2omg);
  fFKRDM.Cz     = ((fZeta * GA)/3) * (1 + ((3*q2 + nomg) / d));

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("FKR", pDEBUG)
     << "FKR params for DMRES = " << resname << " : " << fFKRDM;
#endif

  // Calculate the Rein-Sehgal Helicity Amplitudes

  const RSHelicityAmplModelDMI * hamplmod = 0;
  if(is_dmres) {
      if (is_p) { hamplmod = fHAmplModelDMp;}
      else      { hamplmod = fHAmplModelDMn;}
    }
  assert(hamplmod);

  const RSHelicityAmplDM & hampl = hamplmod->Compute(resonance, fFKRDM);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("DMHAmpl", pDEBUG)
     << "Helicity Amplitudes for DMRES = " << resname << " : " << hampl;
#endif

  // Compute the cross section structure factors
  double gZp4 = TMath::Power(fgZp, 4.0);
  double XoEn = 1.0-(mchi2/E2); //1 - (mchi/E)^2
  double XoPropMass = TMath::Power(q2 - mZprime2, 2); //(q^2 - mZ'^2)^2
  double Xo = gZp4/(XoEn*XoPropMass);
  double sig0 = 0.0625*(Xo/kPi)*(-q2/Q2)*(W/Mnuc);
  double scLR = 0.5*W/Mnuc;
  double scS  = 0.5*(Mnuc/W)*(Q2/(-q2));
  double sigL = 0.0;
  double sigR = 0.0;
  double sigS = 0.0;
  double sigZ = 0.0;

//Including Hadron spin avg and res-rest frame transformations
  sigL = scLR* (hampl.Amp2Plus3() + hampl.Amp2Plus1());
  sigR = scLR* (hampl.Amp2Minus3 () + hampl.Amp2Minus1 ());
  sigS = scS * (hampl.Amp20Plus () + hampl.Amp20Minus());
  sigZ = scS *(hampl.Ampz20Plus () + hampl.Ampz20Minus());

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("DMRES", pDEBUG) << "sig_{0} = " << sig0;
  LOG("DMRES", pDEBUG) << "sig_{L} = " << sigL;
  LOG("DMRES", pDEBUG) << "sig_{R} = " << sigR;
  LOG("DMRES", pDEBUG) << "sig_{S} = " << sigS;
  LOG("DMRES", pDEBUG) << "sig_{Z} = " << sigZ;
#endif

//DM Cross Section Calculation:

//XSec calc:
  double xsec = 0.0;
  if (is_dm) {
    xsec = sig0*(AL*sigL + AR*sigR + AS*sigS + AZ*sigZ);
  }
  else
  if (is_dmbar) {
    xsec = sig0*(AL*sigR + AR*sigL + AS*sigS + AZ*sigZ);
  }
  xsec = TMath::Max(0.,xsec);


//________________________________________________________________________________________
//X-Sec scaling by BreitWigner,Jacobian,free-nucleon,Pauli-Blocking
//________________________________________________________________________________________


  // Check whether the cross section is to be weighted with a
  // Breit-Wigner distribution (default: true)
  double bw = 1.0;
  if(fWghtBW) {
     bw = utils::bwfunc::BreitWignerL(W,LR,MR,WR,NR);

  }
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
     LOG("DMRES", pDEBUG)
       << "BreitWigner(RES=" << resname << ", W=" << W << ") = " << bw;
#endif
  xsec *= bw;


  //Apply given scaling factor
  double xsec_scale = 1.;
  if      (is_dmres) { xsec_scale = fXSecScaleDM; }
  xsec *= xsec_scale;

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("DMRES", pINFO)
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
//________________________________________________________________________________________
//________________________________________________________________________________________


//____________________________________________________________________________
double DMRESPXSec::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool DMRESPXSec::ValidProcess(const Interaction * interaction) const
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
  bool is_dmres = proc_info.IsDarkMatterResonant();
  bool is_dm   = proc_info.IsDarkMatter();

  if (!is_dm && !is_dmres) return false;

  return true;
}
//____________________________________________________________________________
void DMRESPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DMRESPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DMRESPXSec::LoadConfig(void)
{
  // Cross section scaling factor
  this->GetParam( "RES-DM-XSecScale", fXSecScaleDM ) ;

  //FKR Params
  this->GetParam( "RES-Zeta", fZeta ) ;
  this->GetParam( "RES-Omega", fOmega ) ;

  // velocity dependence of interaction: fermion/scalar DM
  this->GetParamDef("velocity-mode", fVelMode, 0 );
  //DM mass
  double ma, mv ;
  this->GetParam( "RES-Ma", ma ) ;
  this->GetParam( "RES-Mv", mv ) ;
  fMa2 = TMath::Power(ma,2);
  fMv2 = TMath::Power(mv,2);
  //DM Charges
  double QchiL, QchiR;
  this->GetParam( "DarkLeftCharge", QchiL ) ;
  this->GetParam( "DarkRightCharge", QchiR ) ;
  this->GetParam( "DarkScalarCharge", fQchiS ) ;
  fQchiV = 0.5*(QchiL + QchiR);
  fQchiA = 0.5*(- QchiL + QchiR);
  // DM Mediator
  this->GetParam("ZpCoupling", fgZp ) ;
  fMedMass = PDGLibrary::Instance()->Find(kPdgMediator)->Mass();

  this->GetParamDef( "BreitWignerWeight", fWghtBW, true ) ;
  this->GetParamDef( "BreitWignerNorm",   fNormBW, true);


  this->GetParam("FermiMomentumTable", fKFTable);
  this->GetParam("RFG-UseParametrization", fUseRFGParametrization);
  this->GetParam("UsePauliBlockingForRES", fUsePauliBlocking);

  // Load all the sub-algorithms needed

  fHAmplModelDMp    = 0;
  fHAmplModelDMn    = 0;

  AlgFactory * algf = AlgFactory::Instance();

  fHAmplModelDMp = dynamic_cast<const RSHelicityAmplModelDMI *> (
      algf->GetAlgorithm("genie::RSHelicityAmplModelDMp","Default")); //DM + p
  fHAmplModelDMn = dynamic_cast<const RSHelicityAmplModelDMI *> (
      algf->GetAlgorithm("genie::RSHelicityAmplModelDMn","Default")); //DM + n

  assert( fHAmplModelDMp );
  assert( fHAmplModelDMn );

  // Use algorithm within a DIS/RES join scheme. If yes get Wcut
  this->GetParam( "UseDRJoinScheme", fUsingDisResJoin ) ;
  fWcut = 999999;
  if(fUsingDisResJoin) {
    this->GetParam( "Wcut", fWcut ) ;
  }

  double thw ;
  this->GetParam( "WeinbergAngle", thw ) ;
  fSin48w = TMath::Power( TMath::Sin(thw), 4 );
  // NeuGEN limits in the allowed resonance phase space:
  // W < min{ Wmin(physical), (res mass) + x * (res width) }
  // It limits the integration area around the peak and avoids the
  // problem with huge xsec increase at low Q2 and high W.
  // In correspondence with Hugh, Rein said that the underlying problem
  // are unphysical assumptions in the model.
  this->GetParamDef( "MaxNWidthForN2Res", fN2ResMaxNWidths, 2.0 ) ;
  this->GetParamDef( "MaxNWidthForN0Res", fN0ResMaxNWidths, 6.0 ) ;
  this->GetParamDef( "MaxNWidthForGNRes", fGnResMaxNWidths, 4.0 ) ;

//______________________________________________________________________________________


  // load the differential cross section integrator
  fXSecIntegrator =
      dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);
}
//____________________________________________________________________________
