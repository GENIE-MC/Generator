//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

 @ Sep 13, 2007 - CA
   Debugged the model in order to be included in the default event generation
   threads in the next physics release (2.0.2). Rather than using Kovalenko's
   expression for the ZR scaling factor, I apply an ad-hoc scaling factor 
   maintaining the relative strength of the QELC channels but lowering their 
   sum to be consistent with recent NOMAD measurement. The default value of
   M0 has been changed from 0.1 to sqrt(0.1) as in M.Bischofberger's (ETHZ)
   PhD thesis (DISS.ETH NO 16034). For more details see GENIE-PUB/2007/006.
*/
//____________________________________________________________________________

#include <TMath.h>
#include <Math/Integrator.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Physics/Charm/XSection/KovalenkoQELCharmPXSec.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"
/////#include "Numerical/IntegratorI.h"
#include "Physics/PartonDistributions/PDF.h"
#include "Physics/PartonDistributions/PDFModelI.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Numerical/GSLUtils.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
KovalenkoQELCharmPXSec::KovalenkoQELCharmPXSec() :
XSecAlgorithmI("genie::KovalenkoQELCharmPXSec")
{

}
//____________________________________________________________________________
KovalenkoQELCharmPXSec::KovalenkoQELCharmPXSec(string config) :
XSecAlgorithmI("genie::KovalenkoQELCharmPXSec", config)
{

}
//____________________________________________________________________________
KovalenkoQELCharmPXSec::~KovalenkoQELCharmPXSec()
{

}
//____________________________________________________________________________
double KovalenkoQELCharmPXSec::XSec(
                  const Interaction * interaction, KinePhaseSpace_t kps) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  //----- get kinematics & init state - compute auxiliary vars
  const Kinematics &   kinematics  = interaction->Kine();
  const InitialState & init_state  = interaction->InitState();
  const Target &       target      = init_state.Tgt();

  //neutrino energy & momentum transfer
  double E   = init_state.ProbeE(kRfHitNucRest);
  double E2  = E * E;
  double Q2  = kinematics.Q2();

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("QELCharmXSec", pDEBUG) << "E = " << E << ", Q2 = " << Q2;
#endif

  //resonance mass & nucleon mass
  double MR    = this->MRes  (interaction);
  double MR2   = TMath::Power(MR,2);
  double Mnuc  = target.HitNucMass();
  double Mnuc2 = TMath::Power(Mnuc,2);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("QELCharmXSec", pDEBUG) << "M(RES) = " << MR;
#endif

  //----- Calculate the differential cross section dxsec/dQ^2
  double Gf        = kGF2 / (2*kPi);
  double vR        = (MR2 - Mnuc2 + Q2) / (2*Mnuc);
  double xiR       = this->xiBar(Q2, Mnuc, vR);
  double vR2       = vR*vR;
  double vR_E      = vR/E;
  double Q2_4E2    = Q2/(4*E2);
  double Q2_2MExiR = Q2/(2*Mnuc*E*xiR);
  double Z         = this->ZR(interaction);
  double D         = this->DR(interaction);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("QELCharmXSec", pDEBUG) 
     << "Z = " << Z << ", D = " << D << ". xiR = " << xiR << ", vR = " << vR;
#endif

  double xsec = Gf*Z*D * (1 - vR_E + Q2_4E2 + Q2_2MExiR) *
                                     TMath::Sqrt(vR2 + Q2) / (vR*xiR);

  //----- The algorithm computes dxsec/dQ2
  //      Check whether variable tranformation is needed
  if(kps!=kPSQ2fE) {
    double J = utils::kinematics::Jacobian(interaction,kPSQ2fE,kps);
    xsec *= J;
  }

  //----- If requested return the free nucleon xsec even for input nuclear tgt
  if( interaction->TestBit(kIAssumeFreeNucleon) ) return xsec;

  //----- Nuclear cross section (simple scaling here)
  int nuc   = target.HitNucPdg();
  int NNucl = (pdg::IsProton(nuc)) ? target.Z() : target.N();
  xsec *= NNucl;

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("QELCharmXSec", pINFO) 
     << "dsigma/dQ2(E=" << E << ", Q2=" << Q2 << ") = " 
                             << xsec / (1E-40*units::cm2) << " x 1E-40 cm^2";
#endif

  return xsec;
}
//____________________________________________________________________________
double KovalenkoQELCharmPXSec::ZR(const Interaction * interaction) const
{
  const XclsTag &      xcls       = interaction->ExclTag();
  const InitialState & init_state = interaction->InitState();

  int pdgc = xcls.CharmHadronPdg();
  bool isP = pdg::IsProton ( init_state.Tgt().HitNucPdg() );
  bool isN = pdg::IsNeutron( init_state.Tgt().HitNucPdg() );

  if      ( pdgc == kPdgLambdaPc && isN ) return fScLambdaP;
  else if ( pdgc == kPdgSigmaPc  && isN ) return fScSigmaP;
  else if ( pdgc == kPdgSigmaPPc && isP ) return fScSigmaPP;
  else                                    abort();
}
//____________________________________________________________________________
double KovalenkoQELCharmPXSec::DR(const Interaction * interaction) const
{
  const InitialState & init_state = interaction -> InitState();

  // Compute PDFs
  PDF pdfs;
  pdfs.SetModel(fPDFModel);   // <-- attach algorithm

  // Compute integration area = [xi_bar_plus, xi_bar_minus]
  const Kinematics & kinematics = interaction->Kine();

  double Q2     = kinematics.Q2();
  double Mnuc   = init_state.Tgt().HitNucMass();
  double Mnuc2  = TMath::Power(Mnuc,2);
  double MR     = this->MRes(interaction);
  double DeltaR = this->ResDM(interaction);

  double vR_minus  = ( TMath::Power(MR-DeltaR,2) - Mnuc2 + Q2 ) / (2*Mnuc);
  double vR_plus   = ( TMath::Power(MR+DeltaR,2) - Mnuc2 + Q2 ) / (2*Mnuc);

  LOG("QELCharmXSec", pDEBUG)
    << "vR = [plus: " << vR_plus << ", minus: " << vR_minus << "]";

  double xi_bar_minus  = 0.999;
  double xi_bar_plus  = this->xiBar(Q2, Mnuc, vR_plus);

  LOG("QELCharmXSec", pDEBUG) 
    << "Integration limits = [" << xi_bar_plus << ", " << xi_bar_minus << "]";

  int pdgc = init_state.Tgt().HitNucPdg();

  ROOT::Math::IBaseFunctionOneDim * integrand = new 
          utils::gsl::wrap::KovQELCharmIntegrand(&pdfs,Q2,pdgc);
  ROOT::Math::IntegrationOneDim::Type ig_type = 
          utils::gsl::Integration1DimTypeFromString("adaptive");
  
  double abstol   = 1;    // We mostly care about relative tolerance
  double reltol   = 1E-4; 
  int    nmaxeval = 100000;
  ROOT::Math::Integrator ig(*integrand,ig_type,abstol,reltol,nmaxeval);
  double D = ig.Integral(xi_bar_plus, xi_bar_minus);

  delete integrand;

  return D;
}
//____________________________________________________________________________
double KovalenkoQELCharmPXSec::xiBar(double Q2, double Mnuc, double v) const
{
  double Mo2 = fMo*fMo;
  double v2  = v *v;
  double xi  = (Q2/Mnuc) / (v + TMath::Sqrt(v2+Q2));
  double xib = xi * ( 1 + (1 + Mo2/(Q2+Mo2))*Mo2/Q2 );
  return xib;
}
//____________________________________________________________________________
double KovalenkoQELCharmPXSec::ResDM(const Interaction * interaction) const
{
// Resonance Delta M obeys the constraint DM(R+/-) <= |M(R+/-) - M(R)|
// where M(R-) <= M(R) <= M(R+) are the masses of the neighboring
// resonances R+, R-.
// Get the values from the algorithm conf. registry, and if they do not exist
// set them to default values (Eq.(20) in Sov.J.Nucl.Phys.52:934 (1990)

  const XclsTag & xcls = interaction->ExclTag();

  int pdgc = xcls.CharmHadronPdg();

  bool isLambda = (pdgc == kPdgLambdaPc);
  bool isSigma  = (pdgc == kPdgSigmaPc || pdgc == kPdgSigmaPPc);

  if      ( isLambda ) return fResDMLambda;
  else if ( isSigma  ) return fResDMSigma;
  else                 abort();

  return 0;
}
//____________________________________________________________________________
double KovalenkoQELCharmPXSec::MRes(const Interaction * interaction) const
{
  const XclsTag & xcls = interaction->ExclTag();

  int pdgc  = xcls.CharmHadronPdg();
  double MR = PDGLibrary::Instance()->Find(pdgc)->Mass();
  return MR;
}
//____________________________________________________________________________
double KovalenkoQELCharmPXSec::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool KovalenkoQELCharmPXSec::ValidProcess(
                                        const Interaction * interaction) const
{
  // Make sure we are dealing with one of the following channels:
  // v + n --> mu- + Lambda_{c}^{+} (2285)
  // v + n --> mu- + Sigma_{c}^{+} (2455)
  // v + p --> mu- + Sigma_{c}^{++} (2455)

  if(interaction->TestBit(kISkipProcessChk)) return true;

  const XclsTag &      xcls       = interaction->ExclTag();
  const InitialState & init_state = interaction->InitState();
  const ProcessInfo &  proc_info  = interaction->ProcInfo();

  bool is_exclusive_charm = (xcls.IsCharmEvent() && !xcls.IsInclusiveCharm());
  if(!is_exclusive_charm) return false;

  if(!proc_info.IsQuasiElastic()) return false;
  if(!proc_info.IsWeak())         return false;

  bool isP = pdg::IsProton ( init_state.Tgt().HitNucPdg() );
  bool isN = pdg::IsNeutron( init_state.Tgt().HitNucPdg() );

  int pdgc = xcls.CharmHadronPdg();

  bool can_handle = (
     (pdgc == kPdgLambdaPc && isN) || /* v + n -> l + #Lambda_{c}^{+} */
     (pdgc == kPdgSigmaPc  && isN) || /* v + n -> l + #Sigma_{c}^{+}  */
     (pdgc == kPdgSigmaPPc && isP)    /* v + p -> l + #Sigma_{c}^{++} */
  );
  return can_handle;
}
//____________________________________________________________________________
bool KovalenkoQELCharmPXSec::ValidKinematics(
                                        const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipKinematicChk)) return true;

  const InitialState & init_state  = interaction->InitState();
  double E = init_state.ProbeE(kRfHitNucRest);

  //resonance, final state primary lepton & nucleon mass
  double MR    = this -> MRes  (interaction);
  double ml    = interaction->FSPrimLepton()->Mass();
  double Mnuc  = init_state.Tgt().HitNucP4Ptr()->M();
  double Mnuc2 = TMath::Power(Mnuc,2);

  //resonance threshold
  double ER = ( TMath::Power(MR+ml,2) - Mnuc2 ) / (2*Mnuc);

  if(E <= ER) return false;

  return true;
}
//____________________________________________________________________________
void KovalenkoQELCharmPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void KovalenkoQELCharmPXSec::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfig();
}
//____________________________________________________________________________
void KovalenkoQELCharmPXSec::LoadConfig(void)
{
  fPDFModel   = 0;
  ///fIntegrator = 0;

  // Get config values or set defaults
  GetParamDef( "Scale-LambdaP", fScLambdaP,  0.8 * 0.0102 ) ;
  GetParamDef( "Scale-SigmaP",  fScSigmaP ,  0.8 * 0.0028 ) ;
  GetParamDef( "Scale-SigmaPP", fScSigmaPP,  0.8 * 0.0159 ) ;
  GetParamDef( "Res-DeltaM-Lambda", fResDMLambda,  0.56 ) ;      //GeV
  GetParamDef( "Res-DeltaM-Sigma",  fResDMSigma,   0.20 ) ;      //GeV
  GetParamDef( "Mo",                fMo,           sqrt(0.1) );  //GeV

  // get PDF model and integrator

  fPDFModel = dynamic_cast<const PDFModelI *>(this->SubAlg("PDF-Set"));
  assert(fPDFModel);

  // load XSec Integrator
  fXSecIntegrator =
      dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);

  // load numerical integrator for integrand in diff x-section calc.
//  fIntegrator = dynamic_cast<const IntegratorI *>(this->SubAlg("Integrator"));
//  assert(fIntegrator);
}
//____________________________________________________________________________
// Auxiliary scalar function for internal integration
//____________________________________________________________________________
utils::gsl::wrap::KovQELCharmIntegrand::KovQELCharmIntegrand(
           PDF * pdf, double Q2, int nucleon_pdgc) :
ROOT::Math::IBaseFunctionOneDim()
{
  fPDF  = pdf;
  fQ2   = TMath::Max(0.3, Q2);
  fPdgC = nucleon_pdgc;
}
//____________________________________________________________________________
utils::gsl::wrap::KovQELCharmIntegrand::~KovQELCharmIntegrand()
{

}
//____________________________________________________________________________
unsigned int utils::gsl::wrap::KovQELCharmIntegrand::NDim(void) const
{
  return 1;
}
//____________________________________________________________________________
double utils::gsl::wrap::KovQELCharmIntegrand::DoEval(double xin) const
{
  if(xin<0 || xin>1) return 0.;

  fPDF->Calculate(xin, fQ2);
  bool isP = pdg::IsProton(fPdgC);
  double f = (isP) ? fPDF->DownValence() : fPDF->UpValence();

  LOG("QELCharmXSec", pDEBUG) 
        << "f(xin = " << xin << ", Q2 = " << fQ2 << ") = " << f; 

  return f;
}
//____________________________________________________________________________
ROOT::Math::IBaseFunctionOneDim * 
  utils::gsl::wrap::KovQELCharmIntegrand::Clone(void) const
{
  return new utils::gsl::wrap::KovQELCharmIntegrand(fPDF, fQ2, fPdgC);
}
//____________________________________________________________________________
