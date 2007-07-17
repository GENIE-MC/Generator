//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory - June 10, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Algorithm/AlgConfigPool.h"
#include "Base/XSecIntegratorI.h"
#include "Charm/KovalenkoQELCharmPXSec.h"
#include "Conventions/Constants.h"
#include "Conventions/RefFrame.h"
#include "Conventions/Units.h"
#include "Messenger/Messenger.h"
#include "Numerical/IntegratorI.h"
#include "PDF/PDF.h"
#include "PDF/PDFModelI.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGLibrary.h"
#include "Utils/KineUtils.h"

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

  //resonance mass & nucleon mass
  double MR    = this->MRes  (interaction);
  double MR2   = TMath::Power(MR,2);
  double Mnuc  = target.HitNucMass();
  double Mnuc2 = TMath::Power(Mnuc,2);

  //----- Calculate the differential cross section dxsec/dQ^2
  double Gf        = kGF2 / (2*kPi);
  double vR        = (MR2 - Mnuc2 + Q2) / (2*Mnuc);
  double xiR       = this->xiBar(interaction, vR);
  double vR2       = vR*vR;
  double vR_E      = vR/E;
  double Q2_4E2    = Q2/(4*E2);
  double Q2_2MExiR = Q2/(2*Mnuc*E*xiR);
  double Z         = this->ZR(interaction);
  double D         = this->DR(interaction);

  LOG("QELCharmXSec", pDEBUG) << "Z = " << Z << ", D = " << D;

  double xsec = Gf*Z*D * (1 - vR_E + Q2_4E2 + Q2_2MExiR) *
                                     TMath::Sqrt(vR2 + Q2) / (vR*xiR);

  LOG("QELCharmXSec", pDEBUG) 
    << "dxsec/dQ^2[QELCharm,FreeN] (E="<< E << ", Q2=" << Q2<< ") = "<< xsec;

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

  return xsec;
}
//____________________________________________________________________________
double KovalenkoQELCharmPXSec::ZR(const Interaction * interaction) const
{
  const InitialState & init_state  = interaction->InitState();

  double Mo2   = fMo*fMo;
  double Mnuc  = init_state.Tgt().HitNucMass();
  double Mnuc2 = TMath::Power(Mnuc,2);
  double MR    = this->MRes(interaction);
  double MR2   = TMath::Power(MR,2.);
  double D0    = this->DR(interaction, true); // D^R(Q^2=0)
  double sumF2 = this->SumF2(interaction);    // FA^2+F1^2

  double Z  = 2*Mo2*fSin8c2 * sumF2 / (D0 * (MR2-Mnuc2));
  return Z;
}
//____________________________________________________________________________
double KovalenkoQELCharmPXSec::DR(
                             const Interaction * interaction, bool norm) const
{
  const InitialState & init_state = interaction -> InitState();

  //----- compute PDFs
  PDF pdfs;
  pdfs.SetModel(fPDFModel);   // <-- attach algorithm

  //----- compute integration area = [xi_bar_plus, xi_bar_minus]
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

  double xi_bar_minus = this->xiBar(interaction, vR_minus);
  double xi_bar_plus  = this->xiBar(interaction, vR_plus);

  LOG("QELCharmXSec", pDEBUG) << "Integration limits = ["
                             << xi_bar_plus << ", " << xi_bar_minus << "]";

  int pdgc = init_state.Tgt().HitNucPdg();

  KovQELCharmIntegrand * func = new KovQELCharmIntegrand(&pdfs,Q2,pdgc,norm);
  func->SetParam(0,"xi_bar",xi_bar_minus, xi_bar_plus);
  double D = fIntegrator->Integrate(*func);

  delete func;
  return D;
}
//____________________________________________________________________________
double KovalenkoQELCharmPXSec::xiBar(
                              const Interaction * interaction, double v) const
{
  const InitialState & init_state = interaction -> InitState();
  const Kinematics &   kinematics = interaction -> Kine();

  double Q2     = kinematics.Q2();
  double Mnuc   = init_state.Tgt().HitNucMass();
  double Mo2    = fMo*fMo;
  double v2     = v *v;

  LOG("QELCharmXSec", pDEBUG)
                     << "Q2 = " << Q2 << ", Mo = " << fMo << ", v = " << v;

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
double KovalenkoQELCharmPXSec::vR_minus(const Interaction * interaction) const
{
  const InitialState & init_state = interaction -> InitState();
  const Kinematics & kinematics = interaction->Kine();

  double Q2  = kinematics.Q2();
  double dR  = this->ResDM(interaction);
  double MR  = MRes(interaction);
  double MN  = init_state.Tgt().HitNucP4Ptr()->M();
  double MN2 = TMath::Power(MN,2);
  double vR  = (TMath::Power(MR-dR,2) - MN2 + Q2) / (2*MN);
  return vR;
}
//____________________________________________________________________________
double KovalenkoQELCharmPXSec::vR_plus(const Interaction * interaction) const
{
  const InitialState & init_state = interaction -> InitState();
  const Kinematics & kinematics = interaction->Kine();

  double Q2  = kinematics.Q2();
  double dR  = this->ResDM(interaction);
  double MR  = MRes(interaction);
  double MN  = init_state.Tgt().HitNucP4Ptr()->M();
  double MN2 = TMath::Power(MN,2);
  double vR  = (TMath::Power(MR+dR,2) - MN2 + Q2) / (2*MN);
  return vR;
}
//____________________________________________________________________________
double KovalenkoQELCharmPXSec::SumF2(const Interaction * interaction) const
{
// Returns F1^2 (Q^2=0) + FA^2 (Q^2 = 0) for the normalization factor.
// Get the values from the algorithm conf. registry, and if they do not exist
// set them to default values I computed using Sov.J.Nucl.Phys.52:934 (1990)

  const XclsTag &      xcls       = interaction->ExclTag();
  const InitialState & init_state = interaction->InitState();

  int pdgc = xcls.CharmHadronPdg();
  bool isP = pdg::IsProton ( init_state.Tgt().HitNucPdg() );
  bool isN = pdg::IsNeutron( init_state.Tgt().HitNucPdg() );

  if      ( pdgc == kPdgLambdaPc && isN ) return fF2LambdaP;
  else if ( pdgc == kPdgSigmaPc  && isN ) return fF2SigmaP;
  else if ( pdgc == kPdgSigmaPPc && isP ) return fF2SigmaPP;
  else                                    abort();

  return 0;
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
  //----- make sure we are dealing with one of the following channels:
  //
  //   v + n --> mu- + Lambda_{c}^{+} (2285)
  //   v + n --> mu- + Sigma_{c}^{+} (2455)
  //   v + p --> mu- + Sigma_{c}^{++} (2455)

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
  fIntegrator = 0;

  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();

  // Get config values or set defaults
  fF2LambdaP   = fConfig->GetDoubleDef("F1^2+FA^2-LambdaP", 2.07);
  fF2SigmaP    = fConfig->GetDoubleDef("F1^2+FA^2-SigmaP",  0.71);
  fF2SigmaPP   = fConfig->GetDoubleDef("F1^2+FA^2-SigmaPP", 1.42);
  fResDMLambda = fConfig->GetDoubleDef("Res-DeltaM-Lambda", 0.56); /*GeV*/
  fResDMSigma  = fConfig->GetDoubleDef("Res-DeltaM-Sigma",  0.20); /*GeV*/

  // 'proper scale of internal nucleon dynamics'.
  // In the original paper Mo = 0.08 +/- 0.02 GeV.
  fMo = fConfig->GetDoubleDef("Mo", 0.1);

  // cabbibo angle
  double thc = fConfig->GetDoubleDef(
                       "CabbiboAngle", gc->GetDouble("CabbiboAngle"));
  fSin8c2 = TMath::Power(TMath::Sin(thc), 2);

  // get PDF model and integrator

  fPDFModel = dynamic_cast<const PDFModelI *>(this->SubAlg("PDF-Set"));
  assert(fPDFModel);

  // load XSec Integrator
  fXSecIntegrator =
      dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);

  // load numerical integrator for integrand in diff x-section calc.
  fIntegrator = dynamic_cast<const IntegratorI *>(this->SubAlg("Integrator"));
  assert(fIntegrator);
}
//____________________________________________________________________________
// Auxiliary scalar function for internal integration
//____________________________________________________________________________
KovQELCharmIntegrand::KovQELCharmIntegrand(
                          PDF * pdf, double Q2, int nucleon_pdgc, bool norm) :
GSFunc(1)
{
  fPDF  = pdf;
  fQ2   = Q2;
  fPdgC = nucleon_pdgc;
  fNorm = norm;
}
//____________________________________________________________________________
KovQELCharmIntegrand::~KovQELCharmIntegrand()
{

}
//____________________________________________________________________________
double KovQELCharmIntegrand::operator() (const vector<double> & param)
{
  double t = param[0];
  bool isP = pdg::IsProton(fPdgC);

  if(t<0 || t>1) return 0.;
  else {
       if(fNorm) fPDF->Calculate(t, 0.);
       else      fPDF->Calculate(t, fQ2);
  }
  double f = (isP) ? ( fPDF->DownValence() + fPDF->DownSea() ):
                     ( fPDF->UpValence()   + fPDF->UpSea()   );
  return f;
}
//____________________________________________________________________________


