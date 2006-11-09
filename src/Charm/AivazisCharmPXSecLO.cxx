//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - June 10, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Algorithm/AlgConfigPool.h"
#include "Base/XSecIntegratorI.h"
#include "Charm/AivazisCharmPXSecLO.h"
#include "Conventions/Constants.h"
#include "Conventions/RefFrame.h"
#include "Conventions/KineVar.h"
#include "Conventions/Units.h"
#include "Messenger/Messenger.h"
#include "PDF/PDF.h"
#include "PDF/PDFModelI.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGCodes.h"
#include "Utils/MathUtils.h"
#include "Utils/KineUtils.h"
#include "Utils/Range1.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
AivazisCharmPXSecLO::AivazisCharmPXSecLO() :
XSecAlgorithmI("genie::AivazisCharmPXSecLO")
{

}
//____________________________________________________________________________
AivazisCharmPXSecLO::AivazisCharmPXSecLO(string config) :
XSecAlgorithmI("genie::AivazisCharmPXSecLO", config)
{

}
//____________________________________________________________________________
AivazisCharmPXSecLO::~AivazisCharmPXSecLO()
{

}
//____________________________________________________________________________
double AivazisCharmPXSecLO::XSec(
                  const Interaction * interaction, KinePhaseSpace_t kps) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  //----- get init-state & kinematical parameters
  const Kinematics &   kinematics = interaction -> Kine();
  const InitialState & init_state = interaction -> InitState();
  const Target &       target     = init_state.Tgt();
  
  //----- get target information (hit nucleon and quark)
  int  nuc  = target.HitNucPdg();
  bool isP  = pdg::IsProton (nuc);
  bool isN  = pdg::IsNeutron(nuc);
  bool qset = target.HitQrkIsSet();
  int  qpdg = (qset) ? target.HitQrkPdg()   : 0;
  bool sea  = (qset) ? target.HitSeaQrk()   : false;
  bool isd  = (qset) ? pdg::IsDQuark (qpdg) : false;
  bool iss  = (qset) ? pdg::IsSQuark (qpdg) : false;

  if (pdg::IsAntiQuark(qpdg)) return 0.; // prevent qbar -> charm

  //----- compute kinematic & auxiliary parameters
  double E           = init_state.ProbeE(kRfHitNucRest);
  double x           = kinematics.x();
  double y           = kinematics.y();
  double x2          = TMath::Power(x,    2);
  double Mnuc        = target.HitNucMass();
  double Mnuc2       = TMath::Power(Mnuc, 2);
  double Q2          = 2*Mnuc*E*x*y;
  double W2          = Mnuc2 + 2*Mnuc*E*y*(1-x);
  double W           = TMath::Max(0., TMath::Sqrt(W2));
  double inverse_eta = 0.5/x + TMath::Sqrt( 0.25/x2 + Mnuc2/Q2 );
  double eta         = 1 / inverse_eta;
  double xi          = eta * (1 + fMc2/Q2);
  double coshpsi     = (2-y)/y; // hyperbolic-cosine(psi)
  double sinh2psi    = TMath::Power(coshpsi, 2) - 1;

  //----- make sure that the mass-corrected x is in physical region
  if(xi<=0 || xi>1) return 0;

  //----- Calculate the PDFs
  PDF pdfs;
  pdfs.SetModel(fPDFModel);  // <-- attach algorithm
  pdfs.Calculate(xi, Q2);    // <-- calculate

  //----- proton pdfs
  double us = pdfs.UpSea()/xi;
  double uv = pdfs.UpValence()/xi;
  double ds = pdfs.DownSea()/xi;
  double dv = pdfs.DownValence()/xi;
  double s  = pdfs.Strange()/xi;

  //----- if the hit nucleon is a neutron, swap u<->d
  double tmp;
  if (isN) { 
    tmp = uv; uv = dv; dv = tmp;
    tmp = us; us = ds; ds = tmp;
  }

  //----- if a hit quark has been set then switch off other contributions
  if(qset) {
    dv = ( isd && !sea) ? dv : 0.;
    ds = ( isd &&  sea) ? ds : 0.;
    s  = ( iss &&  sea) ? s  : 0.;
  }

  //----- calculate the cross section
  double Gw2    = TMath::Power((kGF/kSqrt2)*(1+Q2/kMw2), 2);
  double f1     = TMath::Power( (1+coshpsi)/2, 2);
  double f2     = 0.25 * (fMc2/Q2) * sinh2psi;
  double xsec_0 = 2 * Gw2 * (y*Q2/kPi) * (f1+f2);
  double xsec_d = xsec_0 * fVcd2 * (dv+ds);
  double xsec_s = xsec_0 * fVcs2 * s;
  double xsec   = xsec_d + xsec_s;

  LOG("DISCharmXSec", pDEBUG)
    << "\n dxsec[DISCharm,FreeN]/dxdy (E= " << E
                 << ", x= " << x << ", y= " << y
                         << ", W= " << W << ", Q2 = " << Q2 << ") = " << xsec;

  //----- The algorithm computes d^2xsec/dxdy
  //      Check whether variable tranformation is needed
  if(kps!=kPSxyfE) {
    double J = utils::kinematics::Jacobian(interaction,kPSxyfE,kps);
    xsec *= J;
  }

  //----- If requested return the free nucleon xsec even for input nuclear tgt
  if( interaction->TestBit(kIAssumeFreeNucleon) ) return xsec;

  //----- Nuclear cross section (simple scaling here)
  int NNucl = (isP) ? target.Z() : target.N();
  xsec *= NNucl;

  return xsec;
}
//____________________________________________________________________________
double AivazisCharmPXSecLO::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool AivazisCharmPXSecLO::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const XclsTag &      xcls       = interaction->ExclTag();
  const InitialState & init_state = interaction->InitState();
  const ProcessInfo &  proc_info  = interaction->ProcInfo();

  if(!proc_info.IsDeepInelastic()) return false;
  if(!proc_info.IsWeak())          return false;

  bool is_inclusive_charm = (xcls.IsCharmEvent() && xcls.IsInclusiveCharm());
  if(!is_inclusive_charm) return false;

  int  nu  = init_state.ProbePdg();
  int  nuc = init_state.Tgt().HitNucPdg();

  if (!pdg::IsProton(nuc)  && !pdg::IsNeutron(nuc))     return false;
  if (!pdg::IsNeutrino(nu) && !pdg::IsAntiNeutrino(nu)) return false;

  return true;
}
//____________________________________________________________________________
void AivazisCharmPXSecLO::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void AivazisCharmPXSecLO::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfig();
}
//____________________________________________________________________________
void AivazisCharmPXSecLO::LoadConfig(void)
{
  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();

  // read mc, Vcd, Vcs from config or set defaults
  fMc    = fConfig->GetDoubleDef("Charm-Mass", gc->GetDouble("Charm-Mass"));
  fVcd   = fConfig->GetDoubleDef("CKM-Vcd",    gc->GetDouble("CKM-Vcd"));
  fVcs   = fConfig->GetDoubleDef("CKM-Vcs",    gc->GetDouble("CKM-Vcs"));

  fMc2   = TMath::Power(fMc,  2);
  fVcd2  = TMath::Power(fVcd, 2);
  fVcs2  = TMath::Power(fVcs, 2);

  // load PDF set
  fPDFModel = dynamic_cast<const PDFModelI *> (this->SubAlg("PDF-Set"));
  assert(fPDFModel);

  // load XSec Integrator
  fXSecIntegrator =
      dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);
}
//____________________________________________________________________________
