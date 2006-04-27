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
double AivazisCharmPXSecLO::XSec(const Interaction * interaction) const
{
  LOG("AivazisCharm", pDEBUG) << *fConfig;

  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  //----- get init-state & kinematical parameters
  const Kinematics &   kinematics = interaction -> GetKinematics();
  const InitialState & init_state = interaction -> GetInitialState();

  //----- compute kinematic & auxiliary parameters
  double E           = init_state.GetProbeE(kRfStruckNucAtRest);
  double x           = kinematics.x();
  double y           = kinematics.y();
  double x2          = TMath::Power(x,    2);
  double Mnuc        = init_state.GetTarget().StruckNucleonMass();
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

  int  nuc = init_state.GetTarget().StruckNucleonPDGCode();
  bool isP = pdg::IsProton (nuc);
  bool isN = pdg::IsNeutron(nuc);

  double d = 0;

  if(!isP && !isN) return 0;
  else if(isP) d = pdfs.DownValence() + pdfs.DownSea();
  else         d = pdfs.UpValence()   + pdfs.UpSea();

  double s = pdfs.Strange();

  d /= xi;
  s /= xi;

  //----- Calculate cross section
  double Gw  = (kGF/kSqrt2) * (1 + Q2/kMw2);
  double Gw2 = TMath::Power(Gw, 2);
  double tmp = Gw2 * (y*Q2/kPi) *
                   (TMath::Power((1+coshpsi)/2, 2) + (0.5*fMc2/Q2)*sinh2psi/2);
  double xsec = 0;

  if(fDContributes) {
      double xsec_d = 2 * fVcd2 * d * tmp;
      xsec += xsec_d;
  }
  if(fSContributes) {
      double xsec_s = 2 * fVcs2 * s * tmp;
      xsec += xsec_s;
  }

  LOG("AivazisCharm", pDEBUG)
    << "\n dxsec[DIS-Charm]/dxdy (E= " << E
                 << ", x= " << x << ", y= " << y
                         << ", W= " << W << ", Q2 = " << Q2 << ") = " << xsec;
  return xsec;
}
//____________________________________________________________________________
bool AivazisCharmPXSecLO::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const XclsTag &      xcls       = interaction->GetExclusiveTag();
  const InitialState & init_state = interaction->GetInitialState();
  const ProcessInfo &  proc_info  = interaction->GetProcessInfo();

  if(!proc_info.IsDeepInelastic()) return false;
  if(!proc_info.IsWeak())          return false;

  bool is_inclusive_charm = (xcls.IsCharmEvent() && xcls.IsInclusiveCharm());
  if(!is_inclusive_charm) return false;

  int  nu  = init_state.GetProbePDGCode();
  int  nuc = init_state.GetTarget().StruckNucleonPDGCode();

  if (!pdg::IsProton(nuc)  && !pdg::IsNeutron(nuc))     return false;
  if (!pdg::IsNeutrino(nu) && !pdg::IsAntiNeutrino(nu)) return false;

  return true;
}
//____________________________________________________________________________
bool AivazisCharmPXSecLO::ValidKinematics(
                                        const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipKinematicChk)) return true;

  const Kinematics &   kinematics = interaction -> GetKinematics();
  const InitialState & init_state = interaction -> GetInitialState();

  double E    = init_state.GetProbeE(kRfStruckNucAtRest);
  double x    = kinematics.x();
  double y    = kinematics.y();

  if(x<=0 || x>=1) {
    LOG("AivazisCharm", pDEBUG) << "x is unphysical or at limit, x = " << x;
    return false;
  }
  if(y<=0 || y>=1) {
    LOG("AivazisCharm", pDEBUG) << "x is unphysical or at limit, x = " << y;
    return false;
  }

  double Mnuc  = init_state.GetTarget().StruckNucleonMass();
  double Mnuc2 = TMath::Power(Mnuc, 2);
  double Q2    = 2*Mnuc*E*x*y;
  double W2    = Mnuc2 + 2*Mnuc*E*y*(1-x);
  double W     = TMath::Max(0., TMath::Sqrt(W2));

  //----- Get the physical W and Q2 range and check whether the current W,Q2
  //      pair is allowed
  Range1D_t rW  = utils::kinematics::KineRange(interaction, kKVW);
  Range1D_t rQ2 = utils::kinematics::KineRange(interaction, kKVQ2);

  bool in_range = utils::math::IsWithinLimits(Q2, rQ2)
                                        && utils::math::IsWithinLimits(W, rW);
  if(!in_range) {
    LOG("AivazisCharm", pDEBUG)
        << "\n W: " << "[" << rW.min << ", " << rW.max << "] GeV"
                 << " Q2: "<< "[" << rQ2.min << ", " << rQ2.max << "] GeV^2";
    LOG("AivazisCharm", pDEBUG)
        << "\n (W = " << W << ", Q2 = " << Q2 << " is not in physical range"
        << " - returning 0";
    return false;
  }
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
  fMc    = fConfig->GetDoubleDef("c-quark-mass", gc->GetDouble("Charm-Mass"));
  fVcd   = fConfig->GetDoubleDef("Vcd", gc->GetDouble("CKM-Vcd"));
  fVcs   = fConfig->GetDoubleDef("Vcs", gc->GetDouble("CKM-Vcs"));

  fMc2   = TMath::Power(fMc,  2);
  fVcd2  = TMath::Power(fVcd, 2);
  fVcs2  = TMath::Power(fVcs, 2);

  // check if we compute contributions from both d and s quarks
  fDContributes = fConfig->GetBoolDef("d-contrib-switch", true);
  fSContributes = fConfig->GetBoolDef("s-contrib-switch", true);

  // load PDF set
  fPDFModel = 0;
  fPDFModel = dynamic_cast<const PDFModelI *> (
                           this->SubAlg("pdf-alg-name", "pdf-param-set"));
  assert(fPDFModel);
}
//____________________________________________________________________________
