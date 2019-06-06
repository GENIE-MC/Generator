//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Physics/Charm/XSection/SlowRsclCharmDISPXSecLO.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/PartonDistributions/PDF.h"
#include "Physics/PartonDistributions/PDFModelI.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Utils/Range1.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
SlowRsclCharmDISPXSecLO::SlowRsclCharmDISPXSecLO() :
XSecAlgorithmI("genie::SlowRsclCharmDISPXSecLO")
{

}
//____________________________________________________________________________
SlowRsclCharmDISPXSecLO::SlowRsclCharmDISPXSecLO(string config) :
XSecAlgorithmI("genie::SlowRsclCharmDISPXSecLO", config)
{

}
//____________________________________________________________________________
SlowRsclCharmDISPXSecLO::~SlowRsclCharmDISPXSecLO()
{

}
//____________________________________________________________________________
double SlowRsclCharmDISPXSecLO::XSec(
                  const Interaction * interaction, KinePhaseSpace_t kps) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  if(interaction->ProcInfo().IsWeakNC()) return 0;

  //----- get kinematics & init-state parameters
  const Kinematics &   kinematics = interaction->Kine();
  const InitialState & init_state = interaction->InitState();
  const Target &       target     = init_state.Tgt();

  double Mnuc = target.HitNucMass();
  double E    = init_state.ProbeE(kRfHitNucRest);
  double x    = kinematics.x();
  double y    = kinematics.y();
  double Q2   = 2*Mnuc*E*x*y;

  //----- get target information (hit nucleon and quark)
  int  nu    = init_state.ProbePdg();
  int  nuc   = target.HitNucPdg();
  bool isP   = pdg::IsProton (nuc);
  bool isN   = pdg::IsNeutron(nuc);
  bool qset  = target.HitQrkIsSet();
  int  qpdg  = (qset) ? target.HitQrkPdg()       : 0;
  bool sea   = (qset) ? target.HitSeaQrk()       : false;
  bool isd   = (qset) ? pdg::IsDQuark (qpdg)     : false;
  bool iss   = (qset) ? pdg::IsSQuark (qpdg)     : false;
  bool isdb  = (qset) ? pdg::IsAntiDQuark (qpdg) : false;
  bool issb  = (qset) ? pdg::IsAntiSQuark (qpdg) : false;
  bool isnu  = pdg::IsNeutrino(nu);
  bool isnub = pdg::IsAntiNeutrino(nu);

  //----- compute slow rescaling variable & check its value
  double xi = x * (1 + fMc2/Q2);
  if(xi<=0 || xi>1) return 0.;

  //----- calculate the PDFs
  PDF pdfs;
  pdfs.SetModel(fPDFModel);   // <-- attach algorithm
  pdfs.Calculate(xi, Q2);     // <-- calculate

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
    if(isnub) { bool pass = (isdb||issb)&&sea; if(!pass) return 0; }
    if(isnu)  { bool pass = isd||(iss&&sea);   if(!pass) return 0; }
    dv = (  isd        && !sea) ? dv : 0.;
    ds = ( (isd||isdb) &&  sea) ? ds : 0.;
    s  = ( (iss||issb) &&  sea) ? s  : 0.;
  }
  //----- in case of a \bar{v}+N calculation (quark not set) zero the d/val contribution
  if(isnub) dv=0;

  //----- calculate cross section
  double Gw2    = TMath::Power( (kGF/kSqrt2)*(1+Q2/kMw2), 2);
  double xsec_0 = Gw2 * 2*Q2/(y*kPi) * (y + xi*(1-y)/x);
  double xsec_d = xsec_0 * fVcd2 * (dv+ds);
  double xsec_s = xsec_0 * fVcs2 * s;
  double xsec   = xsec_d + xsec_s;

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  double Mnuc2 = TMath::Power(Mnuc, 2);
  double W2    = Mnuc2 + 2*Mnuc*E*y*(1-x);
  double W     = TMath::Max(0., TMath::Sqrt(W2));
  LOG("DISCharmXSec", pDEBUG)
    << "\n dxsec[DISCharm,FreeN]/dxdy (E= " << E
                 << ", x= " << x << ", y= " << y
                         << ", W= " << W << ", Q2 = " << Q2 << ") = " << xsec;
#endif

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
double SlowRsclCharmDISPXSecLO::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool SlowRsclCharmDISPXSecLO::ValidProcess(
                                        const Interaction * interaction) const
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
void SlowRsclCharmDISPXSecLO::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void SlowRsclCharmDISPXSecLO::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfig();
}
//____________________________________________________________________________
void SlowRsclCharmDISPXSecLO::LoadConfig(void)
{
  // read mc, Vcd, Vcs from config or set defaults
  GetParam( "Charm-Mass", fMc ) ;
  GetParam( "CKM-Vcd",    fVcd ) ;
  GetParam( "CKM-Vcs",    fVcs ) ;

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
