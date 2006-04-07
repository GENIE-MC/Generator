//____________________________________________________________________________
/*!

\class    genie::SlowRsclCharmDISPXSecLO

\brief    Computes, at Leading Order (LO), the differential cross section for
          neutrino charm production using a slow rescaling model.

          The computed cross section is the D2xsec = d^2(xsec) / dy dx \n
          where \n
            \li \c y is the inelasticity, and
            \li \c x is the Bjorken scaling variable \c x

          Is a concrete implementation of the XSecAlgorithmI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  June 10, 2004

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Charm/SlowRsclCharmDISPXSecLO.h"
#include "Conventions/Constants.h"
#include "Conventions/RefFrame.h"
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
double SlowRsclCharmDISPXSecLO::XSec(const Interaction * interaction) const
{
  LOG("SlowRsclCharm", pDEBUG) << *fConfig;

  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  //----- get kinematics & init-state parameters
  const Kinematics &   kinematics = interaction->GetKinematics();
  const InitialState & init_state = interaction->GetInitialState();

  double Mnuc = init_state.GetTarget().StruckNucleonP4()->M();
  double E    = init_state.GetProbeE(kRfStruckNucAtRest);
  double x    = kinematics.x();
  double y    = kinematics.y();

  //----- compute kinematic & auxiliary parameters
  double Mnuc2 = TMath::Power(Mnuc, 2);
  double Q2    = 2*Mnuc*E*x*y;
  double W2    = Mnuc2 + 2*Mnuc*E*y*(1-x);
  double W     = TMath::Max(0., TMath::Sqrt(W2));

  //----- compute slow rescaling variable & check its value
  double xi = x * (1 + fMc2/Q2);
  if(xi<=0 || xi>1) return 0.;

  //----- Calculate the PDFs
  PDF pdfs;
  pdfs.SetModel(fPDFModel);   // <-- attach algorithm
  pdfs.Calculate(xi, Q2);     // <-- calculate

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
  double Gw  = (kGF/kSqrt2) * (1 + Q2/kMw_2);
  double Gw2 = TMath::Power(Gw, 2);
  double tmp = Gw2 * 2*Q2/(y*kPi) * (y + xi*(1-y)/x);

  double xsec = 0;

  if(fDContributes) {
      double xsec_d = fVcd2 * d * tmp;
      xsec += xsec_d;
  }
  if(fSContributes) {
      double xsec_s = fVcs2 * s * tmp;
      xsec += xsec_s;
  }

  LOG("SlowRsclCharm", pDEBUG)
    << "\n dxsec[DIS-Charm]/dxdy (E= " << E
                 << ", x= " << x << ", y= " << y
                         << ", W= " << W << ", Q2 = " << Q2 << ") = " << xsec;
  return xsec;
}
//____________________________________________________________________________
bool SlowRsclCharmDISPXSecLO::ValidProcess(
                                        const Interaction * interaction) const
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
bool SlowRsclCharmDISPXSecLO::ValidKinematics(
                                        const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipKinematicChk)) return true;
  const Kinematics &   kinematics = interaction -> GetKinematics();
  const InitialState & init_state = interaction -> GetInitialState();

  double E    = init_state.GetProbeE(kRfStruckNucAtRest);
  double x    = kinematics.x();
  double y    = kinematics.y();

  if(x<=0 || x>=1) {
    LOG("SlowRsclCharm", pDEBUG) << "x is unphysical or at limit, x = " << x;
    return false;
  }
  if(y<=0 || y>=1) {
    LOG("SlowRsclCharm", pDEBUG) << "x is unphysical or at limit, x = " << y;
    return false;
  }

  double Mnuc  = init_state.GetTarget().StruckNucleonP4()->M();
  double Mnuc2 = TMath::Power(Mnuc, 2);
  double Q2    = 2*Mnuc*E*x*y;
  double W2    = Mnuc2 + 2*Mnuc*E*y*(1-x);
  double W     = TMath::Max(0., TMath::Sqrt(W2));

  //----- Get the physical W and Q2 range and check whether the current W,Q2
  //      pair is allowed
  Range1D_t rW  = utils::kinematics::WRange     (interaction);
  Range1D_t rQ2 = utils::kinematics::Q2Range_xy (interaction);

  bool in_range = utils::math::IsWithinLimits(Q2, rQ2)
                                        && utils::math::IsWithinLimits(W, rW);
  if(!in_range) {
    LOG("SlowRsclCharm", pDEBUG)
        << "\n W: " << "[" << rW.min << ", " << rW.max << "] GeV"
                 << " Q2: "<< "[" << rQ2.min << ", " << rQ2.max << "] GeV^2";
    LOG("SlowRsclCharm", pDEBUG)
        << "\n (W = " << W << ", Q2 = " << Q2 << " is not in physical range"
        << " - returning 0";
    return false;
  }
  return true;
}
//____________________________________________________________________________
void SlowRsclCharmDISPXSecLO::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
  this->LoadSubAlg();
}
//____________________________________________________________________________
void SlowRsclCharmDISPXSecLO::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfigData();
  this->LoadSubAlg();
}
//____________________________________________________________________________
void SlowRsclCharmDISPXSecLO::LoadConfigData(void)
{
  // default charm mass
  double mcdef = PDGLibrary::Instance()->Find(kPdgCQuark)->Mass();

  // read mc, Vcd, Vcs from config or set defaults
  fMc    = fConfig->GetDoubleDef("c-quark-mass", mcdef);
  fVcd   = fConfig->GetDoubleDef("Vcd", kVcd);
  fVcs   = fConfig->GetDoubleDef("Vcs", kVcs);

  fMc2   = TMath::Power(fMc,  2);
  fVcd2  = TMath::Power(fVcd, 2);
  fVcs2  = TMath::Power(fVcs, 2);

  // check if we compute contributions from both d and s quarks
  fDContributes = fConfig->GetBoolDef("d-contrib-switch", true);
  fSContributes = fConfig->GetBoolDef("s-contrib-switch", true);
}
//____________________________________________________________________________
void SlowRsclCharmDISPXSecLO::LoadSubAlg(void)
{
  // load PDF set
  fPDFModel = 0;
  fPDFModel = dynamic_cast<const PDFModelI *> (
                           this->SubAlg("pdf-alg-name", "pdf-param-set"));
  assert(fPDFModel);
}
//____________________________________________________________________________
