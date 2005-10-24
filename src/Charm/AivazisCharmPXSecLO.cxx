//____________________________________________________________________________
/*!

\class    genie::AivazisCharmPXSecLO

\brief    Computes, at Leading Order (LO), the differential cross section for
          neutrino charm production using the \b Aivazis,Olness,Tung model.

          The computed cross section is the D2xsec = d^2(xsec) / dy dx \n
          where \n
            \li \c y is the inelasticity, and
            \li \c x is the Bjorken scaling variable \c x

          Is a concrete implementation of the XSecAlgorithmI interface.

\ref      M.A.G.Aivazis, F.I.Olness and W.K.Tung
          Phys.Rev.D50, 3085-3101 (1994)

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  June 10, 2004

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Charm/AivazisCharmPXSecLO.h"
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

  //----- get scattering & init-state parameters

  const ScatteringParams & sc_params  = interaction -> GetScatteringParams();
  const InitialState &     init_state = interaction -> GetInitialState();

  TLorentzVector * p4 = init_state.GetProbeP4(kRfStruckNucAtRest);

  double Mnuc = kNucleonMass; // or init_state.TargetMass(); ?
  double E    = p4->Energy();
  double x    = sc_params.x();
  double y    = sc_params.y();

  delete p4;

  //----- make sure that x, y are in the physically acceptable region

  if(x<=0 || x>1) return 0;
  if(y<=0 || y>1) return 0;


  //----- get charm mass & CKM elements -- get constants unless they
  //      are overriden in the config registry

  double mc    = this->CharmMass();
  double Vcd   = this->Vcd();
  double Vcs   = this->Vcs();
  double mc2   = TMath::Power(mc,  2);
  double Vcd2  = TMath::Power(Vcd, 2);
  double Vcs2  = TMath::Power(Vcs, 2);

  //----- compute kinematic & auxiliary parameters

  double x2          = TMath::Power(x,    2);
  double Mnuc2       = TMath::Power(Mnuc, 2);
  double Q2          = 2*Mnuc*E*x*y;
  double W2          = Mnuc2 + 2*Mnuc*E*y*(1-x);
  double W           = TMath::Max(0., TMath::Sqrt(W2));
  double inverse_eta = 0.5/x + TMath::Sqrt( 0.25/x2 + Mnuc2/Q2 );
  double eta         = 1 / inverse_eta;
  double xi          = eta * (1 + mc2/Q2);
  double coshpsi     = (2-y)/y; // hyperbolic-cosine(psi)
  double sinh2psi    = TMath::Power(coshpsi, 2) - 1;

  //----- make sure that the mass-corrected x is in physical region

  if(xi<=0 || xi>1) return 0;

  //----- Get the physical W and Q2 range and check whether the current W,Q2
  //      pair is allowed

  Range1D_t rW  = utils::kinematics::WRange     (interaction);
  Range1D_t rQ2 = utils::kinematics::Q2Range_xy (interaction);

  bool in_range = utils::math::IsWithinLimits(Q2, rQ2)
                                        && utils::math::IsWithinLimits(W, rW);

  if(!in_range) {
    LOG("AivazisCharm", pDEBUG)
        << "\n W: " << "[" << rW.min << ", " << rW.max << "] GeV"
                 << " Q2: "<< "[" << rQ2.min << ", " << rQ2.max << "] GeV^2";
    LOG("AivazisCharm", pDEBUG)
        << "\n (W = " << W << ", Q2 = " << Q2 << " is not in physical range"
        << " - returning 0";
    return 0;
  }

  //----- Calculate the PDFs

  const Algorithm * algbase = this->SubAlg("pdf-alg-name", "pdf-param-set");

  const PDFModelI* pdf_model = dynamic_cast<const PDFModelI *>(algbase);

  PDF pdfs;

  pdfs.SetModel(pdf_model);   // <-- attach algorithm
  pdfs.Calculate(xi, Q2);    // <-- calculate

  bool isP = pdg::IsProton ( init_state.GetTarget().StruckNucleonPDGCode() );
  bool isN = pdg::IsNeutron( init_state.GetTarget().StruckNucleonPDGCode() );

  double d = 0;

  if(!isP && !isN) return 0;
  else if(isP) d = pdfs.DownValence() + pdfs.DownSea();
  else         d = pdfs.UpValence()   + pdfs.UpSea();

  double s = pdfs.Strange();

  d /= xi;
  s /= xi;


  //----- Check if we compute contributions from both d and s quarks

  bool d_contributes = true;
  bool s_contributes = true;

  if( fConfig->Exists("d-contrib-switch") )
                      fConfig->Get("d-contrib-switch", d_contributes);
  if( fConfig->Exists("s-contrib-switch") )
                      fConfig->Get("s-contrib-switch", s_contributes);

  //----- Calculate cross section

  double Gw  = (kGF/kSqrt2) * (1 + Q2/kMw_2);
  double Gw2 = TMath::Power(Gw, 2);
  double tmp = Gw2 * (y*Q2/kPi) *
                   (TMath::Power((1+coshpsi)/2, 2) + (0.5*mc2/Q2)*sinh2psi/2);
  double xsec = 0;

  if(d_contributes) {
      double xsec_d = 2 * Vcd2 * d * tmp;
      xsec += xsec_d;
  }
  if(s_contributes) {
      double xsec_s = 2 * Vcs2 * s * tmp;
      xsec += xsec_s;
  }

  LOG("AivazisCharm", pDEBUG)
    << "\n dxsec[DIS-Charm]/dxdy (E= " << E
                 << ", x= " << x << ", y= " << y
                         << ", W= " << W << ", Q2 = " << Q2 << ") = " << xsec;
  return xsec;
}
//____________________________________________________________________________
double AivazisCharmPXSecLO::CharmMass(void) const
{
// Allows default charm mass to be overriden by XML config / config registry

  return (fConfig->Exists("c-quark-mass")) ?
                  fConfig->GetDouble("c-quark-mass") :
                  PDGLibrary::Instance()->Find(kPdgCQuark)->Mass();
}
//____________________________________________________________________________
double AivazisCharmPXSecLO::Vcd(void) const
{
// Allows default CKM Vcd  to be overriden by XML config / config registry

  return (fConfig->Exists("Vcd")) ? fConfig->GetDouble("Vcd") : kVcd;
}
//____________________________________________________________________________
double AivazisCharmPXSecLO::Vcs(void) const
{
// Allows default CKM Vcs to be overriden by XML config / config registry

  return (fConfig->Exists("Vcs")) ? fConfig->GetDouble("Vcs") : kVcs;
}
//____________________________________________________________________________
