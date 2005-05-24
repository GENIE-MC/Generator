//____________________________________________________________________________
/*!

\class    genie::BardinDISPXSec

\brief    Computes, the differential cross section for neutrino DIS including
          radiative corrections according to the \b Bardin-Dokuchaeva model.
          
          The computed xsec is the double differential d^2(xsec) / dy dx \n
          where \n
           \li \c y is the inelasticity, and
           \li \c x is the Bjorken scaling variable \c x

          Is a concrete implementation of the XSecAlgorithmI interface.

\ref      D.Yu.Bardin, V.A.Dokuchaeva, "On the radiative corrections to the
          Neutrino Deep Inelastic Scattering", JINR-E2-86-260, Apr. 1986

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  July 06, 2004

*/
//____________________________________________________________________________

#include <iostream>

#include <TMath.h>

#include "Algorithm/AlgFactory.h"
#include "Conventions/Constants.h"
#include "Conventions/RefFrame.h"
#include "Conventions/Units.h"
#include "Conventions/Utils.h"
#include "Messenger/Messenger.h"
#include "Numerical/UnifGrid.h"
#include "Numerical/IntegratorI.h"
#include "PDF/PDFModelI.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "Radiative/BardinDISPXSec.h"
#include "Utils/MathUtils.h"
#include "Utils/KineLimits.h"
#include "Utils/Range1.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
BardinDISPXSec::BardinDISPXSec() :
XSecAlgorithmI()
{
  fName     = "genie::BardinDISPXSec";
}
//____________________________________________________________________________
BardinDISPXSec::BardinDISPXSec(const char * param_set) :
XSecAlgorithmI(param_set)
{
  fName = "genie::BardinDISPXSec";

  FindConfig();
}
//____________________________________________________________________________
BardinDISPXSec::~BardinDISPXSec()
{

}
//____________________________________________________________________________
double BardinDISPXSec::XSec(const Interaction * interaction) const
{
  LOG("Bardin", pDEBUG) << *fConfig;

  //----- get scattering & init-state parameters

  const ScatteringParams & sc_params = interaction -> GetScatteringParams();
  const InitialState &    init_state = interaction -> GetInitialState();

  TLorentzVector * p4 = init_state.GetProbeP4(kRfStruckNucAtRest);
  
  double Mnuc = kNucleonMass; // or init_state.TargetMass(); ?
  double E    = p4->Energy();
  double x    = sc_params.x();
  double y    = sc_params.y();

  delete p4;
  
  //----- make sure that x, y are in the physically acceptable region

  if(x<=0 || x>1) return 0;
  if(y<=0 || y>1) return 0;

  //----- compute auxiliary & kinematic parameters

  double Mnuc2 = pow(Mnuc, 2);
  double W2    = Mnuc2 + 2*Mnuc*E*y*(1-x);
  double Q2    = S(interaction) * x * y;
  double W     = TMath::Max(0., TMath::Sqrt(W2));

  //----- Get the physical W and Q2 range and check whether the current W,Q2
  //      pair is allowed

  Range1D_t rW  = kine_limits::WRange     (interaction);
  Range1D_t rQ2 = kine_limits::Q2Range_xy (interaction);

  bool in_range = math_utils::IsWithinLimits(Q2, rQ2)
                                        && math_utils::IsWithinLimits(W, rW);

  if(!in_range) {
       LOG("DISXSec", pDEBUG)
             << "\n *** point (W = " << W
                           << ", Q2 = " << Q2 << " is not in physical range";
       LOG("DISXSec", pDEBUG)
             << "\n Physical W range: "
                               << "[" << rW.min << ", " << rW.max << "] GeV";
       LOG("DISXSec", pDEBUG)
             << "\n Physical Q2 range: "
                           << "[" << rQ2.min << ", " << rQ2.max << "] GeV^2";
       return 0;
   }
  
  //----- make sure it can get all critical configuration parameters

  assert(
          fConfig->Exists("pdf-alg-name")           &&
          fConfig->Exists("pdf-param-set")          &&
          fConfig->Exists("init-quark-pdg-code")    &&
          fConfig->Exists("final-quark-pdg-code")   &&
          fConfig->Exists("final-quark-mass")  
        );

  //----- get init & final quarks

  int init_pdgc = fConfig->GetInt("init-quark-pdg-code");
  int fin_pdgc  = fConfig->GetInt("final-quark-pdg-code");

  //----- check if the user wants to override standard the CKM element
  //      before reading it from the constants

  double Vckm = 0;

  if ( ! fConfig->Exists("Vckm") )
        Vckm = Utils::CkmElement(init_pdgc, fin_pdgc);
  else  Vckm = fConfig->GetDouble("Vckm");

  LOG("Bardin", pDEBUG) << "Vckm = " << Vckm;

  double Vckm2 = Vckm*Vckm;
  
  //----- get the PDF set

  string pdf_alg_name, pdf_param_set;

  fConfig->Get("pdf-alg-name",  pdf_alg_name );
  fConfig->Get("pdf-param-set", pdf_param_set);
  
  //----- Get PDF objects & attach PDFModels

  AlgFactory * algf = AlgFactory::Instance();

  const Algorithm* algbase = algf->GetAlgorithm(pdf_alg_name, pdf_param_set);

  assert(algbase);

  const PDFModelI* pdf_model = dynamic_cast<const PDFModelI *> (algbase);

  PDF pdf_x, pdf_xi;

  pdf_x.SetModel  (pdf_model);   // <-- attach algorithm
  pdf_xi.SetModel (pdf_model);   // <-- attach algorithm

  //----- Get init quark PDF at (x,Q2)
  
  pdf_x.Calculate(x, Q2);

  LOG("Bardin", pDEBUG) << pdf_x;
  
  double f_x  = PDFFunc( pdf_x,  init_pdgc )  / x;
  
  //-- Define the range for variable xi
  //
  //   Note: Integration over xi is performed for xi e [x, 1].
  //   During the numerical integration, When we compute values of the
  //   cross section terms at various xi, then :
  //   - xi = x can not be included because of a singularity { 1/(xi-x) term }
  //   - xi = 1 can not be included because PDF calculators complain

  const int    nxi   = 201;
  const double xi_min = x + 0.001;
  const double xi_max = 0.999;
  const double dxi   = (xi_max-xi_min)/(nxi-1);

  UnifGrid grid;

  grid.AddDimension(nxi, xi_min, xi_max);

  FunctionMap term2_vs_xi(grid);
  FunctionMap term3_vs_xi(grid);

  //-- internal loop : integration over xi for the given (x,y,E)
  //   xi e (x,1)

  double term1 =  x * f_x * (1 + (kAem/kPi)*DeltaCCi(interaction));

  for(int ix = 0; ix < nxi; ix++) {

     double xi = xi_min + ix * dxi;

     //----- Get init quark PDF at (xi,Q2)
     
     pdf_xi.Calculate(xi, Q2);

     double f_xi = PDFFunc( pdf_xi, init_pdgc )  / xi;

     double term2_xi = f_xi * PhiCCi(xi, interaction);

     double term3_xi = (xi*f_xi*Ii(xi,interaction)-x*f_x*Ii(x, interaction))/(xi-x);

     term2_vs_xi.AddPoint( term2_xi, ix);
     term3_vs_xi.AddPoint( term3_xi, ix);
  }

  //----- Numerical integration

  //-- get specified integration algorithm from the config. registry
  //   or use Simpson1D if no one else is defined

  string integrator_name;

  if( fConfig->Exists("integrator") )
                          fConfig->Get("integrator", integrator_name );
  else integrator_name = "genie::Simpson1D";

  //-- ask the AlgFactory for the integrator algorithm

  const Algorithm * alg_base_intg = algf->GetAlgorithm(integrator_name);

  const IntegratorI * integrator =
                          dynamic_cast<const IntegratorI *> (alg_base_intg);

  //-- integrate

  double term2 = integrator->Integrate( term2_vs_xi );
  double term3 = integrator->Integrate( term3_vs_xi );

  //-- form the differential cross section
  
  double Gfactor = pow(kGF,2) * S(interaction) * Vckm2 / kPi;

  double d2xsec_dxdy = Gfactor * ( term1 + (kAem/kPi) * (term2 + term3) );

  LOG("Bardin", pDEBUG) << "term1 = " << term1;
  LOG("Bardin", pDEBUG) << "term2 = " << term2;
  LOG("Bardin", pDEBUG) << "term3 = " << term3;

  LOG("Bardin", pINFO)
             << "d2xsec/dxdy (E = " << E << ", x = " << x
                                   << ", y = " << y << ") = " << d2xsec_dxdy;

  assert(d2xsec_dxdy >= 0);
  
  return d2xsec_dxdy;
}
//____________________________________________________________________________
double BardinDISPXSec::PhiCCi(double xi, const Interaction * interaction) const
{
  const ScatteringParams & scp = interaction->GetScatteringParams();

  double x    = scp.GetDouble("x");
  double y    = scp.GetDouble("y");
  int    pdg  = fConfig->GetInt("init-quark-pdg-code");

  double mqf  = fConfig->GetDouble("final-quark-mass");
  double mqi  = InitQuarkMass(xi);
  double mqf2 = pow(mqf,2);
  double mqi2 = pow(mqi,2);

  double ml   = interaction->GetFSPrimaryLepton()->Mass();

  double Q2   = S(interaction) * x * y;
  double t    = tau(xi, interaction);
  double st   = St(xi, interaction);
  double u    = U(xi, interaction);
  double su   = Su(xi, interaction);

  double f    = Utils::QuarkCharge(pdg);

  double ml2  = ml * ml;
  double st2  = st * st;
  double u2   = u  * u;
  double f2   = f  * f;
  double xi2  = xi * xi;
  double x2   = x  * x;
  double y2   = y  * y;

  double term1 = pow(1+f, 2) * (0.25*Q2 /t) * (1 - mqf2/t);

  double term2 = 0.25 + 0.75*y - 0.5*y*(1 + xi/(xi-y*(xi-x)) ) * log(st2/(t*ml2));

  double term3 = f * ( 2 - (y*xi/(xi-y*(xi-x))) * log(st2/(t*ml2)) );

  double term4 = f * ( y*(y-2-y*x/xi) * log(st*u/(t*su)) + (0.5+y)*x/xi );

  double term5 = f2 * ( 2 - (1 + 0.5*x/xi + 0.5*x2/xi2) * log(u2/(t*mqi2)) );

  double term6 = f2 * ( 2 - 0.5*y + 0.25*y2 ) * x/xi;

  double term7 = f2 * (1.5 - 0.5*y - 0.25*y2) * x2/xi2;

  double phi   = term1 + term2 + term3 + term4 + term5 + term6 + term7;

  LOG("Bardin", pDEBUG) << "PhiCC = " << phi;

  return phi;
}
//__________________________________________________________________________
double BardinDISPXSec::DeltaCCi(const Interaction * interaction) const
{
  const ScatteringParams & scp = interaction->GetScatteringParams();

  double x    = scp.GetDouble("x");
  double y    = scp.GetDouble("y");
  int    pdg  = fConfig->GetInt("init-quark-pdg-code");

  double mqf  = fConfig->GetDouble("final-quark-mass");
  double mqi  = InitQuarkMass(interaction);
  double mqf2 = pow(mqf,2);
  double mqi2 = pow(mqi,2);

  double ml   = interaction->GetFSPrimaryLepton()->Mass();
  double ml2  = ml*ml;

  double s       = S(interaction);
  double sq      = Sq(interaction);
  double Iixx    = Ii(x, interaction);
  double sq_mqf2 = sq / mqf2;
  double sq_ml2  = sq / ml2;
  double sq_mZ2  = sq / kMz_2;
  double Q2      = S(interaction) * x * y;

  double f       = Utils::QuarkCharge(pdg);
  double f2      = f*f;

  LOG("Bardin", pDEBUG) << "ml      = " << ml;
  LOG("Bardin", pDEBUG) << "s       = " << s;
  LOG("Bardin", pDEBUG) << "sq      = " << sq;
  LOG("Bardin", pDEBUG) << "Iixx    = " << Iixx;
  LOG("Bardin", pDEBUG) << "sq_Mgf2 = " << sq_mqf2;
  LOG("Bardin", pDEBUG) << "sq_Ml2  = " << sq_ml2;
  LOG("Bardin", pDEBUG) << "sq_MZ2  = " << sq_mZ2;
  LOG("Bardin", pDEBUG) << "Q2      = " << Q2;
  LOG("Bardin", pDEBUG) << "f       = " << f;

  double delta =  Iixx * log(s * y * (1-x) / mqf2) + 0.75 +
                   (1.75 - log(sq_ml2)) * log(sq_mqf2) -
                   0.5 * pow( log(sq_mqf2), 2 ) +
                   pow(kPi, 2) / 2. -
                   1.5 * log(sq_mZ2) +
                   0.75 * log(sq_ml2) +
                   f * (1.75 + (1.5 + 2*log(1-y)) * log(Q2/mqf2) -
                   pow( log(Q2/mqf2), 2) +
                   pow(kPi, 2) / 2. -
                   1.5 * log(sq_mZ2 * (1-y)) -
                   log(y) * log(1-y) ) +
                   f2 * ( -1 + (1.75 - log(Q2/mqi2))*log(Q2/mqf2) -
                   0.5 * pow( log(Q2/mqf2), 2) + 0.75*log(Q2/mqi2));

  LOG("Bardin", pDEBUG) << "delta = " << delta;

  return delta;
}
//__________________________________________________________________________
double BardinDISPXSec::Ii(double xi, const Interaction * interaction) const
{
  const ScatteringParams & scp = interaction->GetScatteringParams();

  double x    = scp.GetDouble("x");
  double y    = scp.GetDouble("y");
  int    pdg  = fConfig->GetInt("init-quark-pdg-code");

  double ml   = interaction->GetFSPrimaryLepton()->Mass();
  double mqi  = InitQuarkMass(xi);

  double st  = St(xi, interaction);
  double u   = U(xi, interaction);
  double su  = Su(xi, interaction);
  double t   = tau(xi, interaction);

  double f    = Utils::QuarkCharge(pdg);

  double ml2  = ml  * ml;
  double mqi2 = mqi * mqi;
  double st2  = st  * st;
  double u2   = u   * u;
  double f2   = f   * f;

  double term1 = ( xi / (xi - y*(xi-x)) ) * log( st2/(t*ml2) ) - 2;

  double term2 = f * (2 * log(st*u/(t*su)) - 2 +
                        ( y*(xi-x)/(xi-y*(xi-x)) ) * log( st2/(t*ml2) ) );

  double term3 = f2 * ( log( u2/(t*mqi2) ) - 2 );

  double I = term1 + term2 + term3;

  return I;
}
//__________________________________________________________________________
double BardinDISPXSec::S(const Interaction * interaction) const
{
  const InitialState & init_state = interaction -> GetInitialState();

  TLorentzVector * p4 = init_state.GetProbeP4(kRfStruckNucAtRest);

  double E = p4->Energy();

  delete p4;
  
  return 2 * kNucleonMass * E;
}
//__________________________________________________________________________
double BardinDISPXSec::U(double xi, const Interaction * interaction) const
{
  double y  = interaction->GetScatteringParams().GetDouble("y");

  return S(interaction) * y * xi;
}
//__________________________________________________________________________
double BardinDISPXSec::tau(double xi, const Interaction * interaction) const
{
  double x    = interaction->GetScatteringParams().GetDouble("x");
  double y    = interaction->GetScatteringParams().GetDouble("y");

  double mqf  = fConfig->GetDouble("final-quark-mass");
  double mqf2 = pow(mqf,2);

  return S(interaction) * y * (xi-x) + mqf2;
}
//__________________________________________________________________________
double BardinDISPXSec::St(double xi, const Interaction * interaction) const
{
  double x  = interaction->GetScatteringParams().GetDouble("x");
  double y  = interaction->GetScatteringParams().GetDouble("y");

  return S(interaction) * (xi - y*(xi-x));
}
//__________________________________________________________________________
double BardinDISPXSec::Su(double xi, const Interaction * interaction) const
{
  double y  = interaction->GetScatteringParams().GetDouble("y");

  return S(interaction) * xi * (1-y);
}
//__________________________________________________________________________
double BardinDISPXSec::Sq(const Interaction * interaction) const
{
  double x  = interaction->GetScatteringParams().GetDouble("x");

  return S(interaction) * x;
}
//__________________________________________________________________________
double BardinDISPXSec::InitQuarkMass(const Interaction * interaction) const
{
  double x  = interaction->GetScatteringParams().GetDouble("x");

  return x * kNucleonMass;
}
//__________________________________________________________________________
double BardinDISPXSec::InitQuarkMass(double xi) const
{
  return xi * kNucleonMass;
}
//__________________________________________________________________________
double BardinDISPXSec::PDFFunc(const PDF & pdf, int pdgc) const
{
  if     ( pdg::IsUQuark(pdgc) ) return (pdf.UpValence()   + pdf.UpSea()  );
  else if( pdg::IsDQuark(pdgc) ) return (pdf.DownValence() + pdf.DownSea());
  else {
     LOG("Bardin", pERROR)
         << "Scattering off quarks with pdg = " << pdgc << "is not handled";
  }
  return 0;  
}
//__________________________________________________________________________

