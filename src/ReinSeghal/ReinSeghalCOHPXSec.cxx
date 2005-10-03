//____________________________________________________________________________
/*!

\class    genie::ReinSeghalCOHPXSec

\brief    Computes the double differential cross section for coherent pi0
          production according to the \b Rein-Seghal model.

          The computed cross section is the d^3 xsec/ dx dy dt

          where \n
            \li \c x : Bjorken x = Q2/2Mv
            \li \c y : Inelasticity y=v/E, v=E-E'

          The t dependence is analytically integrated out.

          Is a concrete implementation of the XSecAlgorithmI interface.

\ref      D.Rein and L.M.Seghal, Coherent pi0 production in neutrino
          reactions, Nucl.Phys.B223:29-144 (1983)

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  March 11, 2005

____________________________________________________________________________*/

#include <TMath.h>

#include "Algorithm/AlgFactory.h"
#include "Conventions/Constants.h"
#include "Conventions/RefFrame.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "ReinSeghal/ReinSeghalCOHPXSec.h"
#include "Utils/XSecUtils.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
ReinSeghalCOHPXSec::ReinSeghalCOHPXSec() :
XSecAlgorithmI()
{
  fName = "genie::ReinSeghalCOHPXSec";
}
//____________________________________________________________________________
ReinSeghalCOHPXSec::ReinSeghalCOHPXSec(const char * param_set) :
XSecAlgorithmI(param_set)
{
  fName = "genie::ReinSeghalCOHPXSec";

  FindConfig();
}
//____________________________________________________________________________
ReinSeghalCOHPXSec::~ReinSeghalCOHPXSec()
{

}
//____________________________________________________________________________
double ReinSeghalCOHPXSec::XSec(const Interaction * interaction) const
{
  LOG("ReinSeghalCoh", pDEBUG) << *fConfig;

  //----- Get scattering & init-state parameters

  const ScatteringParams & sc_params  = interaction -> GetScatteringParams();
  const InitialState &     init_state = interaction -> GetInitialState();

  double E      = init_state.GetProbeE(kRfLab); // neutrino energy
  double Mnuc   = kNucleonMass;
  double Mpi    = kPionMass;
  double x      = sc_params.x(); // bjorken x
  double y      = sc_params.y(); // inelasticity y

  if (E<=Mpi || x<=0 || x>=1 || y<=Mpi/E || y>=1) {

    LOG("ReinSeghalCoh", pINFO)
                 << "d2xsec/dxdy[COH] (x= " << x << ", y="
                                             << y << ", E=" << E << ") = 0";
    return 0.;
  }

  //----- Compute the coherent NC pi0 production d2xsec/dxdy
  //      see page 34 in Nucl.Phys.B223:29-144 (1983)

  double Q2     = 2.*x*y*Mnuc*E;  // momentum transfer Q2>0
  double A      = (double) init_state.GetTarget().A(); // mass number
  double A2     = TMath::Power(A,2.);
  double A_3    = TMath::Power(A,1./3.);
  double Gf     = kGF_2 * Mnuc / (32 * kPi_3);
  double fp     = 0.93 * Mpi; // pion decay constant
  double fp2    = TMath::Power(fp,2.);
  double Epi    = y*E; // pion energy
  double ma     = this->Ma(); // Axial Mass [default can be overriden]
  double ma2    = TMath::Power(ma,2);
  double propg  = TMath::Power(ma2/(ma2+Q2),2.); // propagator term
  double r      = this->ReImPiApl(); // Re/Im Fwd Ampl. [def. can be overriden]
  double r2     = TMath::Power(r,2.);
  double sTot   = utils::xsec::TotalPionNucleonXSec(Epi); // tot. pi+N xsec
  double sTot2  = TMath::Power(sTot,2.);
  double sInel  = utils::xsec::InelasticPionNucleonXSec(Epi); // inel. pi+N xsec
  double Ro     = this->NuclSizeScale(); // nuclear size scale parameter
  double Ro2    = TMath::Power(Ro,2.);

  // effect of pion absorption in the nucleus
  double Fabs   = TMath::Exp( -9.*A_3*sInel / (16.*kPi*Ro2) );

  // the xsec in Nucl.Phys.B223:29-144 (1983) is d^3xsec/dxdydt but the only
  // t-dependent factor is an exp(-bt) so it can be integrated analyticaly
  double Epi2   = TMath::Power(Epi,2.);
  double Mpi2   = kPionMass_2;
  double R      = Ro * A_3; // nuclear radius
  double R2     = TMath::Power(R,2.);
  double b      = 0.33333 * R2;
  double tA     = 1. + Mnuc*x/Epi - 0.5*Mpi2/Epi2;
  double tB     = TMath::Sqrt(1. + 2*Mnuc*x/Epi) * TMath::Sqrt(1.-Mpi2/Epi2);
  double tmin   = 2*Epi2 * (tA-tB);
  double tmax   = 2*Epi2 * (tA+tB);

  double tint   = (TMath::Exp(-b*tmin) - TMath::Exp(-b*tmax))/b; // t integration factor

  double xsec = Gf * fp2 * A2 * E * (1-y) * sTot2 * (1+r2) * propg * Fabs * tint;

  LOG("ReinSeghalCoh", pDEBUG)
      << "\n momentum transfer ................ Q2    = " << Q2
      << "\n mass number ...................... A     = " << A
      << "\n pion energy ...................... Epi   = " << Epi
      << "\n propagator term .................. propg = " << propg
      << "\n Re/Im of fwd pion scat. ampl. .... Re/Im = " << r
      << "\n total pi+N cross section ......... sigT  = " << sTot
      << "\n inelastic pi+N cross section ..... sigI  = " << sInel
      << "\n nuclear size scale ............... Ro    = " << Ro
      << "\n pion absorption factor ........... Fabs  = " << Fabs
      << "\n t integration range .............. [" << tmin << "," << tmax << "]"
      << "\n t integration factor ............. tint  =" << tint;


  LOG("ReinSeghalCoh", pINFO)
                << "d2xsec/dxdy[COH] (x= " << x << ", y="
                                         << y << ", E=" << E << ") = "<< xsec;
  return xsec;
}
//____________________________________________________________________________
double ReinSeghalCOHPXSec::Ma(void) const
{
// allow the default Ma to be overriden by the value at the config registry

  if (fConfig->Exists("Ma")) return fConfig->GetDouble("Ma");
  else                       return kCohMa;
}
//____________________________________________________________________________
double ReinSeghalCOHPXSec::ReImPiApl(void) const
{
// allow the default Re/Im {f_{pi+N}(0)} to be overriden by the value at
// the config registry

  if (fConfig->Exists("Re-Im-Ampl")) return fConfig->GetDouble("Re-Im-Ampl");
  else                               return kCohReImAmpl;
}
//____________________________________________________________________________
double ReinSeghalCOHPXSec::NuclSizeScale(void) const
{
// allow the default nuclear size scale param to be overriden by the value at
// the config registry

  if (fConfig->Exists("Ro")) return fConfig->GetDouble("Ro");
  else                       return kCohR0;
}
//____________________________________________________________________________
