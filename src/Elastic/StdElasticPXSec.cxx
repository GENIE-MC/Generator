//____________________________________________________________________________
/*!

\class    genie::StdElasticPXSec

\brief    Standard differential cross section dxsec/dQ^2 for v+N / vbar+N
          elastic scattering.

          StdElasticPXSec is a concrete implementation of the
          XSecAlgorithmI interface. \n

\ref      L.A.Ahrens et al., Physical Review D, VOL 35,3:785 (1987)

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  Fabruary 15, 2005

*/
//____________________________________________________________________________

#include <TMath.h>

#include "AlgFactory/AlgFactory.h"
#include "Conventions/Constants.h"
#include "Conventions/RefFrame.h"
#include "Elastic/StdElasticPXSec.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGUtils.h"
#include "Utils/KineLimits.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
StdElasticPXSec::StdElasticPXSec() :
XSecAlgorithmI()
{
  fName = "genie::StdElasticPXSec";
}
//____________________________________________________________________________
StdElasticPXSec::StdElasticPXSec(const char * param_set) :
XSecAlgorithmI(param_set)
{
  fName = "genie::StdElasticPXSec";

  FindConfig();
}
//____________________________________________________________________________
StdElasticPXSec::~StdElasticPXSec()
{

}
//____________________________________________________________________________
double StdElasticPXSec::XSec(const Interaction * interaction) const
{
  //----- get initial & final state information

  const InitialState & init_state = interaction -> GetInitialState();

  TLorentzVector * nu_p4 = init_state.GetProbeP4(kRfStruckNucAtRest);

  double E  = nu_p4->Energy();
  
  delete nu_p4;

  double Q2  = interaction->GetScatteringParams().Q2();

  Range1D_t rQ2 = kine_limits::Q2Range_M(interaction);

  if ( Q2 < rQ2.min || Q2 > rQ2.max ) return 0;

  //-- compute cross section
  
  double M   = init_state.GetTarget().StruckNucleonMass();
  double M2  = M*M;
  double M4  = M2*M2;
  double E2  = E*E;

  double sig0 = kGF_2*M2 / (8*kPi*E2);

  double su   = 4*M*E - Q2; // s-u
  double su2  = su*su;      // (s-u)^2
  double qm2  = Q2 / M2;

  double ga   = this->GA(interaction);
  double f1   = this->F1(interaction);
  double f2   = this->F2(interaction);

  double ga2  = ga*ga;
  double f12  = f1*f1;
  double f22  = f2*f2;
  
  double A = qm2 * ( ga2 * (1 + 0.25*qm2)                    +
                     (0.25 * f22*qm2 - f12) * (1 - 0.25*qm2) +
                     f1*f2*qm2
                   );  
  double B = qm2 * ga * (f1+f2);
  double C = 0.25 * (ga2 + f12 + 0.25*f22*qm2);
  
  int sign = 1;
  if( pdg::IsAntiNeutrino(init_state.GetProbePDGCode()) ) sign = -1;

  double dsig_dQ2 = sig0 * (A + sign*B*su/M2 + C*su2/M4);

  LOG("InverseMuDecay", pDEBUG)
      << "dXSec[el]/dQ2 (Ev = " << E << ", Q2 = " << Q2 << ") = " << dsig_dQ2;

  return dsig_dQ2;
}
//____________________________________________________________________________
double StdElasticPXSec::GA(const Interaction * interaction) const
{
// Elastic nucleon form factor GA

  double Q2  = interaction->GetScatteringParams().Q2();
  double M   = interaction->GetInitialState().GetTarget().StruckNucleonMass();
  double M2  = M*M;
  
  double qm2 = Q2 / M2;
  
  double Ga  = 0.5 * kElGa0 / TMath::Power(1+qm2,2);

  return Ga;
}
//____________________________________________________________________________
double StdElasticPXSec::F1(const Interaction * interaction) const
{
// Elastic nucleon form factor F1

  double Q2  = interaction->GetScatteringParams().Q2();
  double M   = interaction->GetInitialState().GetTarget().StruckNucleonMass();
  double M2  = M*M;

  double qm2 = Q2 / M2;
  double d   = TMath::Power( 1.+qm2, 2.);

  double Gv3   = 0.5 * (1+kMuP-kMuN) / d;
  double Gv0   = 1.5 * (1+kMuP+kMuN) / d;
  
  double f2 = this->F2(interaction);

  double f1 = this->Alpha() * Gv3 + this->Gamma() * Gv0 - f2;

  return f1;
}
//____________________________________________________________________________
double StdElasticPXSec::F2(const Interaction * interaction) const
{
// Elastic nucleon form factor F2

  double Q2  = interaction->GetScatteringParams().Q2();
  double M   = interaction->GetInitialState().GetTarget().StruckNucleonMass();
  double M2  = M*M;
  
  double qm2 = Q2 / M2;
  double d   = (1. + 0.25*qm2) * TMath::Power( 1.+qm2, 2.);

  double Fv3   = 0.5 * (kMuP-kMuN) / d;
  double Fv0   = 1.5 * (kMuP+kMuN) / d;

  double f2    = this->Alpha() * Fv3 + this->Gamma() * Fv0;

  return f2;
}
//____________________________________________________________________________
double StdElasticPXSec::Alpha(void) const
{
  return (1.-2.*kSin8w_2);
}  
//____________________________________________________________________________
double StdElasticPXSec::Gamma(void) const
{
  return (-0.66666667*kSin8w_2);
}
//____________________________________________________________________________

