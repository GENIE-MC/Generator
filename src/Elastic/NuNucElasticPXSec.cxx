//____________________________________________________________________________
/*!

\class    genie::NuNucElasticPXSec

\brief    Differential cross section dxsec/dQ^2 for v+N / vbar+N elastic 
          scattering. \n
          NuNucElasticPXSec is a concrete implementation of the
          XSecAlgorithmI interface. \n

\ref      L.A.Ahrens et al., Physical Review D, VOL 35,3:785 (1987)

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  Fabruary 15, 2005

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Algorithm/AlgFactory.h"
#include "Conventions/Constants.h"
#include "Conventions/RefFrame.h"
#include "Elastic/NuNucElasticPXSec.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGUtils.h"
#include "Utils/KineUtils.h"
#include "Utils/MathUtils.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
NuNucElasticPXSec::NuNucElasticPXSec() :
XSecAlgorithmI("genie::NuNucElasticPXSec")
{

}
//____________________________________________________________________________
NuNucElasticPXSec::NuNucElasticPXSec(string config) :
XSecAlgorithmI("genie::NuNucElasticPXSec", config)
{

}
//____________________________________________________________________________
NuNucElasticPXSec::~NuNucElasticPXSec()
{

}
//____________________________________________________________________________
double NuNucElasticPXSec::XSec(const Interaction * interaction) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  const InitialState & init_state = interaction -> GetInitialState();
  const Kinematics &   kinematics = interaction -> GetKinematics();

  double E     = init_state.GetProbeE(kRfStruckNucAtRest);
  double Q2    = kinematics.Q2();
  double M     = init_state.GetTarget().StruckNucleonMass();
  double M2    = M*M;
  double M4    = M2*M2;
  double E2    = E*E;
  double sig0  = kGF_2*M2 / (8*kPi*E2);
  double su    = 4*M*E - Q2; // s-u
  double su2   = su*su;      // (s-u)^2
  double qm2   = Q2 / M2;
  double qmf   = TMath::Power(1+qm2,2);
  double alpha = 1.-2.*kSin8w_2;
  double gamma = -0.66666667*kSin8w_2;
  double df1   = TMath::Power( 1.+qm2, 2.);
  double Gv3   = 0.5 * (1+kMuP-kMuN) / df1;
  double Gv0   = 1.5 * (1+kMuP+kMuN) / df1;
  double df2   = (1. + 0.25*qm2) * TMath::Power( 1.+qm2, 2.);
  double Fv3   = 0.5 * (kMuP-kMuN) / df2;
  double Fv0   = 1.5 * (kMuP+kMuN) / df2;
  double ga    = 0.5 * kElGa0 / qmf;           // El.nucl. form factor GA
  double f2    = alpha*Fv3 + gamma*Fv0;        // El.nucl. form factor F2
  double f1    = alpha*Gv3 + gamma* Gv0 - f2;  // El.nucl. form factor F1
  double ga2   = ga*ga;
  double f12   = f1*f1;
  double f22   = f2*f2;

  double A = qm2 * ( ga2 * (1 + 0.25*qm2)                    +
                     (0.25 * f22*qm2 - f12) * (1 - 0.25*qm2) +
                     f1*f2*qm2
                   );
  double B = qm2 * ga * (f1+f2);
  double C = 0.25 * (ga2 + f12 + 0.25*f22*qm2);

  int sign = 1;
  if( pdg::IsAntiNeutrino(init_state.GetProbePDGCode()) ) sign = -1;

  double dsig_dQ2 = sig0 * (A + sign*B*su/M2 + C*su2/M4);

  LOG("NuNucEl", pDEBUG)
     << "dXSec[el]/dQ2 (Ev = " << E << ", Q2 = " << Q2 << ") = " << dsig_dQ2;

  return dsig_dQ2;
}
//____________________________________________________________________________
bool NuNucElasticPXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  return true;
}
//____________________________________________________________________________
bool NuNucElasticPXSec::ValidKinematics(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipKinematicChk)) return true;

  const Kinematics & kinematics = interaction -> GetKinematics();
  double Q2 = kinematics.Q2();

  Range1D_t rQ2 = utils::kinematics::Q2Range_M(interaction);

  return (utils::math::IsWithinLimits(Q2, rQ2));
}
//____________________________________________________________________________

