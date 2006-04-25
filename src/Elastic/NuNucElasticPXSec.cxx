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

#include "Algorithm/AlgConfigPool.h"
#include "Conventions/RefFrame.h"
#include "Conventions/Constants.h"
#include "Elastic/NuNucElasticPXSec.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGUtils.h"
#include "Utils/KineUtils.h"
#include "Utils/MathUtils.h"
#include "Utils/NuclearUtils.h"

using namespace genie;
using namespace genie::utils;
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
  double M2    = TMath::Power(M, 2);
  double M4    = TMath::Power(M2,2);
  double E2    = TMath::Power(E, 2);

  //----- compute the form factors

  double qmv2  = TMath::Power(1 + Q2/fMv2, 2);
  double qma2  = TMath::Power(1 + Q2/fMa2, 2);
  double tau   = 0.25 * Q2/M2;

  double Gv3   = 0.5 * (1+fMuP-fMuN) / qmv2;
  double Gv0   = 1.5 * (1+fMuP+fMuN) / qmv2;
  double Fv3   = 0.5 * (fMuP-fMuN) / ((1+tau)*qmv2);
  double Fv0   = 1.5 * (fMuP+fMuN) / ((1+tau)*qmv2);
  double ga    = -0.5 * fFa0 * (1+fEta)/ qma2;    // El.nucl. form factor GA
  double f2    = fkAlpha*Fv3 + fkGamma*Fv0;       // El.nucl. form factor F2
  double f1    = fkAlpha*Gv3 + fkGamma* Gv0 - f2; // El.nucl. form factor F1
  double ga2   = TMath::Power(ga,2);
  double f12   = TMath::Power(f1,2);
  double f22   = TMath::Power(f2,2);

  //----- compute the free nucleon cross section

  double sig0  = kGF2*M2 / (8*kPi*E2);
  double su    = 4*M*E - Q2;  // s-u
  double su2   = TMath::Power(su,2);    

  double A = 4*tau*(ga2*(1+tau) - f12*(1-tau) + f22*tau*(1-tau) + 4*tau*f1*f2);
  double B = 4*tau*ga*(f1+f2);
  double C = 0.25*(ga2 + f12 + f22*tau);

  int sign = 1;
  if( pdg::IsAntiNeutrino(init_state.GetProbePDGCode()) ) sign = -1;

  double xsec = sig0 * (A + sign*B*su/M2 + C*su2/M4);

  LOG("NuNucEl", pDEBUG)
    << "dXSec[vN,El]/dQ2 [FreeN](Ev = "<< E<< ", Q2 = "<< Q2 << ") = "<< xsec;

  //----- if requested return the free nucleon xsec even for input nuclear tgt
  if( interaction->TestBit(kIAssumeFreeNucleon) ) return xsec;

  //----- compute nuclear suppression factor
  //      (R(Q2) is adapted from NeuGEN - see comments therein)
  double R = nuclear::NuclQELXSecSuppression("Default", 0.5, interaction);

  //----- number of scattering centers in the target
  const Target & target = init_state.GetTarget();
  int nucpdgc = target.StruckNucleonPDGCode();
  int NNucl = (pdg::IsProton(nucpdgc)) ? target.Z() : target.N();

  LOG("NuNucEl", pDEBUG)
       << "Nuclear suppression factor R(Q2) = " << R << ", NNucl = " << NNucl;

  xsec *= (R*NNucl); // nuclear xsec

  LOG("NuNucEl", pDEBUG)
   << "dXSec[vN,El]/dQ2 [Nuclear](Ev = "<< E<< ", Q2 = "<< Q2<< ") = "<< xsec;

  return xsec;
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
void NuNucElasticPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void NuNucElasticPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void NuNucElasticPXSec::LoadConfig(void)
{
  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();

  // alpha and gamma
  double thw = fConfig->GetDoubleDef(
                          "weinberg-angle", gc->GetDouble("WeinbergAngle"));
  double sin2thw = TMath::Power(TMath::Sin(thw), 2);

  fkAlpha = 1.-2.*sin2thw;
  fkGamma = -0.66666667*sin2thw;

  // eta and Fa(q2=0)
  fEta = fConfig->GetDoubleDef("eta", gc->GetDouble("EL-Axial-Eta"));
  fFa0 = fConfig->GetDoubleDef("Fa0", gc->GetDouble("QEL-FA0"));

  // axial and vector masses
  double ma = fConfig->GetDoubleDef("Ma", gc->GetDouble("QEL-Ma"));
  double mv = fConfig->GetDoubleDef("Mv", gc->GetDouble("QEL-Mv"));

  fMa2 = TMath::Power(ma,2);
  fMv2 = TMath::Power(mv,2);

  // anomalous magnetic moments
  fMuP = fConfig->GetDoubleDef("MuP", gc->GetDouble("AnomMagnMoment-P"));
  fMuN = fConfig->GetDoubleDef("MuN", gc->GetDouble("AnomMagnMoment-N"));
}
//____________________________________________________________________________
