//____________________________________________________________________________
/*
  Copyright (c) 2003-2020, The GENIE Collaboration
  For the full text of the license visit http://copyright.genie-mc.org

  Author: Iker de Icaza <i.de-icaza-astiz \at sussex.ac.uk>
  University of Sussex

  Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
  University of Liverpool & STFC Rutherford Appleton Laboratory
*/
//____________________________________________________________________________


#include <TMath.h>
#include <Math/Integrator.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/Numerical/GSLUtils.h"
#include "Physics/NuclearState/NuclearUtils.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Physics/DarkNeutrino/XSection/BertuzzoDNuCOHPXSec.h"
#include "Physics/DarkNeutrino/XSection/EngelFormFactor.h"
// #include "Physics/Coherent/XSection/NuclDensityMomentIntegrand.h"

using namespace genie;
using namespace genie::utils;
using namespace genie::constants;

//____________________________________________________________________________
BertuzzoDNuCOHPXSec::BertuzzoDNuCOHPXSec() :
XSecAlgorithmI("genie::BertuzzoDNuCOHPXSec")
{

}
//____________________________________________________________________________
BertuzzoDNuCOHPXSec::BertuzzoDNuCOHPXSec(string config) :
XSecAlgorithmI("genie::BertuzzoDNuCOHPXSec", config)
{

}
//____________________________________________________________________________
BertuzzoDNuCOHPXSec::~BertuzzoDNuCOHPXSec()
{

}
//____________________________________________________________________________
double BertuzzoDNuCOHPXSec::XSec(
  const Interaction * interaction, KinePhaseSpace_t kps) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  const InitialState & init_state = interaction -> InitState();
  const Kinematics &   kinematics = interaction -> Kine();
  const Target &       target     = init_state.Tgt();

  // User inputs to the calculation
  const double E  = init_state.ProbeE(kRfLab); // neutrino energy, units: GeV
  const double Q2 = kinematics.Q2(); // momentum transfer, units: GeV^2
  const double DNuE =  kinematics.FSLeptonP4().E(); // E_N is the energy of the dark neutrino
  const unsigned int Z = target.Z(); // number of protons
  const unsigned int N = target.N(); // number of nucleons

  // Target atomic mass number and mass calculated from inputs
  const unsigned int A = Z + N;
  const int target_nucleus_pdgc = pdg::IonPdgCode(A, Z);
  const double M = PDGLibrary::Instance()->Find(target_nucleus_pdgc)->Mass(); // units: GeV
  LOG("DNu", pDEBUG) << "M = " << M << " GeV";

  // LOG("DNu", pDEBUG)
  //   << "Q2 = " << Q2 << " GeV^2, E = " << E << " GeV "
  //   << "--> TA = " << TA << " GeV";

  const double FF = fFF->FormFactor(Q2, target);

  // auxiliary variables
  const double E2  = E * E;
  const double DNuE2  = DNuE * DNuE;
  const double Z2 = Z * Z;
  const double FF2 = FF * FF;

  //TODO DNu: is this one required?
  const double elec2 = 0.3*0.3;

  const double const_factor = .125 * elec2 / kPi;
  const double model_params = fEps2 * fTheta2 * fgD2;

  const double num_fact1 = ( FF2  * fDNuMass) * Z2;
  const double num_fact2 = (DNuE+E+M)*fDNuMass2  - 2.*M*(DNuE2 + M*DNuE + E2 - E*M);
  const double den_fact1 = 1. / (E2*M);
  const double den_fact2 = TMath::Power((fDMediatorMass2 - 2.*DNuE*M + 2*E*M), -2.);

  if(kps==kPSEDNufE) {
    // TODO DNu: BUG xsec is negative
    // const double xsec = const_factor * model_params
      // * num_fact1 * num_fact2 * den_fact1 * den_fact2;
    const double xsec_approx = (2.*const_factor) * model_params *
      (E - DNuE) * (M * fDNuMass *Z2) * (1./E2) * den_fact2;
    // const double cross_section = (elec2 * FF2 * eps2 * theta2 * g_D2 * DNu_mass) *
    //   ( (DNu_energy+E+M)*DNu_mass2  - 2*M*(E2 + M*DNu_energy + E2 - E*M) )*Z2 *
    //   (1./ (8*kPi*E2*M *  TMath::Power((DZ_mass2 - 2*DNu*M + 2*E*M),2)));
    // return xsec;
    return xsec_approx;
  }
  return 0.;
}
//____________________________________________________________________________
double BertuzzoDNuCOHPXSec::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool BertuzzoDNuCOHPXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const ProcessInfo & proc_info = interaction->ProcInfo();
  if(!proc_info.IsCoherentElastic()) return false;

  // TODO DNU: add other requirements

  const InitialState & init_state = interaction->InitState();
  const Target & target = init_state.Tgt();
  if(!target.IsNucleus()) return false;

  return true;
}
//____________________________________________________________________________
bool BertuzzoDNuCOHPXSec::ValidKinematics(const Interaction* interaction) const
{
  if(interaction->TestBit(kISkipKinematicChk)) return true;

  return interaction->PhaseSpace().IsAboveThreshold();

  // Pedro suggests, but it's an overkill:
  // const double M = interaction->InitState().Tgt().Mass();
  // const double DNuEnergy = interaction->Kine().FSLeptonP4().E();
  // const double t1 = 4. * fDNuMass2 * (DNuEnergy - E) * (M + E);
  // const double t2 = 4. * M * (DNuEnergy - E) * (DNuEnergy*(M+2.*E) - M*E);
  // const double t3 = fDNuMass2 * fDNuMass2;
  // if( (t1 - t2 - t3) <= 0. ) return false;
}
//____________________________________________________________________________
void BertuzzoDNuCOHPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void BertuzzoDNuCOHPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void BertuzzoDNuCOHPXSec::LoadConfig(void)
{
  double DKineticMixing = 0.;    // \varepsilon
  this->GetParam("Dark-KineticMixing", DKineticMixing);
  fEps2 = DKineticMixing * DKineticMixing;

  double DTheta = 0.;            // \theta
  this->GetParam("Dark-Theta", DTheta);
  fTheta2 = DTheta * DTheta;

  double DGaugeCoupling = 0.;   // g_D
  this->GetParam("Dark-GaugeCoupling", DGaugeCoupling);
  fgD2 = DGaugeCoupling * DGaugeCoupling;

  fDNuMass = 0.;
  this->GetParam("Dark-NeutrinoMass", fDNuMass);
  fDNuMass2 = fDNuMass * fDNuMass;

  fDMediatorMass = 0.;
  this->GetParam("Dark-MediatorMass", fDMediatorMass);
  fDMediatorMass2 = fDMediatorMass * fDMediatorMass;

  fXSecIntegrator =
      dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);
  fFF = dynamic_cast<const EngelFormFactor *> (this->SubAlg("FormFactor"));
  assert(fFF);
}
//____________________________________________________________________________
