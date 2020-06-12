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

  const InitialState & init_state = interaction -> InitState();
  const Kinematics &   kinematics = interaction -> Kine();
  const Target &       target     = init_state.Tgt();

  // User inputs to the calculation
  const double E  = init_state.ProbeE(kRfLab); // neutrino energy, units: GeV
  const double Q2 = kinematics.Q2(); // momentum transfer, units: GeV^2
  const double DNuEnergy =  kinematics.FSLeptonP4().E(); // E_N is the energy of the dark neutrino
  const unsigned int Z = target.Z(); // number of protons
  const unsigned int N = target.N(); // number of nucleons

  // Target atomic mass number and mass calculated from inputs
  const unsigned int A = Z + N;
  const int target_nucleus_pdgc = pdg::IonPdgCode(A, Z);
  const double M = PDGLibrary::Instance()->Find(target_nucleus_pdgc)->Mass(); // units: GeV
  LOG("DNu", pDEBUG) << "M = " << M << " GeV";

  // Calculation of nuclear recoil kinetic energy computed from input Q2
  // double TA = Q2*E / (2*E*M+Q2); // nuclear recoil kinetic energy

  // LOG("DNu", pDEBUG)
  //   << "Q2 = " << Q2 << " GeV^2, E = " << E << " GeV "
  //   << "--> TA = " << TA << " GeV";

  const double FF = fFF->FormFactor(Q2, target);

  // auxiliary variables
  const double E2  = E*E;
  const double Z2 = Z * Z;
  const double FF2 = FF * FF;

  //TODO DNu: is this one required?
  const double elec2 = 0.3*0.3;

  // double TA2 = TA*TA;
  // double Q4  = Q2*Q2;
  // double Q6  = Q2*Q4;

  if(! this -> ValidKinematics (interaction, DNuEnergy, fDNuMass2, E, M) ) return 0.;

  const double const_factor = .125 * elec2 / kPi;
  const double model_params = fEps2 * fTheta2 * fgD2;

  const double num_fact1 = ( FF2  * fDNuMass) * Z2;
  const double num_fact2 = (DNuEnergy+E+M)*fDNuMass2  - 2*M*(E2 + M*DNuEnergy + E2 - E*M);
  const double den_fact1 = 1. / (E2*M);
  const double den_fact2 = TMath::Power((fDMediatorMass2 - 2.*DNuEnergy*M + 2*E*M), -2.);

  const double xsec = const_factor * model_params
    * num_fact1 * num_fact2 * den_fact1 * den_fact2;
  // const double cross_section = (elec2 * FF2 * eps2 * theta2 * g_D2 * DNu_mass) *
  //   ( (DNu_energy+E+M)*DNu_mass2  - 2*M*(E2 + M*DNu_energy + E2 - E*M) )*Z2 *
  //   (1./ (8*kPi*E2*M *  TMath::Power((DZ_mass2 - 2*DNu*M + 2*E*M),2)));




  // // Calculation of weak charge
  // // double Qw = N - Z*(1-fSin2thw);
  // // Qw^2/4 x-section factor in arXiv:1207.0693v1 not needed here.
  // // 1/4 was absorbed in the constant front factor (below) and Qw^2 factor would
  // // have cancelled with ignored 1/Qw factor in the form factor F.

  // // Calculation of nuclear density moments used for the evaluation
  // // of the neutron form factor
  // double avg_density = this->NuclearDensityMoment(A, 0); // units:: fm^-3
  // double Rn2 = this->NuclearDensityMoment(A, 2) / avg_density; // units: fm^2
  // double Rn4 = this->NuclearDensityMoment(A, 4) / avg_density; // units: fm^4
  // double Rn6 = this->NuclearDensityMoment(A, 6) / avg_density; // units: fm^6

  // LOG("CEvNS", pDEBUG)
  //   << "Nuclear density moments:"
  //   << " <Rn^2> = " << Rn2 << " fm^2,"
  //   << " <Rn^4> = " << Rn4 << " fm^4,"
  //   << " <Rn^6> = " << Rn6 << " fm^6";

  // Rn2 *= TMath::Power(units::fm, 2.); // units: GeV^-2
  // Rn4 *= TMath::Power(units::fm, 4.); // units: GeV^-4
  // Rn6 *= TMath::Power(units::fm, 6.); // units: GeV^-6

  // // Calculation of proton form factor
  // // Form factor is neglected since it is multiplied with a small factor 1-4sin^2(\theta_{w})
  // double Fp = 0; // units: -
  // // Calculation of neutron form factor
  // // Using a Taylor expansion of sin(Qr) and keeping the first three terms (shown to be
  // // sufficient for approximating the full Fn calculation, even for heavy nuclei)
  // double Fn = N * (1 - Q2*Rn2/6. + Q4*Rn4/120. - Q6*Rn6/5040.); // units: -
  // // Overall form factor
  // double F  = (Fn - (1-4*fSin2thw)*Fp); // units: -
  // F = TMath::Max(0.,F);
  // double F2 = F*F; // units: -

  // LOG("CEvNS", pDEBUG)
  //   << "Form factors: Fp = " << Fp << ", Fn = " << Fn << ", F = " << F;

  // // dsig/dTA calculation
  // double const_factor = 0.125*kGF2/kPi; // units: GeV^-4
  // double kinematic_term = M * (2 - 2*TA/E + TA2/E2 - M*TA/E2); // units: GeV
  // kinematic_term = TMath::Max(0., kinematic_term);

  // LOG("CEvNS", pDEBUG)
  //   << "kinematic term: " << kinematic_term;

  // double xsec = const_factor * kinematic_term * F2; // units: GeV^-3 (area/GeV)

  // LOG("CEvNS", pINFO)
  //   << "dsig[vA,CEvNS]/dTA (Ev =  "
  //   << E << " GeV, Q2 = "<< Q2 << " GeV^2; TA = " << TA << " GeV) = "
  //   << xsec/(units::cm2) << " cm^2/GeV";

  // // The algorithm computes dxsec/dTA
  // // Check whether variable tranformation is needed
  // if(kps!=kPSTAfE) {
  //   // TODO: Move the calculation in utils::kinematics
  //   // double J = utils::kinematics::Jacobian(interaction,kPSQ2fE,kps);
  //   double J = 0;
  //   if(kps==kPSQ2fE) {
  //       J = 2*E2*M / TMath::Power(2*E*M+Q2, 2.); // units: GeV^-1
  //   }
  //   xsec *= J; // units: GeV^-4 (area/GeV^2)
  // }

  return xsec;
}
//____________________________________________________________________________
// double BertuzzoDNuCOHPXSec::NuclearDensityMoment(int A, int k) const
// {
//   // Calculate moments of the nuclear density
//   // Inputs:
//   //   - atomic mass number, A
//   //   - integer k specifying required nuclear density moment
//   // Output:
//   //   - nuclear density moment in units of fm^k
//   //
//   // THINGS TO DO:
//   // 1) The calculation can be stored, as it is required only once per nucleus.
//   //    The calculation is very fast so it doesn't matter.

//   ROOT::Math::IBaseFunctionOneDim * integrand = new
//               utils::gsl::wrap::NuclDensityMomentIntegrand(A,k);

//   ROOT::Math::IntegrationOneDim::Type ig_type =
//           utils::gsl::Integration1DimTypeFromString("adaptive");

//   double R0 = utils::nuclear::Radius(A); // units: fm
//   double rmin = 0; // units: fm
//   double rmax = fNuclDensMomentCalc_UpperIntegrationLimit * R0; // units: fm

//   ROOT::Math::Integrator ig(
//     *integrand,ig_type,
//     fNuclDensMomentCalc_AbsoluteTolerance,
//     fNuclDensMomentCalc_RelativeTolerance,
//     fNuclDensMomentCalc_MaxNumOfEvaluations);
//   double moment = 2 * constants::kPi * ig.Integral(rmin, rmax); // units: fm^k

//   delete integrand;

//   return moment;
// }
// //____________________________________________________________________________
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
bool BertuzzoDNuCOHPXSec::ValidKinematics(const Interaction* interaction,
                                          const double DNu_energy,
                                          const double DNu_mass2,
                                          const double E,
                                          const double M) const
{
  const double t1 = 4. * DNu_mass2 * (DNu_energy - E) * (M + E);
  const double t2 = 4. * M * (DNu_energy - E) * (DNu_energy*(M+2.*E) - M*E);
  const double t3 = DNu_mass2 * DNu_mass2;
  if( (t1 - t2 - t3) <= 0. ) return false;

  if(interaction->TestBit(kISkipKinematicChk)) return true;
  return true;
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
