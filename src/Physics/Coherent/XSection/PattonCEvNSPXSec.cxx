//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Author: Costas Andreopoulos <costas.andreopoulos \at cern.ch>
         University of Liverpool & STFC Rutherford Appleton Laboratory

 For the class documentation see the corresponding header file.
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
#include "Physics/Coherent/XSection/PattonCEvNSPXSec.h"
#include "Physics/Coherent/XSection/NuclDensityMomentIntegrand.h"

using namespace genie;
using namespace genie::utils;
using namespace genie::constants;

//____________________________________________________________________________
PattonCEvNSPXSec::PattonCEvNSPXSec() :
XSecAlgorithmI("genie::PattonCEvNSPXSec")
{

}
//____________________________________________________________________________
PattonCEvNSPXSec::PattonCEvNSPXSec(string config) :
XSecAlgorithmI("genie::PattonCEvNSPXSec", config)
{

}
//____________________________________________________________________________
PattonCEvNSPXSec::~PattonCEvNSPXSec()
{

}
//____________________________________________________________________________
double PattonCEvNSPXSec::XSec(
  const Interaction * interaction, KinePhaseSpace_t kps) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  const InitialState & init_state = interaction -> InitState();
  const Kinematics &   kinematics = interaction -> Kine();
  const Target &       target     = init_state.Tgt();

  // User inputs to the calculation
  double E  = init_state.ProbeE(kRfLab); // neutrino energy, units: GeV
  double Q2 = kinematics.Q2(); // momentum transfer, units: GeV^2
  int    Z  = target.Z(); // number of protons
  int    N  = target.N(); // number of nucleons

  // Target atomic mass number and mass calculated from inputs
  int A   = Z+N;
  int target_nucleus_pdgc = pdg::IonPdgCode(A,Z);
  double M = PDGLibrary::Instance()->Find(target_nucleus_pdgc)->Mass(); // units: GeV
  LOG("CEvNS", pDEBUG) << "M = " << M << " GeV";

  // Calculation of nuclear recoil kinetic energy computed from input Q2
  double TA = Q2*E / (2*E*M+Q2); // nuclear recoil kinetic energy

  LOG("CEvNS", pDEBUG)
    << "Q2 = " << Q2 << " GeV^2, E = " << E << " GeV "
    << "--> TA = " << TA << " GeV";

  // auxiliary variables
  double E2  = E*E;
  double TA2 = TA*TA;
  double Q4  = Q2*Q2;
  double Q6  = Q2*Q4;

  // Calculation of weak charge
  // double Qw = N - Z*(1-fSin2thw);
  // Qw^2/4 x-section factor in arXiv:1207.0693v1 not needed here.
  // 1/4 was absorbed in the constant front factor (below) and Qw^2 factor would
  // have cancelled with ignored 1/Qw factor in the form factor F.

  // Calculation of nuclear density moments used for the evaluation
  // of the neutron form factor
  double avg_density = this->NuclearDensityMoment(A, 0); // units:: fm^-3
  double Rn2 = this->NuclearDensityMoment(A, 2) / avg_density; // units: fm^2
  double Rn4 = this->NuclearDensityMoment(A, 4) / avg_density; // units: fm^4
  double Rn6 = this->NuclearDensityMoment(A, 6) / avg_density; // units: fm^6

  LOG("CEvNS", pDEBUG)
    << "Nuclear density moments:"
    << " <Rn^2> = " << Rn2 << " fm^2,"
    << " <Rn^4> = " << Rn4 << " fm^4,"
    << " <Rn^6> = " << Rn6 << " fm^6";

  Rn2 *= TMath::Power(units::fm, 2.); // units: GeV^-2
  Rn4 *= TMath::Power(units::fm, 4.); // units: GeV^-4
  Rn6 *= TMath::Power(units::fm, 6.); // units: GeV^-6

  // Calculation of proton form factor
  // Form factor is neglected since it is multiplied with a small factor 1-4sin^2(\theta_{w})
  double Fp = 0; // units: -
  // Calculation of neutron form factor
  // Using a Taylor expansion of sin(Qr) and keeping the first three terms (shown to be
  // sufficient for approximating the full Fn calculation, even for heavy nuclei)
  double Fn = N * (1 - Q2*Rn2/6. + Q4*Rn4/120. - Q6*Rn6/5040.); // units: -
  // Overall form factor
  double F  = (Fn - (1-4*fSin2thw)*Fp); // units: -
  F = TMath::Max(0.,F);
  double F2 = F*F; // units: -

  LOG("CEvNS", pDEBUG)
    << "Form factors: Fp = " << Fp << ", Fn = " << Fn << ", F = " << F;

  // dsig/dTA calculation
  double const_factor = 0.125*kGF2/kPi; // units: GeV^-4
  double kinematic_term = M * (2 - 2*TA/E + TA2/E2 - M*TA/E2); // units: GeV
  kinematic_term = TMath::Max(0., kinematic_term);

  LOG("CEvNS", pDEBUG)
    << "kinematic term: " << kinematic_term;

  double xsec = const_factor * kinematic_term * F2; // units: GeV^-3 (area/GeV)

  LOG("CEvNS", pINFO)
    << "dsig[vA,CEvNS]/dTA (Ev =  "
    << E << " GeV, Q2 = "<< Q2 << " GeV^2; TA = " << TA << " GeV) = "
    << xsec/(units::cm2) << " cm^2/GeV";

  // The algorithm computes dxsec/dTA
  // Check whether variable tranformation is needed
  if(kps!=kPSTAfE) {
    // TODO: Move the calculation in utils::kinematics
    // double J = utils::kinematics::Jacobian(interaction,kPSQ2fE,kps);
    double J = 0;
    if(kps==kPSQ2fE) {
        J = 2*E2*M / TMath::Power(2*E*M+Q2, 2.); // units: GeV^-1
    }
    xsec *= J; // units: GeV^-4 (area/GeV^2)
  }

  return xsec;
}
//____________________________________________________________________________
double PattonCEvNSPXSec::NuclearDensityMoment(int A, int k) const
{
  // Calculate moments of the nuclear density
  // Inputs:
  //   - atomic mass number, A
  //   - integer k specifying required nuclear density moment
  // Output:
  //   - nuclear density moment in units of fm^k
  //
  // THINGS TO DO:
  // 1) The calculation can be stored, as it is required only once per nucleus. 
  //    The calculation is very fast so it doesn't matter.

  ROOT::Math::IBaseFunctionOneDim * integrand = new
              utils::gsl::wrap::NuclDensityMomentIntegrand(A,k);

  ROOT::Math::IntegrationOneDim::Type ig_type =
          utils::gsl::Integration1DimTypeFromString("adaptive");

  double R0 = utils::nuclear::Radius(A); // units: fm
  double rmin = 0; // units: fm
  double rmax = fNuclDensMomentCalc_UpperIntegrationLimit * R0; // units: fm

  ROOT::Math::Integrator ig(
    *integrand,ig_type,
    fNuclDensMomentCalc_AbsoluteTolerance,
    fNuclDensMomentCalc_RelativeTolerance,
    fNuclDensMomentCalc_MaxNumOfEvaluations);
  double moment = 2 * constants::kPi * ig.Integral(rmin, rmax); // units: fm^k

  delete integrand;

  return moment;
}
//____________________________________________________________________________
double PattonCEvNSPXSec::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool PattonCEvNSPXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const ProcessInfo & proc_info = interaction->ProcInfo();
  if(!proc_info.IsCoherentElastic()) return false;

  const InitialState & init_state = interaction->InitState();
  const Target & target = init_state.Tgt();
  if(!target.IsNucleus()) return false;

  return true;
}
//____________________________________________________________________________
void PattonCEvNSPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void PattonCEvNSPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void PattonCEvNSPXSec::LoadConfig(void)
{
  double thw = 0.;
  this->GetParam("WeinbergAngle", thw);
  fSin2thw = TMath::Power(TMath::Sin(thw), 2.);

  this->GetParamDef(
          "nuclear-density-moment-gsl-upper-limit", 
          fNuclDensMomentCalc_UpperIntegrationLimit, 
          10.);  // in nuclear radii
  this->GetParamDef(
          "nuclear-density-moment-gsl-rel-tol",    
          fNuclDensMomentCalc_RelativeTolerance, 
          1E-3); 
  this->GetParamDef(
          "nuclear-density-moment-gsl-abs-tol",    
          fNuclDensMomentCalc_AbsoluteTolerance, 
          1.); 
  this->GetParamDef(
          "nuclear-density-moment-gsl-max-eval",  
          fNuclDensMomentCalc_MaxNumOfEvaluations, 
          10000); 

  fXSecIntegrator =
      dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);
}
//____________________________________________________________________________
