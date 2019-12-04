//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Steven Gardiner <gardiner \at fnal.gov>
         Fermi National Acclerator Laboratory

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "Math/IFunction.h"
#include "Math/Integrator.h"

#include "Physics/QuasiElastic/XSection/UnifiedQELPXSec.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Physics/QuasiElastic/XSection/QELFormFactors.h"
#include "Physics/QuasiElastic/XSection/QELFormFactorsModelI.h"
#include "Physics/QuasiElastic/XSection/QELUtils.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Physics/NuclearState/NuclearModelI.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;
using namespace genie::utils;

//____________________________________________________________________________
UnifiedQELPXSec::UnifiedQELPXSec() :
XSecAlgorithmI("genie::UnifiedQELPXSec")
{

}
//____________________________________________________________________________
UnifiedQELPXSec::UnifiedQELPXSec(string config) :
XSecAlgorithmI("genie::UnifiedQELPXSec", config)
{

}
//____________________________________________________________________________
UnifiedQELPXSec::~UnifiedQELPXSec()
{

}
//____________________________________________________________________________
double UnifiedQELPXSec::XSec(const Interaction* interaction,
  KinePhaseSpace_t kps) const
{
  if ( !this->ValidProcess(interaction) ) return 0.;
  if ( !this->ValidKinematics(interaction) ) return 0.;

  // Get kinematics and init-state parameters
  const Kinematics&   kinematics = interaction->Kine();
  const InitialState& init_state = interaction->InitState();
  const Target& target = init_state.Tgt();

  TLorentzVector* temp_probeP4 = init_state.GetProbeP4( kRfLab );
  TLorentzVector probeP4 = *temp_probeP4;
  delete temp_probeP4;
  double p_probe = probeP4.P();

  const TLorentzVector& p4Ni = target.HitNucP4();
  double mNi = target.HitNucMass(); // on-shell initial hit nucleon mass
  double pNi = p4Ni.P();
  double E_NiOnShell = std::sqrt(pNi*pNi + mNi*mNi);
  TLorentzVector p4NiOnShell = TLorentzVector(p4Ni.Vect(), E_NiOnShell);
  double epsilon_B = E_NiOnShell - p4Ni.E();

  const TLorentzVector& lepP4 = kinematics.FSLeptonP4();
  double E_lep = lepP4.E();

  const TLorentzVector& p4Nf = kinematics.HadSystP4();
  double E_Nf = p4Nf.E();

  // Start xsec calculation with overall phase space factor and
  // scaling factor from XML
  double xsec = fXSecScale / (64.*kPi*kPi * p_probe * E_NiOnShell
    * E_lep * E_Nf);

  // If we're dealing with a nuclear target, then apply Pauli blocking as
  // needed
  if ( target.IsNucleus() && !interaction->TestBit(kIAssumeFreeNucleon) ) {
    double kF = fPauliBlocker->GetFermiMomentum(target,
      interaction->RecoilNucleonPdg(), target.HitNucPosition());
    if ( p4Nf.P() < kF ) return 0.;
  }

  // Scale cross section by the number of active nucleons
  int hit_nuc_pdg = target.HitNucPdg();
  // Number of active nucleons in the target
  int num_active = pdg::IsProton(hit_nuc_pdg) ? target.Z() : target.N();

  xsec *= num_active;

  // Compute form factors using Q2tilde (the effective Q2 value after
  // binding energy corrections)
  TLorentzVector qP4 = probeP4 - lepP4;
  TLorentzVector qTildeP4 = qP4;
  qTildeP4.SetE( qP4.E() - epsilon_B );

  double Q2 = -1. * qP4.M2();
  double Q2tilde = -1. * qTildeP4.M2();

  // Make sure Q2 is physical
  if ( Q2 <= 0. ) {
    LOG("UnifiedQELPXSec", pWARN) << "Q2 <= 0, returning xsec = 0";
    return 0.;
  }

  // Get the correct couplings and form factors model for the current
  // interaction
  double coupling_factor = 1.;
  const genie::ProcessInfo& proc_info = interaction->ProcInfo();
  if ( proc_info.IsWeakCC() ) {
    coupling_factor = kGF2 * fCos8c2 / 2.;
    fFormFactors.SetModel( fCCFormFactorsModel );
  }
  else if ( proc_info.IsWeakNC() ) {
    coupling_factor = kGF2 / 2.;
    fFormFactors.SetModel( fNCFormFactorsModel );
  }
  else if ( proc_info.IsEM() ) {
    // Lorentz-Heaviside natural units: e^4 = 16 * kPi^2 * kAem2
    coupling_factor = 16. * kPi * kPi * kAem2 / ( Q2*Q2 );
    fFormFactors.SetModel( fEMFormFactorsModel );
  }
  else {
    LOG("UnifiedQELPXSec", pERROR) << "Unrecognized process type encountered"
      << " in genie::CBFSpecFuncQELPXSec::XSec()";
    return 0.;
  }

  // Apply the coupling factor to the differential cross section
  xsec *= coupling_factor;

  // For a bound hit nucleon, the corrected energy transfer qTilde0 needs to be
  // positive to enforce energy conservation
  if ( qTildeP4.E() <= 0. && interaction->InitState().Tgt().IsNucleus()
    && !interaction->TestBit(kIAssumeFreeNucleon) ) return 0.;

  // Set Q2 to Q2tilde while computing form factors
  interaction->KinePtr()->SetQ2( Q2tilde );

  // Evaluate the form factors
  fFormFactors.Calculate( interaction );

  // Now that we've calculated them, store the true Q2 value
  interaction->KinePtr()->SetQ2( Q2 );

  // Compute the tensor contraction
  LeptonTensor L_munu( *interaction );
  FreeNucleonTensor ATilde_munu( *interaction, fFormFactors );

  std::complex<double> contraction = L_munu * ATilde_munu;

  if ( std::abs(contraction.imag()) > kASmallNum ) {
    LOG("UnifiedQELPXSec", pWARN) << "Tensor contraction has nonvanishing"
      << " imaginary part!";
  }

  // Apply the tensor contraction to the cross section
  xsec *= contraction.real();

  // Multiply by the analytic solution of the energy-conserving delta function
  // used by the kPSQELEvGen phase space.
  xsec *= genie::utils::EnergyDeltaFunctionSolutionQEL( *interaction );

  // Check whether variable tranformation is needed
  if ( kps != kPSQELEvGen ) {
    // Compute the appropriate Jacobian for transformation to the requested
    // phase space
    double J = utils::kinematics::Jacobian(interaction, kPSQELEvGen, kps);
    xsec *= J;
  }

  return xsec;
}
//____________________________________________________________________________
double UnifiedQELPXSec::Integral(const Interaction* in) const
{
  // Intended for use with genie::NewQELXSec, which is smart
  // enough to handle free nucleon vs. nuclear targets, different
  // nuclear models, etc.
  return fXSecIntegrator->Integrate(this, in);
}
//____________________________________________________________________________
bool UnifiedQELPXSec::ValidProcess(const Interaction* interaction) const
{
  if ( interaction->TestBit(kISkipProcessChk) ) return true;

  const InitialState& init_state = interaction->InitState();
  const ProcessInfo&  proc_info  = interaction->ProcInfo();

  if ( !proc_info.IsQuasiElastic() ) return false;

  if ( proc_info.IsEM() || proc_info.IsWeakNC() ) return true;

  // For weak CC interactions, check that the hit nucleon "works"
  else if ( proc_info.IsWeakCC() ) {
    int nucleon_pdg = init_state.Tgt().HitNucPdg();
    int probe_pdg = init_state.ProbePdg();

    if ( pdg::IsNeutron(nucleon_pdg) && pdg::IsNeutrino(probe_pdg) ) {
      return true;
    }
    else if ( pdg::IsProton(nucleon_pdg) && pdg::IsAntiNeutrino(probe_pdg) ) {
      return true;
    }
    else return false;
  }

  return false;
}
//____________________________________________________________________________
void UnifiedQELPXSec::Configure(const Registry& config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void UnifiedQELPXSec::Configure(std::string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void UnifiedQELPXSec::LoadConfig(void)
{
  double thc;
  GetParam( "CabibboAngle", thc ) ;
  fCos8c2 = TMath::Power(TMath::Cos(thc), 2);

  // cross section scaling factor
  GetParam( "QEL-CC-XSecScale", fXSecScale ) ;

  // load QEL form factors models
  fCCFormFactorsModel = dynamic_cast<const QELFormFactorsModelI*>(
    this->SubAlg("CCFormFactorsAlg") );
  assert( fCCFormFactorsModel );

  fNCFormFactorsModel = dynamic_cast<const QELFormFactorsModelI*>(
    this->SubAlg("NCFormFactorsAlg") );
  assert( fNCFormFactorsModel );

  fEMFormFactorsModel = dynamic_cast<const QELFormFactorsModelI*>(
    this->SubAlg("EMFormFactorsAlg") );
  assert( fEMFormFactorsModel );

  // Attach CC model for now. This will be updated later.
  fFormFactors.SetModel( fCCFormFactorsModel );

  // load xsec integrator
  fXSecIntegrator = dynamic_cast<const XSecIntegratorI*>(
    this->SubAlg("XSec-Integrator") );
  assert(fXSecIntegrator);

  // get nuclear model
  fNuclModel = dynamic_cast<const NuclearModelI*>(
    this->SubAlg("IntegralNuclearModel") );
  assert(fNuclModel);

  // get the algorithm ID for the PauliBlocker
  RgAlg pauliBlockID;
  GetParamDef( "PauliBlockerAlg", pauliBlockID,
    RgAlg("genie::PauliBlocker", "Default") );
  AlgId pbID = AlgId( pauliBlockID );

  AlgFactory* algf = AlgFactory::Instance();
  fPauliBlocker = dynamic_cast<const PauliBlocker*>(
    algf->GetAlgorithm(pbID) );
  assert( fPauliBlocker );
}
