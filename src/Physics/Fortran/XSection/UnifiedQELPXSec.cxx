//____________________________________________________________________________
/*
 Copyright (c) 2003-2023, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Noah Steinberg <nsteinbe \at fnal.gov>
         Steven Gardiner <gardiner \at fnal.gov>
         Fermi National Acclerator Laboratory

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "Math/IFunction.h"
#include "Math/Integrator.h"

#include "UnifiedQELPXSec.h"
#include "fortran_functions.h"
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
#include "Physics/NuclearState/NuclearUtils.h"

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
  double E_probe = probeP4.E();

  TLorentzVector p4Ni = target.HitNucP4();
  double mNi = target.HitNucMass(); // on-shell initial hit nucleon mass
  double pNi = p4Ni.P();
  double E_NiOnShell = std::sqrt(pNi*pNi + mNi*mNi);
  TLorentzVector p4NiOnShell = TLorentzVector(p4Ni.Vect(), E_NiOnShell);
  double epsilon_B = E_NiOnShell - p4Ni.E();

  TLorentzVector lepP4 = kinematics.FSLeptonP4();
  double E_lep = lepP4.E();
  double P_lep = lepP4.P();

  TLorentzVector p4Nf = kinematics.HadSystP4();
  double E_Nf = p4Nf.E();
  
 
  // Include phase space factors in cross section
  double xsec = fXSecScale / (E_lep * E_probe * E_NiOnShell * E_Nf) / 32. / kPi / kPi;

  // If we're dealing with a nuclear target, then apply Pauli blocking as
  // needed  
  if (fDoPauliBlocking && target.IsNucleus() && !interaction->TestBit(kIAssumeFreeNucleon) ) {
    double kF = fPauliBlocker->GetFermiMomentum(target, interaction->RecoilNucleonPdg(), target.HitNucPosition());
    if ( p4Nf.P() < kF ) {return 0.;}
  }

  // Scale cross section by the number of active nucleons
  int hit_nuc_pdg = target.HitNucPdg();
  // Number of active nucleons in the target
  int num_active = pdg::IsProton(hit_nuc_pdg) ? target.Z() : target.N();

  xsec *= num_active;

  // Do we need to rotate so that \vec{q} || \vec{z}?
  if (fDoqAlongZ) {
    std::vector<TLorentzVector> otherp4 {p4Ni, p4NiOnShell, p4Nf};
    genie::utils::Rotate_qvec_alongZ(probeP4, lepP4, otherp4);
    p4Ni = otherp4[0];
    p4NiOnShell = otherp4[1];
    p4Nf = otherp4[2];
  }

  // Compute form factors using Q2tilde (the effective Q2 value after
  // binding energy corrections)
  TLorentzVector qP4 = probeP4 - lepP4;
  TLorentzVector qTildeP4 = qP4;
  qTildeP4.SetE( qP4.E() - epsilon_B );

  double Q2 = -1. * qP4.M2();
  double Q2tilde = -1. * qTildeP4.M2();

  double Q2min = genie::controls::kMinQ2Limit; // CC/NC limit
  //std::cout << Q2min << "\n";
  if( interaction->ProcInfo().IsEM() ) Q2min = genie::utils::kinematics::electromagnetic::kMinQ2Limit; // EM limit
  // Make sure Q2 is physical
  if ( Q2 < Q2min ) {
    return 0.;
  }

  // Get the correct couplings and form factors model for the current
  // interaction
  double coupling_factor = 1.;
  const genie::ProcessInfo& proc_info = interaction->ProcInfo();
  if ( proc_info.IsWeakCC() ) {
    coupling_factor = kGF2 * fCos8c2;
    fFormFactors.SetModel( fCCFormFactorsModel );
  }

  else if ( proc_info.IsWeakNC() ) {
    coupling_factor = kGF2 * P_lep;
    fFormFactors.SetModel( fNCFormFactorsModel );
  }
  else if ( proc_info.IsEM() ) {
    coupling_factor = 32. * kPi * kPi * kAem2 / ( Q2*Q2 ) ;
    fFormFactors.SetModel( fEMFormFactorsModel );
  }
  else {
    LOG("UnifiedQELPXSec", pERROR) << "Unrecognized process type encountered"
      << " in genie::UnifiedQELPXSec::XSec()";
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

  double f1v = fFormFactors.F1V();
  double f2v = fFormFactors.xiF2V();
  double ffa = fFormFactors.FA();
  double ffp = fFormFactors.Fp();

  // Now that we've calculated them, store the true Q2 value
  interaction->KinePtr()->SetQ2( Q2 );

  // Compute the leptonic tensor
  // Note we have to pass a bool to LeptonTensor to use the SF conventions
  InteractionType_t type = interaction->ProcInfo().InteractionTypeId();
  LeptonTensor L_munu( probeP4, lepP4, init_state.ProbePdg(), type, true);

  // Set up the 4x4 hadronic response tensor
  std::complex<double> HadronTensor[4][4]; 

  // For CC use mean of proton/neutron mass
  // works for NC and EM as well 
  double xmn = ( mNi + interaction->RecoilNucleon()->Mass() ) / 2.;

  // Get energy and momentum transfer values
  double w = qP4.E();
  double wt = qTildeP4.E(); 
  double pNi_x, pNi_y, pNi_z, q_x, q_y, q_z;
  
  pNi_x = p4NiOnShell.X();
  pNi_y = p4NiOnShell.Y();
  pNi_z = p4NiOnShell.Z();
  q_x = qP4.X();
  q_y = qP4.Y();
  q_z = qP4.Z();

  // Compute hadron tensor 
  Get_hadrontensor_from_fortran(fFortranTensorModel, xmn, w, wt, pNi_x, pNi_y, pNi_z, q_x, q_y, q_z, f1v, f2v, ffa, ffp, HadronTensor);
  
  // Convert to a GENIE Rank2LorentzTensor object
  ManualResponseTensor ATilde_munu(HadronTensor);
      
  // Contract hadron and lepton tensors
  std::complex<double> contraction = L_munu * ATilde_munu;

  if ( std::abs(contraction.imag()) > kASmallNum ) {
    LOG("UnifiedQELPXSec", pWARN) << "Tensor contraction has nonvanishing imaginary part!";
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

  // Calculation is only appropriate for complex nuclear targets,
  // not free nucleons.
  if ( !init_state.Tgt().IsNucleus() ) return false;

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

  GetParamDef("FortranTensorModel", fFortranTensorModel, std::string("compute_hadron_tensor_SF"));  

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

  // Decide whether or not it should be used in XSec()
  GetParamDef( "DoPauliBlocking", fDoPauliBlocking, true );
  
  // Decide whether or not we rotate so that q is along z
  GetParamDef( "DoRotate_qAlong_z", fDoqAlongZ, false );

}

