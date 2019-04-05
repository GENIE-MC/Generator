//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Steven Gardiner <gardiner \at fnal.gov>
         Fermi National Accelerator Laboratory
*/
//____________________________________________________________________________

// standard library includes
#include <cmath>
#include <limits>
#include <sstream>

// ROOT includes
#include "TDatabasePDG.h"
#include "TLorentzVector.h"

// GENIE includes
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/KinePhaseSpace.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/Interaction/InitialState.h"
#include "Framework/Interaction/Kinematics.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Physics/QuasiElastic/XSection/QELUtils.h"


// Helper functions local to this source file
namespace {

  TVector3 COMframe2Lab(const genie::InitialState& initialState)
  {
    TLorentzVector* k4 = initialState.GetProbeP4( genie::kRfLab );
    TLorentzVector* p4 = initialState.TgtPtr()->HitNucP4Ptr();
    TLorentzVector totMom = *k4 + *p4;

    TVector3 beta = totMom.BoostVector();

    delete k4;

    return beta;
  }

}

double genie::utils::EnergyDeltaFunctionSolutionQEL(
  const genie::Interaction& inter)
{
  // Get the final-state lepton and struck nucleon 4-momenta in the lab frame
  const TLorentzVector& lepton = inter.Kine().FSLeptonP4();
  const TLorentzVector& outNucleon = inter.Kine().HadSystP4();

  // Get beta for a Lorentz boost from the lab frame to the COM frame and
  // vice-versa
  TLorentzVector* probe = inter.InitStatePtr()->GetProbeP4( kRfLab );
  const TLorentzVector& hit_nucleon = inter.InitStatePtr()->TgtPtr()
    ->HitNucP4();
  TLorentzVector total_p4 = (*probe) + hit_nucleon;
  TVector3 beta_COM_to_lab = total_p4.BoostVector();
  TVector3 beta_lab_to_COM = -beta_COM_to_lab;

  // Square of the Lorentz gamma factor for the boosts
  double gamma2 = std::pow(total_p4.Gamma(), 2);

  // Get the final-state lepton 4-momentum in the COM frame
  TLorentzVector leptonCOM = TLorentzVector( lepton );
  leptonCOM.Boost( beta_lab_to_COM );

  // Angle between COM outgoing lepton and lab-frame velocity of COM frame
  double theta0 = leptonCOM.Angle( beta_COM_to_lab );
  double cos_theta0_squared = std::pow(std::cos( theta0 ), 2);
  //double sin_theta0_squared = 1. - cos_theta0_squared; // sin^2(x) + cos^2(x) = 1

  // Difference in lab frame velocities between lepton and outgoing nucleon
  TVector3 lepton_vel = ( lepton.Vect() * (1. / lepton.E()) );
  TVector3 outNucleon_vel = ( outNucleon.Vect() * (1. / outNucleon.E()) );
  double vel_diff = ( lepton_vel - outNucleon_vel ).Mag();

  // Compute the factor in the kPSQELEvGen phase space that arises from
  // eliminating the energy-conserving delta function
  double factor = std::sqrt( 1. + (1. - cos_theta0_squared)*(gamma2 - 1.) )
    * std::pow(leptonCOM.P(), 2) / vel_diff;

  delete probe;

  return factor;
}

double genie::utils::ComputeFullQELPXSec(genie::Interaction* interaction,
  const genie::NuclearModelI* nucl_model, const genie::XSecAlgorithmI* xsec_model,
  double cos_theta_0, double phi_0, double& Eb,
  genie::QELEvGen_BindingMode_t hitNucleonBindingMode, double min_angle_EM,
  bool bind_nucleon)
{
  // If requested, set the initial hit nucleon 4-momentum to be off-shell
  // according to the binding mode specified in the function call
  if ( bind_nucleon ) {
    genie::utils::BindHitNucleon(*interaction, *nucl_model, Eb,
      hitNucleonBindingMode);
  }

  // Mass of the outgoing lepton
  double lepMass = interaction->FSPrimLepton()->Mass();

  // Look up the (on-shell) mass of the final nucleon
  TDatabasePDG *tb = TDatabasePDG::Instance();
  double mNf = tb->GetParticle( interaction->RecoilNucleonPdg() )->Mass();

  // Mandelstam s for the probe/hit nucleon system
  double s = std::pow( interaction->InitState().CMEnergy(), 2 );

  // Return a differential cross section of zero if we're below threshold (and
  // therefore need to sample a new event)
  if ( std::sqrt(s) < lepMass + mNf ) return 0.;

  double outLeptonEnergy = ( s - mNf*mNf + lepMass*lepMass ) / (2 * std::sqrt(s));

  if (outLeptonEnergy*outLeptonEnergy - lepMass*lepMass < 0.) return 0.;
  double outMomentum = TMath::Sqrt(outLeptonEnergy*outLeptonEnergy - lepMass*lepMass);

  // Compute the boost vector for moving from the COM frame to the
  // lab frame, i.e., the velocity of the COM frame as measured
  // in the lab frame.
  TVector3 beta = COMframe2Lab( interaction->InitState() );

  // FullDifferentialXSec depends on theta_0 and phi_0, the lepton COM
  // frame angles with respect to the direction of the COM frame velocity
  // as measured in the lab frame. To generate the correct dependence
  // here, first set the lepton COM frame angles with respect to +z
  // (via TVector3::SetTheta() and TVector3::SetPhi()).
  TVector3 lepton3Mom(0., 0., outMomentum);
  lepton3Mom.SetTheta( TMath::ACos(cos_theta_0) );
  lepton3Mom.SetPhi( phi_0 );

  // Then rotate the lepton 3-momentum so that the old +z direction now
  // points along the COM frame velocity (beta)
  TVector3 zvec(0., 0., 1.);
  TVector3 rot = ( zvec.Cross(beta) ).Unit();
  double angle = beta.Angle( zvec );

  // Handle the edge case where beta is along -z, so the
  // cross product above vanishes
  if ( beta.Perp() == 0. && beta.Z() < 0. ) {
    rot = TVector3(0., 1., 0.);
    angle = genie::constants::kPi;
  }

  // Rotate if the rotation vector is not 0
  if ( rot.Mag() >= genie::controls::kASmallNum ) {
    lepton3Mom.Rotate(angle, rot);
  }

  // Construct the lepton 4-momentum in the COM frame
  TLorentzVector lepton(lepton3Mom, outLeptonEnergy);

  // The final state nucleon will have an equal and opposite 3-momentum
  // in the COM frame and will be on the mass shell
  TLorentzVector outNucleon(-1*lepton.Px(),-1*lepton.Py(),-1*lepton.Pz(),
    TMath::Sqrt(outMomentum*outMomentum + mNf*mNf));

  // Boost the 4-momenta for both particles into the lab frame
  lepton.Boost(beta);
  outNucleon.Boost(beta);

  // Check if event is at a low angle - if so return 0 and stop wasting time
  if (180 * lepton.Theta() / genie::constants::kPi < min_angle_EM && interaction->ProcInfo().IsEM()) {
    return 0;
  }

  TLorentzVector * nuP4 = interaction->InitState().GetProbeP4( genie::kRfLab );
  TLorentzVector qP4 = *nuP4 - lepton;
  delete nuP4;
  double Q2 = -1 * qP4.Mag2();

  interaction->KinePtr()->SetFSLeptonP4( lepton );
  interaction->KinePtr()->SetHadSystP4( outNucleon );
  interaction->KinePtr()->SetQ2( Q2 );

  // Check the Q2 range. If we're outside of it, don't bother
  // with the rest of the calculation.
  Range1D_t Q2lim = interaction->PhaseSpace().Q2Lim();
  if (Q2 < Q2lim.min || Q2 > Q2lim.max) return 0.;

  // Compute the QE cross section for the current kinematics
  double xsec = xsec_model->XSec(interaction, genie::kPSQELEvGen);

  return xsec;
}

genie::QELEvGen_BindingMode_t genie::utils::StringToQELBindingMode(
  const std::string& binding_mode)
{
  // Translate the string setting the binding mode to the appropriate
  // enum value, or complain if one couldn't be found
  if ( binding_mode == "UseGroundStateRemnant" ) {
    return kUseGroundStateRemnant;
  }
  else if ( binding_mode == "UseNuclearModel" ) {
    return kUseNuclearModel;
  }
  else if ( binding_mode == "OnShell" ) {
    return kOnShell;
  }
  else {
    LOG("QELEvent", pFATAL) << "Unrecognized setting \"" << binding_mode
      << "\" requested in genie::utils::StringToQELBindingMode()";
    gAbortingInErr = true;
    std::exit(1);
  }
}

double genie::utils::CosTheta0Max(const genie::Interaction& interaction) {

  // q0 > 0 only needs to be enforced (indirectly via a calculation of
  // CosTheta0Max) for bound nucleons. The Q2 limits should take care of valid
  // kinematics for free nucleons.
  if ( !interaction.InitState().Tgt().IsNucleus()
    || interaction.TestBit(kIAssumeFreeNucleon) ) return 1.;

  double probe_E_lab = interaction.InitState().ProbeE( genie::kRfLab );

  TVector3 beta = COMframe2Lab( interaction.InitState() );
  double gamma = 1. / std::sqrt(1. - beta.Mag2());

  double sqrt_s = interaction.InitState().CMEnergy();
  double mNf = interaction.RecoilNucleon()->Mass();
  double ml = interaction.FSPrimLepton()->Mass();
  double lepton_E_COM = (sqrt_s*sqrt_s + ml*ml - mNf*mNf) / (2.*sqrt_s);

  // If there isn't enough available energy to create an on-shell
  // final lepton, then don't bother with the rest of the calculation
  // NOTE: C++11 would allow use to use lowest() here instead
  if ( lepton_E_COM <= ml ) return -std::numeric_limits<double>::max();

  double lepton_p_COM = std::sqrt( std::max(lepton_E_COM*lepton_E_COM
    - ml*ml, 0.) );

  // Possibly off-shell initial struck nucleon total energy
  // (BindHitNucleon() should have been called previously if needed)
  const TLorentzVector& p4Ni = interaction.InitState().Tgt().HitNucP4();
  double ENi = p4Ni.E();
  // On-shell mass of initial struck nucleon
  double mNi = interaction.InitState().Tgt().HitNucMass();
  // On-shell initial struck nucleon energy
  double ENi_on_shell = std::sqrt( mNi*mNi + p4Ni.Vect().Mag2() );
  // Energy needed to put initial nucleon on the mass shell
  double epsilon_B = ENi_on_shell - ENi;

  double cos_theta0_max = ( probe_E_lab - gamma*lepton_E_COM - epsilon_B )
    / ( gamma * lepton_p_COM * beta.Mag() );
  return cos_theta0_max;
}

void genie::utils::BindHitNucleon(genie::Interaction& interaction,
  const genie::NuclearModelI& nucl_model, double& Eb,
  genie::QELEvGen_BindingMode_t hitNucleonBindingMode)
{
  genie::Target* tgt = interaction.InitState().TgtPtr();
  TLorentzVector* p4Ni = tgt->HitNucP4Ptr();

  // Initial nucleon 3-momentum (lab frame)
  TVector3 p3Ni = nucl_model.Momentum3();

  // Look up the (on-shell) mass of the initial nucleon
  TDatabasePDG* tb = TDatabasePDG::Instance();
  double mNi = tb->GetParticle( tgt->HitNucPdg() )->Mass();

  // Set the (possibly off-shell) initial nucleon energy based on
  // the selected binding energy mode. Always put the initial nucleon
  // on shell if it is not part of a composite nucleus
  double ENi = 0.;
  if ( tgt->IsNucleus() && hitNucleonBindingMode != genie::kOnShell ) {

    // For a nuclear target with a bound initial struck nucleon, take binding
    // energy effects and Pauli blocking into account when computing QE
    // differential cross sections
    interaction.ResetBit( kIAssumeFreeNucleon );

    // Initial nucleus mass
    double Mi = tgt->Mass();

    // Final nucleus mass
    double Mf = 0.;

    // If we use the removal energy reported by the nuclear
    // model, then it implies a certain value for the final
    // nucleus mass
    if ( hitNucleonBindingMode == genie::kUseNuclearModel ) {
      Eb = nucl_model.RemovalEnergy();
      // This equation is the definition that we assume
      // here for the "removal energy" (Eb) returned by the
      // nuclear model. It matches GENIE's convention for
      // the Bodek/Ritchie Fermi gas model.
      Mf = Mi + Eb - mNi;
    }
    // We can also assume that the final nucleus is in its
    // ground state. In this case, we can just look up its
    // mass from the standard table. This implies a particular
    // binding energy for the hit nucleon.
    else if ( hitNucleonBindingMode == genie::kUseGroundStateRemnant ) {
      // Determine the mass and proton numbers for the remnant nucleus
      int Af = tgt->A() - 1;
      int Zf = tgt->Z();
      if ( genie::pdg::IsProton( tgt->HitNucPdg()) ) --Zf;
      Mf = genie::PDGLibrary::Instance()->Find( genie::pdg::IonPdgCode(Af, Zf) )->Mass();

      // Deduce the binding energy from the final nucleus mass
      Eb = Mf - Mi + mNi;
    }

    // The (lab-frame) off-shell initial nucleon energy is the difference
    // between the lab frame total energies of the initial and remnant nuclei
    ENi = Mi - std::sqrt( Mf*Mf + p3Ni.Mag2() );
  }
  else {
    // Keep the struck nucleon on shell either because
    // hitNucleonBindingMode == kOnShell or because
    // the target is a single nucleon
    ENi = std::sqrt( p3Ni.Mag2() + std::pow(mNi, 2) );
    Eb = 0.;

    // If we're dealing with a nuclear target but using the on-shell
    // binding mode, set the "assume free nucleon" flag. This turns off
    // Pauli blocking and the requirement that q0 > 0 in the QE cross section
    // models (an on-shell nucleon *is* a free nucleon)
    if ( tgt->IsNucleus() ) interaction.SetBit( kIAssumeFreeNucleon );
  }

  // Update the initial nucleon lab-frame 4-momentum in the interaction with
  // its current components
  p4Ni->SetVect( p3Ni );
  p4Ni->SetE( ENi );

}
