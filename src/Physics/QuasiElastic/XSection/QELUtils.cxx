//____________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Steven Gardiner <gardiner \at fnal.gov>
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
  genie::QELEvGen_BindingMode_t hitNucleonBindingMode, double fMinAngleEM, bool bind_nucleon)
{
  // If requested, set the initial hit nucleon 4-momentum to be off-shell
  // according to the binding mode specified in the function call
  if ( bind_nucleon ) {
    genie::utils::BindHitNucleon(*interaction, *nucl_model, Eb,
      hitNucleonBindingMode);
  }

  // A very high-momentum bound nucleon (which is far off the mass shell)
  // can have a momentum greater than its total energy. This leads to numerical
  // issues (NaNs) since the invariant mass of the nucleon becomes imaginary.
  // In such cases, just return zero to avoid trouble.
  if ( interaction->InitState().Tgt().HitNucP4().M() <= 0. ) return 0.;

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

  //LOG("QELEvent", pWARN) << "Outgoing nucleon momentum: Px = " << outNucleon.X() << ", Py = " << outNucleon.Y() << ", Pz = " << outNucleon.Z();

  TLorentzVector * nuP4 = interaction->InitState().GetProbeP4( genie::kRfLab );
  TLorentzVector qP4 = *nuP4 - lepton;
  delete nuP4;
  double Q2 = -1 * qP4.Mag2();

  // Check the Q2 range. If we're outside of it, don't bother
  // with the rest of the calculation.
  Range1D_t Q2lim = interaction->PhaseSpace().Q2Lim();
  if (Q2 < Q2lim.min || Q2 > Q2lim.max) return 0.;

  interaction->KinePtr()->SetFSLeptonP4( lepton );
  interaction->KinePtr()->SetHadSystP4( outNucleon );
  interaction->KinePtr()->SetQ2( Q2 );

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
      if ( nucl_model.ModelType(*tgt) != kNucmSpectralFunc ) {
        Eb = nucl_model.RemovalEnergy();
        // For all nuclear models except SpectralFunc, this equation is the
        // definition that we assume for the "removal energy" (Eb). It matches
        // GENIE's convention for the Bodek/Ritchie Fermi gas model.
        Mf = Mi + Eb - mNi;
      }
      else {
        // The SpectralFunc nuclear model returns a removal energy
        // which includes the kinetic energy of the final-state nucleus.
        // We account for this difference here.
        double E = nucl_model.RemovalEnergy();
        Mf = std::sqrt( std::max(0., std::pow(Mi + E - mNi, 2) - p3Ni.Mag2()) );
        Eb = Mf + mNi - Mi;
      }
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

// Assumes that the probe 4-momentum and struck nucleon PDG code are already
// set in the input interaction.
// TODO: adjust to avoid assuming that the probe is traveling along the +z
// direction
double genie::utils::ComputeTestQELPXSec(genie::Interaction* interaction,
  const genie::XSecAlgorithmI* xsec_model, double omega, double ctl,
  double phiLep, double phiNi, double pNi, double Eremove, KinePhaseSpace_t kps)
{
  // Get the probe 4-momentum
  TLorentzVector* tempProbeP4 = interaction->InitState().GetProbeP4(
    genie::kRfLab );
  TLorentzVector nuP4 = *tempProbeP4;
  delete tempProbeP4;

  // Prepare the outgoing lepton 4-momentum
  double ml = interaction->FSPrimLepton()->Mass();

  double Ev = nuP4.E();
  double El = Ev - omega;
  double pl = std::sqrt( std::max(0., El*El - ml*ml) );
  if ( El < ml ) return 0.;

  double stl = std::sqrt( std::max(0., 1. - ctl*ctl) );
  TLorentzVector lP4( std::cos(phiLep)*stl*pl,
    std::sin(phiLep)*stl*pl, ctl*pl, El );

  // Get the 4-momentum transfer
  TLorentzVector qP4 = nuP4 - lP4;
  
  // 3-momentum transfer
  double qMag = qP4.Vect().Mag();
  // Get the cosine between the initial nucleon 3-momentum and the 3-momentum
  // transfer

  // On-shell mass of initial struck nucleon
  double mNi = interaction->InitState().Tgt().HitNucMass();
  // Off-shell total energy of initial struck nucleon
  //double ENi = mNi - Eremove;
  double ENi = std::sqrt(pNi*pNi + mNi*mNi); 
  // Look up the (on-shell) mass of the final nucleon
  double mNf = interaction->RecoilNucleon()->Mass();
  
  //Get omega_tilde
  double wt = omega - ENi - std::abs(Eremove) + mNi;
  if (wt < 0) return 0.;
  //std::cout << "wt from QEL: " << wt << std::endl; 
  //Get angle between initial and final nucleon momenta
  double costheta_p = (pow(wt + ENi,2) - pNi*pNi - qMag*qMag - mNf*mNf)/(2.*pNi*qMag);
  if (std::abs(costheta_p) > 1.0) return 0.;

  // Final nucleon 3-momentum
  double pNf = std::sqrt(pNi*pNi + qMag*qMag + 2.*qMag*pNi*costheta_p);
  // Final nucleon total energy
  double ENf = std::sqrt(pNf*pNf + mNf*mNf);
  
 // std::cout << "ENf = " << ENf << "pNf = " << pNf << std::endl;
  double pCos = ( pNf*pNf - pNi*pNi - qMag*qMag ) / ( 2.*pNi*qMag );
  if ( std::abs(pCos) > 1. ) return 0.;

  double pSin = std::sqrt( std::max(0., 1. - pCos*pCos) );

  // Set up the initial struck nucleon 3-momentum in a frame where the
  // 3-momentum transfer points along +z
  TVector3 p3Ni( std::cos(phiNi)*pSin*pNi, std::sin(phiNi)*pSin*pNi,
    pCos*pNi );

  // Build the initial struck nucleon 4-momentum
  TLorentzVector p4Ni( p3Ni, ENi );

  // Use it and the 4-momentum transfer to get the final nucleon 4-momentum
  TVector3 p3Nf( std::cos(phiNi)*pSin*pNi, std::sin(phiNi)*pSin*pNi,
    pCos*pNi + qMag );
  TLorentzVector p4Nf(p3Nf, ENf);

  // Set the needed 4-momenta in the input interaction
  interaction->InitState().TgtPtr()->SetHitNucP4( p4Ni );
  interaction->KinePtr()->SetFSLeptonP4( lP4 );
  interaction->KinePtr()->SetHadSystP4( p4Nf );

  // A very high-momentum bound nucleon (which is far off the mass shell)
  // can have a momentum greater than its total energy. This leads to numerical
  // issues (NaNs) since the invariant mass of the nucleon becomes imaginary.
  // In such cases, just return zero to avoid trouble.
  if ( interaction->InitState().Tgt().HitNucP4().M() <= 0. ) return 0.;

  // Mandelstam s for the probe/hit nucleon system
  double s = std::pow( interaction->InitState().CMEnergy(), 2 );

  // Return a differential cross section of zero if we're below threshold
  if ( std::sqrt(s) < ml + mNf ) return 0.;

  // Check that Q^2 is within the allowed range
  double Q2 = -1 * qP4.Mag2();

  //std::cout << "here1? " << std::endl;
  // Check the Q2 range. If we're outside of it, don't bother
  // with the rest of the calculation.
  Range1D_t Q2lim = interaction->PhaseSpace().Q2Lim();
  if ( Q2 < Q2lim.min || Q2 > Q2lim.max ) return 0.;
  //std::cout << "here2? " << std::endl;

  // Set Q^2 in the input interaction
  interaction->KinePtr()->SetQ2( Q2 );

  // Compute the QE cross section for the current kinematics
  double xsec = xsec_model->XSec(interaction, kps);
  //std::cout << "KPS = " << kps << std::endl;
  if (kps == genie::kPSQELEvGen) {
	//std::cout << "Using kPSQELEvGen" << std::endl;
 	 // Adjust some factors for the testing phase space (which differs from
  	// kPSQELEvGen)
  	double delta_sol = EnergyDeltaFunctionSolutionQEL( *interaction );

  	double factor = lP4.E() * lP4.P() * p4Nf.E() * p4Ni.P()
    	/ ( qMag * delta_sol );

  	xsec *= factor;
  }

  return xsec;
}

// Function takes ingoing and outgoing lepton and
// std::vector<> of other 4 momenta
// and rotates the system so that \vec{q} is along \hat{z}
void genie::utils::Rotate_qvec_alongZ(TLorentzVector &probe_leptonP4, 
TLorentzVector &out_leptonP4, std::vector<TLorentzVector> &otherP4) 
{

  TVector3 probeLMom3 = probe_leptonP4.Vect();
  TVector3 outLMom3 = out_leptonP4.Vect();
  
  std::vector<TVector3> otherMom3;
 
  for (auto vect4: otherP4 ) {
	TVector3 vectorPart(vect4.X(), vect4.Y(), vect4.Z());
	otherMom3.push_back(vectorPart);
  }

  TVector3 q3Vec = probeLMom3 - outLMom3;
  TVector3 zvec(0.0, 0.0, 1.0);
  TVector3 rot = ( q3Vec.Cross(zvec) ).Unit(); // Vector to rotate about
  // Angle between the z direction and q
  double angle = zvec.Angle( q3Vec );

  // Handle the edge case where q3Vec is along -z, so the
  // cross product above vanishes
  if ( q3Vec.Perp() == 0. && q3Vec.Z() < 0. ) {
    rot = TVector3(0., 1., 0.);
    angle = genie::constants::kPi;
  }

  // Rotate if the rotation vector is not 0
  if ( rot.Mag() >= genie::controls::kASmallNum ) {

    probeLMom3.Rotate(angle,rot);
    probe_leptonP4.SetVect(probeLMom3);

    outLMom3.Rotate(angle,rot);
    out_leptonP4.SetVect(outLMom3);

    for(int i = 0; i < otherP4.size(); i++) {
    	otherMom3[i].Rotate(angle,rot);
        otherP4[i].SetVect(otherMom3[i]);
    }

  }

} 









