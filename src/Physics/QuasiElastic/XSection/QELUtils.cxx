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

// ROOT includes
#include "TLorentzVector.h"

// GENIE includes
#include "Framework/Interaction/InitialState.h"
#include "Framework/Interaction/Kinematics.h"
#include "Physics/QuasiElastic/XSection/QELUtils.h"

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
