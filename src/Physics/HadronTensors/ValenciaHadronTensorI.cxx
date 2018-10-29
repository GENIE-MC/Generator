// standard library includes
#include <cmath>

// GENIE includes
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/HadronTensors/ValenciaHadronTensorI.h"

double genie::ValenciaHadronTensorI::contraction(
  const Interaction* interaction, double Q_value) const
{
  // Don't do anything if you've been handed a nullptr
  if ( !interaction ) return 0.;

  int probe_pdg  = interaction->InitState().ProbePdg();
  double E_probe = interaction->InitState().ProbeE(kRfLab);
  double m_probe = interaction->InitState().Probe()->Mass();
  double Tl      = interaction->Kine().GetKV(kKVTl);
  double cos_l   = interaction->Kine().GetKV(kKVctl);
  double ml      = interaction->FSPrimLepton()->Mass();

  return contraction(probe_pdg, E_probe, m_probe, Tl, cos_l, ml, Q_value);
}

double genie::ValenciaHadronTensorI::contraction(int probe_pdg, double E_probe,
  double m_probe, double Tl, double cos_l, double ml, double Q_value) const
{
  // Final state lepton total energy
  double El = Tl + ml;

  // Energy transfer (uncorrected)
  double q0 = E_probe - El;

  // The corrected energy transfer takes into account the binding
  // energy of the struck nucleon(s)
  double q0_corrected = q0 - Q_value;

  // Magnitude of the initial state lepton 3-momentum
  double k_initial = std::sqrt( std::max(0., std::pow(E_probe, 2)
    - std::pow(m_probe, 2) ));

  // Magnitude of the final state lepton 3-momentum
  double k_final = std::sqrt( std::max(0., std::pow(Tl, 2) + 2*ml*Tl) );

  // Square of the magnitude of the 3-momentum transfer
  double q_mag2 = std::pow(k_initial, 2) + std::pow(k_final, 2)
    - 2.*k_initial*k_final*cos_l;

  // Square of the magnitude of the 4-momentum transfer
  double q2 = q0_corrected*q0_corrected - q_mag2;

  // Differential cross section
  double xs = dSigma_dT_dCosTheta(probe_pdg, E_probe, m_probe, Tl, cos_l, ml,
    Q_value);

  if ( pdg::IsNeutrino(probe_pdg) || pdg::IsAntiNeutrino(probe_pdg) ) {
    /// See equation (3) of Nieves, Amaro, and Valverde, Phys. Rev. C 70,
    /// 055503 (2004)
    xs *= 2. * genie::constants::kPi * k_initial
      / (k_final * genie::constants::kGF2);
  }
  else if ( probe_pdg == kPdgElectron ) {
    /// See equation (82) of Gil, Nieves, and Oset, Nucl. Phys. A 627, 543
    /// (1997)
    xs *= std::pow(q2 / genie::constants::kAem, 2) * k_initial / k_final
      / (2. * genie::constants::kPi);
  }
  else {
    LOG("TabulatedHadronTensor", pERROR)
      << "Unknown probe PDG code " << probe_pdg
      << "encountered.";
    return 0.;
  }

  return xs;
}
