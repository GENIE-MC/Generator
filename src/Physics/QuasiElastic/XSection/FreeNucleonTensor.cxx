#include "Physics/QuasiElastic/XSection/FreeNucleonTensor.h"
#include "Physics/QuasiElastic/XSection/LeptonTensor.h"
#include "Physics/QuasiElastic/XSection/QELFormFactorsModelI.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGUtils.h"

namespace {

  // Helper function that chooses the correct value of the nucleon mass
  // given an input interaction type
  double GetNucleonMass( genie::InteractionType_t type, int hit_nuc_pdg ) {

    // This variable stores the value of the nucleon mass that should be
    // used when computing the free nucleon tensor.
    double mN = 0.;

    genie::PDGLibrary* pdg_lib = genie::PDGLibrary::Instance();
    double neutron_mass = pdg_lib->Find( genie::kPdgNeutron )->Mass();
    double proton_mass = pdg_lib->Find( genie::kPdgNeutron )->Mass();

    if ( type == genie::kIntWeakCC ) {
      // For charged current interactions, the (on-shell) "nucleon mass" used in
      // the matrix element of the nucleon current is the mean of the neutron and
      // proton masses
      mN = (neutron_mass + proton_mass) / 2.;
    }
    else {
      // For kIntEM or kIntWeakNC, use the (on-shell) mass of the initial
      // struck nucleon.
      assert( genie::pdg::IsNucleon(hit_nuc_pdg) );
      if ( hit_nuc_pdg == genie::kPdgNeutron ) mN = neutron_mass;
      else mN = proton_mass;
    }

    return mN;
  }

}

genie::FreeNucleonTensor::FreeNucleonTensor(
  const genie::Interaction& interaction,
  const genie::QELFormFactors& ff)
{
  // Set nucleon mass appropriately for the interaction type
  genie::InteractionType_t type = interaction.ProcInfo().InteractionTypeId();
  int hit_nuc_pdg = interaction.InitState().Tgt().HitNucPdg();
  fNucleonMass = GetNucleonMass( type, hit_nuc_pdg );

  // Compute form factors
  fF1V = ff.F1V();
  fF2V = ff.xiF2V();
  fF3V = ff.F3V();
   fFA = ff.FA();
   fFP = ff.Fp();
  fF3A = ff.F3A();

  // On-shell initial hit nucleon mass
  double mNi = interaction.InitState().Tgt().HitNucMass();
  // Initial hit nucleon total energy (possibly off-shell)
  double ENi = interaction.InitState().Tgt().HitNucP4().E();
  TVector3 p3Ni = interaction.InitState().Tgt().HitNucP4().Vect();
  double ENi_on_shell = std::sqrt( p3Ni.Mag2() + mNi*mNi );

  // TODO: document. See eq. 3 in Phys. Rev. C 77, 044311 (2008).
  double epsilon_B = ENi_on_shell - ENi;

  TLorentzVector* temp_probeP4 = interaction.InitState().GetProbeP4( kRfLab );
  TLorentzVector probeP4 = *temp_probeP4;
  delete temp_probeP4;

  // Get fqTilde (lab frame). Assume that SetFSLeptonP4() has already been
  // called to store the outgoing lepton's 4-momentum
  const TLorentzVector& p4lep = interaction.Kine().FSLeptonP4();
  fqTilde = probeP4 - p4lep;
  fqTilde.SetE( fqTilde.E() - epsilon_B );

  // Initialize the on-shell initial nucleon 4-momentum
  fPNiOnShell = TLorentzVector( p3Ni, ENi_on_shell );
}

genie::FreeNucleonTensor::FreeNucleonTensor(const TLorentzVector& pNi,
  const TLorentzVector& q, double F1V, double F2V, double F3V, double FA,
  double FP, double F3A, genie::InteractionType_t type, int hit_nuc_pdg)
{
  fNucleonMass = GetNucleonMass( type, hit_nuc_pdg );
  fPNiOnShell = pNi;
  fqTilde = q;

  fF1V = F1V;
  fF2V = F2V;
  fF3V = F3V;
  fFA = FA;
  fFP = FP;
  fF3A = F3A;
}

std::complex<double> genie::FreeNucleonTensor::operator()(
  genie::TensorIndex_t mu, genie::TensorIndex_t nu) const
{
  double q2 = fqTilde.M2();
  double M2 = std::pow(fNucleonMass, 2);
  double g_mu_nu = Metric(mu,nu);
  std::complex<double> imaginary_i(0., 1.);

  double p_mu = fPNiOnShell[mu];
  double p_nu = fPNiOnShell[nu];
  double q_mu = fqTilde[mu];
  double q_nu = fqTilde[nu];

  std::complex<double> result(0., 0.);

  // F1V squared term
  result += std::pow(fF1V, 2) * (p_mu*q_nu + q_mu*p_nu + 2.*p_mu*p_nu
    + 0.5*q2*g_mu_nu);

  // F2V squared term
  result -= std::pow(fF2V, 2) * ( -0.5*q2*g_mu_nu
    + 0.5*(1. + 0.25*q2/M2)*q_mu*q_nu
    + (0.25*q2/M2)*(p_mu*q_nu + q_mu*p_nu + 2.*p_mu*p_nu) );

  // FA squared term
  result += std::pow(fFA, 2) * ( 2*M2*g_mu_nu * (0.25*q2/M2 - 1.)
    + q_mu*p_nu + p_mu*q_nu + 2*p_mu*p_nu );

  // FP squared term
  result -= std::pow(fFP, 2) * 0.5*q2*q_mu*q_nu / M2;

  // F3V squared term
  result += 2. * std::pow(fF3V, 2) * (1. - 0.25*q2/M2) * q_mu
    * q_nu;

  // F3A squared term
  result -= std::pow(fF3A, 2) * 0.5*q2/M2 * (2.*p_mu + q_mu)
    * (2.*p_nu + q_nu);

  // F1V * F2V term
  result += fF1V * fF2V * (q2 * g_mu_nu - q_mu * q_nu);

  // (F1V + F2V) * FA term
  result -= 2.*imaginary_i * (fF1V + fF2V) * fFA
    * this->LeviCivitaProduct(mu, nu, fPNiOnShell, fqTilde);

  // FA * FP term
  result -= 2. * fFA * fFP * q_mu * q_nu;

  // Remaining terms (involving second-class currents)
  result += (0.5*fF2V*fF3V/M2 - q2*fF3A*fFP/M2 - 2.*fF3A*fFA
    + 2.*fF1V*fF3V) * (p_mu*q_nu + q_mu*p_nu + q_mu*q_nu);

  // ** Overall factor **
  // Note that the expression I used to write this code (equation (A1) in the appendix of
  // https://arxiv.org/pdf/1506.02355.pdf) is off by a factor of two from the definition of J^{\mu\nu}
  // (see equation equation 9 therein) given earlier in the paper. We're using FreeNucleonTensor = J^{\mu\nu}
  // as defined in equation 9, so the overall factor that should be used is 2, not 4. You can verify this
  // yourself, e.g., try computing the Fp^2 term from the Dirac traces.
  result *= 2.;

  return result;
}
