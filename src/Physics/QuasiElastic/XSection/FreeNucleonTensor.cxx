#include "Physics/QuasiElastic/XSection/FreeNucleonTensor.h"
#include "Physics/QuasiElastic/XSection/LeptonTensor.h"
#include "Physics/QuasiElastic/XSection/QELFormFactorsModelI.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGLibrary.h"

namespace {

  double curly_Q(const TLorentzVector& probeP4, const TLorentzVector& qTildeP4)
  {
    TVector3 qP3 = qTildeP4.Vect();
    double k_dot_q = probeP4.Vect().Dot( qP3 );
    double qMag = qP3.Mag();
    double Q = ( 2. * k_dot_q / qMag ) - qMag;
    return Q;
  }

  double curly_P(const TLorentzVector& hitNucleonOnShellP4,
    const TLorentzVector& qTildeP4, double epsilon_B)
  {
    double qMag = qTildeP4.P();
    double omega = qTildeP4.E() + epsilon_B;
    double P = ( 2.*hitNucleonOnShellP4.E() + omega ) / (2. * qMag);
    return P;
  }

  double C1(const TLorentzVector& probeP4, const TLorentzVector qTildeP4,
    double epsilon_B, double& Q)
  {
    double omegaTilde = qTildeP4.E();
    double omega = omegaTilde + epsilon_B;
    Q = curly_Q(probeP4, qTildeP4);
    double k_dot_q = probeP4.Vect().Dot( qTildeP4.Vect() );
    double coeff1 = (omega - omegaTilde)
      * ( (Q*Q - omega*omega)*(omega + omegaTilde)
      - 4.*(probeP4.P()*qTildeP4.P()*Q - omega*k_dot_q) );
    return coeff1;
  }

  double C2(const TLorentzVector& probeP4, const TLorentzVector& qTildeP4,
    const TLorentzVector& hitNucleonOnShellP4, double epsilon_B,
    double& coeff1)
  {
    double Q = 0.;
    coeff1 = C1(probeP4, qTildeP4, epsilon_B, Q);
    double P = curly_P(hitNucleonOnShellP4, qTildeP4, epsilon_B);
    double omegaTilde = qTildeP4.E();
    double omega = omegaTilde + epsilon_B;
    double k_dot_p = probeP4.Vect().Dot( hitNucleonOnShellP4.Vect() );
    double p_dot_q = hitNucleonOnShellP4.Vect().Dot( qTildeP4.Vect() );
    double k_dot_q = probeP4.Vect().Dot( qTildeP4.Vect() );
    double coeff2 = coeff1 * std::pow(P, 2);
    coeff2 += 4. * (omega - omegaTilde) * P * Q * ( k_dot_p
      - p_dot_q*k_dot_q / qTildeP4.Vect().Mag2() );
    return coeff2;
  }

}

genie::FreeNucleonTensor::FreeNucleonTensor(
  const genie::Interaction& interaction,
  const genie::QELFormFactors& ff)
{
  fInteractionType = interaction.ProcInfo().InteractionTypeId();

  // Get the (on-shell) mass of the initial hit nucleon
  double mNi = interaction.InitState().Tgt().HitNucMass();

  if ( fInteractionType == kIntWeakCC ) {
    // For charged current interactions, the (on-shell) "nucleon mass" used in
    // the matrix element of the nucleon current is the mean of the neutron and
    // proton masses
    genie::PDGLibrary* pdg_lib = genie::PDGLibrary::Instance();
    double neutron_mass = pdg_lib->Find( genie::kPdgNeutron )->Mass();
    double proton_mass = pdg_lib->Find( genie::kPdgNeutron )->Mass();
    fNucleonMass = (neutron_mass + proton_mass) / 2.;
  }
  else {
    // For kIntEM or kIntWeakNC, use the mass of the initial struck
    // nucleon. Note that Target::HitNucMass() gets its value from the
    // PDGLibrary, so this will correctly be the on-shell mass.
    fNucleonMass = mNi;
  }

  // Compute form factors
  fF1V = ff.F1V();
  fF2V = ff.xiF2V();
  fF3V = ff.F3V();
   fFA = ff.FA();
   fFP = ff.Fp();
  fF3A = ff.F3A();

  double ENi = interaction.InitState().Tgt().HitNucP4().E();
  TVector3 p3Ni = interaction.InitState().Tgt().HitNucP4().Vect();
  double ENi_on_shell = std::sqrt( p3Ni.Mag2() + mNi*mNi );
  // TODO: document. See eq. 3 in Phys. Rev. C 77, 044311 (2008).
  fEpsilonB = ENi_on_shell - ENi;

  TLorentzVector* temp_probeP4 = interaction.InitState().GetProbeP4( kRfLab );
  TLorentzVector probeP4 = *temp_probeP4;
  delete temp_probeP4;

  // Get fqTilde (lab frame). Assume that SetFSLeptonP4() has already been
  // called to store the outgoing lepton's 4-momentum
  const TLorentzVector& p4lep = interaction.Kine().FSLeptonP4();
  fqTilde = probeP4 - p4lep;
  fqTilde.SetE( fqTilde.E() - fEpsilonB );

  // Initialize the on-shell initial nucleon 4-momentum
  fPNiOnShell = TLorentzVector( p3Ni, ENi_on_shell );
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

std::complex<double> genie::FreeNucleonTensor::Contract(
  const Rank2LorentzTensorI& other) const
{
  std::complex<double> result = genie::Rank2LorentzTensorI::Contract( other );

  // If we're dealing with the contraction of the EM nucleon current with the
  // leptonic current, then we need to add a couple of terms to correct for
  // current conservation
  // TODO: describe the trick you're using
  const genie::LeptonTensor* lep = dynamic_cast<const genie::LeptonTensor*>(
    &other );
  if ( lep && fInteractionType == kIntEM ) {

    const TLorentzVector& probeP4 = lep->GetProbeP4();

    double qTilde2 = fqTilde.M2();
    double tau = -qTilde2 / ( 4. * std::pow(fNucleonMass, 2) );

    double coeff1 = 0.;
    double coeff2 = C2(probeP4, fqTilde, fPNiOnShell, fEpsilonB, coeff1);

    double H1t = tau * std::pow(fF1V + fF2V, 2);
    double H2t = std::pow(fF1V, 2) + ( tau * std::pow(fF2V, 2) );

    double correction = (std::pow(fNucleonMass, 2) / qTilde2) * coeff1 * H1t;
    correction += coeff2 * H2t;

    // Note that the paper from which these current-conservation corrections
    // were taken, https://arxiv.org/abs/0711.2031, uses a different
    // normalization convention for the lepton and free nucleon tensors. To
    // convert to ours (that of https://arxiv.org/pdf/1506.02355.pdf), we need
    // to multiply by 16.
    correction *= 16.;

    result += correction;
  }

  return result;
}
