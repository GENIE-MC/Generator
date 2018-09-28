// standard library includes
#include <cmath>
#include <fstream>

// GENIE includes
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/Multinucleon/XSection/TabulatedHadronTensor.h"

namespace {
  /// Enumerated type that represents the format used to read in
  /// grid points from the hadron tensor data file
  enum kHadronTensorGridFlag_t {
    /// A starting value and a step size are used to define a regular grid
    kStartAndStep = 0,
    /// An explicit table of grid points is used to define an irregular grid
    kExplicitValues = 1,
    /// Dummy value used for error checking
    kHadronTensorGridFlag_COUNT = 2
  };

  /// Definition of sqrt() that returns zero if the argument is negative.
  /// Used to prevent spurious NaNs due to numerical roundoff.
  double real_sqrt(double x) {
    if (x < 0.) return 0.;
    else return std::sqrt(x);
  }
}

genie::TabulatedHadronTensor::TabulatedHadronTensor(
  const std::string& table_file_name)
  : fGrid(&fq0Points, &fqmagPoints, &fEntries)
{
  // Read in the table
  std::ifstream in_file( table_file_name );


  // Skip the initial comment line
  std::string dummy;
  std::getline(in_file, dummy);

  /// \todo Add error checks
  std::string type_name;
  int Z, A, num_q0, num_q_mag;

  /// \todo Use type name
  in_file >> Z >> A >> type_name >> num_q0 >> num_q_mag;

  int q0_flag;
  in_file >> q0_flag;
  read1DGridValues(num_q0, q0_flag, in_file, fq0Points);

  int q_mag_flag;
  in_file >> q_mag_flag;
  read1DGridValues(num_q_mag, q_mag_flag, in_file, fqmagPoints);

  set_pdg( genie::pdg::IonPdgCode(A, Z) );

  for (int j = 0; j < num_q0; ++j) {
    for (int k = 0; k < num_q_mag; ++k) {

      fEntries.push_back( TableEntry() );
      TableEntry& entry = fEntries.back();

      in_file >> entry.W00 >> entry.ReW0z >> entry.Wxx
        >> entry.ImWxy >> entry.Wzz;
    }
  }
}

genie::TabulatedHadronTensor::~TabulatedHadronTensor()
{
}

std::complex<double> genie::TabulatedHadronTensor::tt(
  double q0, double q_mag) const
{
  TableEntry entry = fGrid.interpolate(q0, q_mag);
  return std::complex<double>(entry.W00, 0.);
}

std::complex<double> genie::TabulatedHadronTensor::tz(
  double q0, double q_mag) const
{
  TableEntry entry = fGrid.interpolate(q0, q_mag);
  // Currently only the real part of W0z is tabulated
  /// \todo Think about adding the imaginary part even though it's not needed
  /// for the cross section calculation
  return std::complex<double>(entry.ReW0z, 0.);
}

std::complex<double> genie::TabulatedHadronTensor::xx(
  double q0, double q_mag) const
{
  TableEntry entry = fGrid.interpolate(q0, q_mag);
  return std::complex<double>(entry.Wxx, 0.);
}

std::complex<double> genie::TabulatedHadronTensor::xy(
  double q0, double q_mag) const
{
  TableEntry entry = fGrid.interpolate(q0, q_mag);
  // The Wxy element is purely imaginary
  return std::complex<double>(0., entry.ImWxy);
}

std::complex<double> genie::TabulatedHadronTensor::zz(
  double q0, double q_mag) const
{
  TableEntry entry = fGrid.interpolate(q0, q_mag);
  // The Wxy element is purely imaginary
  return std::complex<double>(entry.Wzz, 0.);
}

double genie::TabulatedHadronTensor::W1(double q0,
  double q_mag, double Mi) const
{
  TableEntry entry = fGrid.interpolate(q0, q_mag);
  return W1(q0, q_mag, entry) / Mi;
}

double genie::TabulatedHadronTensor::W2(double q0,
  double q_mag, double Mi) const
{
  TableEntry entry = fGrid.interpolate(q0, q_mag);
  return W2(q0, q_mag, entry) / Mi;
}

double genie::TabulatedHadronTensor::W3(double q0,
  double q_mag, double /*Mi*/) const
{
  TableEntry entry = fGrid.interpolate(q0, q_mag);
  return W3(q0, q_mag, entry);
}

double genie::TabulatedHadronTensor::W4(double q0,
  double q_mag, double Mi) const
{
  TableEntry entry = fGrid.interpolate(q0, q_mag);
  return W4(q0, q_mag, entry) * Mi;
}

double genie::TabulatedHadronTensor::W5(double q0,
  double q_mag, double /*Mi*/) const
{
  TableEntry entry = fGrid.interpolate(q0, q_mag);
  return W5(q0, q_mag, entry);
}

double genie::TabulatedHadronTensor::W6(double /* q0 */,
  double /* q_mag */, double /* Mi */) const
{
  return 0.;
}

void genie::TabulatedHadronTensor::read1DGridValues(int num_points,
  int flag, std::ifstream& in_file, std::vector<double>& vec_to_fill)
{
  vec_to_fill.clear();

  if (flag >= kHadronTensorGridFlag_COUNT) {
    LOG("TabulatedHadronTensor", pERROR)
      << "Invalid hadron tensor grid flag value \"" << flag
      << "\" encountered.";
    return;
  }
  else if (flag == kStartAndStep) {
    double start_val, step;
    in_file >> start_val >> step;
    for (int k = 0; k < num_points; ++k)
      vec_to_fill.push_back( start_val + (k * step) );
  }
  else if (flag == kExplicitValues) {
    double val;
    for (int k = 0; k < num_points; ++k) {
      in_file >> val;
      vec_to_fill.push_back( val );
    }
  }

}

double genie::TabulatedHadronTensor::W1(double /*q0*/,
 double /*q_mag*/, const TableEntry& entry) const
{
  return entry.Wxx / 2.;
}

double genie::TabulatedHadronTensor::W2(double q0, double q_mag,
  const TableEntry& entry) const
{
  double temp_1 = ( std::pow(q0, 2) / std::pow(q_mag, 2) )
    * (entry.Wzz - entry.Wxx);
  double temp_2 = 2. * q0 * entry.ReW0z / q_mag;
  return ( entry.W00 + entry.Wxx + temp_1 - temp_2 ) / 2.;
}

double genie::TabulatedHadronTensor::W3(double /*q0*/, double q_mag,
  const TableEntry& entry) const
{
  return ( entry.ImWxy / q_mag );
}

double genie::TabulatedHadronTensor::W4(double /*q0*/, double q_mag,
  const TableEntry& entry) const
{
  return ( entry.Wzz - entry.Wxx ) / ( 2. * std::pow(q_mag, 2) );
}

double genie::TabulatedHadronTensor::W5(double q0, double q_mag,
  const TableEntry& entry) const
{
  return ( entry.ReW0z - q0 * (entry.Wzz - entry.Wxx) / q_mag ) / q_mag;
}

double genie::TabulatedHadronTensor::W6(double /*q0*/,
  double /*q_mag*/, const TableEntry& /*entry*/) const
{
  // Currently, \Im W^{0z} has not been tabulated, so we can't really
  // evaluate W6. This isn't a huge problem, however, because W6 doesn't
  // contribute to neutrino cross sections.
  return 0.;
}

double genie::TabulatedHadronTensor::dSigma_dT_dCosTheta(
 const Interaction* interaction, double Q_value) const
{
  // Don't do anything if you've been handed a nullptr
  if ( !interaction ) return 0.;

  int probe_pdg     = interaction->InitState().ProbePdg();
  double E_probe    = interaction->InitState().ProbeE(kRfLab);
  double m_probe    = interaction->InitState().Probe()->Mass();
  double Tl      = interaction->Kine().GetKV(kKVTl);
  double cos_l   = interaction->Kine().GetKV(kKVctl);
  double ml      = interaction->FSPrimLepton()->Mass();

  return dSigma_dT_dCosTheta(probe_pdg, E_probe, m_probe, Tl, cos_l, ml,
    Q_value);
}

double genie::TabulatedHadronTensor::dSigma_dT_dCosTheta(int probe_pdg,
  double E_probe, double m_probe, double Tl, double cos_l, double ml,
  double Q_value) const
{
  // dSigma_dT_dCosTheta in GeV^(-3)
  double xsec = 0.;

  /// \todo Add check if in grid. If not, return 0

  // Final state lepton total energy
  double El = Tl + ml;

  // Energy transfer (uncorrected)
  double q0 = E_probe - El;

  // The corrected energy transfer takes into account the binding
  // energy of the struck nucleon(s)
  double q0_corrected = q0 - Q_value;

  // Magnitude of the initial state lepton 3-momentum
  double k_initial = real_sqrt( std::pow(E_probe, 2) - std::pow(m_probe, 2) );

  // Magnitude of the final state lepton 3-momentum
  double k_final = real_sqrt( std::pow(Tl, 2) + 2*ml*Tl );

  // Magnitude of the 3-momentum transfer
  double q_mag2 = std::pow(k_initial, 2) + std::pow(k_final, 2)
    - 2.*k_initial*k_final*cos_l;
  double q_mag = real_sqrt( q_mag2 );

  // Find the appropriate values of the hadron tensor elements for the
  // given combination of q0_corrected and q_mag
  TableEntry entry = fGrid.interpolate(q0_corrected, q_mag);

  // The half-angle formulas come in handy here. See, e.g.,
  // http://mathworld.wolfram.com/Half-AngleFormulas.html
  double s2_half = (1. - cos_l) / 2.; // sin^2(theta/2) = (1 - cos(theta)) / 2
  double c2_half = 1. - s2_half; // cos^2(theta / 2) = 1 - sin^2(theta / 2)

  // Calculate a nonzero cross section only for incident (anti)electrons
  // or (anti)neutrinos
  int abs_probe_pdg = std::abs(probe_pdg);

  /// \todo Add any needed changes for positron projectiles
  if ( probe_pdg == genie::kPdgElectron ) {
    // (e,e') differential cross section

    double q2 = std::pow(q0, 2) - q_mag2;

    double mott = std::pow(
      genie::constants::kAem / (2. * E_probe * s2_half), 2) * c2_half;

    // Longitudinal part
    double l_part = std::pow(q2 / q_mag2, 2) * entry.W00;

    // Transverse part
    double t_part = ( (2. * s2_half / c2_half) - (q2 / q_mag2) ) * entry.Wxx;

    xsec = (2. * genie::constants::kPi) * mott * (l_part + t_part);
  }
  else if ( abs_probe_pdg == genie::kPdgNuE
    || abs_probe_pdg == genie::kPdgNuMu
    || abs_probe_pdg == genie::kPdgNuTau )
  {
    // Simplify the expressions below by pre-computing the structure functions
    double w1 = this->W1(q0, q_mag, entry);
    double w2 = this->W2(q0, q_mag, entry);
    double w4 = this->W4(q0, q_mag, entry);
    double w5 = this->W5(q0, q_mag, entry);

    // Flip the sign of the terms proportional to W3 for antineutrinos
    double w3 = this->W3(q0, q_mag, entry);
    if (probe_pdg < 0) w3 *= -1;

    double part_1 = (2. * w1 * s2_half) + (w2 * c2_half)
      - (w3 * (E_probe + El) * s2_half);

    double part_2 = (w1 * cos_l) - (w2/2. * cos_l)
      + (w3/2. * (El + k_final - (E_probe + El)*cos_l))
      + (w4/2. * (std::pow(ml, 2)*cos_l + 2*El*(El + k_final)*s2_half))
      - (w5/2. * (El + k_final));

    double all_terms = part_1 + std::pow(ml, 2)
      * part_2 / (El * (El + k_final));

    xsec = (2. / genie::constants::kPi) * k_final * El
      * genie::constants::kGF2 * all_terms;
  }

  return xsec;
}
