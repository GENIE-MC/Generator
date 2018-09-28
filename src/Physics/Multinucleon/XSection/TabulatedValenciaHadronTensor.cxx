// standard library includes
#include <fstream>

// GENIE includes
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/Multinucleon/XSection/TabulatedValenciaHadronTensor.h"

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


  // Definition of sqrt() that returns zero if the argument is negative.
  // Used to prevent spurious NaNs due to numerical roundoff.
  double real_sqrt(double x) {
    if (x < 0.) return 0.;
    else return std::sqrt(x);
  }
}

genie::TabulatedValenciaHadronTensor::TabulatedValenciaHadronTensor(
  const std::string& table_file_name)
  : fGrid(&fq0Points, &fqmagPoints, &fEntries)
{
  // Read in the table
  std::ifstream in_file( table_file_name );


  // Skip the initial comment line
  std::string dummy;
  std::getline(in_file, dummy);

  /// \todo Add error checks
  int Z, A, num_q0, num_q_mag;

  in_file >> Z >> A >> num_q0 >> num_q_mag;

  set_pdg( nucleus_pdg(Z, A) );

  double q0, q_mag;

  for (int j = 0; j < num_q0; ++j) {
    for (int k = 0; k < num_q_mag; ++k) {

      in_file >> q0 >> q_mag;

      fq0Points.push_back(q0);
      fqmagPoints.push_back(q_mag);

      fEntries.push_back( TableEntry() );
      TableEntry& entry = fEntries.back();

      in_file >> entry.W00 >> entry.ReW0z >> entry.Wxx
        >> entry.ImWxy >> entry.Wzz;
    }
  }
}

genie::TabulatedValenciaHadronTensor::~TabulatedValenciaHadronTensor()
{
}

std::complex<double> genie::TabulatedValenciaHadronTensor::tt(
  double q0, double q_mag) const
{
  TableEntry entry = fGrid.interpolate(q0, q_mag);
  return std::complex<double>(entry.W00, 0.);
}

std::complex<double> genie::TabulatedValenciaHadronTensor::tz(
  double q0, double q_mag) const
{
  TableEntry entry = fGrid.interpolate(q0, q_mag);
  // Currently only the real part of W0z is tabulated
  /// \todo Think about adding the imaginary part even though it's not needed
  /// for the cross section calculation
  return std::complex<double>(entry.ReW0z, 0.);
}

std::complex<double> genie::TabulatedValenciaHadronTensor::xx(
  double q0, double q_mag) const
{
  TableEntry entry = fGrid.interpolate(q0, q_mag);
  return std::complex<double>(entry.Wxx, 0.);
}

std::complex<double> genie::TabulatedValenciaHadronTensor::xy(
  double q0, double q_mag) const
{
  TableEntry entry = fGrid.interpolate(q0, q_mag);
  // The Wxy element is purely imaginary
  return std::complex<double>(0., entry.ImWxy);
}

std::complex<double> genie::TabulatedValenciaHadronTensor::zz(
  double q0, double q_mag) const
{
  TableEntry entry = fGrid.interpolate(q0, q_mag);
  // The Wxy element is purely imaginary
  return std::complex<double>(entry.Wzz, 0.);
}

double genie::TabulatedValenciaHadronTensor::W1(double q0,
  double q_mag, double Mi) const
{
  TableEntry entry = fGrid.interpolate(q0, q_mag);
  return W1(q0, q_mag, entry) / Mi;
}

double genie::TabulatedValenciaHadronTensor::W2(double q0,
  double q_mag, double Mi) const
{
  TableEntry entry = fGrid.interpolate(q0, q_mag);
  return W2(q0, q_mag, entry) / Mi;
}

double genie::TabulatedValenciaHadronTensor::W3(double q0,
  double q_mag, double /*Mi*/) const
{
  TableEntry entry = fGrid.interpolate(q0, q_mag);
  return W3(q0, q_mag, entry);
}

double genie::TabulatedValenciaHadronTensor::W4(double q0,
  double q_mag, double Mi) const
{
  TableEntry entry = fGrid.interpolate(q0, q_mag);
  return W4(q0, q_mag, entry) * Mi;
}

double genie::TabulatedValenciaHadronTensor::W5(double q0,
  double q_mag, double /*Mi*/) const
{
  TableEntry entry = fGrid.interpolate(q0, q_mag);
  return W5(q0, q_mag, entry);
}

double genie::TabulatedValenciaHadronTensor::W6(double /* q0 */,
  double /* q_mag */, double /* Mi */) const
{
  return 0.;
}

void genie::TabulatedValenciaHadronTensor::read1DGridValues(int num_points,
  int flag, std::ifstream& in_file, std::vector<double>& vec_to_fill)
{
  vec_to_fill.clear();

  if (flag >= kHadronTensorGridFlag_COUNT) {
    LOG("TabulatedValenciaHadronTensor", pERROR)
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

double genie::TabulatedValenciaHadronTensor::W1(double /*q0*/,
 double /*q_mag*/, const TableEntry& entry) const
{
  return entry.Wxx / 2.;
}

double genie::TabulatedValenciaHadronTensor::W2(double q0, double q_mag,
  const TableEntry& entry) const
{
  double temp_1 = ( std::pow(q0, 2) / std::pow(q_mag, 2) )
    * (entry.Wzz - entry.Wxx);
  double temp_2 = 2. * q0 * entry.ReW0z / q_mag;
  return ( entry.W00 + entry.Wxx + temp_1 + temp_2 ) / 2.;
}

double genie::TabulatedValenciaHadronTensor::W3(double /*q0*/, double q_mag,
  const TableEntry& entry) const
{
  return ( entry.ImWxy / q_mag );
}

double genie::TabulatedValenciaHadronTensor::W4(double /*q0*/, double q_mag,
  const TableEntry& entry) const
{
  return ( entry.Wzz - entry.Wxx ) / ( 2. * std::pow(q_mag, 2) );
}

double genie::TabulatedValenciaHadronTensor::W5(double q0, double q_mag,
  const TableEntry& entry) const
{
  return ( entry.ReW0z - q0 * (entry.Wzz - entry.Wxx) / q_mag ) / q_mag;
}

double genie::TabulatedValenciaHadronTensor::W6(double /*q0*/,
  double /*q_mag*/, const TableEntry& /*entry*/) const
{
  // Currently, \Im W^{0z} has not been tabulated, so we can't really
  // evaluate W6. This isn't a huge problem, however, because W6 doesn't
  // contribute to neutrino cross sections.
  return 0.;
}

double genie::TabulatedValenciaHadronTensor::dSigma_dT_dCosTheta(
 const Interaction* interaction, double Q_value) const
{
  // Don't do anything if you've been handed a nullptr
  if ( !interaction ) return 0.;

  int nu_pdg     = interaction->InitState().ProbePdg();
  double E_nu    = interaction->InitState().ProbeE(kRfLab);
  double Tl      = interaction->Kine().GetKV(kKVTl);
  double cos_l   = interaction->Kine().GetKV(kKVctl);
  double ml      = interaction->FSPrimLepton()->Mass();

  return dSigma_dT_dCosTheta(nu_pdg, E_nu, Tl, cos_l, ml, Q_value);
}

double genie::TabulatedValenciaHadronTensor::dSigma_dT_dCosTheta(int nu_pdg,
  double E_nu, double Tl, double cos_l, double ml, double Q_value) const
{
  /// \todo Add check if in grid. If not, return 0

  // Lepton total energy
  double El = Tl + ml;

  // Energy transfer (uncorrected)
  double q0 = E_nu - El;

  // The corrected energy transfer takes into account the binding
  // energy of the struck nucleon(s)
  double q0_corrected = q0 - Q_value;

  // Magnitude of the final state lepton 3-momentum
  double k_prime = real_sqrt( std::pow(Tl, 2) + 2*ml*Tl );

  // Magnitude of the 3-momentum transfer
  double q_mag = real_sqrt( std::pow(E_nu, 2) + std::pow(k_prime, 2)
    - 2.*E_nu*k_prime*cos_l );

  // Find the appropriate values of the hadron tensor elements for the
  // given combination of q0_corrected and q_mag
  TableEntry entry = fGrid.interpolate(q0_corrected, q_mag);

  // The half-angle formulas come in handy here. See, e.g.,
  // http://mathworld.wolfram.com/Half-AngleFormulas.html
  double s2_half = (1. - cos_l) / 2.; // sin^2(theta/2) = (1 - cos(theta)) / 2
  double c2_half = 1. - s2_half; // cos^2(theta / 2) = 1 - sin^2(theta / 2)

  // Simplify the expressions below by pre-computing the structure functions
  double w1 = this->W1(q0, q_mag, entry);
  double w2 = this->W1(q0, q_mag, entry);
  double w4 = this->W1(q0, q_mag, entry);
  double w5 = this->W1(q0, q_mag, entry);

  // Flip the sign of the terms proportional to W3 for antineutrinos
  double w3 = this->W3(q0, q_mag, entry);
  if (nu_pdg < 0) w3 *= -1;

  double part_1 = (2. * w1 * s2_half) + (w2 * c2_half)
    - (w3 * (E_nu + El) * s2_half);

  double part_2 = (w1 * cos_l) - (w2/2. * cos_l)
    + (w3/2. * (El + k_prime - (E_nu + El)*cos_l))
    + (w4/2. * (std::pow(ml, 2)*cos_l + 2*El*(El + k_prime)*s2_half))
    - (w5/2. * (El + k_prime));

  double all_terms = part_1 + std::pow(ml, 2) * part_2 / (El * (El + k_prime));

  double xsec = (2. / genie::constants::kPi) * k_prime * El
    * genie::constants::kGF2 * all_terms;

  // \todo Fix units!
  return xsec;
}
