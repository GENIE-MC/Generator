// standard library includes
#include <fstream>

// GENIE includes
#include "Physics/Multinucleon/XSection/TabulatedValenciaHadronTensor.h"

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
  return W1(q0, q_mag, Mi, entry);
}

double genie::TabulatedValenciaHadronTensor::W2(double q0,
  double q_mag, double Mi) const
{
  TableEntry entry = fGrid.interpolate(q0, q_mag);
  return W2(q0, q_mag, Mi, entry);
}

double genie::TabulatedValenciaHadronTensor::W3(double q0,
  double q_mag, double Mi) const
{
  TableEntry entry = fGrid.interpolate(q0, q_mag);
  return W3(q0, q_mag, Mi, entry);
}

double genie::TabulatedValenciaHadronTensor::W4(double q0,
  double q_mag, double Mi) const
{
  TableEntry entry = fGrid.interpolate(q0, q_mag);
  return W4(q0, q_mag, Mi, entry);
}

double genie::TabulatedValenciaHadronTensor::W5(double q0,
  double q_mag, double Mi) const
{
  TableEntry entry = fGrid.interpolate(q0, q_mag);
  return W5(q0, q_mag, Mi, entry);
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
 double /*q_mag*/, double Mi, const TableEntry& entry) const
{
  return entry.Wxx / (2. * Mi);
}

double genie::TabulatedValenciaHadronTensor::W2(double q0, double q_mag,
  double Mi, const TableEntry& entry) const
{
  double temp_1 = ( std::pow(q0, 2) / std::pow(q_mag, 2) )
    * (entry.Wzz - entry.Wxx);
  double temp_2 = 2. * q0 * entry.ReW0z / q_mag;
  return ( entry.W00 + entry.Wxx + temp_1 + temp_2 ) / (2. * Mi);
}

double genie::TabulatedValenciaHadronTensor::W3(double /*q0*/, double q_mag,
  double /*Mi*/, const TableEntry& entry) const
{
  return ( entry.ImWxy / q_mag );
}

double genie::TabulatedValenciaHadronTensor::W4(double /*q0*/, double q_mag,
  double Mi, const TableEntry& entry) const
{
  return Mi * ( entry.Wzz - entry.Wxx ) / ( 2. * std::pow(q_mag, 2) );
}

double genie::TabulatedValenciaHadronTensor::W5(double q0, double q_mag,
  double /*Mi*/, const TableEntry& entry) const
{
  return ( entry.ReW0z - q0 * (entry.Wzz - entry.Wxx) / q_mag ) / q_mag;
}

double genie::TabulatedValenciaHadronTensor::W6(double /*q0*/,
  double /*q_mag*/, double /*Mi*/, const TableEntry& /*entry*/) const
{
  // Currently, \Im W^{0z} has not been tabulated, so we can't really
  // evaluate W6. This isn't a huge problem, however, because W6 doesn't
  // contribute to neutrino cross sections.
  return 0.;
}
