//____________________________________________________________________________
/*!

\class    genie::TabulatedLabFrameHadronTensor

\brief    Computes the elements and structure functions of the hadron
          tensor \f$W^{\mu\nu}\f$ (using the conventions of the Valencia model)
          using precomputed tables.
          Is a concrete implementation of the HadronTensorI interface.

\author   Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  August 23, 2018

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef TABULATED_VALENCIA_HADRON_TENSOR_H
#define TABULATED_VALENCIA_HADRON_TENSOR_H

// standard library includes
#include <vector>

// GENIE includes
#include "Framework/Numerical/BLI2DNonUnifObjectGrid.h"
#include "Physics/HadronTensors/LabFrameHadronTensorI.h"

namespace genie {

class TabulatedLabFrameHadronTensor : public LabFrameHadronTensorI {

  public:

  TabulatedLabFrameHadronTensor(const std::string& table_file_name);
  virtual ~TabulatedLabFrameHadronTensor();

  // \todo Enable override specifiers when GENIE modernizes to C++11

  virtual std::complex<double> tt(double q0, double q_mag) const /*override*/;

  virtual std::complex<double> tz(double q0, double q_mag) const /*override*/;

  virtual std::complex<double> xx(double q0, double q_mag) const /*override*/;

  virtual std::complex<double> xy(double q0, double q_mag) const /*override*/;

  virtual std::complex<double> zz(double q0, double q_mag) const /*override*/;

  virtual double W1(double q0, double q_mag, double Mi) const /*override*/;
  virtual double W2(double q0, double q_mag, double Mi) const /*override*/;
  virtual double W3(double q0, double q_mag, double Mi) const /*override*/;
  virtual double W4(double q0, double q_mag, double Mi) const /*override*/;
  virtual double W5(double q0, double q_mag, double Mi) const /*override*/;
  virtual double W6(double q0, double q_mag, double Mi) const /*override*/;

  //virtual double contraction(const Interaction* interaction) const /*override*/;

  virtual double dSigma_dT_dCosTheta(const Interaction* interaction,
    double Q_value) const /*override*/;

  virtual double dSigma_dT_dCosTheta(int probe_pdg, double E_probe,
    double m_probe, double Tl, double cos_l, double ml, double Q_value)
    const /*override*/;

  virtual double dSigma_dT_dCosTheta_rosenbluth(const Interaction* interaction,
    double Q_value) const /*override*/;

  virtual double dSigma_dT_dCosTheta_rosenbluth(int probe_pdg, double E_probe,
    double m_probe, double Tl, double cos_l, double ml, double Q_value)
    const /*override*/;

  inline virtual double q0Min() const /*override*/ { return fGrid.x_min(); }
  inline virtual double q0Max() const /*override*/ { return fGrid.x_max(); }
  inline virtual double qMagMin() const /*override*/ { return fGrid.y_min(); }
  inline virtual double qMagMax() const /*override*/ { return fGrid.y_max(); }

  protected:

  /// Helper function that allows this class to handle variations in the
  /// data file format for the 1D \f$q_0\f$ and
  /// \f$\left|\overrightarrow{q}\right|\f$ grids
  /// \param[in] num_points The number of grid points to be read from the file
  /// \param[in] flag A numerical flag describing the grid data format
  /// \param[inout] in_file A reference to the file being read
  /// \param[out] vec_to_fill The vector that will store the grid points
  void read1DGridValues(int num_points, int flag, std::ifstream& in_file,
    std::vector<double>& vec_to_fill);

  class TableEntry {

    public:

    double W00;
    double Wxx;
    double Wzz;
    double ImWxy;
    double ReW0z;

    // Scalar operations
    TableEntry operator*(double d) const
      { return apply_to_elements(d, &TableEntry::multiply); }

    TableEntry operator/(double d) const
      { return apply_to_elements(d, &TableEntry::divide); }

    // TableEntry operations
    TableEntry operator+(const TableEntry& rhs) const
      { return apply_to_elements(rhs, &TableEntry::add); }

    TableEntry operator-(const TableEntry& rhs) const
      { return apply_to_elements(rhs, &TableEntry::subtract); }

    TableEntry operator*(const TableEntry& rhs) const
      { return apply_to_elements(rhs, &TableEntry::multiply); }

    TableEntry operator/(const TableEntry& rhs) const
      { return apply_to_elements(rhs, &TableEntry::divide); }

    bool operator==(const TableEntry& rhs) const {
      bool are_equal = this->W00 == rhs.W00;
      if ( are_equal ) are_equal = this->Wxx == rhs.Wxx;
      if ( are_equal ) are_equal = this->Wzz == rhs.Wzz;
      if ( are_equal ) are_equal = this->ImWxy == rhs.ImWxy;
      if ( are_equal ) are_equal = this->ReW0z == rhs.ReW0z;
      return are_equal;
    }

    bool operator!=(const TableEntry& rhs) const {
      return !operator==(rhs);
    }

    protected:

    TableEntry apply_to_elements(double d,
      double (*my_function)(double, double) ) const
    {
      // Create a new TableEntry object by applying the function to
      // each pair of elements
      TableEntry result;
      result.W00 = (*my_function)( this->W00, d );
      result.Wxx = (*my_function)( this->Wxx, d );
      result.Wzz = (*my_function)( this->Wzz, d );
      result.ImWxy = (*my_function)( this->ImWxy, d );
      result.ReW0z = (*my_function)( this->ReW0z, d );

      return result;
    }

    TableEntry apply_to_elements(const TableEntry& rhs,
      double (*my_function)(double, double) ) const
    {
      // Create a new TableEntry object by applying the function to
      // each pair of elements
      TableEntry result;
      result.W00 = (*my_function)( this->W00, rhs.W00 );
      result.Wxx = (*my_function)( this->Wxx, rhs.Wxx );
      result.Wzz = (*my_function)( this->Wzz, rhs.Wzz );
      result.ImWxy = (*my_function)( this->ImWxy, rhs.ImWxy );
      result.ReW0z = (*my_function)( this->ReW0z, rhs.ReW0z );
      return result;
    }

    private:

    // Dummy functions defined so that they can be applied to each class member
    // by passing a function pointer to apply_to_elements()
    /// \todo Replace these with lambdas (or something similar) when GENIE
    /// updates to C++11. This technique is a bit of a hack.
    inline static double add(double x, double y) { return x + y; }
    inline static double subtract(double x, double y) { return x - y; }
    inline static double multiply(double x, double y) { return x * y; }
    inline static double divide(double x, double y) { return x / y; }
  };

  /// \name Structure function helpers
  /// \brief These helper functions allow multiple structure
  /// function values (e.g., \f$W_1\f$ and \f$W_2\f$) to be computed
  /// without having to perform bilinear interpolation every time.
  /// \details Because the differential cross section
  /// \f$\frac{ d^2\sigma_{\nu l} }
  /// { d\cos(\theta^\prime) dE^\prime_l }\f$ does not depend on the
  /// initial nucleus's mass \f$M_i\f$, the explicit factors of \f$M_i\f$
  /// have been removed from these internally-used functions.
  /// \param[in] Hadronic tensor table entry that has been pre-interpolated
  /// @{
  virtual double W1(double q0, double q_mag, const TableEntry& entry) const;

  virtual double W2(double q0, double q_mag, const TableEntry& entry) const;

  virtual double W3(double q0, double q_mag, const TableEntry& entry) const;

  virtual double W4(double q0, double q_mag, const TableEntry& entry) const;

  virtual double W5(double q0, double q_mag, const TableEntry& entry) const;

  virtual double W6(double q0, double q_mag, const TableEntry& entry) const;
  ///@}

  std::vector<double> fq0Points;
  std::vector<double> fqmagPoints;
  std::vector<TableEntry> fEntries;

  BLI2DNonUnifObjectGrid<TableEntry> fGrid;

}; // class TabulatedLabFrameHadronTensor

}  // genie namespace
#endif
