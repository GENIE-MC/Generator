//____________________________________________________________________________
/*!

\class    genie::TabulatedValenciaHadronTensor

\brief    Computes the elements and structure functions of the Valencia MEC
          model hadron tensor using precomputed tables.
          Is a concrete implementation of the ValenciaHadronTensorI interface.

\author   Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  August 23, 2018

\cpright  Copyright (c) 2003-2018, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#pragma once

#include <deque>
#include <iterator>
#include <vector>

#include "Physics/Multinucleon/XSection/ValenciaHadronTensorI.h"

namespace genie {

template<typename ZObject, typename IndexType = int, typename XType = double,
  typename YType = double> class TabulatedValenciaHadronTensorGrid
{
  public:

  TabulatedValenciaHadronTensorGrid(const std::vector<XType>* X,
    const std::vector<YType>* Y, const std::deque<ZObject>* Z,
    bool extrapolate = false) : fX(X), fY(Y), fZ(Z), fExtrapolate(extrapolate)
  {}

  // Grid boundaries
  inline XType x_min() const { return fX->front(); }
  inline XType x_max() const { return fX->back(); }

  inline YType y_min() const { return fY->front(); }
  inline YType y_max() const { return fY->back(); }

  IndexType index_Z(IndexType ix, IndexType iy) {
    IndexType num_y_points = fY->size();

    return (num_y_points * ix) + iy;
  }

  template<class MemberType, MemberType ZObject::* member>
    MemberType interpolate(double x, double y)
  {
    IndexType ix_lo, ix_hi, iy_lo, iy_hi;

    // For points outside the grid, evaluate the function at the end points
    // unless extrapolation is enabled.
    XType evalx = x;
    if (!fExtrapolate) {
      evalx = std::min(x, x_max());
      evalx = std::max(x, x_min());
    }

    YType evaly = y;
    if (!fExtrapolate) {
      evaly = std::min(y, y_max());
      evaly = std::max(y, y_min());
    }

    get_bound_indices(fX, evalx, ix_lo, ix_hi);
    get_bound_indices(fY, evaly, iy_lo, iy_hi);

    XType x1 = fX->at( ix_lo );
    XType x2 = fX->at( ix_hi );
    YType y1 = fY->at( iy_lo );
    YType y2 = fY->at( iy_hi );

    MemberType z11 = fZ->at( this->index_Z(ix_lo, iy_lo) ).*member;
    MemberType z21 = fZ->at( this->index_Z(ix_hi, iy_lo) ).*member;
    MemberType z12 = fZ->at( this->index_Z(ix_lo, iy_hi) ).*member;
    MemberType z22 = fZ->at( this->index_Z(ix_hi, iy_hi) ).*member;

    MemberType z1  = z11 * (x2-evalx)/(x2-x1) + z21 * (evalx-x1)/(x2-x1);
    MemberType z2  = z12 * (x2-evalx)/(x2-x1) + z22 * (evalx-x1)/(x2-x1);
    MemberType z   = z1  * (y2-evaly)/(y2-y1) + z2  * (evaly-y1)/(y2-y1);

    return z;
  }

protected:

  bool fExtrapolate;

  const std::vector<XType>* fX;
  const std::vector<YType>* fY;
  const std::deque<ZObject>* fZ;

  template <typename Type> bool get_bound_indices(const std::vector<Type>* vec,
    Type val, int& lower_index, int& upper_index)
  {
    /// \todo Check that the vector contains at least two entries

    bool within = true;

    typedef typename std::vector<Type>::const_iterator Iterator;
    Iterator begin = vec->cbegin();
    Iterator end = vec->cend();

    // std::lower_bound returns an iterator to the first element of the
    // container which is not less than the supplied value
    Iterator not_less_point = std::lower_bound(begin, end, val);

    Iterator lower_point;

    // Check whether the requested grid point is within the grid limits
    if (not_less_point == begin) {
      lower_point = begin;
      // first element of vec > val
      if (*begin != val) within = false;
    }
    else if (not_less_point == end) {
      // last element of vec < val
      within = false;
      lower_point = end - 2;
    }
    else {
      // x is within the grid limits
      lower_point = not_less_point - 1;
    }

    lower_index = std::distance(begin, lower_point);
    upper_index = lower_index + 1;

    return within;
  }

};

class TabulatedValenciaHadronTensor : public ValenciaHadronTensorI {

public:

  TabulatedValenciaHadronTensor(const std::string& table_file_name);
    // : fGrid(&fq0Points, &fq3Points, &fEntries)
  virtual ~TabulatedValenciaHadronTensor();

  void Set_q0_q3(double q0, double q3);

  // \todo Enable override specifiers when GENIE modernizes to C++11

  virtual std::complex<double> tt() const /*override*/;

  virtual std::complex<double> tx() const /*override*/;

  virtual std::complex<double> ty() const /*override*/;

  virtual std::complex<double> tz() const /*override*/;

  virtual std::complex<double> xx() const /*override*/;

  virtual std::complex<double> xy() const /*override*/;

  virtual std::complex<double> xz() const /*override*/;

  virtual std::complex<double> yy() const /*override*/;

  virtual std::complex<double> yz() const /*override*/;

  virtual std::complex<double> zz() const /*override*/;

  virtual double W1() const /*override*/;
  virtual double W2() const /*override*/;
  virtual double W3() const /*override*/;
  virtual double W4() const /*override*/;
  virtual double W5() const /*override*/;
  virtual double W6() const /*override*/;

  double fq0;
  double fq3;

  struct TableEntry {
    double W00;
    double Wxx;
    double Wzz;
    double ImWxy;
    double ReW0z;
  };

  std::vector<double> fq0Points;
  std::vector<double> fq3Points;
  std::deque<TableEntry> fEntries;

  TabulatedValenciaHadronTensorGrid<TableEntry> fGrid;
};

}       // genie namespace
