#ifndef BLI2DNONUNIF_OBJECT_GRID_H_
#define BLI2DNONUNIF_OBJECT_GRID_H_

namespace genie {

//____________________________________________________________________________
/*!

\brief    A class template that performs bilinear interpolation on a
          non-uniform grid with an implementation similar to that of
          genie::BLI2DNonUnifGrid.

\details  The main differences between this class template and
          genie::BLI2DNonUnifGrid are

            * Values for the z coordinate can be any arbitrary object Object
              that implements the member functions operator*(double),
              operator*(const Object&), and operator+(const Object&)

            * Rather than C-style arrays, the grid values are accessed via
              pointers to std::vector objects

            * Upper and lower bounds on the grid are found using
              std::lower_bound() rather than a manual linear search

            * The genie::BLI2DNonUnifGrid object does not take ownership of the
              grid vectors, which must be stored elsewhere

\tparam   ZObject Type of the object describing each z coordinate
\tparam   IndexType Type to use when computing indices in the vectors
\tparam   XType Type used to represent x coordinates
\tparam   YType Type used to represent y coordinates

\todo     Think about how to have this class inherit from the other BLI2D
          classes

\author   Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  August 23, 2018

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

template<typename ZObject, typename IndexType = int, typename XType = double,
  typename YType = double> class BLI2DNonUnifObjectGrid
{
  public:

  /// \param[in] X Pointer to a vector of x coordinates
  /// \param[in] Y Pointer to a vector of y coordinates
  /// \param[in] Z Pointer to a vector of z coordinates
  /// \param[in] extrapolate Whether to allow bilinear extrapolation (true)
  /// or to use the grid endpoints (false) when evaluating z values outside
  /// of the grid
  BLI2DNonUnifObjectGrid(const std::vector<XType>* X,
    const std::vector<YType>* Y, const std::vector<ZObject>* Z,
    bool extrapolate = false) : fX(X), fY(Y), fZ(Z), fExtrapolate(extrapolate)
  {}

  /// Retrieve the minimum x value
  inline XType x_min() const { return fX->front(); }

  /// Retrieve the maximum x value
  inline XType x_max() const { return fX->back(); }

  /// Retrieve the minimum y value
  inline YType y_min() const { return fY->front(); }

  /// Retrieve the maximum y value
  inline YType y_max() const { return fY->back(); }

  /// Calculates the index in the vector of z coordinates that
  /// corresponds to a given set of x and y indices
  /// \param[in] ix Index of the desired grid point on the x axis
  /// \param[in] iy Index of the desired grid point on the y axis
  IndexType index_Z(IndexType ix, IndexType iy) const {
    IndexType num_y_points = fY->size();

    return (num_y_points * ix) + iy;
  }

  /// Uses bilinear interpolation to compute the z coordinate (represented
  /// by an object of type ZObject) that corresponds to the given x and y
  /// coordinates
  ZObject interpolate(double x, double y) const
  {
    IndexType ix_lo, ix_hi, iy_lo, iy_hi;

    // For points outside the grid, evaluate the function at the end points
    // unless extrapolation is enabled.
    XType evalx = x;
    if (!fExtrapolate) {
      evalx = std::min(x, x_max());
      evalx = std::max(evalx, x_min());
    }

    YType evaly = y;
    if (!fExtrapolate) {
      evaly = std::min(y, y_max());
      evaly = std::max(evaly, y_min());
    }

    // Find the indices of the grid points on either side of the
    // desired x and y values. If the desired point is outside of
    // the x or y grid limits, get the indices of the two closest
    // grid points to use for possible extrapolation.
    get_bound_indices(fX, evalx, ix_lo, ix_hi);
    get_bound_indices(fY, evaly, iy_lo, iy_hi);

    // Get the x and y values corresponding to the lower (x1, y1) and
    // upper (x2, y2) bounds found previously
    XType x1 = fX->at( ix_lo );
    XType x2 = fX->at( ix_hi );
    YType y1 = fY->at( iy_lo );
    YType y2 = fY->at( iy_hi );

    // Retrieve the z values corresponding to each of the four locations
    // that will be used for the bilinear interpolation
    const ZObject& z11 = fZ->at( this->index_Z(ix_lo, iy_lo) );
    const ZObject& z21 = fZ->at( this->index_Z(ix_hi, iy_lo) );
    const ZObject& z12 = fZ->at( this->index_Z(ix_lo, iy_hi) );
    const ZObject& z22 = fZ->at( this->index_Z(ix_hi, iy_hi) );

    // Perform the interpolation (first y, then x)
    ZObject z1  = z11 * (y2-evaly)/(y2-y1) + z12 * (evaly-y1)/(y2-y1);
    ZObject z2  = z21 * (y2-evaly)/(y2-y1) + z22 * (evaly-y1)/(y2-y1);
    ZObject z   = z1  * (x2-evalx)/(x2-x1) + z2  * (evalx-x1)/(x2-x1);

    return z;
  }

protected:

  const std::vector<XType>* fX; ///< Pointer to the vector of x coordinates
  const std::vector<YType>* fY; ///< Pointer to the vector of y coordinates

  /// Pointer to the vector of z coordinate objects
  const std::vector<ZObject>* fZ;

  /// Whether to allow bilinear extrapolation (true) or to compute z values for
  /// x and coordinates outside of the grid using the grid endpoints (false)
  bool fExtrapolate;

  /// Determines the indices for the two gridpoints surrounding a requested
  /// x or y coordinate. If the x or y coordinate is outside of the grid,
  /// this function returns the two closest grid points.
  /// \param[in] vec A vector of grid point coordinates
  /// \param[in] val The requested x or y coordinate
  /// \param[out] lower_index The index of the closest grid point less than or
  /// equal to the requested value, or the lower of the two nearest grid points
  /// if the value falls outside of the grid
  /// \param[out] upper_index The index of the closest grid point greater than
  /// the requested value, or the higher of the two nearest grid points
  /// if the value falls outside of the grid
  /// \return true if the requested value is within the grid, or false
  /// otherwise
  template <typename Type> bool get_bound_indices(const std::vector<Type>* vec,
    Type val, int& lower_index, int& upper_index) const
  {
    /// \todo Check that the vector contains at least two entries

    bool within = true;

    typedef typename std::vector<Type>::const_iterator Iterator;
    Iterator begin = vec->begin();
    Iterator end = vec->end();

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

}; // template class BLI2DNonUnifObjectGrid

} // namespace genie
#endif
