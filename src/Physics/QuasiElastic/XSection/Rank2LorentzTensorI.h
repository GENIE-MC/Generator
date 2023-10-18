//____________________________________________________________________________
/*!

\class    genie::Rank2LorentzTensor

\brief    Abstract interface for an object that computes the elements of a rank
          2 Lorentz tensor \f$A^{\mu\nu}\f$.

\author   Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  January 17, 2019

\cpright  Copyright (c) 2003-2018, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef RANK_2_LORENTZ_TENSORI_H
#define RANK_2_LORENTZ_TENSORI_H

#include <complex>

#include "TLorentzVector.h"

namespace genie {

// Order is chosen to match TLorentzVector convention
typedef enum ELorentzTensorIndices {
  kTIdx_SpaceX = 0,
  kTIdx_SpaceY = 1,
  kTIdx_SpaceZ = 2,
  kTIdx_Time = 3,
  kTIdx_NumDimensions = 4
} TensorIndex_t;

class Rank2LorentzTensorI {

public:

  inline virtual ~Rank2LorentzTensorI() {}

  /// \name Tensor elements
  /// \brief Functions that return the elements of the tensor.
  /// @{

  /// Retrieves a tensor element corresponding to the given indices
  virtual std::complex<double> operator()(TensorIndex_t mu, TensorIndex_t nu)
    const = 0;

  ///// The tensor element \f$A^{00}\f$
  //virtual std::complex<double> tt() const = 0;

  ///// The tensor element \f$A^{0x}\f$
  //virtual std::complex<double> tx() const = 0;

  ///// The tensor element \f$A^{0y}\f$
  //virtual std::complex<double> ty() const = 0;

  ///// The tensor element \f$A^{0z}\f$
  //virtual std::complex<double> tz() const = 0;

  ///// The tensor element \f$A^{x0} = (A^{0x})^*\f$
  //virtual std::complex<double> xt() const = 0;

  ///// The tensor element \f$A^{xx}\f$
  //virtual std::complex<double> xx() const = 0;

  ///// The tensor element \f$A^{xy}\f$
  //virtual std::complex<double> xy() const = 0;

  ///// The tensor element \f$A^{xz}\f$
  //virtual std::complex<double> xz() const = 0;

  ///// The tensor element \f$A^{y0} = (A^{0y})^*\f$
  //virtual std::complex<double> yt() const = 0;

  ///// The tensor element \f$A^{yx} = (A^{xy})^*\f$
  //virtual std::complex<double> yx() const = 0;

  ///// The tensor element \f$A^{yy}\f$
  //virtual std::complex<double> yy() const = 0;

  ///// The tensor element \f$A^{yz}\f$
  //virtual std::complex<double> yz() const = 0;

  ///// The tensor element \f$A^{z0} = (A^{0z})^*\f$
  //virtual std::complex<double> zt() const = 0;

  ///// The tensor element \f$A^{zx} = (A^{xz})^*\f$
  //virtual std::complex<double> zx() const = 0;

  ///// The tensor element \f$A^{zy} = (A^{yz})^*\f$
  //virtual std::complex<double> zy() const = 0;

  ///// The tensor element \f$A^{zz}\f$
  //virtual std::complex<double> zz() const = 0;

  /// @}

  /// \name Operations
  /// @{
  virtual std::complex<double> Contract(const Rank2LorentzTensorI& other) const;
  inline std::complex<double> operator*(const Rank2LorentzTensorI& other) const
  {
    return this->Contract( other );
  }
  /// @}

  // Returns the corresponding element of the metric tensor
  double Metric(TensorIndex_t mu, TensorIndex_t nu) const {
    // TODO: deal with invalid enum values
    if (mu != nu) return 0.;
    else if (mu == kTIdx_Time) return 1.;
    else return -1.;
  }

  // TODO: document this
  double LeviCivitaProduct(genie::TensorIndex_t mu, genie::TensorIndex_t nu,
    const TLorentzVector& p1, const TLorentzVector& p2) const;

  double LeviCivitaProductSF(genie::TensorIndex_t mu, genie::TensorIndex_t nu,
    const TLorentzVector& p1, const TLorentzVector& p2) const;

}; // class Rank2LorentzTensorI

} // genie namespace
#endif
