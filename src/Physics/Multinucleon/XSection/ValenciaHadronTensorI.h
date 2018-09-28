//____________________________________________________________________________
/*!

\class    genie::ValenciaHadronTensorI

\brief    Abstract interface for an object that computes the elements
          (\f$W^{xx}\f$, \f$W^{0z}\f$, etc.) and structure functions
          (\f$W_1\f$, \f$W_2\f$, etc.) of
          the hadron tensor \f$W^{\mu\nu}\f$ as defined in equation (8)
          and footnote 2 of the MEC Valencia model paper.

          The calculation is carried out in the lab frame (i.e., the frame
          where the target nucleus has initial 4-momentum
          \f$P^\mu = (M_i, \overrightarrow{0})\f$) with the 3-momentum
          transfer \f$\overrightarrow{q}\f$ chosen to lie along the
          z axis, i.e., \f$q = (q^0, |\overrightarrow{q}|
          \overrightarrow{u}_z)\f$.

\ref      J. Nieves, J. E. Amaro, and M. Valverde,
          "Inclusive Quasi-Elastic Charged-Current Neutrino-Nucleus Reactions,"
          Phys. Rev. C 70, 055503 (2004)
          https://doi.org/10.1103/PhysRevC.70.055503
          https://arxiv.org/abs/nucl-th/0408005v3

\author   Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  August 23, 2018

\cpright  Copyright (c) 2003-2018, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#pragma once

#include <complex>

namespace genie {

class ValenciaHadronTensorI {

public:

  virtual ~ValenciaHadronTensorI();

  /// \name Functions that return the elements of the tensor. Since it is
  /// Hermitian, only ten elements are independent. Although a return type of
  /// std::complex<double> is used for all elements, note that hermiticity
  /// implies that the diagonal elements must be real.
  //@{

  /// Returns the tensor element \f$W^{00}\f$
  virtual std::complex<double> tt() const = 0;

  /// Returns the tensor element \f$W^{0x}\f$
  virtual std::complex<double> tx() const = 0;

  /// Returns the tensor element \f$W^{0y}\f$
  virtual std::complex<double> ty() const = 0;

  /// Returns the tensor element \f$W^{0z}\f$
  virtual std::complex<double> tz() const = 0;

  /// Returns the tensor element \f$W^{x0} = (W^{0x})^*\f$
  inline virtual std::complex<double> xt() const
    { return std::conj( this->tx() ); }

  /// Returns the tensor element \f$W^{xx}\f$
  virtual std::complex<double> xx() const = 0;

  /// Returns the tensor element \f$W^{xy}\f$
  virtual std::complex<double> xy() const = 0;

  /// Returns the tensor element \f$W^{xz}\f$
  virtual std::complex<double> xz() const = 0;

  /// Returns the tensor element \f$W^{y0} = (W^{0y})^*\f$
  inline virtual std::complex<double> yt() const
    { return std::conj( this->ty() ); }

  /// Returns the tensor element \f$W^{yx} = (W^{xy})^*\f$
  inline virtual std::complex<double> yx() const
    { return std::conj( this->xy() ); }

  /// Returns the tensor element \f$W^{yy}\f$
  virtual std::complex<double> yy() const = 0;

  /// Returns the tensor element \f$W^{yz}\f$
  virtual std::complex<double> yz() const = 0;

  /// Returns the tensor element \f$W^{z0} = (W^{0z})^*\f$
  inline virtual std::complex<double> zt() const
    { return std::conj( this->tz() ); }

  /// Returns the tensor element \f$W^{zx} = (W^{xz})^*\f$
  inline virtual std::complex<double> zx() const
    { return std::conj( this->xz() ); }

  /// Returns the tensor element \f$W^{zy} = (W^{yz})^*\f$
  inline virtual std::complex<double> zy() const
    { return std::conj( this->yz() ); }

  /// Returns the tensor element \f$W^{zz}\f$
  virtual std::complex<double> zz() const = 0;
  //@}

  /// \name Structure functions
  //@{
  virtual double W1() const = 0;
  virtual double W2() const = 0;
  virtual double W3() const = 0;
  virtual double W4() const = 0;
  virtual double W5() const = 0;
  virtual double W6() const = 0;
  //@}

protected:

  ValenciaHadronTensorI();
};

}       // genie namespace
