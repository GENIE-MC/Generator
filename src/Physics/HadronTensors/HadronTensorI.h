//____________________________________________________________________________
/*!

\class    genie::HadronTensorI

\brief    Abstract interface for an object that computes the elements
          a hadron tensor \f$W^{\mu\nu}\f$. Also computes the contraction
          of the hadron tensor with the lepton tensor
          \f$L_{\mu\nu}W^{\mu\nu}\f$ for one or more kinds of projectile
          (e.g., neutrinos, electrons)

\author   Steven Gardiner <gardiner \at fnal.gov>
          Liang Liu <liangliu \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  August 23, 2018

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE

*/
//____________________________________________________________________________

#ifndef HADRON_TENSORI_H
#define HADRON_TENSORI_H

// standard library includes
#include <complex>

// GENIE includes
#include "Framework/Interaction/Interaction.h"
#include "Framework/ParticleData/PDGUtils.h"

namespace genie {

/// Enumerated type that describes the physics represented by a particular
/// hadron tensor
/// \todo Document enum values
typedef enum HadronTensorType {
  kHT_Undefined = -1,
  kHT_MEC_FullAll,
  kHT_MEC_Fullpn,
  kHT_MEC_DeltaAll,
  kHT_MEC_Deltapn,
  kHT_MEC_EM,
  kHT_MEC_EM_pn,
  kHT_MEC_EM_pp,
  kHT_QE_EM,
  kHT_MEC_FullAll_Param,
  kHT_MEC_FullAll_wImag,
  kHT_QE_Full,
  kHT_MEC_EM_wImag,

  kHT_QE_CRPA_Low,
  kHT_QE_CRPA_Medium,
  kHT_QE_CRPA_High,

  kHT_QE_CRPA_anu_Low,
  kHT_QE_CRPA_anu_Medium,
  kHT_QE_CRPA_anu_High,

  kHT_QE_HF_Low,
  kHT_QE_HF_Medium,
  kHT_QE_HF_High,

  kHT_QE_HF_anu_Low,
  kHT_QE_HF_anu_Medium,
  kHT_QE_HF_anu_High,

  kHT_QE_CRPAPW_Low,
  kHT_QE_CRPAPW_Medium,
  kHT_QE_CRPAPW_High,

  kHT_QE_CRPAPW_anu_Low,
  kHT_QE_CRPAPW_anu_Medium,
  kHT_QE_CRPAPW_anu_High,

  kHT_QE_HFPW_Low,
  kHT_QE_HFPW_Medium,
  kHT_QE_HFPW_High,

  kHT_QE_HFPW_anu_Low,
  kHT_QE_HFPW_anu_Medium,
  kHT_QE_HFPW_anu_High,

  kHT_QE_SuSABlend,
  kHT_QE_SuSABlend_anu,

  kHT_QE_EM_proton,
  kHT_QE_EM_neutron
}
HadronTensorType_t;

class HadronTensorI {

public:

  inline virtual ~HadronTensorI() {}

  /// \name Tensor elements
  /// \brief Functions that return the elements of the tensor.
  /// \param[in] q0 The energy transfer \f$q^0\f$ in the lab frame (GeV)
  /// \param[in] q_mag The magnitude of the 3-momentum transfer
  /// \f$\left|\overrightarrow{q}\right|\f$ in the lab frame (GeV)
  /// \retval std::complex<double> The value of the hadronic tensor element
  /// @{

  /// The tensor element \f$W^{00}\f$
  virtual std::complex<double> tt(double q0, double q_mag) const = 0;

  /// The tensor element \f$W^{0x}\f$
  virtual std::complex<double> tx(double q0, double q_mag) const = 0;

  /// The tensor element \f$W^{0y}\f$
  virtual std::complex<double> ty(double q0, double q_mag) const = 0;

  /// The tensor element \f$W^{0z}\f$
  virtual std::complex<double> tz(double q0, double q_mag) const = 0;

  /// The tensor element \f$W^{x0} = (W^{0x})^*\f$
  virtual std::complex<double> xt(double q0, double q_mag) const = 0;

  /// The tensor element \f$W^{xx}\f$
  virtual std::complex<double> xx(double q0, double q_mag) const = 0;

  /// The tensor element \f$W^{xy}\f$
  virtual std::complex<double> xy(double q0, double q_mag) const = 0;

  /// The tensor element \f$W^{xz}\f$
  virtual std::complex<double> xz(double q0, double q_mag) const = 0;

  /// The tensor element \f$W^{y0} = (W^{0y})^*\f$
  virtual std::complex<double> yt(double q0, double q_mag) const = 0;

  /// The tensor element \f$W^{yx} = (W^{xy})^*\f$
  virtual std::complex<double> yx(double q0, double q_mag) const = 0;

  /// The tensor element \f$W^{yy}\f$
  virtual std::complex<double> yy(double q0, double q_mag) const = 0;

  /// The tensor element \f$W^{yz}\f$
  virtual std::complex<double> yz(double q0, double q_mag) const = 0;

  /// The tensor element \f$W^{z0} = (W^{0z})^*\f$
  virtual std::complex<double> zt(double q0, double q_mag) const = 0;

  /// The tensor element \f$W^{zx} = (W^{xz})^*\f$
  virtual std::complex<double> zx(double q0, double q_mag) const = 0;

  /// The tensor element \f$W^{zy} = (W^{yz})^*\f$
  virtual std::complex<double> zy(double q0, double q_mag) const = 0;

  /// The tensor element \f$W^{zz}\f$
  virtual std::complex<double> zz(double q0, double q_mag) const = 0;
  /// @}

  /// Computes the contraction \f$L_{\mu\nu}W^{\mu\nu}\f$ of the hadron tensor
  /// with the appropriate lepton tensor for a given type of projectile
  /// (e.g., neutrino, electron)
  /// \param[in] interaction An Interaction object storing information about the
  /// initial and final states
  /// \param[in] Q_value The Q-value that should be used to correct
  /// the energy transfer \f$q_0\f$ (GeV)
  /// \returns The tensor contraction \f$L_{\mu\nu}W^{\mu\nu}\f$ (GeV)
  virtual double contraction(const Interaction* interaction,
    double Q_value) const = 0;

  /// PDG code of the target nucleus
  inline int pdg() const { return fTargetPDG; }

  /// Atomic number of the target nucleus
  inline int Z() const { return genie::pdg::IonPdgCodeToZ(fTargetPDG); }

  /// Mass number of the target nucleus
  inline int A() const { return genie::pdg::IonPdgCodeToA(fTargetPDG); }

  /// Set the target nucleus PDG code
  inline void set_pdg(int pdg) { fTargetPDG = pdg; }

  /// The minimum value of the energy transfer \f$q^0\f$ for which this
  /// hadron tensor may be used to compute cross sections
  virtual double q0Min() const = 0;

  /// The maximum value of the energy transfer \f$q^0\f$ for which this
  /// hadron tensor may be used to compute cross sections
  virtual double q0Max() const = 0;

  /// The minimum value of the magnitude of the 3-momentum transfer
  /// \f$\left|\overrightarrow{q}\right|\f$ for which this
  /// hadron tensor may be used to compute cross sections
  virtual double qMagMin() const = 0;

  /// The maximum value of the magnitude of the 3-momentum transfer
  /// \f$\left|\overrightarrow{q}\right|\f$ for which this
  /// hadron tensor may be used to compute cross sections
  virtual double qMagMax() const = 0;

protected:

  inline HadronTensorI(int pdg = 0) : fTargetPDG(pdg) {}

  inline HadronTensorI(int Z, int A)
    : fTargetPDG( genie::pdg::IonPdgCode(A, Z) ) {}

  ///< PDG code for the target nucleus represented by the tensor
  int fTargetPDG;

}; // class HadronTensorI

} // genie namespace
#endif
