//____________________________________________________________________________
/*!

\class    genie::LabFrameHadronTensorI

\brief    Abstract interface for an object that computes the elements
          (\f$W^{xx}\f$, \f$W^{0z}\f$, etc.) and structure functions
          (\f$W_1\f$, \f$W_2\f$, etc.) of
          the hadron tensor \f$W^{\mu\nu}\f$ defined according to the
          conventions of the Valencia model.

\details  For (anti)neutrino projectiles, the hadron tensor \f$W^{\mu\nu}\f$
          is defined as in equations (8) and (9) of

          J. Nieves, J. E. Amaro, and M. Valverde,
          "Inclusive Quasi-Elastic Charged-Current Neutrino-Nucleus Reactions,"
          Phys. Rev. C 70, 055503 (2004)
          https://doi.org/10.1103/PhysRevC.70.055503
          https://arxiv.org/abs/nucl-th/0408005v3

          Note that the associated erratum includes an important correction
          in the formula for the differential cross section:

          J. Nieves, J. E. Amaro, and M. Valverde,
          "Erratum: Inclusive quasielastic charged-current neutrino-nucleus
           reactions [Phys. Rev. C 70, 055503 (2004)],"
          Phys. Rev. C. 72, 019902 (2005)
          http://doi.org/10.1103/PhysRevC.72.019902

          The calculation is carried out in the lab frame (i.e., the frame
          where the target nucleus has initial 4-momentum
          \f$P^\mu = (M_i, \overrightarrow{0})\f$) with the 3-momentum
          transfer \f$\overrightarrow{q}\f$ chosen to lie along the
          z axis, i.e., \f$q = (q^0, |\overrightarrow{q}|
          \overrightarrow{u}_z)\f$. With this choice of frame, the only
          nonzero elements are \f$W^{00}\f$, \f$W^{0z} = (W^{z0})^*\f$,
          \f$W^{xx} = W^{yy}\f$, \f$W^{xy} = (W^{yx})^*\f$, and \f$W^{zz}\f$.

          For an electron projectile, the hadron tensor is defined as in
          equation (83) of

          A. Gil, J. Nieves, and E. Oset,
          "Many-body approach to the inclusive \f$(e,e^\prime)\f$ reaction
          from the quasielastic to the \f$\Delta\f$ excitation region,"
          Nuclear Physics A 627, 543-598 (1997)
          http://doi.org/10.1016/S0375-9474(97)00513-7

          It is evaluated in the same reference frame as the neutrino case.

\author   Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  August 23, 2018

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef VALENCIA_HADRON_TENSORI_H
#define VALENCIA_HADRON_TENSORI_H

// standard library includes
#include <complex>

// GENIE includes
#include "Framework/Interaction/Interaction.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Physics/HadronTensors/HadronTensorI.h"

namespace genie {

class LabFrameHadronTensorI : public HadronTensorI {

public:

  inline virtual ~LabFrameHadronTensorI() {}

  /// \name Tensor elements
  /// \brief Functions that return the elements of the tensor. Since it is
  /// Hermitian, only ten elements are independent. Although a return type of
  /// std::complex<double> is used for all elements, note that hermiticity
  /// implies that the diagonal elements must be real.
  /// \param[in] q0 The energy transfer \f$q^0\f$ in the lab frame (GeV)
  /// \param[in] q_mag The magnitude of the 3-momentum transfer
  /// \f$\left|\overrightarrow{q}\right|\f$ in the lab frame (GeV)
  /// \retval std::complex<double> The value of the hadronic tensor element
  /// @{

  /// The tensor element \f$W^{0x}\f$
  inline virtual std::complex<double> tx(double /*q0*/, double /*q_mag*/) const
    { return std::complex<double>(0., 0.); }

  /// The tensor element \f$W^{0y}\f$
  inline virtual std::complex<double> ty(double /*q0*/, double /*q_mag*/) const
    { return std::complex<double>(0., 0.); }

  /// The tensor element \f$W^{x0} = (W^{0x})^*\f$
  inline virtual std::complex<double> xt(double q0, double q_mag) const
    { return std::conj( this->tx(q0, q_mag) ); }

  /// The tensor element \f$W^{xz}\f$
  inline virtual std::complex<double> xz(double /*q0*/, double /*q_mag*/) const
    { return std::complex<double>(0., 0.); }

  /// The tensor element \f$W^{y0} = (W^{0y})^*\f$
  inline virtual std::complex<double> yt(double q0, double q_mag) const
    { return std::conj( this->ty(q0, q_mag) ); }

  /// The tensor element \f$W^{yx} = (W^{xy})^*\f$
  inline virtual std::complex<double> yx(double q0, double q_mag) const
    { return std::conj( this->xy(q0, q_mag) ); }

  /// The tensor element \f$W^{yy}\f$
  inline virtual std::complex<double> yy(double q0, double q_mag) const
    { return this->xx(q0, q_mag); }

  /// The tensor element \f$W^{yz}\f$
  inline virtual std::complex<double> yz(double /*q0*/, double /*q_mag*/) const
    { return std::complex<double>(0., 0.); }

  /// The tensor element \f$W^{z0} = (W^{0z})^*\f$
  inline virtual std::complex<double> zt(double q0, double q_mag) const
    { return std::conj( this->tz(q0, q_mag) ); }

  /// The tensor element \f$W^{zx} = (W^{xz})^*\f$
  inline virtual std::complex<double> zx(double q0, double q_mag) const
    { return std::conj( this->xz(q0, q_mag) ); }

  /// The tensor element \f$W^{zy} = (W^{yz})^*\f$
  inline virtual std::complex<double> zy(double q0, double q_mag) const
    { return std::conj( this->yz(q0, q_mag) ); }

  /// @}

  /// \name Structure functions
  /// \param[in] q0 The energy transfer \f$q^0\f$ in the lab frame (GeV)
  /// \param[in] q_mag The magnitude of the 3-momentum transfer
  /// \f$\left|\overrightarrow{q}\right|\f$ in the lab frame (GeV)
  /// \param[in] Mi The mass of the target nucleus \f$M_i\f$ (GeV)
  /// @{

  /// The structure function \f$W_1 = \frac{ W^{xx} }{ 2M_i }\f$
  virtual double W1(double q0, double q_mag, double Mi) const = 0;

  /// The structure function \f$W_2 = \frac{ 1 }{ 2M_i }
  /// \left( W^{00} + W^{xx} + \frac{ (q^0)^2 }
  /// { \left|\overrightarrow{q}\right|^2 } ( W^{zz} - W^{xx} )
  /// - 2\frac{ q^0 }{ \left|\overrightarrow{q}\right| }\Re W^{0z}
  /// \right) \f$
  virtual double W2(double q0, double q_mag, double Mi) const = 0;

  /// The structure function \f$ W_3 = -i \frac{ W^{xy} }
  /// { \left|\overrightarrow{q}\right| } \f$
  virtual double W3(double q0, double q_mag, double Mi) const = 0;

  /// The structure function \f$ W_4 = \frac{ M_i }
  /// { 2 \left|\overrightarrow{q}\right|^2 } ( W^{zz} - W^{xx} ) \f$
  virtual double W4(double q0, double q_mag, double Mi) const = 0;

  /// The structure function \f$ W_5 = \frac{ 1 }
  /// { \left|\overrightarrow{q}\right| }
  /// \left( \Re W^{0z} - \frac{ q^0 }{ \left|\overrightarrow{q}\right| }
  /// ( W^{zz} - W^{xx} ) \right) \f$
  virtual double W5(double q0, double q_mag, double Mi) const = 0;

  /// The structure function \f$ W_6 = \frac{ \Im W^{0z} }
  /// { \left|\overrightarrow{q}\right| } \f$
  virtual double W6(double q0, double q_mag, double Mi) const = 0;
  /// @}

  virtual double contraction(const Interaction* interaction,
    double Q_value) const /*override*/;

  /// \copybrief contraction(const Interaction*)
  /// \param[in] probe_pdg The PDG code for the incident projectile
  /// \param[in] E_probe The lab frame energy of the incident projectile (GeV)
  /// \param[in] m_probe The mass of the incident projectile (GeV)
  /// \param[in] Tl The lab frame kinetic energy of the final state lepton (GeV)
  /// \param[in] cos_l The angle between the direction of the incident
  /// neutrino and the final state lepton (radians)
  /// \param[in] ml The mass of the final state lepton (GeV)
  /// \param[in] Q_value The Q-value that should be used to correct
  /// the energy transfer \f$q_0\f$ (GeV)
  /// \returns The tensor contraction \f$L_{\mu\nu}W^{\mu\nu}\f$ (GeV)
  virtual double contraction(int probe_pdg, double E_probe, double m_probe,
    double Tl, double cos_l, double ml, double Q_value) const;

  /// Computes the differential cross section \f$\frac{ d^2\sigma } { dT_\ell
  /// d\cos(\theta^\prime) }\f$ for the reaction represented by this hadron
  /// tensor
  /// \param[in] interaction An Interaction object storing information about
  /// the initial and final states
  /// \param[in] Q_value The Q-value that should be used to correct
  /// the energy transfer \f$q_0\f$ (GeV)
  /// \returns The differential cross section (GeV<sup>-3</sup>)
  virtual double dSigma_dT_dCosTheta(const Interaction* interaction,
    double Q_value) const = 0;

  /// \copybrief dSigma_dT_dCosTheta(const Interaction*, double)
  /// \param[in] probe_pdg The PDG code for the incident projectile
  /// \param[in] E_probe The lab frame energy of the incident projectile (GeV)
  /// \param[in] m_probe The mass of the incident projectile (GeV)
  /// \param[in] Tl The lab frame kinetic energy of the final state lepton (GeV)
  /// \param[in] cos_l The angle between the direction of the incident
  /// neutrino and the final state lepton (radians)
  /// \param[in] ml The mass of the final state lepton (GeV)
  /// \param[in] Q_value The Q-value that should be used to correct
  /// the energy transfer \f$q_0\f$ (GeV)
  /// \returns The differential cross section (GeV<sup>-3</sup>)
  virtual double dSigma_dT_dCosTheta(int probe_pdg, double E_probe,
    double m_probe, double Tl, double cos_l, double ml, double Q_value)
    const /*override*/ = 0;


  //As above but using the rosenbluth formalism: should give the same result!
  virtual double dSigma_dT_dCosTheta_rosenbluth(const Interaction* interaction,
    double Q_value) const = 0;
  virtual double dSigma_dT_dCosTheta_rosenbluth(int probe_pdg, double E_probe,
    double m_probe, double Tl, double cos_l, double ml, double Q_value)
    const /*override*/ = 0;

protected:

  inline LabFrameHadronTensorI(int pdg = 0) : HadronTensorI(pdg) {}

  inline LabFrameHadronTensorI(int Z, int A) : HadronTensorI(Z, A) {}

}; // class LabFrameHadronTensorI

} // genie namespace
#endif
