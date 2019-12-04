//____________________________________________________________________________
/*!

\class    genie::FreeNucleonTensor

\brief    Computes the elements of the free nucleon tensor \f$A^{\mu\nu}\f$

\author   Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  January 17, 2019

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef FREE_NUCLEON_TENSOR_H
#define FREE_NUCLEON_TENSOR_H

// GENIE includes
#include "Framework/Interaction/Interaction.h"
#include "Physics/QuasiElastic/XSection/QELFormFactors.h"
#include "Physics/QuasiElastic/XSection/Rank2LorentzTensorI.h"

namespace genie {

class FreeNucleonTensor : public Rank2LorentzTensorI {

public:

  FreeNucleonTensor(const genie::Interaction& interaction,
    const genie::QELFormFactors& ff);

  inline virtual ~FreeNucleonTensor() {}

  // Override the Contract() method of Rank2LorentzTensorI to add in a
  // correction for conservation of the EM current (see Phys. Rev. C 77, 044311
  // (2008) eqs. (4)-(5) and accompanying text for an explanation)
  virtual std::complex<double> Contract(const Rank2LorentzTensorI& other) const
    /*override*/;

  /// Retrieves a tensor element corresponding to the given indices
  virtual std::complex<double> operator()(TensorIndex_t mu, TensorIndex_t nu)
    const /*override*/;

protected:

  TLorentzVector fPNiOnShell; ///< Initial nucleon 4-momentum (forced on-shell)

  /// Momentum transfer 4-vector (corrected for binding energy effects
  /// according to the de Forest prescription)
  TLorentzVector fqTilde;

  double fNucleonMass; ///< Nucleon mass to use when computing tensor elements

  /// Momentum-dependent binding energy of the initial hit nucleon
  double fEpsilonB;

  /// Interaction type (EM, CC, NC) to assume when computing the tensor
  /// elements
  InteractionType_t fInteractionType;

  /// \name Form factors
  // TODO: reference conventions!
  /// @{
  double fF1V;
  double fF2V;
  double fF3V;
  double fFA;
  double fFP;
  double fF3A;
  /// @}

}; // class FreeNucleonTensor

} // genie namespace
#endif
