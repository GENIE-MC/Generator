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
#include "Framework/Interaction/InteractionType.h"
#include "Physics/QuasiElastic/XSection/QELFormFactors.h"
#include "Physics/QuasiElastic/XSection/Rank2LorentzTensorI.h"

namespace genie {

class FreeNucleonTensor : public Rank2LorentzTensorI {

public:

  /// Constructor for typical uses
  FreeNucleonTensor(const genie::Interaction& interaction,
    const genie::QELFormFactors& ff);

  /// Constructor to facilitate testing
  FreeNucleonTensor(const TLorentzVector& pNi, const TLorentzVector& q,
    double F1V, double F2V, double F3V, double FA, double FP, double F3A,
    InteractionType_t type, int hit_nuc_pdg);

  inline virtual ~FreeNucleonTensor() {}

  /// Retrieves a tensor element corresponding to the given indices
  virtual std::complex<double> operator()(TensorIndex_t mu, TensorIndex_t nu)
    const /*override*/;

protected:

  TLorentzVector fPNiOnShell; ///< Initial nucleon 4-momentum (forced on-shell)

  /// Momentum transfer 4-vector (corrected for binding energy effects
  /// according to the de Forest prescription)
  TLorentzVector fqTilde;

  double fNucleonMass; ///< Nucleon mass to use when computing tensor elements

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
