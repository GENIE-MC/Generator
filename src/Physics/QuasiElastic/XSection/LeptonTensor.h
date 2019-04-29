//____________________________________________________________________________
/*!

\class    genie::LeptonTensor

\brief    Computes the elements of the lepton tensor \f$L_{\mu\nu}\f$

\author   Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  January 17, 2019

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef LEPTON_TENSOR_H
#define LEPTON_TENSOR_H

// GENIE includes
#include "Physics/QuasiElastic/XSection/Rank2LorentzTensorI.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Interaction/InteractionType.h"

namespace genie {

class LeptonTensor : public Rank2LorentzTensorI {

public:

  LeptonTensor(const genie::Interaction& interaction);

  inline virtual ~LeptonTensor() {}

  inline const TLorentzVector& GetProbeP4() const { return fProbeP4; }
  inline const TLorentzVector& GetFSLeptonP4() const { return fFSLepP4; }

  // Override the Contract() method of Rank2LorentzTensorI to add in a
  // correction for conservation of the EM current (see Phys. Rev. C 77, 044311
  // (2008) eqs. (4)-(5) and accompanying text for an explanation) when this
  // tensor is contracted with the FreeNucleonTensor.
  virtual std::complex<double> Contract(const Rank2LorentzTensorI& other) const
    /*override*/;

  virtual std::complex<double> operator()(genie::TensorIndex_t mu,
    genie::TensorIndex_t nu) const /*override*/;

protected:

  std::complex<double> LeviCivitaProduct(genie::TensorIndex_t mu,
    genie::TensorIndex_t nu) const;

  TLorentzVector fProbeP4; ///< Initial lepton 4-momentum
  TLorentzVector fFSLepP4; ///< Final lepton 4-momentum
  int fInitialLeptonPDG; ///< PDG code for the initial lepton
  /// Kind of interaction represented by this tensor (EM, CC, NC, etc.)
  genie::InteractionType_t fInteractionType;
  double fMLep2; ///< Square of the initial lepton mass

}; // class LeptonTensor

} // genie namespace
#endif
