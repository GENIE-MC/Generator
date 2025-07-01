//____________________________________________________________________________
/*!

\class    genie::HybridXSecAlgorithm

\brief    Defines an XSecAlgorithmI that delegates the actual calculation
          to one or more sub-algorithms (each of which is itself an XSecAlgorithmI)
          based on the input interaction. The choice of sub-algorithms is configurable
          via XML.

          The current use case is to allow GENIE to simultaneously use the
          SuSAv2 quasielastic cross section model for complex targets and the
          Llewellyn-Smith model for free nucleons.

          Is a concrete implementation of the XSecAlgorithmI interface.

\author   Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Acclerator Laboratory

\created  November 4, 2019

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _HYBRID_XSEC_ALG_H_
#define _HYBRID_XSEC_ALG_H_

#include "Framework/EventGen/XSecAlgorithmI.h"

namespace genie {

class HybridXSecAlgorithm : public XSecAlgorithmI {

public:

  HybridXSecAlgorithm();
  HybridXSecAlgorithm(string config);
  virtual ~HybridXSecAlgorithm();

  // XSecAlgorithmI interface implementation
  double XSec(const Interaction* i, KinePhaseSpace_t k) const;
  double Integral(const Interaction* i) const;
  bool   ValidProcess(const Interaction* i) const;

  // override the Algorithm::Configure methods to load configuration
  // data to private data members
  void Configure (const Registry & config);
  void Configure (string config);

private:

  /// Load algorithm configuration
  void LoadConfig (void);

  /// Retrieve a pointer to the appropriate cross section algorithm
  /// using the map. If no suitable algorithm was found, return a
  /// null pointer.
  const XSecAlgorithmI* ChooseXSecAlg(const Interaction& interaction) const;

  /// Map specifying the managed cross section algorithms. Keys are strings
  /// generated with Interaction::AsString() (identical to those used for
  /// splines). Values are pointers to the corresponding cross section
  /// algorithms.
  mutable std::map<string, const XSecAlgorithmI*> fXSecAlgMap;

  /// Optional XSecAlgorithmI to use by default
  const XSecAlgorithmI* fDefaultXSecAlg;
};

} // genie namespace
#endif // _HYBRID_XSEC_ALG_H_
