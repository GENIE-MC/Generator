#ifndef _LLEWELLYN_SMITH_QELCC_CROSS_SECTION_H_
#define _LLEWELLYN_SMITH_QELCC_CROSS_SECTION_H_

#include "Framework/EventGen/XSecAlgorithmI.h"

#include "Tools/VMC/Key.h"

class TGraph;

namespace genie {
namespace vmc {

class EvtLibPXSec : public XSecAlgorithmI {

public:
  EvtLibPXSec();
  EvtLibPXSec(string config);
  virtual ~EvtLibPXSec();

  // XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k) const;
  double Integral        (const Interaction * i) const;
  bool   ValidProcess    (const Interaction * i) const;

  // Override the Algorithm::Configure methods to load configuration
  // data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

protected:
  TGraph* GetXSec(const Interaction* in) const;
  void LoadXSecs() const;

  mutable std::map<Key, TGraph*> fXSecs;
};

} // vmc namespace
} // genie namespace

#endif
