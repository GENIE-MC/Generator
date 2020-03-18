#ifndef _LLEWELLYN_SMITH_QELCC_CROSS_SECTION_H_
#define _LLEWELLYN_SMITH_QELCC_CROSS_SECTION_H_

#include "Physics/NuclearState/NuclearModelI.h"
#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Physics/QuasiElastic/XSection/QELFormFactors.h"
#include "Physics/QuasiElastic/XSection/QELUtils.h"
#include "Physics/NuclearState/PauliBlocker.h"

namespace genie {
namespace vmc {

class QELFormFactorsModelI;
class XSecIntegratorI;

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
};

} // vmc namespace
} // genie namespace

#endif
