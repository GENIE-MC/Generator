#ifndef _BY2021CHARM_H
#define _BY2021CHARM_H

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Physics/DeepInelastic/XSection/QPMDISPXSec.h"
namespace genie {

class PDFModelI;
class XSecIntegratorI;

class BYCharmXSec : public XSecAlgorithmI {
public:
  BYCharmXSec()  : XSecAlgorithmI("genie::BYCharmXSec") {};
  BYCharmXSec(string config): XSecAlgorithmI("genie::BYCharmXSec", config) {};
  virtual ~BYCharmXSec() = default;

  //-- XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k) const;
  double Integral        (const Interaction * i) const;
  bool   ValidProcess    (const Interaction * i) const;

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

private:

  const QPMDISPXSec * fDISModel ;
  void LoadConfig();
};
}

#endif