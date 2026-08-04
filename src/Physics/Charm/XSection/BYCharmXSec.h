#ifndef _BYCharmXSec_H
#define _BYCharmXSec_H

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Physics/DeepInelastic/XSection/DISStructureFunc.h"
namespace genie {


class BYCharmXSec : public XSecAlgorithmI {
public:
  BYCharmXSec()  : XSecAlgorithmI("genie::BYCharmXSec") {};
  BYCharmXSec(string config): XSecAlgorithmI("genie::BYCharmXSec", config) {};
  virtual ~BYCharmXSec() = default;

  //-- XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k) const ;
  double Integral        (const Interaction * i) const;
  bool   ValidProcess    (const Interaction * i) const;

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

private:
  bool initilized = false;
  const DISStructureFuncModelI * fDISSFModel ;
  const XSecIntegratorI *        fXSecIntegrator;     ///< diff. xsec integrator
  mutable DISStructureFunc fDISSF;
  void LoadConfig();
  double fCCScale;            ///< cross section scaling factor
  mutable double fF1;
  mutable double fF2;
  mutable double fF3;
  mutable double fF4;
  mutable double fF5;
  mutable double fF6;
  bool fUse2016Corrections;
};
}

#endif