#ifndef _BYCharmXSec_H
#define _BYCharmXSec_H

#include "Framework/EventGen/XSecAlgorithmI.h"
namespace genie {

class QPMDISPXSec;

class BYCharmXSec : public XSecAlgorithmI {
public:
  BYCharmXSec()  : XSecAlgorithmI("genie::BYCharmXSec") {};
  BYCharmXSec(string config): XSecAlgorithmI("genie::BYCharmXSec", config) {};
  virtual ~BYCharmXSec() = default;

  //-- XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k) const ;
  double Integral        (const Interaction * i) const {return this->fDISModel->Integral(i);};
  bool   ValidProcess    (const Interaction * i) const {return this->fDISModel->ValidProcess(i);};
  bool ValidKinematics (const Interaction * i) const {return this->fDISModel->ValidKinematics(i);};

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

private:

  const XSecAlgorithmI * fDISModel ;
  void LoadConfig();
};
}

#endif