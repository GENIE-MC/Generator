//____________________________________________________________________________
/*!

\class    genie::NuclearTargetXSec

\brief    Computes the cross section for a nuclear target from the free
          nucleon cross sections.

          Is a concrete implementation of the XSecAlgorithmI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  June 10, 2004

*/
//____________________________________________________________________________

#ifndef _NUCLEAR_TARGET_XSEC_H_
#define _NUCLEAR_TARGET_XSEC_H_

#include "Base/XSecAlgorithmI.h"

namespace genie {

class NuclearTargetXSec : public XSecAlgorithmI {

public:

  NuclearTargetXSec();
  NuclearTargetXSec(string config);
  virtual ~NuclearTargetXSec();

  //-- XSecAlgorithmI interface implementation
  double XSec            (const Interaction * interaction) const;
  bool   ValidProcess    (const Interaction * interaction) const;
  bool   ValidKinematics (const Interaction * interaction) const;

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

private:

  void LoadSubAlg(void);

  const XSecAlgorithmI * fFreeParticleXSecAlg;
};

}       // genie namespace

#endif  // _NUCLEAR_TARGET_XSEC_H_
