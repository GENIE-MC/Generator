//____________________________________________________________________________
/*!

\class    genie::COHXSec

\brief    Computes the cross section for COH neutrino-nucleus scattering.

          Is a concrete implementation of the XSecAlgorithmI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 04, 2004

*/
//____________________________________________________________________________

#ifndef _COH_XSEC_H_
#define _COH_XSEC_H_

#include "Base/XSecAlgorithmI.h"

namespace genie {

class IntegratorI;

class COHXSec : public XSecAlgorithmI {

public:

  COHXSec();
  COHXSec(string config);
  virtual ~COHXSec();

  //-- XSecAlgorithmI interface implementation
  double XSec            (const Interaction * interaction) const;
  bool   ValidProcess    (const Interaction * interaction) const;
  bool   ValidKinematics (const Interaction * interaction) const;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void LoadConfig (void);

  const XSecAlgorithmI * fPartialXSecAlg;
  const IntegratorI *    fIntegrator;
};

}       // genie namespace
#endif  // _COH_XSEC_H_
