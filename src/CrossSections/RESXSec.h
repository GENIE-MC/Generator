//____________________________________________________________________________
/*!

\class    genie::RESXSec

\brief    Computes the RES Cross Section.\n
          Is a concrete implementation of the XSecAlgorithmI interface.\n

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 04, 2004

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _RES_XSEC_H_
#define _RES_XSEC_H_

#include "Base/XSecAlgorithmI.h"
#include "Utils/Range1.h"

namespace genie {

class IntegratorI;

class RESXSec : public XSecAlgorithmI {

public:
  RESXSec();
  RESXSec(string param_set);
  virtual ~RESXSec();

  //-- XSecAlgorithmI interface implementation
  double XSec            (const Interaction * interaction) const;
  bool   ValidProcess    (const Interaction * interaction) const;
  bool   ValidKinematics (const Interaction * interaction) const;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void      LoadConfig (void);
  Range1D_t WRange     (const Interaction * interaction) const;
  Range1D_t Q2Range    (const Interaction * interaction) const;

  const XSecAlgorithmI * fPartialXSecAlg;
  const IntegratorI *    fIntegrator;
};

}       // genie namespace
#endif  // _RES_XSEC_H_
