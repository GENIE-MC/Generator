//____________________________________________________________________________
/*!

\class    genie::IMDXSec

\brief    Computes the Inverse Muon Decay cross section.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  Fabruary 14, 2005

*/
//____________________________________________________________________________

#ifndef _IMD_XSEC_H_
#define _IMD_XSEC_H_

#include "Base/XSecAlgorithmI.h"
#include "Numerical/IntegratorI.h"

namespace genie {

class IntegratorI;

class IMDXSec : public XSecAlgorithmI {

public:

  IMDXSec();
  IMDXSec(string config);
  virtual ~IMDXSec();

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

  int   fNBins;
  const XSecAlgorithmI * fDiffXSecModel;
  const IntegratorI *    fIntegrator;
};

}       // genie namespace

#endif  // _IMD_XSEC_H_
