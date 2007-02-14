//____________________________________________________________________________
/*!

\class    genie::IMDXSec

\brief    Computes the Inverse Muon Decay cross section.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  Fabruary 14, 2005

\cpright  Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _IMD_XSEC_H_
#define _IMD_XSEC_H_

#include "Base/XSecIntegratorI.h"

namespace genie {

class IMDXSec : public XSecIntegratorI {

public:
  IMDXSec();
  IMDXSec(string config);
  virtual ~IMDXSec();

  //! XSecIntegratorI interface implementation
  double Integrate(const XSecAlgorithmI * model, const Interaction * i) const;

  //! overload the Algorithm::Configure() methods to load private data
  //! members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void LoadConfig (void);
};

}       // genie namespace
#endif  // _IMD_XSEC_H_
