//____________________________________________________________________________
/*!

\class    genie::COHXSec

\brief    Computes the cross section for COH neutrino-nucleus scattering.\n
          Is a concrete implementation of the XSecIntegratorI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  May 04, 2004

\cpright  Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _COH_XSEC_H_
#define _COH_XSEC_H_

#include "Base/XSecIntegratorI.h"

namespace genie {

class COHXSec : public XSecIntegratorI {
public:
  COHXSec();
  COHXSec(string config);
  virtual ~COHXSec();

  //! XSecIntegratorI interface implementation
  double Integrate(const XSecAlgorithmI * model, const Interaction * i) const;

  //! overload the Algorithm::Configure() methods to load private data
  //!  members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void LoadConfig (void);
};

}       // genie namespace
#endif  // _COH_XSEC_H_
