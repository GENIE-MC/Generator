//____________________________________________________________________________
/*!

\class    genie::ASKXSec

\brief    Computes the cross section for COH neutrino-nucleus pi production.\n
          Is a concrete implementation of the XSecIntegratorI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  May 04, 2004

\cpright  Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _ASK_XSEC_H_
#define _ASK_XSEC_H_

#include "Base/XSecIntegratorI.h"

#include <Math/Integrator.h>

namespace genie {

class ASKXSec : public XSecIntegratorI {
public:
  ASKXSec();
  ASKXSec(string config);
  virtual ~ASKXSec();

  // XSecIntegratorI interface implementation
  double Integrate(const XSecAlgorithmI * model, const Interaction * i) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

protected:
  bool fSplitIntegral;

private:
  void LoadConfig (void);
};

}       // genie namespace
#endif  // _COH_XSEC_H_
