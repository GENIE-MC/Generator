//____________________________________________________________________________
/*!
*/
//____________________________________________________________________________

#ifndef _HE_NUELECTRON_PXSEC_H_
#define _HE_NUELECTRON_PXSEC_H_

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Physics/HELepton/XSection/Born.h"

namespace genie {

class XSecIntegratorI;

class HENuElectronPXSec : public XSecAlgorithmI {

public:
  HENuElectronPXSec ();
  HENuElectronPXSec (string config);
  virtual ~HENuElectronPXSec ();

  // XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k) const;
  double Integral        (const Interaction * i) const;
  bool   ValidProcess    (const Interaction * i) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void   LoadConfig (void);

  const XSecIntegratorI *        fXSecIntegrator;     ///< diff. xsec integrator

  Born * born;

};

}       // genie namespace

#endif  // _HE_NUELECTRON_PXSEC_H_
