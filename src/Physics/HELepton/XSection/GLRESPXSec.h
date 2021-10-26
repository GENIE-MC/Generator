//____________________________________________________________________________
/*!
*/
//____________________________________________________________________________

#ifndef _GLASHOW_RESONANCE_PXSEC_H_
#define _GLASHOW_RESONANCE_PXSEC_H_

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Physics/HELepton/XSection/Born.h"

namespace genie {

class XSecIntegratorI;

class GLRESPXSec : public XSecAlgorithmI {

public:
  GLRESPXSec ();
  GLRESPXSec (string config);
  virtual ~GLRESPXSec ();

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

  double fWmin;            ///< Minimum value of W

  Born * born;

};

}       // genie namespace

#endif  // _GLASHOW_RESONANCE_XSEC_H_
