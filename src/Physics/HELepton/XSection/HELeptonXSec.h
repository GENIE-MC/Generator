//____________________________________________________________________________
/*!
*/
//____________________________________________________________________________

#ifndef _HE_LEPTON_XSEC_H_
#define _HE_LEPTON_XSEC_H_

#include "Physics/XSectionIntegration/XSecIntegratorI.h"

namespace genie {

class XSecAlgorithmI;
class Interaction;

class HELeptonXSec : public XSecIntegratorI {

public:
  HELeptonXSec();
  HELeptonXSec(string config);
  virtual ~HELeptonXSec();

  //! XSecIntegratorI interface implementation
  double Integrate(const XSecAlgorithmI * model, const Interaction * i) const;

  //! Overload the Algorithm::Configure() methods to load private data
  //! members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void LoadConfig (void);
};

}       // genie namespace
#endif  // _HE_LEPTON_XSEC_H_
