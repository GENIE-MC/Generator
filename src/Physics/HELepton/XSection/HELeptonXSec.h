//____________________________________________________________________________
/*!

\class    genie::HELeptonXSec

\brief    Total cross section integrator for neutrino-electron

\author   Alfonso Garcia <aagarciasoto \at km3net.de>
          IFIC & Harvard University

\ref      Phys. Rev. D 100, 091301 (2019)

\created  Dec 8, 2021

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
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
