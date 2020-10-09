//____________________________________________________________________________
/*!

\class    genie::GLRESXSec

\brief    nubar + e- scattering glashow resonance. Integrates the loaded
          differential cross section model. An analytical cross section
          model also exists, so you cal also use that if you do not apply
          any kinematical cuts.

          The cross section algorithm handles:
             - nuebar + e- -> nuebar   + e-   [CC + NC + interference]
             - nuebar + e- -> numubar  + mu-  [CC]
             - nuebar + e- -> nutaubar + tau- [CC]
             - nuebar + e- -> hadrons         [CC]

          Is a concrete implementation of the XSecIntegratorI interface. \n

\author   Alfonso Garcia <alfonsog \at nikhef.nl>
          NIKHEF (Amsterdam)

\created  November 8, 2019

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _GLASHOW_RESONANCE_XSEC_H_
#define _GLASHOW_RESONANCE_XSEC_H_

#include "Physics/XSectionIntegration/XSecIntegratorI.h"

namespace genie {

class GLRESXSec : public XSecIntegratorI {

public:
  GLRESXSec();
  GLRESXSec(string config);
  virtual ~GLRESXSec();

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
#endif  // _GLASHOW_RESONANCE_XSEC_H_
