//____________________________________________________________________________
/*!

\class    genie::NuElectronXSec

\brief    nu/nubar + e- scattering cross section. Integrates the loaded
          differential cross section model. An analytical cross section
          model also exists, so you cal also use that if you do not apply
          any kinematical cuts.

          The cross section algorithm handles:
             - nue/nuebar + e- -> nue/nuebar + e- [CC + NC + interference]
             - numu/nutau + e- -> numu/nutau + e- [NC]
             - numubar/nutaubar + e- -> numubar/nutaubar + e- [NC]
             - numu/nutau + e- -> l- + nu_e [CC]

          NuElectronXSec is a concrete implementation of the XSecAlgorithmI
          interface. \n

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  February 10, 2006

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _NU_ELECTRON_XSEC_H_
#define _NU_ELECTRON_XSEC_H_

#include "Base/XSecAlgorithmI.h"

namespace genie {

class IntegratorI;

class NuElectronXSec : public XSecAlgorithmI {

public:
  NuElectronXSec();
  NuElectronXSec(string config);
  virtual ~NuElectronXSec();

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

  const XSecAlgorithmI * fDiffXSecModel;
  const IntegratorI *    fIntegrator;
};

}       // genie namespace
#endif  // _NU_ELECTRON_XSEC_H_
