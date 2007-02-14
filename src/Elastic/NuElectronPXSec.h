//____________________________________________________________________________
/*!

\class    genie::NuElectronPXSec

\brief    nu/nubar + e- scattering differential cross section \n
          The cross section algorithm handles:
             - nue/nuebar + e- -> nue/nuebar + e- [CC + NC + interference]
             - numu/nutau + e- -> numu/nutau + e- [NC]
             - numubar/nutaubar + e- -> numubar/nutaubar + e- [NC]
             - numu/nutau + e- -> l- + nu_e [CC]

          Is a concrete implementation of the XSecAlgorithmI interface. \n

\ref      W.Greiner and B.Muller, Gauge Theory of Weak  Interactions, Springer
          F.Halzen and A.Martin, Quarks and Leptons, J.Wiley

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  February 10, 2006

\cpright  Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _NU_ELECTRON_PARTIAL_XSEC_H_
#define _NU_ELECTRON_PARTIAL_XSEC_H_

#include "Base/XSecAlgorithmI.h"

namespace genie {

class IntegratorI;
class XSecIntegratorI;

class NuElectronPXSec : public XSecAlgorithmI {

public:
  NuElectronPXSec();
  NuElectronPXSec(string config);
  virtual ~NuElectronPXSec();

  //-- XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k) const;
  double Integral        (const Interaction * i) const;
  bool   ValidProcess    (const Interaction * i) const;
  bool   ValidKinematics (const Interaction * i) const;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void LoadConfig (void);

  const XSecIntegratorI * fXSecIntegrator;

  double fCv;
  double fCa;
};

}       // genie namespace
#endif  // _NU_ELECTRON_PARTIAL_XSEC_H_
