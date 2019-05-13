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

\ref      W.J.Marciano and Z.Parsa, Neutrino-electron scattering theory,
          J.Phys.G: Nucl.Part.Phys. 29 (2003) 2629-2645

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  February 10, 2006

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _NU_ELECTRON_PARTIAL_XSEC_H_
#define _NU_ELECTRON_PARTIAL_XSEC_H_

#include "Framework/EventGen/XSecAlgorithmI.h"

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

  double fSin28w; // sin^2(theta-weinberg)
  double fSin48w; 
};

}       // genie namespace
#endif  // _NU_ELECTRON_PARTIAL_XSEC_H_
