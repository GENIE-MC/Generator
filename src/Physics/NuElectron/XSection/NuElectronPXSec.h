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

\author   Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
          University of Liverpool & STFC Rutherford Appleton Laboratory

          Changes required to implement the Electron Velocity module
          were installed by Brinden Carlson (Univ. of Florida)

\created  February 10, 2006

\cpright  Copyright (c) 2003-2023, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _NU_ELECTRON_PARTIAL_XSEC_H_
#define _NU_ELECTRON_PARTIAL_XSEC_H_

#include "Physics/NuElectron/XSection/PXSecOnElectron.h"

namespace genie {

class NuElectronPXSec : public PXSecOnElectron {

public:
  NuElectronPXSec();
  NuElectronPXSec(string config);
  virtual ~NuElectronPXSec();

  //-- PXSecOnElectron interface implementation
  double XSec (const Interaction * i, KinePhaseSpace_t k) const override;

protected:
  void LoadConfig (void) override;
private:

  double fSin28w; // sin^2(theta-weinberg)
  double fSin48w;
};

}       // genie namespace
#endif  // _NU_ELECTRON_PARTIAL_XSEC_H_
