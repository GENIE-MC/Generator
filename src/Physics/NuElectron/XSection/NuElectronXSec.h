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

          Is a concrete implementation of the XSecIntegratorI interface. \n

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  February 10, 2006

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _NU_ELECTRON_XSEC_H_
#define _NU_ELECTRON_XSEC_H_

#include "Physics/XSectionIntegration/XSecIntegratorI.h"

namespace genie {

class NuElectronXSec : public XSecIntegratorI {

public:
  NuElectronXSec();
  NuElectronXSec(string config);
  virtual ~NuElectronXSec();

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
#endif  // _NU_ELECTRON_XSEC_H_
