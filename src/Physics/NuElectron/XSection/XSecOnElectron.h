//____________________________________________________________________________
/*!

\class    genie::XSecOnElectron

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

\author   Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
          University of Liverpool & STFC Rutherford Appleton Laboratory

          B. Carlson change name to reflect that cross section applies
          to all nu+e interactions

\created  February 10, 2006

\cpright  Copyright (c) 2003-2023, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _X_Sec_On_Electron_H_
#define _X_Sec_On_Electron_H_

#include "Physics/XSectionIntegration/XSecIntegratorI.h"

namespace genie {

class XSecOnElectron : public XSecIntegratorI {

public:
  XSecOnElectron();
  XSecOnElectron(string config);
  virtual ~XSecOnElectron();

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
#endif  // _X_Sec_On_Electron_H_
