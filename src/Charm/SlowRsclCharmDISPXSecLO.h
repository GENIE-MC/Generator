//____________________________________________________________________________
/*!

\class    genie::SlowRsclCharmDISPXSecLO

\brief    Computes, at Leading Order (LO), the differential cross section for
          neutrino charm production using a slow rescaling model.

          The computed cross section is the D2xsec = d^2(xsec) / dy dx \n
          where \n
            \li \c y is the inelasticity, and
            \li \c x is the Bjorken scaling variable \c x

          Is a concrete implementation of the XSecAlgorithmI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  March 17, 2005

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _SLOW_RESCALING_CHARM_DIS_PARTIAL_XSEC_LO_H_
#define _SLOW_RESCALING_CHARM_DIS_PARTIAL_XSEC_LO_H_

#include "Base/XSecAlgorithmI.h"

namespace genie {

class PDFModelI;

class SlowRsclCharmDISPXSecLO : public XSecAlgorithmI {

public:

  SlowRsclCharmDISPXSecLO();
  SlowRsclCharmDISPXSecLO(string config);
  virtual ~SlowRsclCharmDISPXSecLO();

  //-- XSecAlgorithmI interface implementation
  double XSec            (const Interaction * interaction) const;
  bool   ValidProcess    (const Interaction * interaction) const;
  bool   ValidKinematics (const Interaction * interaction) const;

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

private:

  void LoadConfig (void);

  const PDFModelI * fPDFModel;

  bool   fDContributes;
  bool   fSContributes;
  double fMc;
  double fVcd;
  double fVcs;
  double fMc2;
  double fVcd2;
  double fVcs2;
};

}       // genie namespace

#endif  // _SLOW_RESCALING_CHARM_DIS_PARTIAL_XSEC_LO_H_
