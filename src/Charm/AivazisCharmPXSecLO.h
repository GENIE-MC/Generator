//____________________________________________________________________________
/*!

\class    genie::AivazisCharmPXSecLO

\brief    Computes, at Leading Order (LO), the differential cross section for
          neutrino charm production using the \b Aivazis,Olness,Tung model.
          Is a concrete implementation of the XSecAlgorithmI interface.

\ref      M.A.G.Aivazis, F.I.Olness and W.K.Tung
          Phys.Rev.D50, 3085-3101 (1994)

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  June 10, 2004

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _AIVAZIS_CHARM_PARTIAL_XSEC_LO_H_
#define _AIVAZIS_CHARM_PARTIAL_XSEC_LO_H_

#include "Base/XSecAlgorithmI.h"

namespace genie {

class PDFModelI;

class AivazisCharmPXSecLO : public XSecAlgorithmI {

public:
  AivazisCharmPXSecLO();
  AivazisCharmPXSecLO(string config);
  virtual ~AivazisCharmPXSecLO();

  //-- XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k) const;
  bool   ValidProcess    (const Interaction * i) const;
  bool   ValidKinematics (const Interaction * i) const;

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

private:

  void LoadConfig(void);

  const PDFModelI* fPDFModel;

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

#endif  // _AIVAZIS_CHARM_PARTIAL_XSEC_LO_H_
