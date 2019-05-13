//____________________________________________________________________________
/*!

\class    genie::AivazisCharmPXSecLO

\brief    Computes, at Leading Order (LO), the differential cross section for
          neutrino charm production using the \b Aivazis,Olness,Tung model.
          Is a concrete implementation of the XSecAlgorithmI interface.

\ref      M.A.G.Aivazis, F.I.Olness and W.K.Tung
          Phys.Rev.D50, 3085-3101 (1994)

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  June 10, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _AIVAZIS_CHARM_PARTIAL_XSEC_LO_H_
#define _AIVAZIS_CHARM_PARTIAL_XSEC_LO_H_

#include "Framework/EventGen/XSecAlgorithmI.h"

namespace genie {

class PDFModelI;
class XSecIntegratorI;

class AivazisCharmPXSecLO : public XSecAlgorithmI {

public:
  AivazisCharmPXSecLO();
  AivazisCharmPXSecLO(string config);
  virtual ~AivazisCharmPXSecLO();

  //-- XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k) const;
  double Integral        (const Interaction * i) const;
  bool   ValidProcess    (const Interaction * i) const;

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

private:

  void LoadConfig(void);

  const PDFModelI *       fPDFModel;
  const XSecIntegratorI * fXSecIntegrator;

  //bool   fDContributes;
  //bool   fSContributes;
  double fMc;
  double fVcd;
  double fVcs;
  double fMc2;
  double fVcd2;
  double fVcs2;
};

}       // genie namespace
#endif  // _AIVAZIS_CHARM_PARTIAL_XSEC_LO_H_
