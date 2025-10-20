//____________________________________________________________________________
/*!

\class    genie::PhotonRESPXSec

\brief    Differential cross section for trident production

\author   Alfonso Garcia <aagarciasoto \at km3net.de>
          IFIC & Harvard University

\ref      Phys. Rev. D 100, 091301 (2019)

\created  Dec 8, 2021

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _PHOTON_RES_PXSEC_H_
#define _PHOTON_RES_PXSEC_H_

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Physics/HELepton/XSection/Born.h"

namespace genie {

class XSecIntegratorI;

class PhotonRESPXSec : public XSecAlgorithmI {

public:
  PhotonRESPXSec ();
  PhotonRESPXSec (string config);
  virtual ~PhotonRESPXSec ();

  // XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k) const;
  double Integral        (const Interaction * i) const;
  bool   ValidProcess    (const Interaction * i) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void   LoadConfig (void);

  const XSecIntegratorI *        fXSecIntegrator;     ///< diff. xsec integrator

  double fWmin;            ///< Minimum value of W

  double fQ2PDFmin;
  double fxPDFmin;

  Born * born;


};

}       // genie namespace

#endif  // _PHOTON_RES_PXSEC_H_
