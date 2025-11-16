//____________________________________________________________________________
/*!

\class    genie::PhotonCOHPXSec

\brief    Differential cross section for W boson production

\author   Alfonso Garcia <aagarciasoto \at km3net.de>
          IFIC & Harvard University

\ref      Phys. Rev. D 100, 091301 (2019)

\created  Dec 8, 2021

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _PHOTON_COH_PXSEC_H_
#define _PHOTON_COH_PXSEC_H_

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Physics/HELepton/XSection/Born.h"

namespace genie {

class XSecIntegratorI;

class PhotonCOHPXSec : public XSecAlgorithmI {

public:
  PhotonCOHPXSec ();
  PhotonCOHPXSec (string config);
  virtual ~PhotonCOHPXSec ();

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
  double F2_Q       (double Q, double r0) const; //EM Nuclear Form-factor

  
  const XSecIntegratorI *        fXSecIntegrator;     ///< diff. xsec integrator

  Born * born;

};

}       // genie namespace

#endif  // _PHOTON_COH_PXSEC_H_