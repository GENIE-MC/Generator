//____________________________________________________________________________
/*!

\class    genie::SlowRsclCharmDISPXSecLO

\brief    Computes, at Leading Order (LO), the differential cross section for
          neutrino charm production using a slow rescaling model.
          Is a concrete implementation of the XSecAlgorithmI interface.

\author   Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
          University of Liverpool & STFC Rutherford Appleton Laboratory

\created  March 17, 2005

\cpright  Copyright (c) 2003-2020, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _SLOW_RESCALING_CHARM_DIS_PARTIAL_XSEC_LO_H_
#define _SLOW_RESCALING_CHARM_DIS_PARTIAL_XSEC_LO_H_

#include "Framework/EventGen/XSecAlgorithmI.h"

namespace genie {

class PDFModelI;
class XSecIntegrator;

class SlowRsclCharmDISPXSecLO : public XSecAlgorithmI {

public:
  SlowRsclCharmDISPXSecLO();
  SlowRsclCharmDISPXSecLO(string config);
  virtual ~SlowRsclCharmDISPXSecLO();

  //-- XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k) const;
  double Integral        (const Interaction * i) const;
  bool   ValidProcess    (const Interaction * i) const;

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

private:

  void LoadConfig (void);

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
#endif  // _SLOW_RESCALING_CHARM_DIS_PARTIAL_XSEC_LO_H_
