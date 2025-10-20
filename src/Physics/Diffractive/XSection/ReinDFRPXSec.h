//____________________________________________________________________________
/*!

\class    genie::ReinDFRPXSec

\brief    Neutrino diffractive pion production cross section.

\ref      D.Rein, Nucl.Phys.B278(1986) 61-77

\author   Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  Feb 17th, 2009

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          
*/
//____________________________________________________________________________

#ifndef _REIN_DFRC_PXSEC_H_
#define _REIN_DFRC_PXSEC_H_

#include <string>

#include "Framework/EventGen/XSecAlgorithmI.h"

namespace genie {

class XSecIntegratorI;

class ReinDFRPXSec : public XSecAlgorithmI {

public:
  ReinDFRPXSec();
  ReinDFRPXSec(const std::string & config);
  virtual ~ReinDFRPXSec();

  // XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k) const;
  double Integral        (const Interaction * i) const;
  bool   ValidProcess    (const Interaction * i) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void LoadConfig (void);

  double fMa;      ///< axial mass
  double fBeta;    ///< b in dsig{piN}/dt = dsig0{piN}/dt * exp(-b(t-tmin)), b ~ 0.333 (nucleon_size)^2

  const XSecIntegratorI * fXSecIntegrator;
};

}       // genie namespace
#endif  // _REIN_DFRC_PXSEC_H_
