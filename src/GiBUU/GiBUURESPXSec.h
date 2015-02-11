//____________________________________________________________________________
/*!

\class    genie::GiBUURESPXSec

\brief    Computes the double differential resonance neutrino-production 
          cross section according to the GiBUU model.
          Is a concrete implementation of the XSecAlgorithmI interface.

\ref      T.Leitner, O.Buss, U.Mosel, L.Alvarez-Ruso, 
          Phys. Rev. C 79, 034601 (2009).

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  Jun 03, 2009

\cpright  Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _GIBUU_RES_PXSEC_H_
#define _GIBUU_RES_PXSEC_H_

#include "Base/XSecAlgorithmI.h"

namespace genie {

class XSecIntegratorI;

class GiBUURESPXSec : public XSecAlgorithmI {

public:
  GiBUURESPXSec();
  GiBUURESPXSec(string config);
  virtual ~GiBUURESPXSec();

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

  // configuration data

  double   fMa2;               ///< (axial mass)^2
  double   fMv2;               ///< (vector mass)^2
  bool     fUsingDisResJoin;   ///< use a DIS/RES joining scheme?
  double   fWcut;              ///< apply DIS/RES joining scheme < Wcut

  const XSecIntegratorI * fXSecIntegrator;
};

}       // genie namespace
#endif  // _GIBUU_RES_PXSEC_H_
