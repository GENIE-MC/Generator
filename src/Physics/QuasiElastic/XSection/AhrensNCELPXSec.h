//____________________________________________________________________________
/*!

\class    genie::AhrensNCELPXSec

\brief    Differential cross section for v+N / vbar+N elastic scattering. \n
          Is a concrete implementation of the XSecAlgorithmI interface. \n

\ref      R.E.Hendrick and L.Li, Phys.Rev.D 19:779 (1979)
          L.A.Ahrens et al., Phys.Rev.D 35:785 (1987)

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  Fabruary 15, 2005

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _AHRENS_NCEL_CROSS_SECTION_H_
#define _AHRENS_NCEL_CROSS_SECTION_H_

#include "Framework/EventGen/XSecAlgorithmI.h"

namespace genie {

class XSecIntegratorI;

class AhrensNCELPXSec : public XSecAlgorithmI {

public:
  AhrensNCELPXSec();
  AhrensNCELPXSec(string config);
  virtual ~AhrensNCELPXSec();

  // XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k) const;
  double Integral        (const Interaction * i) const;
  bool   ValidProcess    (const Interaction * i) const;

  // Override the Algorithm::Configure methods to load configuration
  // data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

private:
  void LoadConfig(void);

  const XSecIntegratorI * fXSecIntegrator;

  double fkAlpha;
  double fkGamma; 
  double fEta;
  double fFa0;
  double fMa2;
  double fMv2;
  double fMuP;
  double fMuN;
};

}       // genie namespace

#endif  
