//____________________________________________________________________________
/*!

\class    genie::NuNucElasticPXSec

\brief    Differential cross section for v+N / vbar+N elastic scattering. \n
          Is a concrete implementation of the XSecAlgorithmI interface. \n

\ref      R.E.Hendrick and L.Li, Phys.Rev.D 19:779 (1979)
          L.A.Ahrens et al., Phys.Rev.D 35:785 (1987)

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  Fabruary 15, 2005

\cpright  Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _NU_NUC_ELASTIC_PARTIAL_XSEC_H_
#define _NU_NUC_ELASTIC_PARTIAL_XSEC_H_

#include "Base/XSecAlgorithmI.h"

namespace genie {

class XSecIntegratorI;

class NuNucElasticPXSec : public XSecAlgorithmI {

public:
  NuNucElasticPXSec();
  NuNucElasticPXSec(string config);
  virtual ~NuNucElasticPXSec();

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
#endif  // _NU_NUC_ELASTIC_PARTIAL_XSEC_H_
