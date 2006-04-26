//____________________________________________________________________________
/*!

\class    genie::NuNucElasticPXSec

\brief    Differential cross section dxsec/dQ^2 for v+N / vbar+N elastic 
          scattering. \n
          NuNucElasticPXSec is a concrete implementation of the
          XSecAlgorithmI interface. \n

\ref      L.A.Ahrens et al., Physical Review D, VOL 35,3:785 (1987)

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  Fabruary 15, 2005

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _NU_NUC_ELASTIC_PARTIAL_XSEC_H_
#define _NU_NUC_ELASTIC_PARTIAL_XSEC_H_

#include "Base/XSecAlgorithmI.h"

namespace genie {

class NuNucElasticPXSec : public XSecAlgorithmI {

public:
  NuNucElasticPXSec();
  NuNucElasticPXSec(string config);
  virtual ~NuNucElasticPXSec();

  //-- XSecAlgorithmI interface implementation
  double XSec            (const Interaction * interaction) const;
  bool   ValidProcess    (const Interaction * interaction) const;
  bool   ValidKinematics (const Interaction * interaction) const;

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

private:

  void LoadConfig(void);

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
