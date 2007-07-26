//____________________________________________________________________________
/*!

\class    genie::QELPXSec

\brief    Computes neutrino-nucleon(nucleus) QEL differential cross sections
          Is a concrete implementation of the XSecAlgorithmI interface. \n

\ref      C.H.Llewellyn Smith, Physics Reports (Sect. C of Physics Letters) 3,
          No. 5  (1972) 261-379

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  May 05, 2004

\cpright  Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _QEL_PARTIAL_XSEC_H_
#define _QEL_PARTIAL_XSEC_H_

#include "Base/XSecAlgorithmI.h"
#include "Base/QELFormFactors.h"

namespace genie {

class QELFormFactorsModelI;
class XSecIntegratorI;

class QELPXSec : public XSecAlgorithmI {

public:
  QELPXSec();
  QELPXSec(string config);
  virtual ~QELPXSec();

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

  mutable QELFormFactors fFormFactors;

  const QELFormFactorsModelI * fFormFactorsModel;
  const XSecIntegratorI *      fXSecIntegrator;

  double fCos8c2; ///< cos^2(cabbibo angle)
};

}       // genie namespace
#endif  // _QEL_PARTIAL_XSEC_H_
