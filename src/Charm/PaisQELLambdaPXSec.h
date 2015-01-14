//____________________________________________________________________________
/*!

\class    genie::PaisQELLambdaPXSec

\brief    

\ref      

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  June 10, 2004

\cpright  Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _PAIS_QEL_LAMBDA_PARTIAL_XSEC_H_
#define _PAIS_QEL_LAMBDA_PARTIAL_XSEC_H_

#include "Base/XSecAlgorithmI.h"
#include "Base/QELFormFactors.h"

namespace genie {

class QELFormFactorsModelI;
class XSecIntegratorI;

class PaisQELLambdaPXSec : public XSecAlgorithmI {

public:
  PaisQELLambdaPXSec();
  PaisQELLambdaPXSec(string config);
  virtual ~PaisQELLambdaPXSec();

  //-- XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k) const;
  double Integral        (const Interaction * i) const;
  bool   ValidProcess    (const Interaction * i) const;
  bool   ValidKinematics (const Interaction * i) const;

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

private:
  void  LoadConfig (void);
  double MHyperon(const Interaction * interaction) const;

  mutable QELFormFactors          fFormFactors;
  const   QELFormFactorsModelI *  fFormFactorsModel; 
  const   XSecIntegratorI *       fXSecIntegrator;
  double                          fSin8c2; 

};

} // genie namespace
#endif  // _PAIS_QEL_LAMBDA_PARTIAL_XSEC_H_
