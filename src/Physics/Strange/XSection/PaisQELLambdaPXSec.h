//____________________________________________________________________________
/*!

\class    genie::PaisQELLambdaPXSec

\brief    Implementation of the quasi-elastic scattering formula for 
          production of particles with different masses than the target.

\ref      Weak Interactions at High Energies
          A. Pais, Annals of Physics 63, 361-392 (1971)    
          Implemented here is equation 2.34 of the Pais paper, but ignoring
          lepton mass terms.  This equation is given also as Equation 3.37 in 
          the Llewellyn-Smith paper, though this paper uses slightly 
          different notation than that used in the Pais paper, and introduces
          a small error in the kinematic coefficient of the w2 term.  The 
          notation here by and large follows that of the Llewelyn-Smith paper. 

\author   Hugh Gallagher
          Tufts University

\created  June 10, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _PAIS_QEL_LAMBDA_PARTIAL_XSEC_H_
#define _PAIS_QEL_LAMBDA_PARTIAL_XSEC_H_

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Physics/QuasiElastic/XSection/QELFormFactors.h"

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
