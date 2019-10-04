//____________________________________________________________________________
/*!

\class    genie::XSecIntegratorI

\brief    Cross Section Integrator Interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 03, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _XSEC_INTEGRATOR_I_H_
#define _XSEC_INTEGRATOR_I_H_

#include "Framework/Algorithm/Algorithm.h"
#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/Interaction/Interaction.h"

namespace genie {

class IntegratorI;

 class XSecIntegratorI : public Algorithm {

public:
  virtual ~XSecIntegratorI();

  virtual double Integrate(const XSecAlgorithmI * model, 
                           const Interaction * interaction 
                       /*, const KPhaseSpaceCut * cut=0*/) const= 0;
protected:
  XSecIntegratorI();
  XSecIntegratorI(string name);
  XSecIntegratorI(string name, string config);

  const IntegratorI * fIntegrator; ///< GENIE numerical integrator 

  string fGSLIntgType;                     ///< name of GSL numerical integrator
  double fGSLRelTol;                       ///< required relative tolerance (error)
  int    fGSLMaxEval;                      ///< GSL max evaluations
  int    fGSLMinEval;                      ///< GSL min evaluations. Ignored by some integrators.
  unsigned int fGSLMaxSizeOfSubintervals;  ///< GSL maximum number of sub-intervals for 1D integrator
  unsigned int fGSLRule;                   ///< GSL Gauss-Kronrod integration rule (only for GSL 1D adaptive type)

};

}       // genie namespace
#endif  // _XSEC_INTEGRATOR_I_H_
