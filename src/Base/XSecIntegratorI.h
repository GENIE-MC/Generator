//____________________________________________________________________________
/*!

\class    genie::XSecIntegratorI

\brief    Cross Section Integrator Interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004

\cpright  Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _XSEC_INTEGRATOR_I_H_
#define _XSEC_INTEGRATOR_I_H_

#include "Algorithm/Algorithm.h"
#include "Base/XSecAlgorithmI.h"
#include "Interaction/Interaction.h"

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

  const IntegratorI * fIntegrator; ///< numerical integrator used
};

}       // genie namespace
#endif  // _XSEC_INTEGRATOR_I_H_
