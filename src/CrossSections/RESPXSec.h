//____________________________________________________________________________
/*!

\class    genie::RESPXSec

\brief    Computes single differential RES cross section.

          It can be configured to compute either
            \li dxsec / dQ2  where \c Q2 is the momentum transfer, or
            \li dxsec / dW   where \c W is the hadronic invariant mass

          This is a model-independent algorithm. It merely integrates the
          specified double differential cross section.
            
          Is a concrete implementation of the XSecAlgorithmI interface.
            
\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  March 03, 2005

*/
//____________________________________________________________________________

#ifndef _RES_PXSEC_H_
#define _RES_PXSEC_H_

#include "Base/XSecAlgorithmI.h"

namespace genie {

class IntegratorI;

class RESPXSec : public XSecAlgorithmI {

public:

  RESPXSec();
  RESPXSec(const char * param_set);
  virtual ~RESPXSec();

  //-- XSecAlgorithmI interface implementation

  double XSec (const Interaction * interaction) const;

public:

  const IntegratorI * Integrator(void) const;
};

}       // genie namespace
#endif  // _RES_PXSEC_H_
