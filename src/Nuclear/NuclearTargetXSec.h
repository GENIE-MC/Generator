//____________________________________________________________________________
/*!

\class    genie::NuclearTargetXSec

\brief    Computes the cross section for a nuclear target from the free
          nucleon cross sections.

          Is a concrete implementation of the XSecAlgorithmI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  June 10, 2004

*/
//____________________________________________________________________________

#ifndef _NUCLEAR_TARGET_XSEC_H_
#define _NUCLEAR_TARGET_XSEC_H_

#include "Base/XSecAlgorithmI.h"

namespace genie {

class NuclearTargetXSec : public XSecAlgorithmI {

public:

  NuclearTargetXSec();
  NuclearTargetXSec(const char * param_set);
  virtual ~NuclearTargetXSec();

  //-- XSecAlgorithmI interface implementation
  
  double XSec (const Interaction * interaction) const;
};

}       // genie namespace

#endif  // _NUCLEAR_TARGET_XSEC_H_
