//____________________________________________________________________________
/*!

\class    genie::RESXSec

\brief    Computes the RES Cross Section.

          Is a concrete implementation of the XSecAlgorithmI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 04, 2004

*/
//____________________________________________________________________________

#ifndef _RES_XSEC_H_
#define _RES_XSEC_H_

#include "Base/XSecAlgorithmI.h"
#include "Utils/Range1.h"

namespace genie {

class RESXSec : public XSecAlgorithmI {

public:

  RESXSec();
  RESXSec(string param_set);
  virtual ~RESXSec();

  //-- XSecAlgorithmI interface implementation
  double XSec (const Interaction * interaction) const;

private:

  Range1D_t WRange  (const Interaction * interaction) const;
  Range1D_t Q2Range (const Interaction * interaction) const;
};

}       // genie namespace

#endif  // _RES_XSEC_H_
