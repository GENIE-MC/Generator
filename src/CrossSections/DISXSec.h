//____________________________________________________________________________
/*!

\class    genie::DISXSec

\brief    Computes the DIS Cross Section.

          Is a concrete implementation of the XSecAlgorithmI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 04, 2004

*/
//____________________________________________________________________________

#ifndef _DIS_XSEC_H_
#define _DIS_XSEC_H_

#include "Base/XSecAlgorithmI.h"

namespace genie {

class DISXSec : public XSecAlgorithmI {

public:

  DISXSec();
  DISXSec(string config);
  virtual ~DISXSec();

  //-- XSecAlgorithmI interface implementation

  double XSec (const Interaction * interaction) const;

private:

  bool IsWithinIntegrationRange(const Interaction * interaction) const;
};

}       // genie namespace

#endif  // _DIS_XSEC_H_
