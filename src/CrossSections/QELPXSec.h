//____________________________________________________________________________
/*!

\class    genie::QELPXSec

\brief    Computes the differential Quasi Elastic cross section dxsec/dq^2.

          Is a concrete implementation of the XSecAlgorithmI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 05, 2004

*/
//____________________________________________________________________________

#ifndef _QEL_PARTIAL_XSEC_H_
#define _QEL_PARTIAL_XSEC_H_

#include "Base/XSecAlgorithmI.h"

namespace genie {

class QELPXSec : public XSecAlgorithmI {

public:

  QELPXSec();
  QELPXSec(const char * param_set);
  virtual ~QELPXSec();

  //-- XSecAlgorithmI interface implementation
  
  double XSec (const Interaction * interaction) const;
};

}       // genie namespace

#endif  // _QEL_PARTIAL_XSEC_H_
