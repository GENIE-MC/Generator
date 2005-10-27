//____________________________________________________________________________
/*!

\class    genie::QELXSec

\brief    Computes the Quasi Elastic (QEL) cross section.

          Is a concrete implementation of the XSecAlgorithmI interface. \n

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 04, 2004

*/
//____________________________________________________________________________

#ifndef _QEL_XSEC_H_
#define _QEL_XSEC_H_

#include "Base/XSecAlgorithmI.h"


namespace genie {

class QELXSec : public XSecAlgorithmI {

public:

  QELXSec();
  QELXSec(string config);
  virtual ~QELXSec();

  //-- XSecAlgorithmI interface implementation
  double XSec (const Interaction * interaction) const;
};

}       // genie namespace

#endif  // _QEL_XSEC_H_
