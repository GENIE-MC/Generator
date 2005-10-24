//____________________________________________________________________________
/*!

\class    genie::DISPartonModelPXSec

\brief    Massless Parton Model DIS Partial (d^2xsec/dxdy) Cross Section

\ref

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 05, 2004

*/
//____________________________________________________________________________

#ifndef _DIS_PARTON_MODEL_PARTIAL_XSEC_H_
#define _DIS_PARTON_MODEL_PARTIAL_XSEC_H_

#include "Base/XSecAlgorithmI.h"

namespace genie {

class DISPartonModelPXSec : public XSecAlgorithmI {

public:

  DISPartonModelPXSec();
  DISPartonModelPXSec(string config);
  virtual ~DISPartonModelPXSec();

  //-- XSecAlgorithmI interface implementation

  double XSec (const Interaction * interaction) const;
};

}       // genie namespace

#endif  // _DIS_PARTON_MODEL_PARTIAL_XSEC_H_
