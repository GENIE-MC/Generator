//____________________________________________________________________________
/*!

\class    genie::DISPXSec

\brief    Computes single differential DIS cross section.

          It can be configured to compute either
            \li dxsec / dx  where \c x is the Bjorken scaling variable, or
            \li dxsec / dy  where \c y is the Inelasticity

          Is a concrete implementation of the XSecAlgorithmI interface.
            
\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  July 01, 2004

*/
//____________________________________________________________________________

#ifndef _DIS_PXSEC_H_
#define _DIS_PXSEC_H_

#include "Base/XSecAlgorithmI.h"

namespace genie {

class DISPXSec : public XSecAlgorithmI {

public:

  DISPXSec();
  DISPXSec(const char * param_set);
  virtual ~DISPXSec();

  //-- XSecAlgorithmI interface implementation

  double XSec (const Interaction * interaction) const;
};

}       // genie namespace

#endif  // _DIS_PXSEC_H_
