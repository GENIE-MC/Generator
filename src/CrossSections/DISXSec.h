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
#include "Numerical/IntegratorI.h"

namespace genie {

class DISXSec : public XSecAlgorithmI {

public:

  DISXSec();
  DISXSec(const char * param_set);
  virtual ~DISXSec();

  //-- XSecAlgorithmI interface implementation

  double XSec (const Interaction * interaction) const;

private:

  const XSecAlgorithmI * PartialXSecAlgorithm (void) const;
  const IntegratorI *    Integrator           (void) const;

  bool IsWithinIntegrationRange(const Interaction * interaction) const;

  int    NLogX (void) const;
  int    NLogY (void) const;
  double Xmin  (void) const;
  double Xmax  (void) const;
  double Ymin  (void) const;
  double Ymax  (void) const;
};

}       // genie namespace

#endif  // _DIS_XSEC_H_
