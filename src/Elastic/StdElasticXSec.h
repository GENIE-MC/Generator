//____________________________________________________________________________
/*!

\class    genie::StdElasticXSec

\brief    Standard v+N / vbar+N elastic scattering cross section.

          StdElasticPXSec is a concrete implementation of the
          XSecAlgorithmI interface. \n

\ref      L.A.Ahrens et al., Physical Review D, VOL 35,3:785 (1987)

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  Fabruary 15, 2005

*/
//____________________________________________________________________________

#ifndef _STD_ELASTIC_XSEC_H_
#define _STD_ELASTIC_XSEC_H_

#include "Base/XSecAlgorithmI.h"

namespace genie {

class IntegratorI;

class StdElasticXSec : public XSecAlgorithmI {

public:

  StdElasticXSec();
  StdElasticXSec(string config);
  virtual ~StdElasticXSec();

  //-- XSecAlgorithmI interface implementation

  double XSec (const Interaction * interaction) const;

private:

  const XSecAlgorithmI * DiffXSec   (void) const;
  const IntegratorI *    Integrator (void) const;
};

}       // genie namespace

#endif  // _STD_ELASTIC_XSEC_H_
