//____________________________________________________________________________
/*!

\class    genie::StdElasticPXSec

\brief    Standard differential cross section dxsec/dQ^2 for v+N / vbar+N
          elastic scattering.        
          
          StdElasticPXSec is a concrete implementation of the
          XSecAlgorithmI interface. \n

\ref      L.A.Ahrens et al., Physical Review D, VOL 35,3:785 (1987)

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  Fabruary 15, 2005

*/
//____________________________________________________________________________

#ifndef _STD_ELASTIC_PARTIAL_XSEC_H_
#define _STD_ELASTIC_PARTIAL_XSEC_H_

#include "Base/XSecAlgorithmI.h"

namespace genie {

class StdElasticPXSec : public XSecAlgorithmI {

public:

  StdElasticPXSec();
  StdElasticPXSec(const char * param_set);
  virtual ~StdElasticPXSec();

  //-- XSecAlgorithmI interface implementation

  double XSec (const Interaction * interaction) const;

private:

  // symbols follow the notation in L.A.Ahrens et al. paper
  double GA    (const Interaction * interaction) const;
  double F1    (const Interaction * interaction) const;
  double F2    (const Interaction * interaction) const;
  double Alpha (void) const;
  double Gamma (void) const;
};

}       // genie namespace

#endif  // _STD_ELASTIC_PARTIAL_XSEC_H_
