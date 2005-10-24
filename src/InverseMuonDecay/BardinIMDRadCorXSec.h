//____________________________________________________________________________
/*!

\class    genie::BardinIMDRadCorXSec

\brief    Computes the Inverse Muon Decay cross section using the Bardin -
          Dokuchaeva model which includes all 1-loop radiative corrections. \n

          This algorithm merely integrates the Bardin differential IMD cross
          section. The specific differential cross section algorithm is
          specified in this algorithm's XML config file.

          The exact 'type' of the cross section depends on the specified
          differential IMD cross section algorithm. It can be a 'trully'
          inclusive IMD cross section or a cross section where part of the
          brem cross section contribution is subtracted
          (for futher details, see the documentation of the Bardin-Dokuchaeva
          model diffential IMD cross section algorithms and the Bardin paper,
          cited below).

          BardinIMDRadCorXSec is a concrete implementation of the
          XSecAlgorithmI interface. \n

\ref      D.Yu.Bardin and V.A.Dokuchaeva, Nucl.Phys.B287:839 (1987)

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  Fabruary 14, 2005

*/
//____________________________________________________________________________

#ifndef _BARDIN_IMD_RADIATIVE_CORRECTIONS_XSEC_H_
#define _BARDIN_IMD_RADIATIVE_CORRECTIONS_XSEC_H_

#include "Base/XSecAlgorithmI.h"
#include "Numerical/IntegratorI.h"

namespace genie {

class BardinIMDRadCorXSec : public XSecAlgorithmI {

public:

  BardinIMDRadCorXSec();
  BardinIMDRadCorXSec(string config);
  virtual ~BardinIMDRadCorXSec();

  //-- XSecAlgorithmI interface implementation
  double XSec (const Interaction * interaction) const;
};

}       // genie namespace

#endif  // _BARDIN_IMD_RADIATIVE_CORRECTIONS_XSEC_H_
