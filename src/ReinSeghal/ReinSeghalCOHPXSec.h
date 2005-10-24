//____________________________________________________________________________
/*!

\class    genie::ReinSeghalCOHPXSec

\brief    Computes the double differential cross section for coherent pi0
          production according to the \b Rein-Seghal model.

          The computed cross section is the d^3 xsec/ dx dy dt

          where \n
            \li \c x : Bjorken x = Q2/2Mv
            \li \c y : Inelasticity y=v/E, v=E-E'

          The t dependence is analytically integrated out.

          Is a concrete implementation of the XSecAlgorithmI interface.

\ref      D.Rein and L.M.Seghal, Coherent pi0 production in neutrino
          reactions, Nucl.Phys.B223:29-144 (1983)

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  March 11, 2005

____________________________________________________________________________*/

#ifndef _REIN_SEGHAL_COH_PARTIAL_XSEC_H_
#define _REIN_SEGHAL_COH_PARTIAL_XSEC_H_

#include "Base/XSecAlgorithmI.h"

namespace genie {

class ReinSeghalCOHPXSec : public XSecAlgorithmI {

public:

  ReinSeghalCOHPXSec();
  ReinSeghalCOHPXSec(string config);
  virtual ~ReinSeghalCOHPXSec();

  //-- XSecAlgorithmI interface implementation

  double XSec (const Interaction * interaction) const;

private:

  //-- methods for overiding default parameters by setting them in the
  //   configuration registry
  double Ma           (void) const; ///< override default axial mass
  double ReImPiApl    (void) const; ///< override default Re/Im Fwd pi scat. ampl.
  double NuclSizeScale(void) const; ///< override default nuclear size scale param
};

}       // genie namespace

#endif  // _REIN_SEGHAL_COH_PARTIAL_XSEC_H_
