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
  double XSec            (const Interaction * interaction) const;
  bool   ValidProcess    (const Interaction * interaction) const;
  bool   ValidKinematics (const Interaction * interaction) const;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void LoadConfigData(void);

  //-- private data members loaded from config Registry or set to defaults
  double fMa;   ///< axial mass
  double fReIm; ///< Re/Im {forward pion scattering amplitude}
  double fRo;   ///< nuclear size scale parameter
};

}       // genie namespace

#endif  // _REIN_SEGHAL_COH_PARTIAL_XSEC_H_
