//____________________________________________________________________________
/*!

\class    genie::RSSingleRESPXSec

\brief    Computes the (+/- breit-wigner weighted) differential RES xsec
          d^2 xsec/ dQ^2 dW for a *single* Baryon Resonance.

          The actual Rein-Seghal model is coded at ReinSeghalPartialXSec.

          This algorithm merely decides which of the ReinSeghalPartialXSec
          algorithm instances to run depending on the initial state:
          1) CC (n+p), 2) NC (n), 3) NC (p).

          Concrete implementation of the XSecAlgorithmI interface.

\ref      D.Rein and L.M.Seghal, Neutrino Excitation of Baryon Resonances
          and Single Pion Production, Ann.Phys.133, 79 (1981)

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  November 22, 2004

*/
//____________________________________________________________________________

#ifndef _RS_SINGLE_RESONANCE_PARTIAL_XSEC_H_
#define _RS_SINGLE_RESONANCE_PARTIAL_XSEC_H_

#include "Base/XSecAlgorithmI.h"

namespace genie {

class RSSingleRESPXSec : public XSecAlgorithmI {

public:

  RSSingleRESPXSec();
  RSSingleRESPXSec(string config);
  virtual ~RSSingleRESPXSec();

  //-- XSecAlgorithmI interface implementation
  double XSec            (const Interaction * interaction) const;
  bool   ValidProcess    (const Interaction * interaction) const;
  bool   ValidKinematics (const Interaction * interaction) const;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void LoadSubAlg(void);

  const XSecAlgorithmI * fRSCC;
  const XSecAlgorithmI * fRSNCp;
  const XSecAlgorithmI * fRSNCn;
};

}       // genie namespace

#endif  // _RS_SINGLE_RESONANCE_PARTIAL_XSEC_H_
