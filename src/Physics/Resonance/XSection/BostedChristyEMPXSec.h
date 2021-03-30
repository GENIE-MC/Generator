//____________________________________________________________________________
/*!

\class    genie::BostedChristyEMPXSec

\brief    An empirical fit of inclusive inelastic electron-proton cross sections in the
          kinematic range of four-momentum transfer 0 < Q2 < 8 GeV^2 and final-state invariant 
          mass 1.1 < W <3.1 GeV constrained by the high-precision longitudinal and 
          transverse separated cross section measurements from Jefferson Lab Hall C. 
          


\author   Igor Kakorin <kakorin@jinr.ru>
          Joint Institute for Nuclear Research

\ref      M.E. Christy, P.E.Bosted, "Empirical fit to precision inclusive 
          electron-proton cross sections in the resonance region", PRC 81 (2010) 055213

\created  April 3, 2021

\cpright  Copyright (c) 2003-2021, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________
#ifndef _BOOSTED_CHRISTY_EM_PXSEC_H_
#define _BOOSTED_CHRISTY_EM_PXSEC_H_

#include "Framework/EventGen/XSecAlgorithmI.h"

namespace genie {

class BostedChristyEMPXSec : public XSecAlgorithmI {

public:
  BostedChristyEMPXSec();
  BostedChristyEMPXSec(string config);
  virtual ~BostedChristyEMPXSec();

  // implement the XSecAlgorithmI interface
  double XSec         (const Interaction * i, KinePhaseSpace_t k) const;
  double Integral     (const Interaction * i) const;
  bool   ValidProcess (const Interaction * i) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void LoadConfig (void);
  double sigmaR(int, double, double) const;
  double sigmaNR(int, double, double) const;

  double fbr[7][3];            ///< branching ratios of resonances
  
  int fang[7];                 ///< meson angular momentum
  
  double fmass[7];             ///< meson mass
  
  double fwidth[7];            ///< meson width
  
  double frescoefT[7][4];      ///< tunable parameters from Table III for resonance \sigma_T
  double frescoefL[7][3];      ///< tunable parameters from Table III for resonance \sigma_L
  
  double fnrcoefT[2][5];       ///< tunable parameters from Table III for resonance \sigma_T
  double fnrcoefL[5];       ///< tunable parameters from Table III for resonance \sigma_L
  
  double fXSecScaleEM;         ///< external EM xsec scaling factor

  const XSecIntegratorI * fXSecIntegrator;
};

}       // genie namespace




#endif  // _BOOSTED_CHRISTY_EM_PXSEC_H_
