//____________________________________________________________________________
/*!

\class    genie::ReinSeghalSPPXSec

\brief    Computes the cross section for an exclusive 1pi reaction through
          resonance neutrinoproduction according to the Rein-Seghal model.

          This algorithm produces in principle what you could also get from 
          the genie::RESXSec algorithm (RES cross section integrator) by 
          specifying the genie::ReinSeghalSPPPXSec as the differential 
          (d2xsec/fQ2dW) cross section model. However, ReinSeghalSPPXSec
          offers a faster alternative. Before computing any SPP cross section
          this algorithm computes and caches splines for resonance neutrino-
          production cross sections. This improves the speed since it is 
          reducing the number of calculations (the generic algorithm needs to
          recompute resonance production xsec for every exclusive channel).

          In this algorithm we follow the non-coherent approach: we sum
          the weighted resonance production cross sections rather than the
          resonance production amplitudes.

          Is a concrete implementation of the XSecAlgorithmI interface.\n

\ref      D.Rein and L.M.Seghal, Neutrino Excitation of Baryon Resonances
          and Single Pion Production, Ann.Phys.133, 79 (1981)

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  March 09, 2006

*/
//____________________________________________________________________________

#ifndef _REIN_SEGHAL_SPP_XSEC_H_
#define _REIN_SEGHAL_SPP_XSEC_H_

#include "Base/XSecAlgorithmI.h"
#include "BaryonResonance/BaryonResList.h"
#include "Utils/Range1.h"

namespace genie {

class IntegratorI;

class ReinSeghalSPPXSec : public XSecAlgorithmI {

public:
  ReinSeghalSPPXSec();
  ReinSeghalSPPXSec(string param_set);
  virtual ~ReinSeghalSPPXSec();

  //-- XSecAlgorithmI interface implementation
  double XSec            (const Interaction * interaction) const;
  bool   ValidProcess    (const Interaction * interaction) const;
  bool   ValidKinematics (const Interaction * interaction) const;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void CacheResExcitationXSec (const Interaction * interaction) const;

  void LoadConfig      (void);
  void GetResonanceList(void);

  Range1D_t WRange (const Interaction * interaction) const;
  Range1D_t Q2Range(const Interaction * interaction) const;

  double fWminCut;
  double fWmaxCut;
  double fQ2minCut;
  double fQ2maxCut;
  double fEMax;

  const XSecAlgorithmI * fSingleResXSecModel;
  const IntegratorI *    fIntegrator;
  BaryonResList          fResList;
};

}       // genie namespace
#endif  // _REIN_SEGHAL_SPP_XSEC_H_

