//____________________________________________________________________________
/*!

\class    genie::ReinSeghalRESXSecWithCache

\brief    An ABC that caches resonance neutrinoproduction cross sections on free 
          nucleons according to the Rein-Seghal model. This significantly speeds 
          the cross section calculation for multiple nuclear targets (eg at the
          spline construction phase)

\ref      D.Rein and L.M.Seghal, Neutrino Excitation of Baryon Resonances
          and Single Pion Production, Ann.Phys.133, 79 (1981)

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  March 09, 2006

*/
//____________________________________________________________________________

#ifndef _REIN_SEGHAL_RES_XSEC_WITH_CACHE_H_
#define _REIN_SEGHAL_RES_XSEC_WITH_CACHE_H_

#include "Base/XSecAlgorithmI.h"
#include "BaryonResonance/BaryonResList.h"
#include "BaryonResonance/BaryonResonance.h"
#include "Utils/Range1.h"

namespace genie {

class IntegratorI;

class ReinSeghalRESXSecWithCache : public XSecAlgorithmI {

protected:
  ReinSeghalRESXSecWithCache();
  ReinSeghalRESXSecWithCache(string name);
  ReinSeghalRESXSecWithCache(string name, string config);
  virtual ~ReinSeghalRESXSecWithCache();

  // Don't implement the XSecAlgorithmI interface - leave it for the concrete
  // subclasses. Just define utility methods and data

  void   CacheResExcitationXSec (const Interaction * interaction) const;
  string CacheBranchName(Resonance_t r, InteractionType_t it, int nu, int nuc) const;

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

#endif  // _REIN_SEGHAL_RES_XSEC_WITH_CACHE_H_

