//____________________________________________________________________________
/*!

\class    genie::ReinSehgalRESXSecWithCache

\brief    An ABC that caches resonance neutrinoproduction cross sections on free 
          nucleons according to the Rein-Sehgal model. This significantly speeds 
          the cross section calculation for multiple nuclear targets (eg at the
          spline construction phase)

\ref      D.Rein and L.M.Sehgal, Neutrino Excitation of Baryon Resonances
          and Single Pion Production, Ann.Phys.133, 79 (1981)

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  March 09, 2006

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _REIN_SEHGAL_RES_XSEC_WITH_CACHE_H_
#define _REIN_SEHGAL_RES_XSEC_WITH_CACHE_H_

#include "Framework/ParticleData/BaryonResList.h"
#include "Framework/ParticleData/BaryonResonance.h"
#include "Framework/Utils/Range1.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"

namespace genie {

class ReinSehgalRESXSecWithCache : public XSecIntegratorI {

protected:
  ReinSehgalRESXSecWithCache();
  ReinSehgalRESXSecWithCache(string name);
  ReinSehgalRESXSecWithCache(string name, string config);
  virtual ~ReinSehgalRESXSecWithCache();

  // Don't implement the XSecIntegratorI interface - leave it for the concrete
  // subclasses. Just define utility methods and data
  void   CacheResExcitationXSec (const Interaction * interaction) const;
  string CacheBranchName(Resonance_t r, InteractionType_t it, int nu, int nuc) const;

  bool   fUsingDisResJoin;
  double fWcut;
  double fEMax;

  mutable const XSecAlgorithmI * fSingleResXSecModel;
  BaryonResList fResList;
};

}       // genie namespace
#endif  // _REIN_SEHGAL_RES_XSEC_WITH_CACHE_H_

