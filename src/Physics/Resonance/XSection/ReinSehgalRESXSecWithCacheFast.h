//____________________________________________________________________________
/*!

\class    genie::ReinSehgalRESXSecWithCacheFast

\brief    Class that caches resonance neutrinoproduction cross sections on free
          nucleons according to the Rein-Sehgal model. This significantly speeds
          the cross section calculation for multiple nuclear targets (eg at the
          spline construction phase). This class integrates cross sections faster,
          than ReinSehgalRESXSecWithCache because of integration area transformation.

\ref      D.Rein and L.M.Sehgal, Neutrino Excitation of Baryon Resonances
          and Single Pion Production, Ann.Phys.133, 79 (1981)

\author   Igor Kakorin <kakorin@jinr.ru>
          Joint Institute for Nuclear Research

          based on code of
          Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  March 01, 2017

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _REIN_SEHGAL_RES_XSEC_WITH_CACHE_FAST_H_
#define _REIN_SEHGAL_RES_XSEC_WITH_CACHE_FAST_H_

#include <Math/IFunction.h>
#include <Math/IntegratorMultiDim.h>

#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Framework/ParticleData//BaryonResList.h"
#include "Framework/ParticleData/BaryonResonance.h"
#include "Framework/Utils/Range1.h"
#include "Framework/Interaction/KPhaseSpace.h"

namespace genie {

class ReinSehgalRESXSecWithCacheFast : public XSecIntegratorI {

protected:
  ReinSehgalRESXSecWithCacheFast();
  ReinSehgalRESXSecWithCacheFast(string name);
  ReinSehgalRESXSecWithCacheFast(string name, string config);
  virtual ~ReinSehgalRESXSecWithCacheFast();

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


class XSecAlgorithmI;
class Interaction;

namespace utils {
namespace gsl   {



//.....................................................................................
//
// genie::utils::gsl::d2XSecRESFast_dWQ2_E
// A 2-D cross section function: d2xsec/dWdQ2 = f(W,Q2)|(fixed E)
//
class d2XSecRESFast_dWQ2_E: public ROOT::Math::IBaseFunctionMultiDim
{
public:
  d2XSecRESFast_dWQ2_E(const XSecAlgorithmI * m, const Interaction * i);
 ~d2XSecRESFast_dWQ2_E();

  // ROOT::Math::IBaseFunctionMultiDim interface
  unsigned int                        NDim   (void)               const;
  double                              DoEval (const double * xin) const;
  ROOT::Math::IBaseFunctionMultiDim * Clone  (void)               const;

private:
  const XSecAlgorithmI * fModel;
  const Interaction *    fInteraction;
  double fWmin;
  double fWmax;
  bool isfWcutLessfWmin;
  KPhaseSpace * kps;
};

} // gsl   namespace
} // utils namespace
} // genie namespace

#endif  // _REIN_SEHGAL_RES_XSEC_WITH_CACHE_H_
