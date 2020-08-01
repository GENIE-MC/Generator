//____________________________________________________________________________
/*!

\class    genie::MKSPPXSecWithCache

\brief    Class that caches neutrino resonance SPP cross sections on free 
          nucleons according to the MK model. This significantly speeds 
          the cross section calculation for multiple nuclear targets (eg at the
          spline construction phase, but only for case without Pauli-blocking). 
          This class integrates cross sections faster,
           
Computes the cross section for an neutrino resonance SPP reaction
         according to the MK model.


\authors  Igor Kakorin <kakorin@jinr.ru>, Joint Institute for Nuclear Research \n
          Konstantin Kuzmin <kkuzmin@theor.jinr.ru >,  Joint Institute for Nuclear Research \n
          Vadim Naumov <vnaumov@theor.jinr.ru >,  Joint Institute for Nuclear Research \n
          based on code of Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  November 12, 2019

\cpright  Copyright (c) 2003-2017, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _MK_SPP_XSEC_WITH_CACHE_H_
#define _MK_SPP_XSEC_WITH_CACHE_H_

#include <Math/IFunction.h>
#include <Math/IntegratorMultiDim.h>

#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Framework/ParticleData//BaryonResList.h"
#include "Framework/ParticleData/BaryonResonance.h"
#include "Framework/Utils/Range1.h"
#include "Framework/Interaction/KPhaseSpace.h"
#include "Framework/Interaction/SppChannel.h"

namespace genie {

class MKSPPXSecWithCache : public XSecIntegratorI {

protected:
  MKSPPXSecWithCache();
  MKSPPXSecWithCache(string name);
  MKSPPXSecWithCache(string name, string config);
  virtual ~MKSPPXSecWithCache();

  // Don't implement the XSecIntegratorI interface - leave it for the concrete
  // subclasses. Just define utility methods and data
  void   CacheResExcitationXSec (const Interaction * interaction) const;
  string CacheBranchName(SppChannel_t spp_channel, InteractionType_t it, int nu) const;

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
// genie::utils::gsl::d3XSecMK_dWQ2CosTheta_E
// A 3-D cross section function: d3xsec/dWdQ2dCosTheta = f(W, Q2, CosTheta)|(fixed E)
//
class d3XSecMK_dWQ2CosTheta_E: public ROOT::Math::IBaseFunctionMultiDim
{
public:
  d3XSecMK_dWQ2CosTheta_E(const XSecAlgorithmI * m, const Interaction * i);
 ~d3XSecMK_dWQ2CosTheta_E();

  // ROOT::Math::IBaseFunctionMultiDim interface
  unsigned int                        NDim   (void)               const;
  double                              DoEval (const double * xin) const;
  ROOT::Math::IBaseFunctionMultiDim * Clone  (void)               const;

private:
  const XSecAlgorithmI * fModel;
  Interaction *    fInteraction;
  Range1D_t Wl;
  bool isZero;
  KPhaseSpace * kps;
};

} // gsl   namespace
} // utils namespace
} // genie namespace

#endif  // _MK_SPP_XSEC_WITH_H_



