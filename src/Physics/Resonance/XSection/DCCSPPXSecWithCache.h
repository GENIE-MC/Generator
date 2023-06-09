//____________________________________________________________________________
/*!

\class    genie::DCCSPPXSecWithCache

\brief    Class that caches neutrino resonance SPP cross sections on free 
          nucleons according to the DCC model. This significantly speeds 
          the cross section calculation for multiple nuclear targets (eg at the
          spline construction phase, but only for case without Pauli-blocking). 
          This class integrates cross sections faster,
           
Computes the cross section for an neutrino resonance SPP reaction
         according to the DCC model.


\authors  Igor Kakorin <kakorin@jinr.ru>, Joint Institute for Nuclear Research \n
          based on code of Costas Andreopoulos <costas.andreopoulos@stfc.ac.ukk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  Apr 12, 2023

\cpright  Copyright (c) 2003-2022, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _DCC_SPP_XSEC_WITH_CACHE_H_
#define _DCC_SPP_XSEC_WITH_CACHE_H_

#include <Math/IFunction.h>
#include <Math/IntegratorMultiDim.h>

#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Framework/ParticleData//BaryonResList.h"
#include "Framework/ParticleData/BaryonResonance.h"
#include "Framework/Utils/Range1.h"
#include "Framework/Interaction/KPhaseSpace.h"
#include "Framework/Interaction/SppChannel.h"

namespace genie {

class DCCSPPXSecWithCache : public XSecIntegratorI {

protected:
  DCCSPPXSecWithCache();
  DCCSPPXSecWithCache(string name);
  DCCSPPXSecWithCache(string name, string config);
  virtual ~DCCSPPXSecWithCache();

  // Don't implement the XSecIntegratorI interface - leave it for the concrete
  // subclasses. Just define utility methods and data
  void   CacheResExcitationXSec (const Interaction * interaction) const;
  string CacheBranchName(SppChannel_t spp_channel, InteractionType_t it, int nu, int helicity) const;
  string ProbeAsString (int probe_pdg, int probe_helicity) const;

  bool   fUsingDisResJoin;
  double fWcut;
  double fEMax;
  double fQ2minForMasslessLepton;   ///<  Minimal boundary of integration by Q2 for massless lepton in EM processes

  mutable const XSecAlgorithmI * fSinglePionProductionXSecModel;
};


class XSecAlgorithmI;
class Interaction;

namespace utils {
namespace gsl   {



//.....................................................................................
//
// genie::utils::gsl::d2XSecDCC_dWQ2_E
// A 2-D cross section function: d2xsec/dWdQ2 = f(W, Q2)|(fixed E)
//
class d2XSecSPP_dWQ2_E: public ROOT::Math::IBaseFunctionMultiDim
{
public:
  d2XSecSPP_dWQ2_E(const XSecAlgorithmI * m, const Interaction * i, double wcut, bool massless, double Q2min);
 ~d2XSecSPP_dWQ2_E();

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
  double fWcut;
  bool   fMasslessLepton;           ///< Is charged lepton massless in EM process?
  double fQ2minForMasslessLepton;   ///<  Minimal boundary of integration by Q2 for massless lepton in EM processes
};

} // gsl   namespace
} // utils namespace
} // genie namespace

#endif  // _DCC_SPP_XSEC_WITH_H_



