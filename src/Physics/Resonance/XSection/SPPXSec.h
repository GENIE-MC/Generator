//____________________________________________________________________________
/*!

\class    genie::SPPXSec

\brief    Computes the cross section for an neutrino resonance SPP reaction.

          Before computing any SPP cross section this algorithm 
          computes and caches splines for neutrino resonance SPP  cross 
          sections. This improves the speed of the GENIE spline construction 
          phase if splines for multiple nuclear targets are to be computed
          (for case without Pauli-blocking).

          Is a concrete implementation of the XSecAlgorithmI interface.\n


\authors  Igor Kakorin <kakorin@jinr.ru>, Joint Institute for Nuclear Research \n
          Vadim Naumov <vnaumov@theor.jinr.ru >,  Joint Institute for Nuclear Research \n
          based on code of Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  November 12, 2019

\cpright  Copyright (c) 2003-2023, GENIE Neutrino MC Generator Collaboration   
          For the full text of the license visit http://copyright.genie-mc.org                         
          or see $GENIE/LICENSE 
*/
//____________________________________________________________________________

#ifndef _SPP_XSEC_H_
#define _SPP_XSEC_H_

#include "Physics/Resonance/XSection/SPPXSecWithCache.h"

namespace genie {

class SPPXSec : public SPPXSecWithCache {

public:
  SPPXSec();
  SPPXSec(string param_set);
  virtual ~SPPXSec();

  // XSecIntegratorI interface implementation
  double Integrate(const XSecAlgorithmI * model, const Interaction * i) const;

  // Overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void LoadConfig(void);
  
  bool fUsePauliBlocking;      ///< account for Pauli blocking?
};

}       // genie namespace
#endif  // _SPP_XSEC_H_

