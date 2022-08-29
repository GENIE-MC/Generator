//____________________________________________________________________________
/*!

\class    genie::DCCSPPXSec

\brief    Computes the cross section for an neutrino resonance SPP
          reaction according to the DCC model.

          Before computing any SPP cross section this algorithm 
          computes and caches splines for neutrino resonance SPP  cross 
          sections. This improves the speed of the GENIE spline construction 
          phase if splines for multiple nuclear targets are to be computed
          (for case without Pauli-blocking).

          Is a concrete implementation of the XSecAlgorithmI interface.\n


\authors  Igor Kakorin <kakorin@jinr.ru>, Joint Institute for Nuclear Research \n
          based on code of Costas Andreopoulos <costas.andreopoulos@stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  Septemper 06, 2022

\cpright  Copyright (c) 2003-2022, GENIE Neutrino MC Generator Collaboration   
          For the full text of the license visit http://copyright.genie-mc.org                         
          or see $GENIE/LICENSE 
*/
//____________________________________________________________________________

#ifndef _DCC_SPP_XSEC_H_
#define _DCC_SPP_XSEC_H_

#include "Physics/Resonance/XSection/DCCSPPXSecWithCache.h"

namespace genie {

class DCCSPPXSec : public DCCSPPXSecWithCache {

public:
  DCCSPPXSec();
  DCCSPPXSec(string param_set);
  virtual ~DCCSPPXSec();

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
#endif  // _DCC_SPP_XSEC_H_

