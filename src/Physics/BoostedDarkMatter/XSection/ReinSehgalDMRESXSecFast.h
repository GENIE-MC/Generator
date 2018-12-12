//____________________________________________________________________________
/*!

\class    genie::ReinSehgalDMRESXSecFast

\brief    Computes the cross section for an exclusive 1pi reaction through
          resonance neutrinoproduction according to the Rein-Sehgal model.

          This algorithm produces in principle what you could also get from 
          the genie::DMRESXSec algorithm (RES cross section integrator) by 
          specifying the genie::ReinSehgalRESPXSec as the differential 
          cross section model. However, ReinSehgalDMRESXSecFast offers a faster 
          alternative. Before computing any RES cross section this algorithm 
          computes and caches splines for resonance neutrino-production cross 
          sections. This improves the speed of the GENIE spline construction 
          phase if splines for multiple nuclear targets are to be computed.
          Also this class integrates cross sections faster, than 
          ReinSehgalDMRESXSec because of integration area transformation.

          Is a concrete implementation of the XSecAlgorithmI interface.\n

\ref      D.Rein and L.M.Sehgal, Neutrino Excitation of Baryon Resonances
          and Single Pion Production, Ann.Phys.133, 79 (1981)

\author   Igor Kakorin <kakorin@jinr.ru>
          Joint Institute for Nuclear Research - March 01, 2017 \n
          based on code of Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  March 01, 2017

\cpright  Copyright (c) 2003-2017, GENIE Neutrino MC Generator Collaboration   
          For the full text of the license visit http://copyright.genie-mc.org                         
          or see $GENIE/LICENSE 
*/
//____________________________________________________________________________

#ifndef _REIN_SEHGAL_DMRES_XSEC_FAST_H_
#define _REIN_SEHGAL_DMRES_XSEC_FAST_H_

#include "Physics/BoostedDarkMatter/XSection/ReinSehgalDMRESXSecWithCacheFast.h"

namespace genie {

class ReinSehgalDMRESXSecFast : public ReinSehgalDMRESXSecWithCacheFast {

public:
  ReinSehgalDMRESXSecFast();
  ReinSehgalDMRESXSecFast(string param_set);
  virtual ~ReinSehgalDMRESXSecFast();

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
#endif  // _REIN_SEHGAL_DMRES_XSEC_H_

