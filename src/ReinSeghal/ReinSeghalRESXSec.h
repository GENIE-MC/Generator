//____________________________________________________________________________
/*!

\class    genie::ReinSeghalRESXSec

\brief    Computes the cross section for an exclusive 1pi reaction through
          resonance neutrinoproduction according to the Rein-Seghal model.

          This algorithm produces in principle what you could also get from 
          the genie::RESXSec algorithm (RES cross section integrator) by 
          specifying the genie::ReinSeghalRESPXSec as the differential 
          cross section model. However, ReinSeghalRESXSec offers a faster 
          alternative. Before computing any RES cross section this algorithm 
          computes and caches splines for resonance neutrino-production cross 
          sections. This improves the speed of the GENIE spline construction 
          phase if splines for multiple nuclear targets are to be computed.

          Is a concrete implementation of the XSecAlgorithmI interface.\n

\ref      D.Rein and L.M.Seghal, Neutrino Excitation of Baryon Resonances
          and Single Pion Production, Ann.Phys.133, 79 (1981)

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  March 09, 2006

\cpright  Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration   
          All rights reserved.                         
          For the licensing terms see $GENIE/USER_LICENSE. 
*/
//____________________________________________________________________________

#ifndef _REIN_SEGHAL_RES_XSEC_H_
#define _REIN_SEGHAL_RES_XSEC_H_

#include "ReinSeghal/ReinSeghalRESXSecWithCache.h"

namespace genie {

class ReinSeghalRESXSec : public ReinSeghalRESXSecWithCache {

public:
  ReinSeghalRESXSec();
  ReinSeghalRESXSec(string param_set);
  virtual ~ReinSeghalRESXSec();

  //-- XSecIntegratorI interface implementation
  double Integrate(const XSecAlgorithmI * model, const Interaction * i) const;

/*
  double XSec            (const Interaction * i, KinePhaseSpace_t k=kPSfE) const;
  bool   ValidProcess    (const Interaction * i) const;
  bool   ValidKinematics (const Interaction * i) const;
*/
  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void LoadConfig(void);
};

}       // genie namespace

#endif  // _REIN_SEGHAL_RES_XSEC_H_

