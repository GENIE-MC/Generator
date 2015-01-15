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

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  March 09, 2006

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration   
          For the full text of the license visit http://copyright.genie-mc.org                         
          or see $GENIE/LICENSE 
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

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void LoadConfig(void);
};

}       // genie namespace
#endif  // _REIN_SEGHAL_RES_XSEC_H_

