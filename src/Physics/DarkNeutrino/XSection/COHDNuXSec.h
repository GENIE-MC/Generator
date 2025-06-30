//____________________________________________________________________________
/*!

\class    genie::COHDNuXSec

\brief    Computes the cross section for coherent dark neutrino scattering.\n
          Is a concrete implementation of the XSecIntegratorI interface.

\author   Author: Iker de Icaza <i.de-icaza-astiz \at sussex.ac.uk>
          University of Sussex

          Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  June 12, 2020

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _COHERENT_DNu_XSEC_H_
#define _COHERENT_DNu_XSEC_H_

#include "Physics/XSectionIntegration/XSecIntegratorI.h"

namespace genie {

class COHDNuXSec : public XSecIntegratorI {

public:
  COHDNuXSec();
  COHDNuXSec(string config);
  virtual ~COHDNuXSec();

  // XSecIntegratorI interface implementation
  double Integrate(const XSecAlgorithmI * model, const Interaction * i) const;

  // Overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

 private:

  void LoadConfig (void);

  double fGSLAbsTol ;
  double fDNuMass;
};

} // genie namespace

#endif  // _COHERENT_DNu_XSEC_H_
