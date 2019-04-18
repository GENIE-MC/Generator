//____________________________________________________________________________
/*!

\class    genie::MartiniEricsonChanfrayMarteauMECPXSec2016

\brief    Computes the Martini, Ericson, Chanfray and Marteau MEC model 
          differential cross section.
          Uses precomputed hadon tensor tables.
          Is a concrete implementation of the XSecAlgorithmI interface. 

\author   Sara Bolognesi <sara.bolognesi@cea.fr>
          CEA Saclay

          Marco Martini
          CEA Saclay

\ref      M. Martini, M. Ericson, G. Chanfray, J. Marteau.
          Neutrino and antineutrino quasielastic interactions with nuclei
          Phys.Rev. C81 (2010) 045502

\created  Mar 30, 2016

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _MARTINI_ERICSON_CHANFRAY_MARTEAU_MEC_PXSEC_2016_H_
#define _MARTINI_ERICSON_CHANFRAY_MARTEAU_MEC_PXSEC_2016_H_

#include <vector>

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Physics/Multinucleon/XSection/MECHadronTensor.h"

using std::vector;

namespace genie {

class XSecIntegratorI;

class MartiniEricsonChanfrayMarteauMECPXSec2016 : public XSecAlgorithmI {

public:
  MartiniEricsonChanfrayMarteauMECPXSec2016();
  MartiniEricsonChanfrayMarteauMECPXSec2016(string config);
  virtual ~MartiniEricsonChanfrayMarteauMECPXSec2016();

  // XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k) const;
  double Integral        (const Interaction * i) const;
  bool   ValidProcess    (const Interaction * i) const;

  // override the Algorithm::Configure methods to load configuration
  // data to private data members
  void Configure (const Registry & config);
  void Configure (string config);

private:

  // Load algorithm configuration
  void LoadConfig (void);

  const XSecIntegratorI *  fXSecIntegrator; // Numerical integrator (GSL)

};
  
}       // genie namespace
#endif  // _MARTINI_ERICSON_CHANFRAY_MARTEAU_MEC_PXSEC_2016_H_
