//____________________________________________________________________________
/*!

\class    genie::IMDAnnihilationPXSec

\brief    nuebar + e- -> mu- + numubar [CC]
          scattering differential cross section \n

          Is a concrete implementation of the XSecAlgorithmI interface. \n

\ref      W.J.Marciano and Z.Parsa, Neutrino-electron scattering theory,
          J.Phys.G: Nucl.Part.Phys. 29 (2003) 2629-2645

\author   Rosen Matev (r.matev@gmail.com)

\created  October 3, 2011

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _IMD_ANNIHILATION_PXSEC_H_
#define _IMD_ANNIHILATION_PXSEC_H_

#include "Framework/EventGen/XSecAlgorithmI.h"

namespace genie {

class IntegratorI;
class XSecIntegratorI;

class IMDAnnihilationPXSec : public XSecAlgorithmI {

public:
  IMDAnnihilationPXSec();
  IMDAnnihilationPXSec(string config);
  virtual ~IMDAnnihilationPXSec();

  //-- XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k) const;
  double Integral        (const Interaction * i) const;
  bool   ValidProcess    (const Interaction * i) const;
  bool   ValidKinematics (const Interaction * i) const;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void LoadConfig (void);

  const XSecIntegratorI * fXSecIntegrator;

};

}       // genie namespace
#endif  // _IMD_ANNIHILATION_PXSEC_H_
