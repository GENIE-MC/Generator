//____________________________________________________________________________
/*!

\class    genie::DfrcPXSec

\brief    Neutrino diffractive scattering cross section.

\ref      

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  

\cpright  Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _DIFFRACTIVE_PXSEC_H_
#define _DIFFRACTIVE_PXSEC_H_

#include "Base/XSecAlgorithmI.h"

namespace genie {

class DfrcPXSec : public XSecAlgorithmI {

public:
  DfrcPXSec();
  DfrcPXSec(string config);
  virtual ~DfrcPXSec();

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

};

}       // genie namespace
#endif  // _DIFFRACTIVE_PXSEC_H_
