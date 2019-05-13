//____________________________________________________________________________
/*!

\class    genie::KLVOxygenIBDPXSec

\brief    An implementation of the neutrino - Oxygen16 cross section.

\ref      E. Kolbe, K. Langanke, P. Vogel ; Phys. Rev. D66, 013007, 2002

\author   Corey Reed <cjreed \at nikhef.nl>
          Nikhef

\created  January 27, 2010

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _KLV_QUASIELASTIC_NU_OXYGEN_XSEC_H_
#define _KLV_QUASIELASTIC_NU_OXYGEN_XSEC_H_

#include "Framework/EventGen/XSecAlgorithmI.h"

class TSpline3;

namespace genie {

class KLVOxygenIBDPXSec : public XSecAlgorithmI {

public:
  static const double kO16NubarThr; //  energy threshold for 16O + nu_e_bar
  static const double kO16NuMinE;   //  minimum energy for 16O + nu_e
  static const double kMaxE;        //  maximum neutrino energy for this xsec model
  
  KLVOxygenIBDPXSec();
  KLVOxygenIBDPXSec(string config);
  virtual ~KLVOxygenIBDPXSec();

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
  
  void MakeAntiNuESpline(void);
  void MakeNuESpline(void);

  TSpline3* fXsplNue; //! a spline around the 16O+nu_e xsec points listed in the reference paper
  TSpline3* fXsplNuebar; //! a spline around the 16O+nu_e_bar xsec points listed in the reference paper

public:
  ClassDef(KLVOxygenIBDPXSec, 1) // Oxygen16 - (anti)neutrino cross section estimator based on a Kolbe/Langanke/Vogel paper
};
}       // genie namespace
#endif  // _KLV_QUASIELASTIC_NU_OXYGEN_XSEC_H_
