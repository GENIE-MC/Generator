//____________________________________________________________________________
/*!

\class    genie::GlashowResXSec

\brief    Nuebar cross section at the Glashow resonance (nuebar + e- -> W-).\n
          Is a concrete implementation of the XSecAlgorithmI interface. \n

\ref      T.K.Gaisser, F.Halzen and T.Stanev, Physics Reports 258:173 (1995)

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  May 04, 2005

\cpright  Copyright (c) 2003-2008, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _GLASHOW_RES_XSEC_H_
#define _GLASHOW_RES_XSEC_H_

#include "Base/XSecAlgorithmI.h"

namespace genie {

class IntegratorI;

class GlashowResXSec : public XSecAlgorithmI {

public:
  GlashowResXSec();
  GlashowResXSec(string config);
  virtual ~GlashowResXSec();

  //-- XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k=kPSfE) const;
  double Integral        (const Interaction * i) const;
  bool   ValidProcess    (const Interaction * i) const;
  bool   ValidKinematics (const Interaction * i) const;
};

}       // genie namespace

#endif  // _GLASHOW_RES_XSEC_H_
