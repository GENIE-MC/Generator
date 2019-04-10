//____________________________________________________________________________
/*!

\class    genie::GLRESPXSec

\brief    Nuebar cross section at the Glashow resonance (nuebar + e- -> W-).
          Is a concrete implementation of the XSecAlgorithmI interface. 

\ref      T.K.Gaisser, F.Halzen and T.Stanev, Physics Reports 258:173 (1995)

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 04, 2005

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _GLASHOW_RESONANCE_PXSEC_H_
#define _GLASHOW_RESONANCE_PXSEC_H_

#include "Framework/EventGen/XSecAlgorithmI.h"

namespace genie {

class GLRESPXSec : public XSecAlgorithmI {

public:
  GLRESPXSec ();
  GLRESPXSec (string config);
  virtual ~GLRESPXSec ();

  // XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k=kPSfE) const;
  double Integral        (const Interaction * i) const;
  bool   ValidProcess    (const Interaction * i) const;
  bool   ValidKinematics (const Interaction * i) const;
};

}       // genie namespace

#endif  // _GLASHOW_RESONANCE_XSEC_H_
