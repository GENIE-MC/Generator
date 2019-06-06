//____________________________________________________________________________
/*!

\class    genie::H3AMNuGammaPXSec

\brief    An anomaly-mediated neutrino-photon interaction cross section model
          Is a concrete implementation of the XSecAlgorithmI interface. 

\ref      J.A.Harvey, C.T.Hill and R.J.Hill, PRL99, 261601 (2007)

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  February 15, 2008

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _H3_ANOMALY_MEDIATED_NUGAMMA_PXSEC_H_
#define _H3_ANOMALY_MEDIATED_NUGAMMA_PXSEC_H_

#include "Framework/EventGen/XSecAlgorithmI.h"

namespace genie {

class H3AMNuGammaPXSec : public XSecAlgorithmI {

public:
  H3AMNuGammaPXSec();
  H3AMNuGammaPXSec(string config);
  virtual ~H3AMNuGammaPXSec();

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

  double fGw; 
};

}       // genie namespace
#endif  // _H3_ANOMALY_MEDIATED_NUGAMMA_PXSEC_H_
