//____________________________________________________________________________
/*!

\class    genie::AhrensDMELPXSec

\brief    Differential cross section for DM+N elastic scattering. \n
          Is a concrete implementation of the XSecAlgorithmI interface. \n

\ref      R.E.Hendrick and L.Li, Phys.Rev.D 19:779 (1979)
          L.A.Ahrens et al., Phys.Rev.D 35:785 (1987)

\author   Joshua Berger <jberger \at physics.wisc.edu>
          University of Wisconsin-Madison

          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  September 4, 2017

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _AHRENS_DMEL_CROSS_SECTION_H_
#define _AHRENS_DMEL_CROSS_SECTION_H_

#include "Framework/EventGen/XSecAlgorithmI.h"

namespace genie {

class XSecIntegratorI;

class AhrensDMELPXSec : public XSecAlgorithmI {

public:
  AhrensDMELPXSec();
  AhrensDMELPXSec(string config);
  virtual ~AhrensDMELPXSec();

  // XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k) const;
  double Integral        (const Interaction * i) const;
  bool   ValidProcess    (const Interaction * i) const;

  // Override the Algorithm::Configure methods to load configuration
  // data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

private:
  void LoadConfig(void);

  const XSecIntegratorI * fXSecIntegrator;

  double fQchiV;
  double fQchiA;
  double fQchiS;
  double fQuV;
  double fQuA;
  double fQdV;
  double fQdA;
  double fQsV;
  double fQsA;
  double fMa2;
  double fMv2;
  double fMp2;
  double fMpi2;
  double fMeta2;
  double fMuP;
  double fMuN;
  double fDelu;
  double fDeld;
  double fDels;
  double fDeluP;
  double fDeldP;
  double fDelsP;
  int    fVelMode;
  double fMedMass;
  double fgZp;
};

}       // genie namespace

#endif  
