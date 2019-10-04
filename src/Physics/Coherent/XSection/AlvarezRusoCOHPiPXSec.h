//____________________________________________________________________________
/*!

\class    genie::AlvarezRusoCOHPiPXSec

\brief    Implementation of the Alvarez-Ruso coherent pion production model

          Is a concrete implementation of the XSecAlgorithmI interface.

\ref      

\author   Steve Dennis
          University of Warwick, Rutherford Appleton Laboratory

\created  October 5, 2012

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _ALVAREZ_RUSO_COH_XSEC_H_
#define _ALVAREZ_RUSO_COH_XSEC_H_

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Physics/Coherent/XSection/AlvarezRusoCOHPiPDXSec.h"

namespace genie {

class XSecIntegratorI; 
class Interaction;

class AlvarezRusoCOHPiPXSec : public XSecAlgorithmI {

public:
  AlvarezRusoCOHPiPXSec();
  AlvarezRusoCOHPiPXSec(string config);
  virtual ~AlvarezRusoCOHPiPXSec();

  //-- XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k) const;
  double Integral        (const Interaction * i) const;
  bool   ValidProcess    (const Interaction * i) const;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void LoadConfig(void);

  //-- private data members loaded from config Registry or set to defaults

  const XSecIntegratorI * fXSecIntegrator;
  
  mutable alvarezruso::AlvarezRusoCOHPiPDXSec * fMultidiff;
  mutable const Interaction * fLastInteraction;
  //Parameters
  //bool fUseLookupTable;
  //double fa4;
  //double fa5;
  //double fb4;
  //double fb5;
  //double ffPi;
  //double ffStar;
  //double fMa;
  //double fRo;
};

}      // genie namespace

#endif  // _REIN_SEGHAL_COHPI_PXSEC_H_
