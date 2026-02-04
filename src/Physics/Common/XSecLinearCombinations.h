//____________________________________________________________________________
/*!

\class    genie::XSecLinearCombinations

\brief    Computes the xsec as a linear combination of different XSecSlgorithmI
          See GENIE docdb 252

\author   Code contributed by J.Tena Vidal and M.Roda

\created  Mar 12, 2020

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _XSEC_LINEAR_COMBINATIONS_H_
#define _XSEC_LINEAR_COMBINATIONS_H_

#include "Framework/EventGen/XSecAlgorithmI.h"

namespace genie {

class XSecIntegratorI;

class XSecLinearCombinations : public XSecAlgorithmI {

public:
  XSecLinearCombinations();
  XSecLinearCombinations(string config);
  virtual ~XSecLinearCombinations();

  // XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k) const;
  double Integral        (const Interaction * i) const;
  bool   ValidProcess    (const Interaction * i) const;

  // override the Algorithm::Configure methods to load configuration
  // data to private data members
  void Configure (const Registry & config);
  void Configure (string config);

 protected:

  // Load algorithm configuration
  void LoadConfig (void);

 private: 

  std::vector<const XSecAlgorithmI*> fXSections ;
  std::vector<double> fLinearCoefficients ;

};

}       // genie namespace
#endif  // _XSEC_LINEAR_COMBINATIONS_H_
