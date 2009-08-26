//____________________________________________________________________________
/*!

\class    genie::MECPXSec

\brief    Computes the MEC differential cross section.
          Is a concrete implementation of the XSecAlgorithmI interface. \n

\ref      

\author   

\created  May 05, 2009

\cpright  Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _MEC_PXSEC_H_
#define _MEC_PXSEC_H_

#include "Base/XSecAlgorithmI.h"

namespace genie {

class MECPXSec : public XSecAlgorithmI {

public:
  MECPXSec();
  MECPXSec(string config);
  virtual ~MECPXSec();

  //-- XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k) const;
  double Integral        (const Interaction * i) const;
  bool   ValidProcess    (const Interaction * i) const;

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

private:
  void LoadConfig (void);

  double fMq2d;   ///< toy model param: `mass' in dipole (Q2 - dependence) form factor (GeV)
  double fMass;   ///< toy model param: peak  of W distribution (GeV)
  double fWidth;  ///< toy model param: width of W distribution (GeV)
  double fEc;     ///< toy model param: low energy cutoff (GeV) 
  double fNorm;   ///< toy model param: MEC rate / nucleon (in 1E-38cm^2)
};

}       // genie namespace
#endif  // _MEC_PARTIAL_XSEC_H_
