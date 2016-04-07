//____________________________________________________________________________
/*!

\class    genie::MECTensorPXSec

\brief    Computes the MEC differential cross section.
          Is a concrete implementation of the XSecAlgorithmI interface. \n

\ref      

\author   

\created  

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________
#ifndef _MEC_TENSOR_PXSEC_H_
#define _MEC_TENSOR_PXSEC_H_

#include "Base/XSecAlgorithmI.h"

namespace genie {

class MECTensorPXSec : public XSecAlgorithmI {


public:
  MECTensorPXSec();
  MECTensorPXSec(string config);
  virtual ~MECTensorPXSec();

  // XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k) const;
  double Integral        (const Interaction * i) const;
  bool   ValidProcess    (const Interaction * i) const;

  // override the Algorithm::Configure methods to load configuration
  // data to private data members
  void Configure (const Registry & config);

private:

  void LoadConfig (void);
  
  //int fTensorModel; // mec model to use for tensors 1=nieves 2=martini

  };
  
}       // genie namespace
#endif  // _MEC_PARTIAL_XSEC_H_
