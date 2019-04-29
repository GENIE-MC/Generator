//____________________________________________________________________________
/*!

\class    genie::RosenbluthPXSec

\brief    Differential cross section for charged lepton elastic scattering. \n
          Is a concrete implementation of the XSecAlgorithmI interface. \n

\ref      See for example: 
          R.Bradford, A.Bodek, H.Budd, J.Arrington, Nucl.Phys.B159 (2006) 127
          
\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  Sep 15, 2009

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _ROSENBLUTH_CROSS_SECTION_H_
#define _ROSENBLUTH_CROSS_SECTION_H_

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Physics/QuasiElastic/XSection/ELFormFactors.h"

namespace genie {

class XSecIntegratorI;

class RosenbluthPXSec : public XSecAlgorithmI {

public:
  RosenbluthPXSec();
  RosenbluthPXSec(string config);
  virtual ~RosenbluthPXSec();

  // XSecAlgorithmI interface implementation
  double XSec         (const Interaction * i, KinePhaseSpace_t k) const;
  double Integral     (const Interaction * i) const;
  bool   ValidProcess (const Interaction * i) const;

  // override the Algorithm::Configure methods to load configuration
  // data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

private:

  void LoadConfig(void);

  const   XSecIntegratorI *     fXSecIntegrator;
  const   ELFormFactorsModelI * fElFFModel;
  mutable ELFormFactors         fELFF;
  bool fCleanUpfElFFModel;
};

}       // genie namespace

#endif  
