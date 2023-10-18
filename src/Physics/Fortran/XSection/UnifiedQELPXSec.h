//____________________________________________________________________________
/*!

\class    genie::UnifiedQELPXSec

\brief    Computes quasielastic neutrino-nucleus differential cross sections
          using a contraction of leptonic and hadronic tensors. Intended for
          use with a spectral function nuclear model.
          Is a concrete implementation of the XSecAlgorithmI interface. \n

\author   Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Acclerator Laboratory

\created  March 25, 2019

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _CBF_SPECTRAL_FUNC_CROSS_SECTION_H_
#define _CBF_SPECTRAL_FUNC_CROSS_SECTION_H_

#include <string>
#include <complex>

#include "Math/IFunction.h"

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Physics/QuasiElastic/XSection/LeptonTensor.h"
#include "Physics/QuasiElastic/XSection/QELFormFactors.h"
#include "Physics/NuclearState/NuclearModelI.h"
#include "Physics/NuclearState/PauliBlocker.h"
#include "Physics/QuasiElastic/XSection/ManualResponseTensor.h"
#include "TVector3.h"
namespace genie {

class QELFormFactorsModelI;
class QELFormFactors;
class XSecIntegratorI;

class UnifiedQELPXSec : public XSecAlgorithmI {

public:

  UnifiedQELPXSec();
  UnifiedQELPXSec(std::string config);
  virtual ~UnifiedQELPXSec();

  // XSecAlgorithmI interface implementation
  double XSec            (const Interaction* i, KinePhaseSpace_t k) const;
  double Integral        (const Interaction* i) const;
  bool ValidProcess      (const Interaction* i) const;

  // Override the Algorithm::Configure methods to load configuration
  // data to private data members
  void Configure (const Registry& config);
  void Configure (std::string param_set);

private:
  void LoadConfig (void);
  
  const QELFormFactorsModelI* fCCFormFactorsModel;
  const QELFormFactorsModelI* fNCFormFactorsModel;
  const QELFormFactorsModelI* fEMFormFactorsModel;
  mutable QELFormFactors fFormFactors;

  const NuclearModelI* fNuclModel;
  const PauliBlocker* fPauliBlocker;
  const XSecIntegratorI* fXSecIntegrator;
  double fCos8c2; ///< cos^2(cabibbo angle)
  double fXSecScale; ///< external xsec scaling factor
  bool fDoPauliBlocking;
  bool fDoqAlongZ;
  std::string fFortranTensorModel;
};

} // genie namespace

#endif
