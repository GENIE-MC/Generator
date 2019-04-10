//____________________________________________________________________________
/*!

\class    genie::LwlynSmithQELCCPXSec

\brief    Computes neutrino-nucleon(nucleus) QELCC differential cross section
          Is a concrete implementation of the XSecAlgorithmI interface. \n

\ref      C.H.Llewellyn Smith, Physics Reports (Sect. C of Physics Letters) 3,
          No. 5  (1972) 261-379

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 05, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _LLEWELLYN_SMITH_QELCC_CROSS_SECTION_H_
#define _LLEWELLYN_SMITH_QELCC_CROSS_SECTION_H_

#include "Physics/NuclearState/NuclearModelI.h"
#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Physics/QuasiElastic/XSection/QELFormFactors.h"
#include "Physics/QuasiElastic/XSection/QELUtils.h"
#include "Physics/NuclearState/PauliBlocker.h"

namespace genie {

class QELFormFactorsModelI;
class XSecIntegratorI;

class LwlynSmithQELCCPXSec : public XSecAlgorithmI {

public:
  LwlynSmithQELCCPXSec();
  LwlynSmithQELCCPXSec(string config);
  virtual ~LwlynSmithQELCCPXSec();

  // XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k) const;
  double Integral        (const Interaction * i) const;
  bool   ValidProcess    (const Interaction * i) const;

  // Override the Algorithm::Configure methods to load configuration
  // data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

private:
  double FullDifferentialXSec(const Interaction * i) const;

  void LoadConfig (void);

  mutable QELFormFactors       fFormFactors;      ///<
  const QELFormFactorsModelI * fFormFactorsModel; ///<
  const XSecIntegratorI *      fXSecIntegrator;   ///<
  double                       fCos8c2;           ///< cos^2(cabibbo angle)

  double                       fXSecScale;        ///< external xsec scaling factor

  // Variables for integrating
  const NuclearModelI *        fNuclModel;
  bool   fLFG;                         ///< If the nuclear model is lfg alway average over nucleons
  bool   fDoAvgOverNucleonMomentum;    ///< Average cross section over hit nucleon monentum?
  double fEnergyCutOff;                ///< Average only for energies below this cutoff defining
                                       ///< the region where nuclear modeling details do matter

  /// Enum specifying the method to use when calculating the binding energy of
  /// the initial hit nucleon during spline generation
  QELEvGen_BindingMode_t fIntegralNucleonBindingMode;

  /// Whether to apply Pauli blocking in FullDifferentialXSec
  bool fDoPauliBlocking;
  /// The PauliBlocker instance to use to apply that correction
  const genie::PauliBlocker* fPauliBlocker;
};

}       // genie namespace

#endif
