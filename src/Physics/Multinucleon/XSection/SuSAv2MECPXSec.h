//____________________________________________________________________________
/*!

\class    genie::SuSAv2MECPXSec

\brief    Computes the SuSAv2-MEC model differential cross section.
          Uses precomputed hadron tensor tables.
          Is a concrete implementation of the XSecAlgorithmI interface.

\author   Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Acclerator Laboratory

          Stephen Dolan <stephen.joseph.dolan \at cern.ch>
          European Organization for Nuclear Research (CERN)

\ref      G.D. Megias et al., "Meson-exchange currents and quasielastic
          predictions for charged-current neutrino-12C scattering in the
          superscaling approach," PRD 91 (2015) 073004

\created  November 2, 2018

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _SUSAV2_MEC_PXSEC_H_
#define _SUSAV2_MEC_PXSEC_H_

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Physics/HadronTensors/HadronTensorModelI.h"
#include "Physics/Common/XSecScaleI.h"
#include "Physics/Common/QvalueShifter.h"

namespace genie {

class XSecIntegratorI;

class SuSAv2MECPXSec : public XSecAlgorithmI {

public:

  SuSAv2MECPXSec();
  SuSAv2MECPXSec(string config);
  virtual ~SuSAv2MECPXSec();

  // XSecAlgorithmI interface implementation
  double XSec(const Interaction* i, KinePhaseSpace_t k) const;
  double Integral(const Interaction* i) const;
  bool   ValidProcess(const Interaction* i) const;

  // override the Algorithm::Configure methods to load configuration
  // data to private data members
  void Configure (const Registry & config);
  void Configure (string config);

  // Method specifically for evaluating np/pp pair probabilities
  double PairRatio(const Interaction* i,
    const std::string& final_state_ratio = "pnFraction") const;

private:

  /// Load algorithm configuration
  void LoadConfig (void);

  // Calculate Qvalue Shift for susa:
  double Qvalue(const Interaction & interaction ) const ;

  /// External scaling factor for this cross section
  double fXSecCCScale;
  double fXSecNCScale;
  double fXSecEMScale;

  const genie::HadronTensorModelI* fHadronTensorModel;

  // Fermi momentum table used for scaling
  string fKFTable;

  // Binding energies:
  double fEbHe;
  double fEbLi;
  double fEbC;
  double fEbO;
  double fEbMg;
  double fEbAr;
  double fEbCa;
  double fEbFe;
  double fEbNi;
  double fEbSn;
  double fEbAu;
  double fEbPb;

  /// GSL numerical integrator
  const XSecIntegratorI*  fXSecIntegrator;

  const XSecScaleI * fMECScaleAlg ; // Optional algorithm to scale the xsec as a function of W
  const QvalueShifter * fQvalueShifter ; // Optional algorithm to retrieve the qvalue shift for a given target

};

} // genie namespace
#endif // _SUSAV2_MEC_PXSEC_H_
