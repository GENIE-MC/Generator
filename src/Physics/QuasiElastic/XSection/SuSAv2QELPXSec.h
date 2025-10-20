//____________________________________________________________________________
/*!

\class    genie::SuSAv2QELPXSec

\brief    Computes the SuSAv2-QE model differential cross section.
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

#ifndef _SUSAV2_QE_PXSEC_H_
#define _SUSAV2_QE_PXSEC_H_

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Physics/HadronTensors/HadronTensorI.h"
#include "Physics/HadronTensors/HadronTensorModelI.h"
#include "Physics/Common/QvalueShifter.h"

namespace genie {

class XSecIntegratorI;

class SuSAv2QELPXSec : public XSecAlgorithmI {

public:

  SuSAv2QELPXSec();
  SuSAv2QELPXSec(string config);
  virtual ~SuSAv2QELPXSec();

  // XSecAlgorithmI interface implementation
  double XSec(const Interaction* i, KinePhaseSpace_t k) const;
  double Integral(const Interaction* i) const;
  bool   ValidProcess(const Interaction* i) const;

  // override the Algorithm::Configure methods to load configuration
  // data to private data members
  void Configure (const Registry & config);
  void Configure (string config);

private:

  /// Load algorithm configuration
  void LoadConfig (void);

  // Apply scaling of the first kind to the xsec
  // (A-scaling with tuning based on Fermi momentum)
  double XSecScaling(double xsec, const Interaction* i, int target_pdg, int tensor_pdg, bool need_to_scale) const;


  /// External scaling factor for this cross section
  double fXSecCCScale;
  double fXSecNCScale;
  double fXSecEMScale;

  /// Blending start/end q0 value for the combined models
  int blendMode;
  double q0BlendStart;
  double q0BlendEnd;
  double qBlendDel;
  double qBlendRef;

  // Which model to run
  enum modelType {
    kMd_Undefined = 0,
    kMd_SuSAv2 = 1,
    kMd_CRPA = 2,
    kMd_HF = 3,
    kMd_CRPASuSAv2Hybrid = 4,
    kMd_HFSuSAv2Hybrid = 5,
    kMd_CRPAPW = 6,
    kMd_HFPW = 7,
    kMd_CRPAPWSuSAv2Hybrid = 8,
    kMd_HFPWSuSAv2Hybrid = 9,
    kMd_SuSAv2Blend = 10
  };

  modelType modelConfig;

  const HadronTensorModelI* fHadronTensorModel;

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

  /// Alternate cross section model for free nucleon targets
  const XSecAlgorithmI* fFreeNucleonXSecAlg;

  const QvalueShifter * fQvalueShifter ; // Gives the option to retrieve a qvalue shift for a given target
};

} // genie namespace
#endif // _SUSAV2_QE_PXSEC_H_
