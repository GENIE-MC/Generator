//____________________________________________________________________________
/*!

\class    genie::MartiniQELPXSec

\brief    Computes the Martini-QE model differential cross section.
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

\cpright  Copyright (c) 2003-2023, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _MARTINI_QE_PXSEC_H_
#define _MARTINI_QE_PXSEC_H_

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Physics/HadronTensors/HadronTensorI.h"
#include "Physics/HadronTensors/HadronTensorModelI.h"
#include "Physics/Common/QvalueShifter.h"

namespace genie {

class XSecIntegratorI;

class MartiniQELPXSec : public XSecAlgorithmI {

public:

  MartiniQELPXSec();
  MartiniQELPXSec(string config);
  virtual ~MartiniQELPXSec();

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

  /// External scaling factor for this cross section
  double fXSecScale;

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
#endif // _MARTINI_QE_PXSEC_H_
