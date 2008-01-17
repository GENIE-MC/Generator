//____________________________________________________________________________
/*!

\class    genie::BardinDISPXSec

\brief    Computes, the differential cross section for neutrino DIS including
          radiative corrections according to the \b Bardin-Dokuchaeva model.
          Is a concrete implementation of the XSecAlgorithmI interface.

\ref      D.Yu.Bardin, V.A.Dokuchaeva, "On the radiative corrections to the
          Neutrino Deep Inelastic Scattering", JINR-E2-86-260, Apr. 1986

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  July 06, 2004

\cpright  Copyright (c) 2003-2008, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _BARDIN_DIS_PARTIAL_XSEC_H_
#define _BARDIN_DIS_PARTIAL_XSEC_H_

#include "Base/XSecAlgorithmI.h"
#include "PDF/PDF.h"

namespace genie {

class PDFModelI;
class XSecIntegratorI;

class BardinDISPXSec : public XSecAlgorithmI {

public:
  BardinDISPXSec();
  BardinDISPXSec(string config);
  virtual ~BardinDISPXSec();

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

  //-- Auxiliary methods used for the Bardin-Dokuchaeva xsec calculation.
  //   xi is the running variable in the integration for computing
  //   the differential cross section d^2(xsec) / dy dx
  double PhiCCi   (double xi, const Interaction * interaction) const;
  double Ii       (double xi, const Interaction * interaction) const;
  double U        (double xi, const Interaction * interaction) const;
  double tau      (double xi, const Interaction * interaction) const;
  double St       (double xi, const Interaction * interaction) const;
  double Su       (double xi, const Interaction * interaction) const;
  double S        (const Interaction * interaction) const;
  double DeltaCCi (const Interaction * interaction) const;
  double Sq       (const Interaction * interaction) const;
  double PDFFunc  (const PDF & pdf, int pgdc) const;

  double fMqf;
  double fVud;
  double fVud2;

  const PDFModelI *       fPDFModel;
  const XSecIntegratorI * fXSecIntegrator;
};

}       // genie namespace
#endif  // _BARDIN_DIS_PARTIAL_XSEC_H_
