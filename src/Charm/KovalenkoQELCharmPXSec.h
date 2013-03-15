//____________________________________________________________________________
/*!

\class    genie::KovalenkoQELCharmPXSec

\brief    Computes the QEL Charm Production Differential Cross Section
          using \b Kovalenko's duality model approach.
          It models the differential cross sections for: \n
             \li v + n \rightarrow mu- + Lambda_{c}^{+} (2285)
             \li v + n \rightarrow mu- + Sigma_{c}^{+}  (2455)
             \li v + p \rightarrow mu- + Sigma_{c}^{++} (2455)
          Is a concrete implementation of the XSecAlgorithmI interface.

\ref      S.G.Kovalenko, Sov.J.Nucl.Phys.52:934 (1990)

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  June 10, 2004

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _KOVALENKO_QEL_CHARM_PARTIAL_XSEC_H_
#define _KOVALENKO_QEL_CHARM_PARTIAL_XSEC_H_

#include "Base/XSecAlgorithmI.h"
#include "Numerical/GSFunc.h"

namespace genie {

class PDF;
class PDFModelI;
class IntegratorI;
class XSecIntegratorI;

class KovalenkoQELCharmPXSec : public XSecAlgorithmI {

public:
  KovalenkoQELCharmPXSec();
  KovalenkoQELCharmPXSec(string config);
  virtual ~KovalenkoQELCharmPXSec();

  //-- XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k) const;
  double Integral        (const Interaction * i) const;
  bool   ValidProcess    (const Interaction * i) const;
  bool   ValidKinematics (const Interaction * i) const;

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

private:
  void  LoadConfig (void);

  double ZR    (const Interaction * interaction)  const;
  double DR    (const Interaction * interaction)  const;
  double MRes  (const Interaction * interaction)  const;
  double ResDM (const Interaction * interaction)  const;
  double xiBar (double Q2, double Mnuc, double v) const;

  const PDFModelI *       fPDFModel;
  const IntegratorI *     fIntegrator;
  const XSecIntegratorI * fXSecIntegrator;

  double fMo;
  double fScLambdaP;
  double fScSigmaP;
  double fScSigmaPP;
  double fResDMLambda;
  double fResDMSigma;
};

} // genie namespace

//____________________________________________________________________________
/*!
\class    genie::KovQELCharmIntegrand

\brief    Auxiliary scalar function for the internal integration in Kovalenko
          QEL charm production cross section algorithm

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  February 20, 2006
*/
//____________________________________________________________________________

namespace genie {

class KovQELCharmIntegrand : public GSFunc
{
public:
  KovQELCharmIntegrand(PDF * pdf, double Q2, int nucleon_pdgc);
  ~KovQELCharmIntegrand();

  double operator () (const vector<double> & x);

private:
  PDF *  fPDF;
  double fQ2;
  int    fPdgC;
};

} // genie namespace

#endif  // _KOVALENKO_QEL_CHARM_PARTIAL_XSEC_H_


