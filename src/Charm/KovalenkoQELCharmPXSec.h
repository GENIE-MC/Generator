//____________________________________________________________________________
/*!

\class    genie::KovalenkoQELCharmPXSec

\brief    Computes the QEL Charm Production Differential Cross Section
          using \b Kovalenko's duality model approach.

          The computed differential cross section is the Dxsec = dxsec/dQ^2
          where \n
            \li \c Q2 is the momentum transfer.

          It models the differential cross sections for: \n
             \li v + n \rightarrow mu- + Lambda_{c}^{+} (2285)
             \li v + n \rightarrow mu- + Sigma_{c}^{+}  (2455)
             \li v + p \rightarrow mu- + Sigma_{c}^{++} (2455)

          Is a concrete implementation of the XSecAlgorithmI interface.

\ref      S.G.Kovalenko, Sov.J.Nucl.Phys.52:934 (1990)

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  June 10, 2004

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

class KovalenkoQELCharmPXSec : public XSecAlgorithmI {

public:
  KovalenkoQELCharmPXSec();
  KovalenkoQELCharmPXSec(string config);
  virtual ~KovalenkoQELCharmPXSec();

  //-- XSecAlgorithmI interface implementation
  double XSec            (const Interaction * interaction) const;
  bool   ValidProcess    (const Interaction * interaction) const;
  bool   ValidKinematics (const Interaction * interaction) const;

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

private:

  void  LoadConfigData (void);
  void  LoadSubAlg     (void);

  double ZR       (const Interaction * interaction) const;
  double DR       (const Interaction * interaction, bool norm = false) const;
  double MRes     (const Interaction * interaction) const;
  double vR_minus (const Interaction * interaction) const;
  double vR_plus  (const Interaction * interaction) const;
  double SumF2    (const Interaction * interaction) const;
  double ResDM    (const Interaction * interaction) const;
  double xiBar    (const Interaction * interaction, double v) const;

  const PDFModelI *   fPDFModel;
  const IntegratorI * fIntegrator;

  double fQ2min;
  double fQ2max;
  double fMo;
  double fF2LambdaP;
  double fF2SigmaP;
  double fF2SigmaPP;
  double fResDMLambda;
  double fResDMSigma;
};

//____________________________________________________________________________
/*!
\class    genie::KovQELCharmIntegrand

\brief    Auxiliary scalar function for the internal integration in Kovalenko
          QEL charm production cross section algorithm

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  February 20, 2006
*/
//____________________________________________________________________________

class KovQELCharmIntegrand : public GSFunc
{
public:
  KovQELCharmIntegrand(PDF * pdf, double Q2, int nucleon_pdgc, bool norm);
  ~KovQELCharmIntegrand();

  double operator () (const vector<double> & x);

private:
  PDF *  fPDF;
  double fQ2;
  int    fPdgC;
  bool   fNorm;
};

}       // genie namespace

#endif  // _KOVALENKO_QEL_CHARM_PARTIAL_XSEC_H_


