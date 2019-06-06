//____________________________________________________________________________
/*!

\class    genie::AlamSimoAtharVacasSKPXSec2014

\brief    Differential cross section model for single kaon production.

\ref      Weak Kaon Production off the Nucleon
          M. Rafi Alam, I. Ruiz Simo, M. Sajjad Athar and M.J. Vicente Vacas 
          Physical Review D82 (2010) 033001; arXiv:1004.5484 [hep-ph]

\author   Chris Marshall and Martti Nirkko

\created  2014-02-14

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _ALAM_SIMO_ATHAR_VACAS_SINGLE_KAON_PXSEC_2014_H_
#define _ALAM_SIMO_ATHAR_VACAS_SINGLE_KAON_PXSEC_2014_H_

#include "Framework/EventGen/XSecAlgorithmI.h"

namespace genie {

class XSecIntegratorI;

class AlamSimoAtharVacasSKPXSec2014 : public XSecAlgorithmI {

public:
  AlamSimoAtharVacasSKPXSec2014();
  AlamSimoAtharVacasSKPXSec2014(string config);
  virtual ~AlamSimoAtharVacasSKPXSec2014();

  // XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k) const;
  double Integral        (const Interaction * i) const;
  bool   ValidProcess    (const Interaction * i) const;

  // Override the Algorithm::Configure methods to load configuration
  // data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

private:

  void LoadConfig(void);

  const XSecIntegratorI * fXSecIntegrator;  ///< cross section integrator
  
  // Calculate matrix elements
  double Amatrix_NN(double theta, double phikq) const;
  double Amatrix_NP(double theta, double phikq) const;
  double Amatrix_PP(double theta, double phikq) const;
  
  // Physics parameters set globally
  // The names of these parameters in the code match the convention in the original FORTRAN code
  // They are used in the matrix element calculations which are thousands of lines long
  // We try to change as little as possible, so keep these names
  double amLam, amEta, Vus, fpi, d, f, g, amup, amun, Fm1, Fm2;
  
  // Interaction parameters set locally
  // These are set event-by-event, and used in the matrix element calculation which is thousands 
  // of lines long. The names are kept from the original FORTRAN code
  // They are mutable because the XSec routine must be const in GENIE, but they are needed in the 
  // matrix element calculations which are separate functions
  mutable int leptonPDG, reactionType;
  mutable double aml, amSig, amk, ampi, am;
  mutable double Enu, Ekaon, pkvec;
  
  // Output calculated by cross-section function
  // essentially returned by reference by the matrix element calculations
  // mutable because XSec routine must be const in GENIE
  mutable double Elep, alepvec, aqvec, angkq, aq0;

};

}       // genie namespace
#endif  // _ALAM_SIMO_ATHAR_VACAS_SINGLE_KAON_PXSEC_2014_H_
