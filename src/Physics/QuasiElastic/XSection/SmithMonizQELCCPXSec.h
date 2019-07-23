//____________________________________________________________________________
/*!

\class    genie::SmithMonizQELCCPXSec

\brief    Computes neutrino-nucleon(nucleus) QELCC differential cross section.
          Is a concrete implementation of the XSecAlgorithmI interface. 

\ref      [1] R.A.Smith and E.J.Moniz, Nuclear Physics  B43, (1972) 605-622 \n
          [2] K.S. Kuzmin, V.V. Lyubushkin, V.A.Naumov Eur. Phys. J. C54, (2008) 517-538

\author   Igor Kakorin <kakorin@jinr.ru>, Joint Institute for Nuclear Research \n
          adapted from  fortran code provided by \n
          Konstantin Kuzmin <kkuzmin@theor.jinr.ru>, Joint Institute for Nuclear Research \n
          Vladimir Lyubushkin, Joint Institute for Nuclear Research \n
          Vadim Naumov <vnaumov@theor.jinr.ru>, Joint Institute for Nuclear Research  \n
          based on code of \n
          Costas Andreopoulos <costas.andreopoulos@stfc.ac.uk>, University of Liverpool & STFC Rutherford Appleton Lab
          

\created  May 05, 2017

\cpright  Copyright (c) 2003-2017, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _SMITH_MONITZ_QELCC_CROSS_SECTION_H_
#define _SMITH_MONITZ_QELCC_CROSS_SECTION_H_

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Physics/QuasiElastic/XSection/QELFormFactors.h"
#include "Physics/QuasiElastic/XSection/SmithMonizUtils.h"


namespace genie {

class QELFormFactorsModelI;
class XSecIntegratorI;

class SmithMonizQELCCPXSec : public XSecAlgorithmI {

public:
  SmithMonizQELCCPXSec();
  SmithMonizQELCCPXSec(string config);
  virtual ~SmithMonizQELCCPXSec();

  // XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t kps) const;
  double Integral        (const Interaction * i) const;
  bool   ValidProcess    (const Interaction * i) const;

  // Override the Algorithm::Configure methods to load configuration
  // data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

private:
  mutable SmithMonizUtils * sm_utils;
  
  void   LoadConfig (void);
  double d3sQES_dQ2dvdkF_SM (const Interaction * interaction) const;
  double dsQES_dQ2_SM(const Interaction * interaction) const;
  double d2sQES_dQ2dv_SM(const Interaction * i) const;
  
  double                       fXSecScale;        ///< external xsec scaling factor
  mutable QELFormFactors       fFormFactors;      
  const QELFormFactorsModelI * fFormFactorsModel; 
  const XSecIntegratorI *      fXSecIntegrator;   
  double                       fVud2;             ///< |Vud|^2(square of magnitude ud-element of CKM-matrix)
  mutable int                          fn_NT;
  mutable double                       fQ2;
  mutable double                       fv;
  mutable double                       fE_nu;
  mutable double                       fE_lep;
  mutable double                       fmm_ini;
  mutable double                       fmm_fin;
  mutable double                       fm_tar;
  mutable double                       fmm_tar;
  mutable double                       fk1;
  mutable double                       fk2;
  mutable double                       fk7;
  mutable double                       fqv;
  mutable double                       fqqv;
  mutable double                       fcosT_k;
  mutable double                       fF_V;
  mutable double                       fF_M;
  mutable double                       fF_A;
  mutable double                       fF_P;
  mutable double                       fFF_V;
  mutable double                       fFF_M;
  mutable double                       fFF_A;
  mutable double                       fW_1;
  mutable double                       fW_2;
  mutable double                       fW_3;
  mutable double                       fW_4;
  mutable double                       fW_5;
   
  
};


} // genie namespace

#endif  //_SMITH_MONITZ_QELCC_CROSS_SECTION_H_

