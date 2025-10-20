//____________________________________________________________________________
/*!

\class    genie::NormXSec

\brief    Normalization channel. Its main property is a constant cross section 
          per nucleon over the whole energy range.

\ref      [1] GENIE docdb 297


\author   Igor Kakorin <kakorin@jinr.ru>
          Joint Institute for Nuclear Research

\created  May 16, 2022

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org      
    

*/
//____________________________________________________________________________

#ifndef _NORM_CROSS_SECTION_H_
#define _NORM_CROSS_SECTION_H_
#include "Framework/EventGen/XSecAlgorithmI.h"

namespace genie {


class NormXSec : public XSecAlgorithmI {

public:
  NormXSec();
  NormXSec(string config);
  virtual ~NormXSec();

  // XSecAlgorithmI interface implementation
  double XSec             (const Interaction * i, KinePhaseSpace_t k) const;
  double Integral         (const Interaction * i) const;
  bool   ValidProcess     (const Interaction * i) const;
  bool   ValidKinematics  (const Interaction * i) const;

  // Override the Algorithm::Configure methods to load configuration
  // data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

private:
  void LoadConfig (void);
  double fNormScale;
};

}       // genie namespace

#endif
