//____________________________________________________________________________
/*!

\class    genie::COHElasticPXSec

\brief    Differential cross section for v+As coherent elastic scattering.
          Is a concrete implementation of the XSecAlgorithmI interface. 

\ref      K.Patton, J.Engel, G.McLaughlin, and N.Schunck, arXiv:1207.0693v1

\author   Costas Andreopoulos <costas.andreopoulos \at cern.ch>
          University of Liverpool & STFC Rutherford Appleton Laboratory

\created  July 12, 2019

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _COHERENT_ELASTIC_PXSEC_H_
#define _COHERENT_ELASTIC_PXSEC_H_

#include "Framework/EventGen/XSecAlgorithmI.h"

namespace genie {

class XSecIntegratorI;

class COHElasticPXSec : public XSecAlgorithmI {

public:
  COHElasticPXSec();
  COHElasticPXSec(string config);
  virtual ~COHElasticPXSec();

  // XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k) const;
  double Integral        (const Interaction * i) const;
  bool   ValidProcess    (const Interaction * i) const;

  // Override the Algorithm::Configure methods to load configuration
  // data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

private:

  double NuclearDensityMoment(int A, int k) const;

  void   LoadConfig(void);

  const XSecIntegratorI * fXSecIntegrator;  ///< cross section integrator
  double fSin2thw;                          ///< sin^2(weinberg angle)
};

}       // genie namespace

#endif  // _COHERENT_ELASTIC_PXSEC_H_
