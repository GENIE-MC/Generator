//____________________________________________________________________________
/*!

\class    genie::COHElasticXSec

\brief    Computes the cross section for coherent elastic scattering.\n
          Is a concrete implementation of the XSecIntegratorI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at cern.ch>
          University of Liverpool & STFC Rutherford Appleton Laboratory

\created  July 12, 2019

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _COH_ELASTIC_XSEC_H_
#define _COH_ELASTIC_XSEC_H_

#include "Physics/XSectionIntegration/XSecIntegratorI.h"

namespace genie {

class COHElasticXSec : public XSecIntegratorI {

public:
  COHElasticXSec();
  COHElasticXSec(string config);
  virtual ~COHElasticXSec();

  // XSecIntegratorI interface implementation
  double Integrate(const XSecAlgorithmI * model, const Interaction * i) const;

  // Overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

 private:

   void LoadConfig (void);

//      double fQ2Min;  ///< lower bound of integration for Q^2 in Berger-Sehgal Model
//      double fQ2Max;  ///< upper bound of integration for Q^2 in Berger-Sehgal Model
//      double fTMax;   ///< upper bound for t = (q - p_pi)^2
};

}       // genie namespace

#endif  // _COH_ELASTIC_XSEC_H_
