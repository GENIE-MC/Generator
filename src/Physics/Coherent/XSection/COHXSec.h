//____________________________________________________________________________
/*!

\class    genie::COHXSec

\brief    Computes the cross section for COH neutrino-nucleus pi production.\n
          Is a concrete implementation of the XSecIntegratorI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 04, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _COH_XSEC_H_
#define _COH_XSEC_H_

#include "Physics/XSectionIntegration/XSecIntegratorI.h"

namespace genie {

  class COHXSec : public XSecIntegratorI {

    public:
      COHXSec();
      COHXSec(string config);
      virtual ~COHXSec();

      // XSecIntegratorI interface implementation
      double Integrate(const XSecAlgorithmI * model, const Interaction * i) const;

      // Overload the Algorithm::Configure() methods to load private data
      // members from configuration options
      void Configure(const Registry & config);
      void Configure(string config);

    private:
      void LoadConfig (void);

      double fQ2Min;  ///< lower bound of integration for Q^2 in Berger-Sehgal Model
      double fQ2Max;  ///< upper bound of integration for Q^2 in Berger-Sehgal Model
      double fTMax;   ///< upper bound for t = (q - p_pi)^2
  };

}       // genie namespace
#endif  // _COH_XSEC_H_
