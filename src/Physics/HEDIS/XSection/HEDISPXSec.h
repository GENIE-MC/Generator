//____________________________________________________________________________
/*!

\class    genie::HEDISPXSec

\brief    Computes the double differential Cross Section for HEDIS. \n
          Is a concrete implementation of the XSecAlgorithmI interface.

\author   Alfonso Garcia <alfonsog \at nikhef.nl>
          NIKHEF

\created  August 28, 2019

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _HEDIS_PXSEC_H_
#define _HEDIS_PXSEC_H_

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Physics/HEDIS/XSection/HEDISStrucFunc.h"

namespace genie {

  class XSecIntegratorI;

  class HEDISPXSec : public XSecAlgorithmI {

    public:
      HEDISPXSec();
      HEDISPXSec(string config);
      virtual ~HEDISPXSec();

      // XSecAlgorithmI interface implementation
      double XSec            (const Interaction * i, KinePhaseSpace_t k) const;
      double Integral        (const Interaction * i) const;
      bool   ValidProcess    (const Interaction * i) const;

      // overload the Algorithm::Configure() methods to load private data
      // members from configuration options
      void Configure(const Registry & config);
      void Configure(string config);

    private:
      void   LoadConfig (void);
      double ds_dxdy      (SF_xQ2 sf, double x, double y) const;
      double ds_dxdy_mass (SF_xQ2 sf, double x, double y, double e, double mt, double ml2) const;

      const XSecIntegratorI *        fXSecIntegrator;     ///< diff. xsec integrator
    
      double fWmin;            ///< Minimum value of W
      bool fMassTerms;         ///< Account for second order effects in DDxsec
      SF_info fSFinfo;         ///< Information used to computed SF
  };

}       // genie namespace

#endif  // _HEDIS_PXSEC_H_
