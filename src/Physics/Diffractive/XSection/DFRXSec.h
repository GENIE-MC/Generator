//____________________________________________________________________________
/*!

\class    genie::DFRXSec

\brief    Computes the cross section for DFR neutrino-nucleus pi production.\n
          Is a concrete implementation of the XSecIntegratorI interface.

\author   Jeremy Wolcott <jeremy.wolcott \at tufts.edu>
          Tufts University

\created  May 13, 2016

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________
#ifndef _DFR_XSEC_H_
#define _DFR_XSEC_H_

#include <string>

#include "Physics/XSectionIntegration/XSecIntegratorI.h"

namespace genie
{
  class Interaction;
  class Registry;
  class XSecAlgorithmI;
}

namespace genie
{

  class DFRXSec : public XSecIntegratorI
  {
    public:
      DFRXSec ();
      DFRXSec (std::string config);
      virtual ~DFRXSec ();

      // XSecIntegratorI interface implementation
      double Integrate(const XSecAlgorithmI * model, const Interaction * i) const;

      // Overload the Algorithm::Configure() methods to load private data
      // members from configuration options
      void Configure(const Registry & config);
      void Configure(std::string config);

    private:
      void LoadConfig (void);

      double fTMax;   ///< upper bound for t = (q - p_pi)^2
  };

} /* namespace genie */

#endif /* _DFR_XSEC_H_ */
