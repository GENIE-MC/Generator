//____________________________________________________________________________
/*!

\class    genie::HEDISXSec

\brief    Computes the HEDIS Cross Section. \n
          Is a concrete implementation of the XSecIntegratorI interface. \n

\author   Alfonso Garcia <alfonsog \at nikhef.nl>
          NIKHEF

\created  August 28, 2019

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _HEDIS_XSEC_H_
#define _HEDIS_XSEC_H_

#include "Physics/XSectionIntegration/XSecIntegratorI.h"

#include <vector>

namespace genie {

class HEDISXSec : public XSecIntegratorI {

public:
  HEDISXSec();
  HEDISXSec(string config);
  virtual ~HEDISXSec();

  //! XSecIntegratorI interface implementation
  double Integrate(const XSecAlgorithmI * model, const Interaction * i) const;

  //! Overload the Algorithm::Configure() methods to load private data
  //! members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void   LoadConfig (void);

  double fSFXmin;   ///< minimum value of x for which SF tables are computed
  double fSFQ2min;  ///< minimum value of Q2 for which SF tables are computed
  double fSFQ2max;  ///< maximum value of Q2 for which SF tables are computed

};

}       // genie namespace
#endif  // _DIS_XSEC_H_
