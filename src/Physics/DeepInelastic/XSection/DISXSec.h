//____________________________________________________________________________
/*!

\class    genie::DISXSec

\brief    Computes the DIS Cross Section. \n
          Is a concrete implementation of the XSecIntegratorI interface. \n

\author   Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  May 04, 2004

\cpright  Copyright (c) 2003-2024, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _DIS_XSEC_H_
#define _DIS_XSEC_H_

#include "Physics/XSectionIntegration/XSecIntegratorI.h"

namespace genie {

class DISXSec : public XSecIntegratorI {

public:
  DISXSec();
  DISXSec(string config);
  virtual ~DISXSec();

  //! XSecIntegratorI interface implementation
  double Integrate(const XSecAlgorithmI * model, const Interaction * i) const;

  //! Overload the Algorithm::Configure() methods to load private data
  //! members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void   LoadConfig (void);

  void   CacheFreeNucleonXSec(const XSecAlgorithmI * model, const Interaction * in) const;
  string CacheBranchName     (const XSecAlgorithmI * model, const Interaction * in) const;

  double fVldEmin;
  double fVldEmax;
};

}       // genie namespace
#endif  // _DIS_XSEC_H_
