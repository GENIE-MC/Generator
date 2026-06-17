//____________________________________________________________________________
/*!

\class    genie::DMRESXSec

\brief    Computes the DM RES Cross Section.\n
          Is a concrete implementation of the XSecIntegratorI interface.\n

\author   Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
          University of Liverpool & STFC Rutherford Appleton Laboratory

          updated for DM RES by Zach Orr, Colorado State University

\created  May 04, 2004

\cpright  Copyright (c) 2003-2023, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _DMRES_XSEC_H_
#define _DMRES_XSEC_H_

#include "Physics/XSectionIntegration/XSecIntegratorI.h"

namespace genie {

class DMRESXSec : public XSecIntegratorI {

public:
  DMRESXSec();
  DMRESXSec(string param_set);
  virtual ~DMRESXSec();

  //! XSecIntegratorI interface implementation
  double Integrate(const XSecAlgorithmI * model, const Interaction * i) const;

  //! Overload the Algorithm::Configure() methods to load private data
  //! members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void LoadConfig (void);
};

}       // genie namespace
#endif  // _DMRES_XSEC_H_
