//____________________________________________________________________________
/*!

\class    genie::DMELXSec

\brief    Computes the Elastic dark matter (DMEL) cross section. \n
          Is a concrete implementation of the XSecIntegratorI interface. \n

\author   Joshua Berger <jberger \at physics.wisc.edu>
          University of Wisconsin-Madison
          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  September 4, 2017

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _DMEL_XSEC_H_
#define _DMEL_XSEC_H_

#include "Physics/XSectionIntegration/XSecIntegratorI.h"

namespace genie {

class NuclearModelI;

class DMELXSec : public XSecIntegratorI {

public:
  DMELXSec();
  DMELXSec(string config);
  virtual ~DMELXSec();

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
#endif  // _DMEL_XSEC_H_
