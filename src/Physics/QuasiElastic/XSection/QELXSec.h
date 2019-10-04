//____________________________________________________________________________
/*!

\class    genie::QELXSec

\brief    Computes the Quasi Elastic (QEL) cross section. \n
          Is a concrete implementation of the XSecIntegratorI interface. \n

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 04, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _QEL_XSEC_H_
#define _QEL_XSEC_H_

#include "Physics/XSectionIntegration/XSecIntegratorI.h"

namespace genie {

class NuclearModelI;

class QELXSec : public XSecIntegratorI {

public:
  QELXSec();
  QELXSec(string config);
  virtual ~QELXSec();

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
#endif  // _QEL_XSEC_H_
